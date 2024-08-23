/**
* This file is part of gnss_comm.
*
* Copyright (C) 2021 Aerial Robotics Group, Hong Kong University of Science and Technology
* Author: CAO Shaozu (shaozu.cao@gmail.com)
*
* gnss_comm is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* gnss_comm is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with gnss_comm. If not, see <http://www.gnu.org/licenses/>.
*/

#include "gnss_spp.hpp"
#include "gnss_utility.hpp"
#include <glog/logging.h>

#define CUT_OFF_DEGREE 15.0
#define MIN_SRN 30.0

namespace gnss_comm
{
    void filter_L1(const std::vector<ObsPtr> &obs, const std::vector<EphemBasePtr> &ephems, 
        std::vector<ObsPtr> &L1_obs, std::vector<EphemBasePtr> &L1_ephems)
    {
        L1_obs.clear();
        L1_ephems.clear();
        for (size_t i = 0; i < obs.size(); ++i)
        {
            const ObsPtr &this_obs = obs[i];
            // check system
            const uint32_t sys = satsys(this_obs->sat, NULL);
            if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_BDS && sys != SYS_GAL)
                continue;
            // check signal frequency
            const double obs_freq = L1_freq(this_obs, NULL);
            if (obs_freq < 0)   continue;

            L1_obs.push_back(this_obs);
            L1_ephems.push_back(ephems[i]);
        }
    }

    std::vector<SatStatePtr> sat_states(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems)
    {
        std::vector<SatStatePtr> all_sv_states;
        const uint32_t num_obs = obs.size();
        for (size_t i = 0; i < num_obs; ++i)
        {
            SatStatePtr sat_state(new SatState());
            all_sv_states.push_back(sat_state);

            ObsPtr this_obs = obs[i];
            const uint32_t sat = this_obs->sat;
            const uint32_t sys = satsys(sat, NULL);
            int l1_idx = -1;
            L1_freq(this_obs, &l1_idx);
            if (l1_idx < 0)   continue;
            double tof = this_obs->psr[l1_idx] / LIGHT_SPEED;

            gtime_t sv_tx = time_add(this_obs->time, -tof);
            double svdt = 0, svddt = 0;
            Eigen::Vector3d sv_pos = Eigen::Vector3d::Zero();
            Eigen::Vector3d sv_vel = Eigen::Vector3d::Zero();
            if (sys == SYS_GLO)
            {
                GloEphemPtr glo_ephem = std::dynamic_pointer_cast<GloEphem>(ephems[i]);
                svdt = geph2svdt(sv_tx, glo_ephem);
                sv_tx = time_add(sv_tx, -svdt);
                sv_pos = geph2pos(sv_tx, glo_ephem, &svdt);
                sv_vel = geph2vel(sv_tx, glo_ephem, &svddt);
            }
            else
            {
                EphemPtr ephem = std::dynamic_pointer_cast<Ephem>(ephems[i]);
                svdt = eph2svdt(sv_tx, ephem);
                sv_tx = time_add(sv_tx, -svdt);
                sv_pos = eph2pos(sv_tx, ephem, &svdt);
                sv_vel = eph2vel(sv_tx, ephem, &svddt);
                sat_state->tgd = ephem->tgd[0];
            }
            sat_state->sat_id = sat;
            sat_state->ttx    = sv_tx;
            sat_state->pos    = sv_pos;
            sat_state->vel    = sv_vel;
            sat_state->dt     = svdt;
            sat_state->ddt    = svddt;
        }
        return all_sv_states;
    }

    void psr_res(const Eigen::Matrix<double, 7, 1> &rcv_state, const std::vector<ObsPtr> &obs, 
        const std::vector<SatStatePtr> &all_sv_states, const std::map<uint32_t, std::vector<double>> &iono_params, bool useDualFreq,
        Eigen::VectorXd &res, Eigen::MatrixXd &J, std::vector<Eigen::Vector2d> &atmos_delay, 
        std::vector<Eigen::Vector2d> &all_sv_azel)
    {
        const uint32_t num_sv = all_sv_states.size();
        // clear output
        res = Eigen::VectorXd::Zero(num_sv);
        J = Eigen::MatrixXd::Zero(num_sv, 7);
        atmos_delay.resize(num_sv);
        all_sv_azel.resize(num_sv);

        for (uint32_t i = 0; i < num_sv; ++i)
        {
            int l1_idx = -1;
            L1_freq(obs[i], &l1_idx);
            if (l1_idx < 0)   continue;

            const SatStatePtr &sat_state = all_sv_states[i];
            uint32_t this_sys = satsys(sat_state->sat_id, NULL);
            Eigen::Vector3d sv_pos = sat_state->pos;

            Eigen::Vector3d rv2sv = sv_pos - rcv_state.topLeftCorner<3,1>();
            Eigen::Vector3d unit_rv2sv = rv2sv.normalized();
            double sagnac_term = EARTH_OMG_GPS*(sv_pos(0)*rcv_state(1,0)- sv_pos(1)*rcv_state(0,0))/LIGHT_SPEED;
            
            double ion_delay=0, ion_delay_f1=0, ion_delay_f2=0, tro_delay=0;
            double azel[2] = {0, M_PI/2.0};
            if (rcv_state.topLeftCorner<3,1>().norm() > 0)
            {
                sat_azel(rcv_state.topLeftCorner<3,1>(), sv_pos, azel);
                Eigen::Vector3d rcv_lla = ecef2geo(rcv_state.head<3>());
                // use satellite signal transmit time instead
                tro_delay = calculate_trop_delay(sat_state->ttx, rcv_lla, azel);
                const auto& iono = iono_params.count(this_sys) ? iono_params.at(this_sys) : std::vector<double>(8,0.0);
                ion_delay_f1 = calculate_ion_delay(sat_state->ttx, iono, rcv_lla, azel);
                if(obs[i]->freqs.size() >= 2 && std::abs(obs[i]->freqs[0] - obs[i]->freqs[1]) > 1){
                    double delay = sagnac_term + rcv_state(3+sys2idx.at(this_sys)) - sat_state->dt*LIGHT_SPEED + tro_delay + ion_delay + sat_state->tgd*LIGHT_SPEED;
                    const double f02 = obs[i]->freqs[0]*obs[i]->freqs[0];
                    const double f12 = obs[i]->freqs[1]*obs[i]->freqs[1];
                    double psr0 = obs[i]->psr[0] - delay;
                    double psr1 = obs[i]->psr[1] - delay;
                    double psr_f2 = (f02*psr0 - f12*psr1) / (f02 - f12);
                    ion_delay_f2 = obs[i]->psr[l1_idx] - delay - psr_f2;
                }else{
                    ion_delay_f2 = ion_delay_f1;
                }
                ion_delay = useDualFreq ? (ion_delay_f2 < 0.0 ? ion_delay_f1 : ion_delay_f2) : ion_delay_f1;
            }

            double psr_estimated = rv2sv.norm() + sagnac_term + rcv_state(3+sys2idx.at(this_sys)) -
             sat_state->dt*LIGHT_SPEED + tro_delay + ion_delay + sat_state->tgd*LIGHT_SPEED;

            J.block(i, 0, 1, 3) = -unit_rv2sv.transpose();
            J(i, 3+sys2idx.at(this_sys)) = 1.0;
            res(i) = psr_estimated - obs[i]->psr[l1_idx];
            // if(std::abs(res(i)) < 1e-2){
            //     std::cout << std::setprecision(8) << std::fixed << "res err = 0 = " << res(i) << ", psr_es=" <<  psr_estimated << ", psr_ms=" << obs[i]->psr[l1_idx] << std::endl;
            // }

            atmos_delay[i] = Eigen::Vector2d(ion_delay_f1, ion_delay_f2);
            all_sv_azel[i] = Eigen::Vector2d(azel[0], azel[1]);
        }
    }

    std::vector<size_t> filterGoodMeas(const std::vector<ObsPtr> &valid_obs, const std::vector<Eigen::Vector2d>& all_sv_azel, 
        const Eigen::MatrixXd& G, const Eigen::VectorXd& b, const double thresh_psr_residual, const std::set<int> excludeSats){
        std::vector<size_t> good_idx;
        assert(valid_obs.size() == G.rows());
        for (uint32_t i = 0; i < valid_obs.size(); ++i)
        {
            if(G.row(i).norm() <= 0)   continue;                         // res not computed
            if(excludeSats.count(valid_obs[i]->sat) > 0) continue;       // 没有选择的卫星
            if(all_sv_azel[i].y() < CUT_OFF_DEGREE/180.0*M_PI) continue; // 高度角太小
            if(std::abs(b[i]) > thresh_psr_residual) continue;           // 伪距残差太大
            int count_low_srn = std::count_if(valid_obs[i]->CN0.begin(), valid_obs[i]->CN0.end(), [](double srn){return srn <= MIN_SRN;});
            if(count_low_srn > 0) continue;                              // SRN太小
            double diffSrn = 0.0; if(valid_obs[i]->CN0.size() >= 2) diffSrn = std::abs(valid_obs[i]->CN0[0] - valid_obs[i]->CN0[1]);
            if(diffSrn > 10.0) continue;                                 // 双频信噪比差异太大
            
            good_idx.push_back(i);
        }
        return good_idx;
    }

    Eigen::Matrix<double, 7, 1> psr_pos(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems, const std::map<uint32_t, std::vector<double>> &iono_params, std::vector<size_t>& good_idx, bool isFinePose,
        bool useDualFreq, const Eigen::Matrix<double, 7, 1>& last_state, const std::set<int> excludeSats, const std::string& ins_name)
    {
        Eigen::Matrix<double, 7, 1> result = Eigen::Matrix<double, 7, 1>::Zero();

        std::vector<ObsPtr> valid_obs;
        std::vector<EphemBasePtr> valid_ephems;
        filter_L1(obs, ephems, valid_obs, valid_ephems);
        if (valid_obs.size() < 6)
        {
            LOG(ERROR) << "[gnss_comm::psr_pos] GNSS observation not enough.\n";
            return result;
        }

        std::vector<SatStatePtr> all_sat_states = sat_states(valid_obs, valid_ephems);
        // for(int i=0; i<all_sat_states.size(); ++i){
        //     printf("%.3f sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f tgd=%.3f\n",
        //       time2sec(obs[0]->time),valid_obs[i]->sat,all_sat_states[i]->pos[0],all_sat_states[i]->pos[1],all_sat_states[i]->pos[2],
        //       all_sat_states[i]->dt*1E9,all_sat_states[i]->tgd*1E9);
        // }

        Eigen::VectorXd b;
        Eigen::MatrixXd G;
        std::vector<Eigen::Vector2d> atmos_delay;
        std::vector<Eigen::Vector2d> all_sv_azel;

        Eigen::Matrix<double, 7, 1> xyzt, Qxx; xyzt.setZero(); Qxx.setZero();
        bool flagExtractLargeRes = false;
        std::vector<size_t> last_good_idx; for(int i=0; i<valid_obs.size(); ++i) last_good_idx.emplace_back(i);
        if(valid_ecef(last_state.head<3>())) {
            psr_res(last_state, valid_obs, all_sat_states, iono_params, useDualFreq, b, G, atmos_delay, all_sv_azel);
            std::vector<size_t> inliners = filterGoodMeas(valid_obs, all_sv_azel, G, b, 50, excludeSats);
            if(inliners.size() >= 10){
                flagExtractLargeRes = true;
                xyzt = last_state;
                last_good_idx = inliners;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        double dx_norm = 1.0, ave_res=0.0, std_res=0.0, thesh_psr_res=0.0;
        auto huber_weight = [&std_res, &ave_res](double r){ double w = std::max(std_res,0.5)/std::abs(r); return 1.0; return w > 1.0 ? 1.0 : w;};
        uint32_t num_iter = 0;
        while(num_iter < MAX_ITER_PVT && dx_norm > EPSILON_PVT)
        {
            psr_res(xyzt, valid_obs, all_sat_states, iono_params, useDualFreq, b, G, atmos_delay, all_sv_azel);
            {
                // 统计 residual 的 ave、std、thresh.
                Eigen::VectorXd b_good(last_good_idx.size()); for(int i=0; i<last_good_idx.size(); ++i) b_good[i] = b[last_good_idx[i]];
                double sum = 0; for(int i=0; i<b_good.size(); ++i) sum+=(b_good[i]-ave_res)*(b_good[i]-ave_res); 
                ave_res = b_good.mean();
                std_res = std::sqrt(sum/(b_good.size() * (b_good.size()-1)));
                thesh_psr_res = std::max(5.0*std_res + ave_res, isFinePose ? 5.0 : 500.0);
            }

            good_idx = filterGoodMeas(valid_obs, all_sv_azel, G, b, thesh_psr_res, excludeSats);
            if(good_idx.size() < 6) {
                LOG(ERROR) << "[gnss_comm::psr_pos] GNSS inliner observation not enough.\n";
                xyzt.setZero(); break;
            }

            int sys_mask[4] = {0, 0, 0, 0};
            for (uint32_t i = 0; i < good_idx.size(); ++i)
            {
                const uint32_t obs_sys = satsys(valid_obs[good_idx[i]]->sat, NULL);
                sys_mask[sys2idx.at(obs_sys)] = 1;
            }
            uint32_t num_extra_constraint = 4;
            for (uint32_t k = 0; k < 4; ++k)    num_extra_constraint -= sys_mask[k];
            LOG_IF(FATAL, num_extra_constraint >= 4) << "[gnss_comm::psr_pos] too many extra-clock constraints.\n";

            const uint32_t good_num = good_idx.size();
            Eigen::MatrixXd good_G(good_num+num_extra_constraint, 7);
            Eigen::VectorXd good_b(good_num+num_extra_constraint);
            Eigen::MatrixXd good_W(good_num+num_extra_constraint, good_num+num_extra_constraint);
            good_W.setZero();
            for (uint32_t i = 0; i < good_num; ++i)
            {
                const uint32_t origin_idx = good_idx[i];
                const ObsPtr &this_obs = valid_obs[origin_idx];
                good_G.row(i) = G.row(origin_idx);
                good_b(i) = b(origin_idx);
                // compose weight
                const double sin_el = sin(all_sv_azel[origin_idx].y());
                double weight = sin_el*sin_el;
                weight *= pow(huber_weight(good_b(i)),2);
                good_W(i, i) = weight;
            }
            uint32_t tmp_count = good_num;
            // add extra pseudo measurement to contraint unobservable clock bias
            for (size_t k = 0; k < 4; ++k)
            {
                if (!sys_mask[k])
                {
                    good_G.row(tmp_count).setZero();
                    good_G(tmp_count, k+3) = 1.0;
                    good_b(tmp_count) = 0;
                    good_W(tmp_count, tmp_count) = 1000;       // large weight
                    ++tmp_count;
                }
            }

            // ready for solving
            Eigen::MatrixXd Q = (good_G.transpose()*good_W*good_G).inverse();
            Eigen::VectorXd dx = -Q * good_G.transpose() * good_W * good_b;
            Qxx = Q.diagonal();

            if(!dx.hasNaN()) {
                dx_norm = dx.norm();
                xyzt += dx;
            }else{
                std::cout << " dx has nan. it=" << num_iter << std::endl;
                std::cout << "G = " << good_G << std::endl;
                std::cout << "b = " << good_b.transpose() <<  std::endl;
                std::cout << "Q = " << Q.transpose() <<  std::endl;
                std::cout << "xyzt = " << xyzt.transpose() << std::endl;
                std::cout << "last xyzt = " << last_state.transpose() << std::endl;
                break;
            }
            // printf("t=%.3f it=%d, nv=%d good=%d pos=[%.3f, %.3f, %.3f] init=[%.3f, %.3f, %.3f]\n", time2sec(obs[0]->time),
            //     num_iter, valid_obs.size(), good_num, xyzt[0], xyzt[1], xyzt[2], last_state[0], last_state[1], last_state[2]);

            if(dx.block<3,1>(0,0).norm() < 1.0) flagExtractLargeRes = true;
            last_good_idx = good_idx;
            ++num_iter;
        }
        if (num_iter == MAX_ITER_PVT)
        {
            LOG(WARNING) << "[gnss_comm::psr_pos] XYZT solver reached maximum iterations.\n";
        }
        auto t2 = std::chrono::high_resolution_clock::now();

        bool isFine = true;
        {
            isFine &= std::abs(ave_res) < 1.0 && std_res < 2.0 && std::abs(ave_res) > 1e-6 && std_res > 1e-6 && Qxx.head<3>().norm() < 100.0 && num_iter < MAX_ITER_PVT; // 误差层是OK的
            isFine &= flagExtractLargeRes && valid_ecef(xyzt.head<3>()); // 结果层是合理的结果
            isFine &= good_idx.size() >= (useDualFreq ? 8 : 10);         // 要有一定的内点个数
            if(!useDualFreq && (1.0 * good_idx.size()/valid_obs.size() < 0.7) && Qxx.head<3>().norm() > 49) isFine = false;   // 当非双频 && 内点比例小 && 误差大时，结果不可靠。
        }

        if(false){
            for(int i=0; i<valid_obs.size(); ++i){
                bool bgood = false; for(auto& g : good_idx) {if(g==i) {bgood = true; break;} }
                // if(obs[i]->freqs.size()==1) continue;
                int l1_idx = -1; L1_freq(obs[i], &l1_idx);
                // printf("%s, it=%d, sat=%d, %s, good=%d, azel=[%.1f,%.1f], res=[%.3f,%.3f], ion=[%.2f,%.2f], srn=[%.1f,%.1f], freq=[%d,%.2f,%.2f], psr=[%.3f,%.3f,%.3f], rs=[%.3f,%.3f,%.3f], ave_res=%.2f, std_res=%.3f\n", 
                //         ins_name.c_str(), num_iter, valid_obs[i]->sat, sat2str(valid_obs[i]->sat).c_str(), bgood, 
                //         all_sv_azel[i][0]*R2D, all_sv_azel[i][1]*R2D, b[i], b[i]-atmos_delay[i][1]+atmos_delay[i][0], atmos_delay[i][0], atmos_delay[i][1],
                //         valid_obs[i]->CN0[l1_idx], (obs[i]->freqs.size()==1 ? obs[i]->CN0[l1_idx] : obs[i]->CN0[l1_idx==0 ? 1 : 0]),
                //         obs[i]->freqs.size(), obs[i]->freqs[0], (obs[i]->freqs.size()==1 ? obs[i]->freqs[0] : obs[i]->freqs[1]),
                //         obs[i]->psr[0], (obs[i]->psr.size()==1 ? obs[i]->psr[0] : obs[i]->psr[1]), (obs[i]->psr.size()==1 ? 0.0 : obs[i]->psr[0] - obs[i]->psr[1]),
                //         all_sat_states[i]->pos[0], all_sat_states[i]->pos[1], all_sat_states[i]->pos[2], ave_res, std_res); 
            }

            double ave_srn = 0.0, ave_srn_good = 0.0, ave_srn_remove = 0.0;
            for(auto& m : obs){ int l1_idx = -1; L1_freq(m, &l1_idx); ave_srn+= m->CN0[l1_idx];}
            for(auto& i : good_idx) {int l1_idx = -1; L1_freq(obs[i], &l1_idx); ave_srn_good+= obs[i]->CN0[l1_idx];}
            ave_srn_remove = (ave_srn - ave_srn_good)/(obs.size() - good_idx.size()); ave_srn /= obs.size(); ave_srn_good /= good_idx.size();

            int count_fine_dual = std::count_if(good_idx.begin(), good_idx.end(), [&atmos_delay](size_t idx){return atmos_delay[idx][1] > 0;});
            double time_all = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() * 0.001;
            printf("%s, %d, %.3f, %.3f, %.3f, %d, %d, %.3f, %d, %d, %d, %d, %.2f, %.3f, %.3f, srn, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f,,,,\n", 
                ins_name.c_str(), useDualFreq, xyzt[0], xyzt[1], xyzt[2], isFine, num_iter, time_all,
                valid_obs.size(), good_idx.size(), count_fine_dual,
                valid_obs.size()-good_idx.size(), 100.0*good_idx.size()/valid_obs.size(),
                ave_res, std_res, ave_srn, ave_srn_good, ave_srn_remove, Qxx[0], Qxx[1], Qxx[2]);
        }
        
        if(isFine) result = xyzt;
        else result.setZero();
        return result;
    }

    void dopp_res(const Eigen::Matrix<double, 4, 1> &rcv_state, const Eigen::Vector3d &rcv_ecef,
                  const std::vector<ObsPtr> &obs, const std::vector<SatStatePtr> &all_sv_states, 
                  Eigen::VectorXd &res, Eigen::MatrixXd &J)
    {
        const uint32_t num_sv = all_sv_states.size();
        //clear output
        res = Eigen::VectorXd::Zero(num_sv);
        J = Eigen::MatrixXd::Zero(num_sv, 4);

        for (uint32_t i = 0; i < num_sv; ++i)
        {
            const SatStatePtr &sat_state = all_sv_states[i];
            Eigen::Vector3d unit_rv2sv = (sat_state->pos - rcv_ecef).normalized();
            double sagnac_term = EARTH_OMG_GPS/LIGHT_SPEED*(
                sat_state->vel(0)*rcv_ecef(1)+ sat_state->pos(0)*rcv_state(1,0) - 
                sat_state->vel(1)*rcv_ecef(0) - sat_state->pos(1)*rcv_state(0,0));
            double dopp_estimated = (sat_state->vel - rcv_state.topLeftCorner<3, 1>()).dot(unit_rv2sv) + 
                    rcv_state(3, 0) + sagnac_term - sat_state->ddt*LIGHT_SPEED;
            int l1_idx = -1;
            const double obs_freq = L1_freq(obs[i], &l1_idx);
            if (obs_freq < 0)   continue;
            const double wavelength = LIGHT_SPEED / obs_freq;
            res(i) = dopp_estimated + obs[i]->dopp[l1_idx]*wavelength;
            J.block(i, 0, 1, 3) = -1.0 * unit_rv2sv.transpose();
            J(i, 3) = 1.0;
        }
    }

    Eigen::Matrix<double, 4, 1> dopp_vel(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems, Eigen::Vector3d &ref_ecef)
    {
        // LOG(INFO) << "try to solve Doppler velocity.";
        Eigen::Matrix<double, 4, 1> result;
        result.setZero();

        if (ref_ecef.norm() == 0)
        {
            // reference point not given, try calculate using pseudorange
            auto inliners = std::vector<size_t>();
            auto psr_result = psr_pos(obs, ephems, std::map<uint32_t, std::vector<double>>(), inliners);
            if (psr_result.head<3>().norm() != 0)
            {
                ref_ecef = psr_result.head<3>();
            }
            else
            {
                LOG(ERROR) << "[gnss_comm::dopp_vel] Unable to initialize reference position for Doppler calculation.\n";
                return result;
            }
        }

        std::vector<ObsPtr> valid_obs;
        std::vector<EphemBasePtr> valid_ephems;
        filter_L1(obs, ephems, valid_obs, valid_ephems);

        if (valid_obs.size() < 4)
        {
            LOG(ERROR) << "[gnss_comm::dopp_vel] GNSS observation not enough for velocity calculation.\n";
            return result;
        }

        std::vector<SatStatePtr> all_sat_states = sat_states(valid_obs, valid_ephems);
        // compute azel for all satellite
        std::vector<Eigen::Vector2d> all_sv_azel;
        for (uint32_t i = 0; i < all_sat_states.size(); ++i)
        {
            double azel[2] = {0, 0};
            sat_azel(ref_ecef, all_sat_states[i]->pos, azel);
            all_sv_azel.emplace_back(azel[0], azel[1]);
        }

        Eigen::Matrix<double, 4, 1> xyzt_dot;
        xyzt_dot.setZero();
        double dx_norm = 1.0;
        uint32_t num_iter = 0;
        while(num_iter < MAX_ITER_PVT && dx_norm > EPSILON_PVT)
        {
            Eigen::MatrixXd G;
            Eigen::VectorXd b;
            dopp_res(xyzt_dot, ref_ecef, valid_obs, all_sat_states, b, G);
            
            std::vector<uint32_t> good_idx;
            for (uint32_t i = 0; i < valid_obs.size(); ++i)
            {
                if (G.row(i).norm() <= 0)   continue;       // res not computed
                if (all_sv_azel[i].y() > CUT_OFF_DEGREE/180.0*M_PI)  
                    good_idx.push_back(i);
            }
            const uint32_t good_num = good_idx.size();
            Eigen::MatrixXd good_G = Eigen::MatrixXd::Zero(good_num, 4);
            Eigen::VectorXd good_b = Eigen::VectorXd::Zero(good_num);
            Eigen::MatrixXd good_W = Eigen::MatrixXd::Zero(good_num, good_num);
            for (uint32_t i = 0; i < good_num; ++i)
            {
                const uint32_t origin_idx = good_idx[i];
                good_G.row(i) = G.row(origin_idx);
                good_b(i) = b(origin_idx);
                const double sin_el = sin(all_sv_azel[origin_idx].y());
                double weight = sin_el*sin_el;
                const ObsPtr &this_obs = valid_obs[origin_idx];
                int l1_idx = -1;
                L1_freq(this_obs, &l1_idx);
                LOG_IF(FATAL, l1_idx < 0) << "[gnss_comm::dopp_vel] no L1 observation found.\n";
                if (this_obs->dopp_std[l1_idx] > 0)
                    weight /= (this_obs->dopp_std[l1_idx]/0.256);

                const uint32_t obs_sys = satsys(this_obs->sat, NULL);
                if (obs_sys == SYS_GPS || obs_sys == SYS_BDS)
                    weight /= valid_ephems[origin_idx]->ura-1;
                else if (obs_sys == SYS_GAL)
                    weight /= valid_ephems[origin_idx]->ura-2;
                else if (obs_sys == SYS_GLO)
                    weight /= 2;
                good_W(i, i) = weight;
            }

            // solve
            Eigen::VectorXd dx = -(good_G.transpose()*good_W*good_G).inverse() * good_G.transpose() * good_W * good_b;
            dx_norm = dx.norm();
            xyzt_dot += dx;

            // dx_norm = 0;
            // LOG(INFO) << "cov is \n" << (G.transpose()*W*G).inverse();
            ++num_iter;
        }

        if (num_iter == MAX_ITER_PVT)
            LOG(WARNING) << "[gnss_comm::dopp_vel] XYZT solver reached maximum iterations.\n";
        
        result = xyzt_dot;
        return result;
    }
}