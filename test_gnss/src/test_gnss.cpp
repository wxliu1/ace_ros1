
// created by wxliu on 2024-8-23

#include <stdio.h>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <sensor_msgs/NavSatFix.h>
#include <gnss_comm/gnss_spp.hpp>

using namespace gnss_comm;
using std::queue;
using std::cerr;

queue<std::vector<ObsPtr>> gnss_meas_buf;

std::condition_variable con;

std::mutex m_buf;

double time_diff_gnss_local;
bool time_diff_valid;
double latest_gnss_time;

bool GNSS_ENABLE = true;
double GNSS_PSR_STD_THRES = 5.0;
double GNSS_DOPP_STD_THRES = 5.0;
uint32_t GNSS_TRACK_NUM_THRES =5;
double GNSS_ELEVATION_THRES=15;

std::vector<double> latest_gnss_iono_params;
Eigen::Matrix<double, 8, 1> latest_spp_status;
std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem;
std::map<uint32_t, std::map<double, size_t>> sat2time_index;
std::map<uint32_t, uint32_t> sat_track_status;

void initGnss()
{
    constexpr double GNSS_LOCAL_TIME_DIFF = 18.0;
    std::vector<double> GNSS_IONO_DEFAULT_PARAMS = { 0.1118E-07,  0.2235E-07, -0.4172E-06,  0.6557E-06, 0.1249E+06, -0.4424E+06,  0.1507E+07, -0.2621E+06};
    
    latest_gnss_time = -1;
    time_diff_valid = true;
    time_diff_gnss_local = GNSS_LOCAL_TIME_DIFF;

    // sat2ephem.clear();
    // sat2time_index.clear();
    sat_track_status.clear();
    latest_gnss_iono_params.clear();
    std::copy(GNSS_IONO_DEFAULT_PARAMS.begin(), GNSS_IONO_DEFAULT_PARAMS.end(), 
        std::back_inserter(latest_gnss_iono_params));
}

void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg)
{
    EphemPtr ephem_ptr = msg2ephem(ephem_msg);
    
    double toe = time2sec(ephem_ptr->toe);
    // if a new ephemeris comes
    if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
    {
        sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);
        sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
    }
}

void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg)
{
    std::vector<ObsPtr> gnss_meas = msg2meas(meas_msg);

    latest_gnss_time = time2sec(gnss_meas[0]->time);

    // cerr << "gnss ts is " << std::setprecision(20) << time2sec(gnss_meas[0]->time) << endl;
    if (!time_diff_valid)   return;

    m_buf.lock();
    gnss_meas_buf.push(std::move(gnss_meas));
    m_buf.unlock();
    con.notify_one();
}

// receive rtk for groudtruth.
void gnss_lla_callback(const sensor_msgs::NavSatFix::ConstPtr &nav_sat_msg)
{
    double t = nav_sat_msg->header.stamp.toSec();
    Eigen::Vector3d gt_geo, gt_ecef;
    gt_geo.x() = nav_sat_msg->latitude;
    gt_geo.y() = nav_sat_msg->longitude;
    gt_geo.z() = nav_sat_msg->altitude;
    if (gt_geo.hasNaN()) return;

    // convert lla to postion in ecef coordinate
    Eigen::Vector3d pos_ecef = gnss_comm::geo2ecef(gt_geo);
    if (pos_ecef.hasNaN() || pos_ecef.norm() < 1e-3)
        return;
    // ecef_rtks.emplace_back(t, pos_ecef);
    //  TODO: output rtk to file.
}

void processGNSS(const std::vector<ObsPtr> &gnss_meas)
{
    std::vector<ObsPtr> valid_meas;
    std::vector<EphemBasePtr> valid_ephems;
    for (auto obs : gnss_meas)
    {
        // filter according to system
        uint32_t sys = satsys(obs->sat, NULL);
        if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
            continue;

        // if not got cooresponding ephemeris yet
        if (sat2ephem.count(obs->sat) == 0)
            continue;
        
        if (obs->freqs.empty())    continue;       // no valid signal measurement
        int freq_idx = -1;
        L1_freq(obs, &freq_idx);
        if (freq_idx < 0)   continue;              // no L1 observation

        int count_low_srn = std::count_if(obs->CN0.begin(), obs->CN0.end(), [](double srn){return srn <= 30;});
        if(count_low_srn > 0 ) { continue; } 
        
        double obs_time = time2sec(obs->time);
        std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
        double ephem_time = EPH_VALID_SECONDS;
        size_t ephem_index = -1;
        for (auto ti : time2index)
        {
            if (std::abs(ti.first - obs_time) < ephem_time)
            {
                ephem_time = std::abs(ti.first - obs_time);
                ephem_index = ti.second;
            }
        }
        if (ephem_time >= EPH_VALID_SECONDS)
        {
            cerr << "ephemeris not valid anymore\n";
            continue;
        }
        const EphemBasePtr &best_ephem = sat2ephem.at(obs->sat).at(ephem_index);

        // filter by tracking status
        LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
        if (obs->psr_std[freq_idx]  > GNSS_PSR_STD_THRES ||
            obs->dopp_std[freq_idx] > GNSS_DOPP_STD_THRES)
        {
            sat_track_status[obs->sat] = 0;
            continue;
        }
        else
        {
            if (sat_track_status.count(obs->sat) == 0)
                sat_track_status[obs->sat] = 0;
            ++ sat_track_status[obs->sat];
        }
        if (sat_track_status[obs->sat] < GNSS_TRACK_NUM_THRES)
            continue;           // not being tracked for enough epochs

        // filter by elevation angle
        /*if (gnss_ready)
        {
            Eigen::Vector3d sat_ecef;
            if (sys == SYS_GLO)
                sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
            else
                sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
            double azel[2] = {0, M_PI/2.0};
            sat_azel(ecef_pos, sat_ecef, azel);
            if (azel[1] < GNSS_ELEVATION_THRES*M_PI/180.0)
                continue;                
        }*/
        valid_meas.push_back(obs);
        valid_ephems.push_back(best_ephem);
    }

    double curTime_gnss = valid_meas.empty() ? 0.0 : gnss_comm::time2sec(valid_meas.front()->time); 
    std::vector<size_t> inliners;
    auto spp_status = gnss_comm::psr_pos(valid_meas, valid_ephems, gnss_comm::getDefaultIonos(), inliners, true, false, latest_spp_status.block<7,1>(0,0));
    // gnss_comm::drawTopSatView(valid_meas, valid_ephems, inliners, ecef_pos, true, GNSS_MODULE);
    // gnss_comm::drawSrnView(valid_meas, inliners, true, GNSS_MODULE);
    if(!gnss_comm::valid_ecef(spp_status.block<3,1>(0,0))){
        valid_meas.clear();
        valid_ephems.clear();
        std::cout << "processGNSS psr_pos fail, clear.";
    }else{
        /*
        Eigen::Vector3d spp_local = ecef2local(spp_status.block<3,1>(0,0));
        Eigen::Vector3d last_spp_local = ecef2local(latest_spp_status.block<3,1>(0,0));
        bool bad_gnss = gnss_comm::valid_ecef(latest_spp_status.block<3,1>(0,0)) &&                                        // 前一状态合法
                curTime_gnss - latest_spp_status(7,0) > 0.0 && curTime_gnss - latest_spp_status(7,0) < 0.35 &&             // 前一状态在合理范围内
                ((spp_local - last_spp_local).block<2,1>(0,0).norm() > 5.0 || (spp_local - last_spp_local).norm() > 10.0); // 前一状态 距离偏移太大

        if(bad_gnss){
            valid_meas.clear();
            valid_ephems.clear();
            std::cout << "processGNSS psr_pos fine, but big gap with last, clear.";
        }else{
            for(int i=0; i<inliners.size(); ++i){
                if(inliners[i] == i) continue;
                valid_meas[i] = valid_meas[inliners[i]];
                valid_ephems[i] = valid_ephems[inliners[i]];
            }
            valid_meas.resize(inliners.size());
            valid_ephems.resize(inliners.size());
        }
        */
    }
    latest_spp_status << spp_status, curTime_gnss;

    if(!valid_meas.empty()) {
        // ecef_spps.emplace_back(time2sec(valid_meas.front()->time), spp_status.block<3,1>(0,0));
        // TODO: output spp to file.
    }
    
    /*if(valid_meas.size() > 6){
        add_opti_gnss(valid_meas, valid_ephems, spp_status.block<3,1>(0,0));
    }*/
    
    /*gnss_meas_buf[frame_count] = valid_meas;
    gnss_ephem_buf[frame_count] = valid_ephems;



    if(!valid_meas.empty()){
        std::cout << "[processGNSS] add buf = " << frame_count << "-" << valid_meas.size() << ",   [";
        for(int i=0; i<frame_count; ++i) std::cout << gnss_meas_buf[i].size() << ",";
        std::cout << "]" << std::endl;
    }else{
        std::cout << "[processGNSS] add buf empty" << std::endl;
    }*/
}

bool
getMeasurements(std::vector<ObsPtr> &gnss_msg)
{
    gnss_msg = gnss_meas_buf.front();
    gnss_meas_buf.pop();

    return true;
}

void process(void* pParam) {
    while(true)
    {
        std::vector<ObsPtr> gnss_msg;

        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
                 {
                    return getMeasurements(gnss_msg);
                 });
        lk.unlock();

        if (GNSS_ENABLE && !gnss_msg.empty())
            processGNSS(gnss_msg);
    }
}

int main(int argc, char* argv[])
{
    ros::init(argc, argv, "test_gnss");
    ros::NodeHandle n("~");

    initGnss();

    ros::Subscriber sub_ephem, sub_glo_ephem, sub_gnss_meas, sub_gnss_iono_params, sub_lla;

    std::string GNSS_EPHEM_TOPIC = "/FM_F20L/ephem";
    std::string GNSS_GLO_EPHEM_TOPIC;
    std::string GNSS_MEAS_TOPIC = "/FM_F20L/range_meas";
    std::string GNSS_IONO_PARAMS_TOPIC;

    sub_ephem = n.subscribe(GNSS_EPHEM_TOPIC, 100, gnss_ephem_callback);
    // sub_glo_ephem = n.subscribe(GNSS_GLO_EPHEM_TOPIC, 100, gnss_glo_ephem_callback);
    sub_gnss_meas = n.subscribe(GNSS_MEAS_TOPIC, 100, gnss_meas_callback);
    // sub_gnss_iono_params = n.subscribe(GNSS_IONO_PARAMS_TOPIC, 100, gnss_iono_params_callback);
    sub_lla = n.subscribe("/FM_F20L/receiver_lla", 100, gnss_lla_callback);

    std::thread measurement_process{process, nullptr};
    measurement_process.join();
    ros::spin();

    return 0;
}