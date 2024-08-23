# This message contains one-epoch measurements from one single satellite

GnssTimeMsg time    # measurement time
uint32 sat          # satellite ID (define by `sat_no` function in `gnss_utility.hpp`)
float64[] freqs     # 观测频率 [Hz]
float64[] CN0       # 载噪密度比（信号强度） [dB-Hz]
uint8[] LLI         # lost-lock indicator (1=lost)
uint8[] code        # 通道代码
float64[] psr       # 伪距测量 [m]
float64[] psr_std   # 伪距标准差 [m]
float64[] cp        # 载波相位测量 [cycle]
float64[] cp_std    # 载波相位标准差 [cycle]
float64[] dopp      # Doppler measurement [Hz]
float64[] dopp_std  # 多普勒标准差 [Hz]
# 跟踪状态。bit_0：psr 有效，bit_1：cp 有效，
# bit_2：半 cp 有效，bit_3：半 cp 减去
uint8[]   status    