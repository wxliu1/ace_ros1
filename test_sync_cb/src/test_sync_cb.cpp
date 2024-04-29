
#include <iostream>

#include <sensor_msgs/Image.h>
#include <sensor_msgs/PointCloud2.h>
#include <message_filters/subscriber.h>
// #include <message_filters/sync_policies/approximate_time.h>
//#include <message_filters/sync_policies/exact_time.h>
#include <message_filters/time_synchronizer.h>
// #include <cv_bridge/cv_bridge.h>


// namespace cb = cv_bridge;
namespace sm = sensor_msgs;
// namespace gm = geometry_msgs;
namespace mf = message_filters;

struct TestSyncCb {
public:
    mf::Subscriber<sm::Image> sub_image0_;
    mf::Subscriber<sm::PointCloud2> sub_point2_;

    using SyncStereo = mf::TimeSynchronizer<sensor_msgs::Image, sensor_msgs::PointCloud2>;
    std::optional<SyncStereo> sync_stereo_;

    explicit TestSyncCb(/*const */ros::NodeHandle& pnh)
    {
        sub_image0_.subscribe(pnh, "/my_img", 5);
        sub_point2_.subscribe(pnh, "/my_point", 5);

        sync_stereo_.emplace(sub_image0_, sub_point2_, 5);
        sync_stereo_->registerCallback(&TestSyncCb::StereoCb, this);
    }

    void StereoCb(const sm::ImageConstPtr& image0_ptr, const sm::PointCloud2ConstPtr& point2_ptr)
    {
        static int i = 0;
        i++;
        std::cout << "i=" << i << std::endl;
    }
};


int main(int argc, char* argv[])
{
    ros::init(argc, argv, "node_test_sync_cb");  // node name

    ros::MultiThreadedSpinner spinner(6);  // use 6 threads

    ros::NodeHandle n("~");

    TestSyncCb test(n);


    // ros::spin();
    spinner.spin();

    return 0;
}