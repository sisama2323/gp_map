/*
Some Notes:
1. Only publish the point that are marked as occupied otherwise it is hard to see what is going on
2. Global variable Sig_mat is empty after passing into callback function, need to fix this
*/

// ros
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <message_filters/subscriber.h>
#include <tf/transform_listener.h>

// pcl
#include <pcl/PCLPointCloud2.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/filters/voxel_grid.h>

// eigen
#include <Eigen/Sparse>

// stl
#include <cmath>
#include <vector>
#include <limits>
#include <std_msgs/String.h>
#include <set>
#include <random>
#include <chrono>
#include <iostream>

// local
#include "markerarray_pub.h"
#include "InformationFilterOccupancyMap.h"


void QuatToRot( double qx, double qy, double qz, double qw,
               Eigen::Matrix<double,3,3>& R )
{
  R << 1.0-2.0*(qy*qy+qz*qz),
	     2.0*qx*qy-2.0*qw*qz,
	     2.0*qw*qy+2.0*qx*qz,
  	   2.0*qx*qy+2.0*qw*qz,
	     1.0-2.0*(qx*qx+qz*qz),
	     2.0*qy*qz-2.0*qw*qx,
  	   2.0*qx*qz-2.0*qw*qy,
	     2.0*qy*qz+2.0*qw*qx,
	     1.0-2.0*(qx*qx+qy*qy);
}


// save time
void WriteVectorToFile(const std::vector<double>& vector,std::string filename)
{
    std::ofstream ofs(filename,std::ios::out | std::ofstream::binary);
    std::ostream_iterator<double> osi{ofs, "\n"};
    std::copy(vector.begin(),vector.end(),osi);
}


class RunIfomFromBag{
    // parameter
    tf::Vector3 last_position;
    tf::Quaternion last_orientation;
    ros::NodeHandle nh_;
    double delta_;   // delta for information vector p(y=1)
    double thresh_;  // threshold for posterier
    //  map range in meter lower
    std::vector<double> devmin_;
    //  map range in meter upper
    std::vector<double> devmax_;
    size_t num_iter_;
    double free_step_size_;
    //  kernel sigma
    std::vector<double> kern_stdev_;
    // for tum bag
    std::string frame_id_;
    bool first;
    double position_change_thresh_;
    double orientation_change_thresh_;
    double sensor_max_range_;          // Maximum range of the sensor
    std::vector<double> map_mu;
    // generate random number
    std::default_random_engine random_engine;
    std::exponential_distribution<> exponential_distribution;
    //  map resolution
    std::vector<double> res_;
    ifom::InformationFilterOccupancyMap info_filter_map;
    //  tum bag store transformation in tf but simulation get pose directly from gazebo
    tf::TransformListener listener;
    // information vecoter map
    gpoctomap::MarkerArrayPub m_pub;
    // true posteior map
    gpoctomap::MarkerArrayPub final_map_pub;
    // information matrix map
    gpoctomap::MarkerArrayPub cov_map_pub;
    bool save_file_;
    std::vector<double> time_file;
    std::string file_name_;
    bool debug_;

public:
    // init
    RunIfomFromBag(ros::NodeHandle nh,
                   const double& sensor_max_range,          // Maximum range of the sensor
                   const std::vector<double>& devmin,
                   const std::vector<double>& devmax,
                   const std::vector<double>& resolution,
                   const int& iter,
                   const std::vector<double>& kern_stdev,
                   const double& noise_cov,   // noise covariacne
                   const double& thresh,  // threshold for posterier
                   const double& delta,   // delta for information vector p(y=1)
                   const double& free_step_size,
                   const double& kernel_radius,       // radius of pseudo point
                   const std::string& frame_id,
                   const bool& debug,
                   const bool& save_file,
                   const std::string& file_name)
        : exponential_distribution(0.4), info_filter_map(devmin, devmax, resolution, kern_stdev, noise_cov, kernel_radius),
            m_pub(nh, "/occ_map", 0.1f), final_map_pub(nh, "/full_map", 0.1f), cov_map_pub(nh, "/cov_map", 0.1f)
        {
        nh_ = nh;
        delta_ = delta;   // delta for information vector p(y=1)
        thresh_ = thresh;  // threshold for posterier
        devmin_ = devmin;  //  map range in meter lower
        devmax_ = devmax;  //  map range in meter upper
        num_iter_ = iter;
        free_step_size_ = free_step_size;  //  kernel sigma
        kern_stdev_ = kern_stdev;
        frame_id_ = frame_id;
        first = true;
        position_change_thresh_ = 0.1;
        orientation_change_thresh_ = 0.2;
        sensor_max_range_ = sensor_max_range;  // Maximum range of the sensor
        res_ = resolution;   //  map resolution

        ROS_INFO_STREAM("Initialization finished");

        // subscriber
        ros::Subscriber point_cloud_sub = nh_.subscribe<sensor_msgs::PointCloud2>("cloud", 1, &RunIfomFromBag::PointCloudCallback, this);
        // recover whole map signal
        ros::Subscriber sub = nh_.subscribe<std_msgs::String>("/chatter", 1000, &RunIfomFromBag::ChatterCallback, this);

        save_file_ = save_file;
        file_name_ = file_name;
        debug_ = debug;

        ros::spin();
    };


    // callback function
    void PointCloudCallback(const sensor_msgs::PointCloud2ConstPtr &cloud_msg_ptr)
    {
    //  laser pose
      std::vector<double> translation;
      Eigen::Matrix3d rotation;

      tf::StampedTransform transform;
      try {
          listener.waitForTransform(frame_id_, cloud_msg_ptr->header.frame_id, cloud_msg_ptr->header.stamp, ros::Duration(3.0));
          listener.lookupTransform(frame_id_, cloud_msg_ptr->header.frame_id, cloud_msg_ptr->header.stamp, transform);
        } catch (tf::TransformException ex) {
          ROS_ERROR("%s", ex.what());
          return;
      }
      tf::Vector3 position = transform.getOrigin();
      tf::Quaternion orientation = transform.getRotation();
      if (first || orientation.angleShortestPath(last_orientation) > orientation_change_thresh_ ||
            position.distance(last_position) > position_change_thresh_) {

          first = false;
          last_position = position;
          last_orientation = orientation;

          QuatToRot((double) orientation.x(),
                    (double) orientation.y(),
                    (double) orientation.z(),
                    (double) orientation.w(),
                    rotation);

          translation = {(double) position.x(),
                         (double) position.y(),
                         (double) position.z()};

          // Transform the points to the world frame and call bresenham
          // use PLC to convert point
          pcl::PCLPointCloud2 pcl_cloud;
          pcl_conversions::toPCL(*cloud_msg_ptr, pcl_cloud);
          pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_msg(new pcl::PointCloud<pcl::PointXYZ>());
          pcl::fromPCLPointCloud2(pcl_cloud, *pcl_msg);

          auto proc_time = tic();

          double sx = translation[0];
          double sy = translation[1];
          double sz = translation[2];

          //  start loop the points
          for (const auto &point : pcl_msg->points) {

              if ((!std::isnan(point.x)) && (!std::isnan(point.y)) && (!std::isnan(point.z))) {
                  // ROS_INFO_STREAM(point);

                  // point in camera frame
                  Eigen::Vector3d pt_c;
                  pt_c << point.x, point.y, point.z;

                  // point in world frame
                  Eigen::Vector3d pt_w = rotation * pt_c;
                  pt_w(0) += translation[0];
                  pt_w(1) += translation[1];
                  pt_w(2) += translation[2];

                  double ex = pt_w(0);
                  double ey = pt_w(1);
                  double ez = pt_w(2);

                  // Remove all the points out of the sensor max range.
                  double ray_length = sqrt(pow((sx - ex), 2) + pow((sy - ey), 2) + pow((sz - ez), 2));

                  if ( ray_length * ray_length < pow(sensor_max_range_, 2)) {
                      // add occupied point to observation, make sure this point is within the size of map
                      if ((pt_w(0) < devmax_[0]) && (pt_w(1) < devmax_[1]) && (pt_w(2) < devmax_[2]) &&
                          (pt_w(0) > devmin_[0]) && (pt_w(1) > devmin_[1]) && (pt_w(2) > devmin_[2])) {

                          // occupied point
                          info_filter_map.addObservation(pt_w, kern_stdev_, true, debug_);
                      }
                  } else {
                      // out of range point mark as free
                      if ((pt_w(0) < devmax_[0]) && (pt_w(1) < devmax_[1]) && (pt_w(2) < devmax_[2]) &&
                          (pt_w(0) > devmin_[0]) && (pt_w(1) > devmin_[1]) && (pt_w(2) > devmin_[2])) {

                          // free points
                          info_filter_map.addObservation(pt_w, kern_stdev_, false, debug_);
                      }
                  }

                  // clear free points
                  Eigen::Vector3d rb_p;
                  // find line in 3D    p = p0 + vt
                  Eigen::Vector3d v = rb_p - pt_w;
                  Eigen::Vector3d p0 = rb_p;

                  // determine which direction to travel
                  double factor;
                  if (sx < ex) {
                      factor = 1;
                  } else {
                      factor = -1;
                  }
                  // n number of free point needed based on the length of the ray
                  int n = static_cast<int> (floor(ray_length / free_step_size_) - 1);
                  double fx = sx;
                  double fy = sy;
                  double fz = sz;
                  auto pf_time = tic();

                  // travel in that direction
                  for (int i=0; i<=n; ++i){
                      if (factor * fx < factor * ex) {
                          if ((fx < devmax_[0]) && (fy < devmax_[1]) && (fz < devmax_[2]) &&
                              (fx > devmin_[0]) && (fy > devmin_[1]) && (fz > devmin_[2])) {

                              rb_p << fx, fy, fz;
                              // free point
                              info_filter_map.addObservation(rb_p, kern_stdev_, false, debug_);
                          }
                      }
                      // new point
                      fx = sx + exponential_distribution(random_engine);    // generate a random step size based on poisson distribution
                      double r = (fx - p0(0)) / v(0);
                      fy = p0(1) + v(1) * r;
                      fz = p0(2) + v(2) * r;
                  }
              }
          }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
          if (save_file_){
              time_file.push_back(toc(proc_time));
          }
          // publish information vector build map for visualization purpose
          ROS_INFO_STREAM("One Cloud finished in " << toc(proc_time) << " sec!");

          m_pub.clear();
          cov_map_pub.clear();

          std::vector<double> observation = info_filter_map.alpha;       // observation value + 1 occupied -1 free
          std::vector<double> num = info_filter_map.obs_num_;                // number of observation lay in this grid
          std::vector<std::map<int, double>> informat = info_filter_map.beta;

          // threshold information matrix directly
          for (int idx=0; idx<observation.size(); idx++) {

            if ( observation[idx] > delta_ ){
              std::vector<int> pt_idx = info_filter_map.converter.ind2subv_rowmajor(idx);

              float xx = pt_idx[0] * res_[0] + devmin_[0];
              float yy = pt_idx[1] * res_[1] + devmin_[1];
              float zz = pt_idx[2] * res_[2] + devmin_[2];
              m_pub.insert_point3d(xx, yy, zz, devmin_[2], devmax_[2], res_[0]);

              double cov = informat[idx][idx];

              cov_map_pub.insert_color_point3d(xx, yy, zz, 0, 150, cov, res_[0]);
            }  // end if loop
          }  // end publish map loop

          m_pub.publish();
          cov_map_pub.publish();
      }
    };


    // subscribe to keyboard, if say recovermap then this will start the map recover stage
    void ChatterCallback(const std_msgs::String::ConstPtr& msg)
    {
       std::string key_input = msg->data.c_str();

       ROS_INFO("! Got Message !");

       std::string first_string = key_input.substr(0, 10);

       auto const pos=key_input.find_last_of('p');

       if(first_string.compare("RecoverMap") == 0){

         auto recover_time = tic();
         ROS_INFO("Recover Map Now... / ...");

         if (map_mu.size() < 1){
             auto t1 = tic();
            map_mu = info_filter_map.latentMean(num_iter_);
            if (save_file_){
                time_file.push_back(toc(t1));
            }
            ROS_INFO_STREAM("map_mu " << toc(t1) << " sec!");
         }

         ROS_INFO_STREAM("Map Recovering computation done in " << toc(recover_time) << " sec!");

         // should only publish points that are occupied
         final_map_pub.clear();

         for (int idx = 0; idx < map_mu.size(); ++idx ){

           if (map_mu[idx] > thresh_){

             if (debug_) {
               ROS_INFO_STREAM("map_mu[idx] " << map_mu[idx]);
             }

             std::vector<int> pt_idx;
             pt_idx = info_filter_map.converter.ind2subv_rowmajor(idx);

             final_map_pub.insert_point3d(pt_idx[0]*res_[0]+devmin_[0],
                                          pt_idx[1]*res_[1]+devmin_[1],
                                          pt_idx[2]*res_[2]+devmin_[2],
                                          devmin_[2], devmax_[2], res_[2]);
           }
         }

         final_map_pub.publish();

         ROS_INFO("Finish Publishing Map");

         if (save_file_){
             WriteVectorToFile(time_file, file_name_);
             ROS_INFO("Finish Save Map");
         }
       }
    };
};


// main function
int main( int argc, char** argv)
{
  // initialize ros
  ros::init(argc, argv, "ifom");
  ros::NodeHandle nh;
  // get parameter
  double sensor_max_range;          // Maximum range of the sensor
  double xmin;
  double ymin;
  double zmin;
  double xmax;
  double ymax;
  double zmax;
  double xres;
  double yres;
  double zres;
  int iter;
  double kernel_sigma_x;   // kernel x width
  double kernel_sigma_y;   // kernel y width
  double kernel_sigma_z;   // kernel z width
  double noise_cov;   // noise covariacne
  double thresh;  // threshold for posterier
  double delta;   // delta for information vector p(y=1)
  double free_step_size;
  double kernel_radius;       // radius of pseudo point
  std::string frame_id("/map");
  bool debug = false;
  bool save_file;
  std::string file_name;

  nh.param<double>("sensor_max_range", sensor_max_range, 10);
  nh.param<double>("xmin", xmin, -10);
  nh.param<double>("ymin", ymin, -10);
  nh.param<double>("zmin", zmin, -10);
  nh.param<double>("xmax", xmax, 10);
  nh.param<double>("ymax", ymax, 10);
  nh.param<double>("zmax", zmax, 10);
  nh.param<double>("xres", xres, 1);
  nh.param<double>("yres", yres, 1);
  nh.param<double>("zres", zres, 1);
  nh.param<int>("iter", iter, 1);
  nh.param<double>("kernel_sigma_x", kernel_sigma_x, 0.75);
  nh.param<double>("kernel_sigma_y", kernel_sigma_y, 0.75);
  nh.param<double>("kernel_sigma_z", kernel_sigma_z, 0.75);
  nh.param<double>("noise_cov", noise_cov, 0.5);
  nh.param<double>("thresh", thresh, 0.5);
  nh.param<double>("delta", delta, 0.1);
  nh.param<double>("free_step_size", free_step_size, 0.1);
  nh.param<double>("kernel_radius", kernel_radius, 0.05);
  nh.param<std::string>("frame_id", frame_id, frame_id);
  nh.param<bool>("debug", debug, debug);
  nh.param<bool>("save_file", save_file, false);
  nh.param<std::string>("file_name", file_name, "/home/parallels/map_ws/plot/1.txt");

  std::vector<double> devmin = {xmin, ymin, zmin};
  std::vector<double> devmax = {xmax, ymax, zmax};
  std::vector<double> resolution = {xres, yres, zres};
  std::vector<double> kern_stdev = {kernel_sigma_x, kernel_sigma_y, kernel_sigma_z};

  // initialize class
  RunIfomFromBag RunIfomFromBag(nh, sensor_max_range,
                                devmin,
                                devmax,
                                resolution,
                                iter,
                                kern_stdev,
                                noise_cov, thresh, delta, free_step_size,
                                kernel_radius,
                                frame_id, debug, save_file, file_name);

  return 0;
}


