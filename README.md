### Dependencies
1. catkin: `sudo apt install python-catkin-tools`
2. octomap and teleop: `sudo apt install ros-melodic-octomap-ros ros-melodic-teleop-twist-joy`
3. catkin simple: `git clone https://github.com/catkin/catkin_simple.git`
4. download data in the gp_map/data folder:
  + Simulation data: https://drive.google.com/file/d/1J0Kzdzu7LQ28e7qOVWhR3hqLXAY_GXNo/view?usp=sharing
  + TUM data: https://vision.in.tum.de/data/datasets/rgbd-dataset
     - `wget https://vision.in.tum.de/rgbd/dataset/freiburg3/rgbd_dataset_freiburg3_cabinet.bag`
     - `wget https://vision.in.tum.de/rgbd/dataset/freiburg3/rgbd_dataset_freiburg3_large_cabinet.bag`
     - `wget https://vision.in.tum.de/rgbd/dataset/freiburg1/rgbd_dataset_freiburg1_teddy.bag`


### Demo (Simulation):
1. Run bag: `roslaunch gp_map play_tum_bag.launch rosbag_path:=$(rospack find gp_map)/data/simulation_2`
2. Run ifom: `roslaunch gp_map mapping_3d_simulation.launch`

### Demo (TUM):
1. Run bag: `roslaunch gp_map play_tum_bag.launch rosbag_path:=$(rospack find gp_map)/data/rgbd_dataset_freiburg3_cabinet`
2. Run ifom: `roslaunch gp_map mapping_3d.launch`

### Package Useage
1. config: rviz config for tum and simulation
2. gp_map: utils, how to communicate with rosbag and read laser and pose messages

### Run from a ros bag
1.  Run bag:
    *  Simulation Bag:  
        `roslaunch gp_map play_tum_bag.launch rosbag_path:=/home/parallels/map2_ws/src/ros_bag/simulation_2`
    *  TUM Bag:  
        `roslaunch gp_map play_tum_bag.launch rosbag_path:=/home/parallels/map2_ws/src/ros_bag/rgbd_dataset_freiburg3_large_cabinet`
2.  Run IFOM:
    *  For simulation:  
     `roslaunch gp_map mapping_3d.launch param_file:=$(rospack find gp_map)/param/simulation.yaml`
    *  For TUM Bag:  
     `roslaunch gp_map mapping_3d_simulation.launch param_file:=$(rospack find gp_map)/param/cabinet.yaml`
3.  When you want to recover the map (compute the posterior) type `RecoverMap` in the same terminal of mapping_3d.launch

### Time and accuracy
1. Time will be print out for every cloud
2. accuracy  
    `roslaunch gp_map accuracy.launch`

### Parameter
1.  mapping_3d: 
    *  `save_file`: whether to save time file
    *  `file_name`: file location
    *  `debug`: will print all sorts of stuff for debug use, use as mean and covariance when update the model
    *  `sensor_max_range`: laser/depth maximum range, anything outside this range will be considered as free
    *  `xmin`: map size, xmin
    *  `ymin`: map size, ymin
    *  `zmin`: map size, zmin
    *  `xmax`: map size, xmax
    *  `ymax`: map size, ymax
    *  `zmax`: map size, zmax
    *  `xres`: x axis grid resolution
    *  `yres`: y axis grid resolution
    *  `zres`: z axis grid resolution
    *  `iter`: conjugate gradient iteration
    *  `frame_id`: laser point frame id 
    *  `kernel_sigma_x`: kernel width in x direction
    *  `kernel_sigma_y`: kernel width in y direction
    *  `kernel_sigma_z`: kernel width in z direction
    *  `free_step_size`: sample free points along a ray of laser, how far away should you sample free points
    *  `kernel_radius`: how many closest points to keep
    *  `thresh`: threshold for posterior to determine a cell is occupied
    *  `delta`: delta for information vector p(y=1), used for visualization before recover the map
    *  `noise_cov`: observation noise
2.  play_tum_bag.launch:
    *  `cloud_input_ns`: input cloud topic name
    *  `rosbag_path`: rosbag path
    *  `rosbag_delay`: bag delay in second
    *  `speed`: play bag speed
