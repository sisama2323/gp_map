#!/usr/bin/env python

import numpy as np
from visualization_msgs.msg import MarkerArray

import rospy
import matplotlib.pyplot as plt


class Accuracy(object):
    def __init__(self, xmin, ymin, zmin, resolution):
        self.truth = None
        self.resolution = resolution
        self.xmin = xmin 
        self.ymin = ymin 
        self.zmin = zmin

        ############################################# initialize subscriber ###########################################################
        rospy.init_node('accuracy', anonymous=True)
        self.map1_sub = rospy.Subscriber("map1_topic", MarkerArray, self.truth_call_back, queue_size = 10)
        self.map2_sub = rospy.Subscriber("truth_map", MarkerArray, self.callback2, queue_size = 100)

    def truth_call_back(self, map_msg):
        # map1
        print("Get truth map")
        self.truth = map_msg
    
    def callback2(self, map_msg):
        print("get map")
        # ground truth
        if self.truth is not None:
            truth_marker = self.truth.markers
            marker2 = map_msg.markers

            point2 = []
            for m in marker2:
                if m.points:
                    for p_tmp in m.points:
                        point2.append([p_tmp.x, p_tmp.y, p_tmp.z])
            point2 = np.floor((np.array(point2)) * self.resolution - np.array([xmin, ymin, zmin])).astype(int)

            truth_point = []
            for m in truth_marker:
                if m.points:
                    for p_tmp in m.points:
                        truth_point.append([p_tmp.x, p_tmp.y, p_tmp.z])
            truth_point = (np.array(truth_point) - np.array([xmin, ymin, zmin]) / self.resolution).astype(int) 
            
            x_min = min(np.min(truth_point[:, 0]), np.min(point2[:, 0]))
            x_max = max(np.max(truth_point[:, 0]), np.max(point2[:, 0]))

            y_min = min(np.min(truth_point[:, 1]), np.min(point2[:, 1]))
            y_max = max(np.max(truth_point[:, 1]), np.max(point2[:, 1]))

            z_min = min(np.min(truth_point[:, 2]), np.min(point2[:, 2]))
            z_max = max(np.max(truth_point[:, 2]), np.max(point2[:, 2]))

            truth_point = truth_point - np.array([x_min, y_min, z_min])

            point2 = point2 - np.array([x_min, y_min, z_min])

            grid_map_t = np.zeros((x_max - x_min + 1, y_max - y_min + 1, z_max - z_min + 1))
            grid_map_t[truth_point[:, 0], truth_point[:, 1], truth_point[:, 2]] += 1

            grid_map = np.zeros((x_max - x_min + 1, y_max - y_min + 1, z_max - z_min + 1))
            grid_map[point2[:, 0], point2[:, 1], point2[:, 2]] += 1

            error = np.sum(abs(grid_map_t - grid_map)) / (grid_map_t.ravel().shape[0])

            print("error rate for if is: ", error)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(point2[:, 0], point2[:, 1], point2[:, 2], 'b.', alpha = 0.4)
            ax.plot(truth_point[:, 0], truth_point[:, 1], truth_point[:, 2], 'r.', alpha = 0.2)
            plt.show()


if __name__ == '__main__':
    try:
        xmin = rospy.get_param('xmin')
        ymin = rospy.get_param('ymin')
        zmin = rospy.get_param('zmin')
        resolution = rospy.get_param('resolution')
        accuracy = Accuracy(xmin, ymin, zmin, resolution)
        rospy.spin()
    except KeyboardInterrupt:
        print("Shutting down ROS")



