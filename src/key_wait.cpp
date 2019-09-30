#include <ros/ros.h>
#include <std_msgs/String.h>

int main(int argc, char** argv)
{
  ros::init(argc, argv, "gp_map");

  ros::NodeHandle nh;
  ros::Publisher p = nh.advertise<std_msgs::String> ("chatter", 1);
  while(true)
  {
    std::string inputString;
    std::cout << "Type [RecoverMap] to start map recovery" << std::endl;
    std::getline(std::cin, inputString);
    std_msgs::String msg;

    std::string first_string = inputString.substr(0, 10);

    if(first_string.compare("RecoverMap") == 0)
    {
      //send a request to the node serving out the messages
      //print out recieved messages.
      msg.data = inputString;
      p.publish(msg);
    } else {
       std::cout << "Don't understand what you want to do. Please type [RecoverMap]" << std::endl; 
    }

    ros::spinOnce();
  }

  return 0;
}