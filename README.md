# IER
Code for balancing slack liner
There are 3 seperate MATLAB codes provided. The description of each is given below.


simulate_slackline.m:
A script that calls the equations of motion of a slackliner. This script is the one that needs to be run. Uncomment line 108 in the script to run the developed controller and uncomment line 110 to run the benchmark controller. 
Please try running only one controller at a time.

equations_of_motion.m:
A function that contains the equations of motion of the slackliner and the implementation of the developed controller. The simulate_slackline.m script calls this function so please make sure this function is in the same folder as simulate_slackline.m.

equations_of_motion_benchmark.m:
A function that contains the equations of motion of the slackliner and an exemplary controller for benchmarking (note that the code is hidden). The simulate_slackline.m script calls this function so please make sure this function is in the same folder as simulate_slackline.m.
