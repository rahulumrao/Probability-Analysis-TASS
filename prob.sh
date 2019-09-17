#!/bin/bash

for i in -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4
#for i in -0.4 
do
echo "=======================================UMBRELLA $i============================="
cd umb_${i}

../WS_analysis_3cvs_ProjectTo_US_TAMD.x -T 300 -tmin 5000 -grid -1.0 3.0 0.02 0.0 7.0 0.02 3.0 6.0 0.02 -pfrqMD 10
#../mtd_analysis_us_test.x -T 300 -tmin 5000 -grid -1.0 3.0 0.02 -pfrqMD 10

mv Pu.dat PROB_${i}
cp PROB_${i} ../../PROB_2D
#cp PROB_${i} ../../PROB
cd ../
done

#-0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4
