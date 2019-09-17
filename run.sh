#!/bin/bash
gfortran Probability_analysis.F90 -o  Probability_analysis.x
#./mtd_analysis_us_test.x -T 300 -tmin 5000 -grid -1.0 3.0 0.02 -pfrqMD 10
#./Probability_analysis.x -T 300 -tmin 5000 -ncv 2 -grid -1.0 3.0 0.02 0.0 6.0 0.02  -pfrqMD 10
./Probability_analysis.x ,\
 -T 300                  ,\
 -tmin 5000              ,\
 -ncv 4                  ,\
 -grid -1.0 3.0 0.02 0.0 7.0 0.02 3.0 6.0 0.02 3.0 6.0 0.02 ,\
 -pfrqMD 10
