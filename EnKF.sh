#!/bin/bash

echo "Numerical illustration (Chapter 6.3.4, Asch et al. 2016 Data Assimilation)"

python src/enkf_example.py

# echo "Stochastic Ensemble Kalman Filter (EnKF)"
# echo "Spin up for one year, and save state variables on the next year"
# echo "(1) Generate the true data: x_true"

# python src/nature_run.py

# mkdir -p Data
# if compgen -G "*.npy" > /dev/null; then
#     echo "Move the true data to /Data directory"
#     mv *.npy Data/
# fi

# echo "(2) Generate the observation data: x_obs"

# python src/make_observation.py

# if compgen -G "*.npy" > /dev/null; then
#     echo "Move the observation data to /Data directory"
#     mv *.npy Data/
# fi

# echo "(3) Algorithm of the EnKF according to Asch et al. 2016"

# if compgen -G "*.npy" > /dev/null; then
#     rm *.npy
# fi

# METHOD='m1 m2 m3'
# MEM='100'
# for M in ${MEM}; do
#     for MET in ${METHOD}; do
#     python src/enkf.py -m ${M} -t ${MET}
# done; done

# INFL='0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.20'
# for i in ${INFL}; do
#     python src/ekf.py -i ${i} -m ${m} -s 20
# done

# echo "(4) Plot the figures of EKF"
# python src/ekf_plot.py -m ${m}

# if compgen -G "*.png" > /dev/null; then
#     echo "Move figures to /Figure"
#     mv *.png Figure/
# fi

# if compgen -G "*.npy" > /dev/null; then
#     mv ekf_*.npy Data/
# fi

# echo "Finished EnKF of Lorenz 96 model"
