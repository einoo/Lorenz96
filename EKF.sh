#!/bin/bash

# echo "Extended Kalman Filter"
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

echo "(3) Algorithm of the EKF according to Asch et al. 2016"

if compgen -G "*.npy" > /dev/null; then
    rm *.npy
fi

# METHOD='m1 m2'
METHOD='m3'
INFL='0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.20'
for m in ${METHOD}; do
    for i in ${INFL}; do
        python src/ekf.py -i ${i} -m ${m} -s 20
    done

    echo "(4) Plot the figures of EKF"
    python src/ekf_plot.py -m ${m}
done

if compgen -G "*.png" > /dev/null; then
    echo "Move figures to /Figure"
    mv *.png Figure/
fi

if compgen -G "*.npy" > /dev/null; then
    mv ekf_*.npy Data/
fi

echo "Finished EKF of Lorenz 96 model"
