#!/bin/bash

echo "3D variational method"

# echo "Spin up for one year, and save state variables on the next year"
# echo "(1) Generate the true data: x_true"

# METHOD='ivp ode rk4'
# for i in ${METHOD}; do
#     python src/nature_run.py -m ${i}
# done
# echo "Use the rk4 to conduct the data assimilation course"
# cp -ar true_rk4.npy true.npy

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

echo "(3) Algorithm of the 3D VAR according to Asch et al. 2016"
if [ -e 3dvar.txt ]; then
    rm 3dvar.txt
fi

INFL='0.05 0.06 0.07 0.08 0.09 0.1 0.12 0.14 0.16 0.18 0.20 0.22 0.225 0.23 0.235 0.24 0.25 0.255 0.26 0.265 0.27 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58 0.60'
for i in ${INFL}; do
    python src/3dvar.py -i ${i}
done

echo "(4) Plot the figures of 3D VAR method"
python src/3dvar_plot.py

if compgen -G "*.png" > /dev/null; then
    echo "Move figures to /Figure"
    mv *.png Figure/
fi

mv 3dvar.txt Data/

echo "Finished 3D VAR method of Lorenz 96 model"
