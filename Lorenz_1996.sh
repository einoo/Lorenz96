#!/bin/bash

echo "Reproduce the error propogation by Lorenz 96 model"

python src/L1996.py

if compgen -G "*.png" > /dev/null; then
    echo "Move figures to /Figure"
    mv *.png Figure/
fi

echo "Finished reproduction of Lorenz 1996 Figure 2"
