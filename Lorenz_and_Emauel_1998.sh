#!/bin/bash

echo "Reproduce Lorenz 96 model with 40 variables"

python src/LE1998.py

if compgen -G "*.png" > /dev/null; then
    echo "Move figures to /Figure"
    mv *.png Figure/
fi

echo "Finished reproduction of Lorenz and Emauel 1998 Figure 1"
