#!/bin/bash

echo -e "\nCompilar"
echo -e "\n"

nvcc 1_varianza.cu -o 1-varianza
nvcc 2_topografia.cu -o 2-topografia
nvcc 3_matricesRGB.cu -o 3-matricesRGB
gcc topografiaC.c -o ctopo