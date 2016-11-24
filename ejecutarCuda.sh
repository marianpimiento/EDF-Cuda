#!/bin/bash

echo -e "\n"
echo -e "\t\tUNIVERSIDAD INDUSTRIAL DE SANTANDER"
echo -e "\t\t\tPROYECTO DE GRADO:"
echo -e "\t\b\bPROCESAMIENTO Y VISUALIZACION EN ARQUITECTURAS EMBARCADAS"
echo -e " \t\t\t\b\b\b\bDE COMPUTO DE ALTO DESEMPENIO"
echo -e "\n"
echo -e "\t\t\tMARIA ANDREA PIMIENTO"
echo -e "\n"

./1-varianza $1 $2
./2-topografia $1 $2
./3-matricesRGB $1 $2
