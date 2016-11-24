#!/bin/bash
for ((i=$1;i<=$2;i+=1)); do
imagen=`ls | grep Image | sort -V | head -n 1`
mv $imagen $i
done
