#!/bin/bash
for ((i=$1;i<=$2;i+=1)); do
mkdir -p $i
done
