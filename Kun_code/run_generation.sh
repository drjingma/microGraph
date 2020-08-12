#!/bin/sh

for p in 20 50 100 150 200
  do
  for n in 100 200 500
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null2 100 none &
  done
done
