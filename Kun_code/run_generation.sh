#!/bin/sh

for p in 200
  do
  for n in 150 250 
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null1.1 200 zinegbin none 0 newref2p &
  done
done
