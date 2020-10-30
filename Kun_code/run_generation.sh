#!/bin/sh



for p in 127
  do
  for n in 100 200 300 500 
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null1.1 200 zinegbin none 0 submission &
  done
done
