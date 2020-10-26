#!/bin/sh


for p in 200
  do
  for n in 100 200 300 500 
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none erdos_renyi 1000 submission &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none erdos_renyi 100 submission &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null2 200 none none 0 submission &
  done
done

for p in 127
  do
  for n in 100 200 300 500 
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin erdos_renyi 1000 submission &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin erdos_renyi 100 submission &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null1.1 200 zinegbin none 0 submission &
  done
done
