#!/bin/sh

for p in 127
  do
  for n in 100 150 200 289
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null1.1 200 zinegbin none 0 &
  done
done

for p in 100 200
  do
  for n in 100 200 500
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p null2 200 none none 0 &
  done
done

for p in 127
  do
  for n in 100 150 200 289
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin erdos_renyi 100 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin erdos_renyi 1000 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin chain_large 0 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt1 200 zinegbin chain_small 0 &
  done
done

for p in 100 200
  do
  for n in 100 200 500
    do
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none erdos_renyi 100 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none erdos_renyi 1000 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none chain_large 0 &
    Rscript ~/Desktop/micro_net/Kun_code/run_data_example.R $n $p alt2 200 none chain_small 0 &
  done
done
