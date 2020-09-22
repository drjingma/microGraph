#!/bin/sh

for n in 100 200 500
  do
  for i in {40..50..1} 
    do
    for j in 1 2 3 4 5 6 7
      do
      qsub -cwd -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 200 null2 100 $i $j rmvnorm_mu04
    done
  done
done
