#!/bin/sh
for p in 127
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 289
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin erdos_renyi 100 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin erdos_renyi 1000 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin chain_large 0 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin chain_small 0 $i $j
      done
	done
  done
done