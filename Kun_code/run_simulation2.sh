#!/bin/sh

for p in 20 50 100 150 200
  do
  for i in {1..20..1} 
    do
    for j in 1 2 3 4 5 6 7
      do
	  for n in 100 200 500
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt2 100 $i $j 2pnetwork
      done
	done
  done
done


