#!/bin/sh



for p in 200
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 200 500
	    do
        qsub -cwd -q shojaie-bigmem.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p null2 500 $i $j rmvnorm_mu04
      done
	done
  done
done


