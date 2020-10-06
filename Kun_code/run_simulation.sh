#!/bin/sh
for p in 200 
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100
	    do
        qsub -cwd -q shojaie-bigmem.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p null2 200 none none 0 $i $j
      done
	done
  done
done


