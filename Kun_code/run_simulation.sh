#!/bin/sh
for p in 200
  do
  for i in {1..20..1} 
    do
    for j in {1..7..1}
      do
	  for n in 150 250 
	    do
        qsub -cwd -q shojaie-bigmem.q -e iotrash/ -o iotrash/ run_data_example_part2.csh newref2p $n $p null1.1 200 zinegbin none 0 $i $j 
      done
	done
  done
done
