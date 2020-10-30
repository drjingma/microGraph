#!/bin/sh



for p in 127
  do
  for i in {101..150..1} 
    do
    for j in {2..7..1}
      do
	  for n in 100 200 300 500
	    do
		qsub -cwd  -e iotrash/ -o iotrash/ run_data_example_part_vary.csh submission $n $p null1.1 200 zinegbin none 0 $i $j & 
      done
	done
  done
done
