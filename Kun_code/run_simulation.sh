#!/bin/sh



for p in 200
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 200 300 500
	    do
		qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh submission $n $p null2 200 none none 0 $i $j & 
      done
	done
  done
done

