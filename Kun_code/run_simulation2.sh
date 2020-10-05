#!/bin/sh
## this is for back-up: these have been run
for p in 127
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 150 200 289
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p null1.1 200 zinegbin none 0 $i $j
      done
	done
  done
done

for p in 200 100
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 200 500
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p null2 200 none none 0 $i $j
      done
	done
  done
done

for p in 127
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 150 200 289
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin erdos_renyi 100 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin erdos_renyi 1000 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin chain_large 0 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt1 200 zinegbin chain_small 0 $i $j
      done
	done
  done
done

for p in 100 200
  do
  for i in {1..200..1} 
    do
    for j in {1..7..1}
      do
	  for n in 100 200 500
	    do
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt2 200 none erdos_renyi 100 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt2 200 none erdos_renyi 1000 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt2 200 none chain_large 0 $i $j
        qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n $p alt2 200 none chain_small 0 $i $j
      done
	done
  done
done


