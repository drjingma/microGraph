#!/bin/sh



for n in 100 200 500
  do
  for i in {1..20..1} 
    do
    for j in 1 2 3 4 5 6 7
      do
      qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 200 null2 100 $i $j
    done
  done
done

for n in 100 200 500
  do
  for i in {1..20..1} 
    do
    for j in 1 2 3 4 5 6 7
      do
      qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 200 alt2 100 $i $j
    done
  done
done

for n in 289 100 120 150 200
  do
  for i in {1..20..1} 
    do
    for j in 1 2 3 4 5 6 7
      do
      qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 127 null1 100 $i $j pois
    done
  done
done