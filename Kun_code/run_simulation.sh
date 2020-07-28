#!/bin/sh
for n in 100 200 500
  do
  for i in {1..10..1} 
    do
    for j in 7
      do
      qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 200 null2 100 $i $j
    done
  done
done
for n in 100 200 500
  do
  for i in {1..10..1} 
    do
    for j in 7
      do
      qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh $n 200 alt1 100 $i $j
    done
  done
done
for i in {1..10..1} 
  do
  for j in 7
    do
    qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh 289 127 null1 100 $i $j
  done
done
for i in {1..10..1} 
  do
  for j in 7
    do
    qsub -cwd -q shojaie-normal.q -e iotrash/ -o iotrash/ run_data_example_part2.csh 289 127 alt2 100 $i $j
  done
done


