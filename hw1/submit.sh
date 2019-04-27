#!/bin/sh

g++ backpack.cpp --std=c++17 -O2

for i in 1 2 3 4 5 6 7 8 9 10
do
  echo "Processing test $i"
  ./a.out "data/$i.public" "results/$i.result"
  ./submit.py "mikhailov_nikita_m" "knapsack" "$i.public" "results/$i.result"
done
