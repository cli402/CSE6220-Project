#!/bin/bash

n=8
l=40
d=2

echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt

d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt

d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt


d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt


d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt


d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt


d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt


d=`expr $d + 2`
echo "$d"
./generate_input $n $l $d ./testFile/Data.txt
mpirun -np 1 ./findmotifs 3 ./testFile/Data.txt ./testFile/Sol.txt



