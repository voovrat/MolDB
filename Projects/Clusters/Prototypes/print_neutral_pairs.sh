#!/bin/bash

for ((i=0;i<7;i++))
do
   for ((j=$i+1;j<7;j++))
   do
       echo c${i}c${j}

   done
 
done

for ((i=0;i<7;i++))
do
   echo c$i
done
