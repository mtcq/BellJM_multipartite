#!/bin/bash
[ -e variables_faacets ] && rm variables_faacets
for i in {1..46}
do
   sed -n 5p $i.yaml>aux1
   sed -n 6p $i.yaml>aux2
   paste aux1 aux2>aux3
   awk -F "[" '{print "["$2}' aux3>aux4
   echo "S(:,$i) = " >aux5
   echo ';' > aux6
   paste aux5 aux4 aux6 >> variables_faacets
done
rm aux1 aux2 aux3 aux4 aux5 aux6
