#!/bin/bash

tail -n +2 $1 > ent.tmp #input parentspec file is $1
awk '{print $1,$2,$3,$4,$6-'$2'}' ent.tmp > ent2.tmp #threshold is $2

min_above_thres=`awk 'BEGIN{a=1000}{if ($5<0+a && $5>0) a=$5} END{print a}' ent2.tmp`;
max_below_thres=`awk 'BEGIN{a=-1000}{if ($5>0+a && $5<0) a=$5} END{print a}' ent2.tmp`;

#echo $min_above_thres
#echo $max_below_thres
echo $min_above_thres - $max_below_thres | bc
rm ent.tmp
rm ent2.tmp


