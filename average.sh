#!/usr/bin/bash

awk 'BEGIN{FS = "\t";}
	{sub(/\t*$/, "", $0);}
    {sum = 0;
	for(i = 1;i <= NF; i++)
	{
		sum += $i;
	}
	average = sum/NF;
  print average;}'
