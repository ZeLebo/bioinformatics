#!/bin/bash

# Get the filename from the command line argument
filename=$1

# Use grep and sed to extract the percent value
percent=$(grep -o -P '\d+\.\d+%' $filename | sed 's/%//')

# trim the first one
first=$(echo $percent | cut -f1 -d' ')

echo "Percent is $first%"

result="BAD"

if (( $(awk 'BEGIN {print ("'$first'" > "90.0")}') )); then 
	result="OK ✅✅✅"
else
	result="BAD ❌❌❌"
fi

echo $result
