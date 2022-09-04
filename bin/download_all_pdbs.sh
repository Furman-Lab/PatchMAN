#!/bin/bash

for a in `ls *matches`;
 do head -50 "$a" | gawk '{print $2}' | cut -d '/' -f 10 | cut -d '.' -f 1 >> all_pdbs;
  tail -50 "$a" | gawk '{print $2}' | cut -d '/' -f 10 | cut -d '.' -f 1 >> all_pdbs;
 done

download_pdbs.py $1
