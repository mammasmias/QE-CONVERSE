#!/bin/sh

echo "Si shielding:"
grep Chemical *_Si*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 837.9127 + tr/3.0}'

echo

echo "O shielding:"
grep Chemical *_O*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 271.0758 +tr/3.0}'

echo

