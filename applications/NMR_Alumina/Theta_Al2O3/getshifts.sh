#!/bin/sh

echo "Al shielding:"
grep Chemical *_Al*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 766.8170 + tr/3.0}'


echo

