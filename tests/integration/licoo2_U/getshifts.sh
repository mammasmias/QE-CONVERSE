#!/bin/sh


echo "Co shielding:"
grep Chemical *_null*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 2018.3630  +tr/3.0}'

echo
grep Chemical *_Co*out | awk ' { tr=$5; getline; tr+=$6; getline; tr+=$7; print 2018.3630  +tr/3.0}'

