#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
export LC_ALL=C

QE_SOURCE=/home/sfioccola/Desktop/q-e-qe-7.2
DEPENDS="${QE_SOURCE}/include ${QE_SOURCE}/iotk/src
         ${QE_SOURCE}/Modules ${QE_SOURCE}/PW/src
         ${QE_SOURCE}/FFTXlib/src ${QE_SOURCE}/UtilXlib
         ${QE_SOURCE}/LAXlib  ${QE_SOURCE}/upflib
         ${QE_SOURCE}/XClib  ${QE_SOURCE}/KS_Solvers
	 ${QE_SOURCE}/dft-d3 ${QE_SOURCE}/cmake
	 ${QE_SOURCE}/external/devxlib/src"

cd src

${QE_SOURCE}/install/moduledep.sh $DEPENDS > make.depend
${QE_SOURCE}/install/includedep.sh $DEPENDS >> make.depend

# handle special cases
sed -i '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d' make.depend
sed -i '/@iso_c_binding@/d;/@ifcore@/d' make.depend

# handle FoX
sed -i "s|@fox_dom@|${QE_SOURCE}/FoX/finclude/fox_dom.mod|" make.depend
sed -i "s|@fox_wxml@|${QE_SOURCE}/FoX/finclude/fox_wxml.mod|" make.depend

    
# check for missing dependencies 
if grep @ make.depend
then
  echo WARNING: dependencies not found in directory 'src'
  exit 1
fi
