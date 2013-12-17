#!/bin/sh

MKLROOT=${MKLROOT:-"/opt/intel/composer_xe_2013_sp1.0.080/mkl"}
echo Using MKLROOT "$MKLROOT".

PROJHOME=`dirname $BASH_SOURCE[0]`

if [ ! -d "$PROJHOME/lib/mkl/" ]; then
  echo "Copying MKL library files into project."
  mkdir -p "$PROJHOME/lib/mkl/";
  cp "${MKLROOT}/lib/intel64/libmkl_intel_lp64.a" \
     "${MKLROOT}/lib/intel64/libmkl_intel_thread.a" \
     "${MKLROOT}/lib/intel64/libmkl_core.a" \
     "${MKLROOT}/../compiler/lib/intel64/libiomp5.a" \
     "$PROJHOME/lib/mkl/"
else
  echo "MKL libraries already present in project."
fi

echo "Compiling linear algebra wrapper functions."
gcc -m64 -w -I"${MKLROOT}/include" -o lib/linear_algebra.o -c linear_algebra.c 

