#!/bin/sh

PROJHOME=`dirname $BASH_SOURCE[0]`
MKLROOT=${MKLROOT:-"/opt/intel/composerxe/mkl"}
MKL_INTEL_LP64="${MKLROOT}/lib/libmkl_intel_lp64.a" 
MKL_INTEL_THREAD="${MKLROOT}/lib/libmkl_intel_thread.a"
MKL_CORE="${MKLROOT}/lib/libmkl_core.a"
MKL_IOMP5="${MKLROOT}/../compiler/lib/libiomp5.a"
echo Using MKLROOT "$MKLROOT".

if [ ! -d "$PROJHOME/lib/mkl/" ]; then
  echo "Copying MKL library files into project."
  mkdir -p "$PROJHOME/lib/mkl/";
  cp "$MKL_INTEL_LP64" "$MKL_INTEL_THREAD" "$MKL_CORE" "$MKL_IOMP5" "$PROJHOME/lib/mkl/"
else
  echo "MKL libraries already present in project."
fi

echo "Compiling linear algebra wrapper functions."
gcc -m64 -w -I"${MKLROOT}/include" -c linear_algebra.c -o lib/linear_algebra.o 

