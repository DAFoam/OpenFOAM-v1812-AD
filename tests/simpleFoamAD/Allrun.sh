#!/usr/bin/env bash

if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit 1
fi

cd simpleFoamMVStateProductReverse && wclean && wmake 2> log && cd - || exit 1
cd simpleFoamMVPointProductReverse && wclean && wmake 2> log && cd - || exit 1
cd run && cp refs/* . && cd - || exit 1
cd run && mpirun -np 4 simpleFoamMVStateProductReverse -parallel && cd - || exit 1
cd run && mpirun -np 4 simpleFoamMVPointProductReverse -parallel && cd - || exit 1
cd run && python checkDerivs.py 4 state && cd - || exit 1
cd run && python checkDerivs.py 4 point && cd - || exit 1
