#!/usr/bin/env bash

if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit 1
fi

cd simpleFoamMVStateProductReverse && wclean && rm log && cd - || exit 1
cd simpleFoamMVPointProductReverse && wclean && rm log && cd - || exit 1
cd run && rm *.txt && cd - || exit 1
