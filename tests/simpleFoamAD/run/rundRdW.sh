#!/bin/bash

for n in `seq 0 2 100`; do

    rm verify_*
    dRdWSimpleFoam -proc 0 -cell $n -var U
    echo "U", $n
    python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt

    rm verify_*
    dRdWSimpleFoam -proc 0 -cell $n -var p
    echo "p", $n
    python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt

    rm verify_*
    dRdWSimpleFoam -proc 0 -cell $n -var phi
    echo "phi", $n
    python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt

done

for i in `seq 0 1 3`; do
    for n in `seq 0 2 50`; do

        rm verify_*
        mpirun -np 4 dRdWSimpleFoam -proc $i -cell $n -var U -parallel
        echo "U", $n, $i
        python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_1_CellI_${n}_FD.txt verify_ProcI_1_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_2_CellI_${n}_FD.txt verify_ProcI_2_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_3_CellI_${n}_FD.txt verify_ProcI_3_CellI_${n}_AD.txt

        rm verify_*
        mpirun -np 4 dRdWSimpleFoam -proc $i -cell $n -var p -parallel
        echo "p", $n, $i
        python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_1_CellI_${n}_FD.txt verify_ProcI_1_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_2_CellI_${n}_FD.txt verify_ProcI_2_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_3_CellI_${n}_FD.txt verify_ProcI_3_CellI_${n}_AD.txt

        rm verify_*
        mpirun -np 4 dRdWSimpleFoam -proc $i -cell $n -var phi -parallel
        echo "phi", $n, $i
        python checkDerivs.py verify_ProcI_0_CellI_${n}_FD.txt verify_ProcI_0_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_1_CellI_${n}_FD.txt verify_ProcI_1_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_2_CellI_${n}_FD.txt verify_ProcI_2_CellI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_3_CellI_${n}_FD.txt verify_ProcI_3_CellI_${n}_AD.txt

    done
done
