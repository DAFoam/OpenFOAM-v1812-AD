#!/bin/bash

for n in `seq 0 2 100`; do

    rm verify_*
    dRdXSimpleFoam -proc 0 -point $n -eps 0.00001
    echo "point", $n
    python checkDerivs.py verify_ProcI_0_PointI_${n}_FD.txt verify_ProcI_0_PointI_${n}_AD.txt

done

for i in `seq 0 1 3`; do
    for n in `seq 0 2 100`; do

        rm verify_*
        mpirun -np 4 dRdXSimpleFoam -proc $i -point $n -eps 0.00001 -parallel
        echo "point", $n, "proc", $i
        python checkDerivs.py verify_ProcI_0_PointI_${n}_FD.txt verify_ProcI_0_PointI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_1_PointI_${n}_FD.txt verify_ProcI_1_PointI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_2_PointI_${n}_FD.txt verify_ProcI_2_PointI_${n}_AD.txt
        python checkDerivs.py verify_ProcI_3_PointI_${n}_FD.txt verify_ProcI_3_PointI_${n}_AD.txt

    done
done


