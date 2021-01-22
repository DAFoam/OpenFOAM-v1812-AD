for n in `seq 16 1 16`; do

  for i in `seq 0 1 2`; do
    mpirun -np 4 simpleFoamPointPartDerivForward -parallel -point $n -comp $i -proc 1 -mode FD 
  done

done
