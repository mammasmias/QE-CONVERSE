#!/bin/sh
NAME="h2o+_gtensor"

for DIR in 1 2 3 
do
cat > ${NAME}_${DIR}.in << EOF

&input_qeconverse
        prefix = 'h2o+'
        outdir = './scratch/'
        diagonalization = 'david'
        verbosity = 'high'
        q_gipaw = 0.01
        dudk_method = 'covariant'
        diago_thr_init = 1d-4
        conv_threshold = 1d-8
        mixing_beta = 0.5
        lambda_so(${DIR}) = 1.0
/
EOF

mpirun -np 6 "your_path"/bin/qe-converse.x  < ${NAME}_${DIR}.in > ${NAME}_${DIR}.out
echo ${NAME}_${DIR}

done

