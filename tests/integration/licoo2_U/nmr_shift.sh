#!/bin/bash
NAME="nmr"
molecule=LiCoO2_U_SF
natoms=1
atoms=("Co")
cart=("null" "x" "y" "z")


for DIR in 1 2 3
do
	base=${molecule}_${atoms[$i]}$i${cart[$DIR]}
	cat > ${base}.in << EOF



&input_qeconverse
        prefix = 'licoo2'
        outdir = './scratchU/'
        diagonalization = 'david'
        verbosity = 'high'
        q_gipaw = 0.01
        dudk_method = 'covariant'
        mixing_beta = 0.5
        lambda_so(1) = 0.0
        m_0(${DIR}) = 1.0
        m_0_atom = 4
/
EOF

mpirun -np 6 /home/sfioccola/Desktop/QE-CONVERSE_DFT+U/QE-CONVERSE/bin/qe-converse.x < ${base}.in > ${base}.out

done

