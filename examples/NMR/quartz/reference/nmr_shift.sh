#!/bin/bash
NAME="nmr"
molecule=quartz
natoms=4
atoms=("null" "Si" "Si" "Si" "O")
cart=("null" "x" "y" "z")


for i in `seq 1 $natoms`
do
for DIR in 1 2 3
do
	base=${molecule}_${atoms[$i]}$i${cart[$DIR]}
	cat > ${base}.in << EOF



&input_qeconverse
        prefix = 'quartz'
        outdir = './scratch/'
        diagonalization = 'david'
        verbosity = 'high'
        q_gipaw = 0.01
        dudk_method = 'covariant'
        mixing_beta = 0.5
        lambda_so(1) = 0.0
        m_0(${DIR}) = 1.0
        m_0_atom = ${i}
/
EOF

mpirun -np 6 "YOUR_QUANTUM_ESPRESSO"/bin/qe-converse.x < ${base}.in > ${base}.out

done
done

