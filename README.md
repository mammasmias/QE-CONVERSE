# QE-CONVERSE
This is the official Git repository of the QE-CONVERSE code for Quantum-Espresso. 
This current version of the code is compatible with the version 7.2 of Quantum-Espresso package.


## Features
the QE-CONVERSE implement a non-perturbative approach (converse) to compute the orbital magnetization in isolated and periodic systems. The calculation of orbital magnetization allows ab-initio computation of macroscopic properties like the Nuclear Magnetic Resonance (NMR) chemical shifts and the Electronic Paramagnetic Resonance (EPR) g tensor.

* NMR shielding tensors
* EPR g-tensor
* It works only with Norm-conserving pseudopotential with GIPAW reconstruction (https://sites.google.com/site/dceresoli/pseudopotentials)
* LDA and GGA functionals
* isolated and periodic systems


## Authors and contributors
S. Fioccola, L. Giacomazzi, D. Ceresoli, N. Richard, A. Hemeryck, L. Martin-Samos


## Build instructions:
1. the Quantum-Espresso package version 7.2 must be previously installed (https://gitlab.com/QEF/q-e/-/releases/qe-7.2). To take advantage of the enhancements in linear algebra operations, the configuration with scaLAPACK package or ELPA library is suggested.
2. git clone https://github.com/SFioccola/QE-CONVERSE.git 
3. cd QE-CONVERSE
4. Copy the source files and the Makefile from /src/ dictory into /PP/ directory of Quantum-Espresso.
5. In the main directory of QE-7.2 type ```make pp```. You should find the binary file ```qe-converse.x``` in the /bin/ directory. 

## How to use it:
$/QE-7.2/bin/qe-converse.x < input > output
