# QE-CONVERSE
This is the official Git repository of the QE-CONVERSE code for Quantum-Espresso.  
This current version of the code is compatible with the version 7.5 of Quantum-Espresso package.


## Features
The QE-CONVERSE implement a non-perturbative approach (converse) to compute the orbital magnetization in isolated and periodic systems. The calculation of orbital magnetization allows ab-initio computation of macroscopic properties like the Nuclear Magnetic Resonance (NMR) chemical shifts and the Electronic Paramagnetic Resonance (EPR) g tensor.

* NMR shielding tensors
* EPR g-tensor
* Electric field gradient (EFG) tensors and NMR/NQR quadrupolar parameters (Cq, eta, nu_Q), via the standalone `qe-efg.x` driver (a ground-state property computed directly from the SCF density)
* It works only with Norm-conserving pseudopotential with GIPAW reconstruction (https://sites.google.com/site/dceresoli/pseudopotentials)
* LDA and GGA functionals
* isolated and periodic systems

## Build instructions:
1. the Quantum-Espresso package version 7.5 must be previously installed (https://gitlab.com/QEF/q-e/-/releases/qe-7.5). To take advantage of the enhancements in linear algebra operations, the configuration with scaLAPACK package or ELPA library is suggested.
2. ```git clone https://github.com/mammasmias/QE-CONVERSE``` 
3. ```cd QE-CONVERSE```
4. ```chmod +x configure```
5. ```./configure --with-qe-source="QE folder containing make.inc"```
6. ```make```
7. ```make test # optionally run tests```
## How to use it:
run ```./qe-converse.x``` for NMR shifts / EPR g-tensor, or ```./qe-efg.x``` for the electric field gradient.

## Directory contents

```/src/```: Contains the source code and the Makefile.
```/doc/```: Contains the User's manual.
```/examples/```: Contains two directories:
 ```/EPR/``` about a EPR g tensor and  ```/NMR/``` about a NMR chemical shift calculation. Inside each one there's a Tutorial.wiki file that explain how to perform the calculation step-by-step.
```/benchmarking/``` : Contains the input files used to benchmark the code. Concerns the EPR g tensor calculation of diatomic paramagnetic radicals.
```/applications/```: Contains the EPR and NMR calculations.
```/tests/```: Contains integration and unit tests.

