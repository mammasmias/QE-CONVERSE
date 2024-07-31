In this directory the input files used for the calculation of EPR g tensor in substitutional Nitrogen in Silicon (paragraph 7.1) and 27Al NMR chemical shift in Alumina (paragraph 7.2) are provided:

In the /EPR_NSi/ directory:

```slurm``` file is an example of script file can be used to compute the SCF calculation with the ```scf_NSi_512off.in``` input file.

```slurm_gtensor``` is and example of script file can be used to run the code qe-converse.x and compute the DeltaG tensor matrix.

NOTE: no output files are provided for this application


In the /NMR_Alumina/ directory:

In the /Alpha_Al2O3/ subdirectory, the input file for the NMR shielding tensor calculation of 27Al in Alpha-Al2O3 is provided.

the script ```getshift.sh``` extrapolate the 27Al shielding tensor which results equal to 531.665 ppm.

In the /Theta_Al2O3/ subdirectory, the input file for the NMR shielding tensor cal
culation of 27Al(oct) and 27Al(tet) is provided.

the script ```getshift.sh``` extrapolate the 27Al shielding tensor:

27Al(oct) = 537.97 ppm
27Al(tet) = 470.69 ppm

To extrapolate the NMR chemical shift:

delta=sigma_ref - sigma_iso


sigma_ref(Alpha-Al2O3) = 531.665

sigma 27Al(oct)  = 531.665 -470.69  = 60.97 ppm
sigma 27Al(tet)  = 531.665 -537.97  = -6.30 ppm


