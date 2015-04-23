
This folder includes:

1) The standard DISORT version 2.0 available at Wiscombe's NASA ftp site: ftp://climate1.gsfc.nasa.gov/wiscombe/ (README_disort.txt is the readme from the DISORT team that accompanies the DISORT code). Minor modifications have been made to DISORT (see Notes below). 

2) Additional code written by P. Rowe and S. Neshyba includes a (Fortran) driver for DISORT, disort_driver_mat to call DISORTv2beta, and the MAKEFILE needed to compile it. After compiling as described below, this driver can be called from Matlab or Octave via a system call. It takes as input a name list that is created from the runDisort_mat codes and outputs a file that can then be loaded in by Matlab or Octave (planned upgrades include a direct call to DISORT). 


Compiling disort_driver_mat:
1) Copy these files to your directory. 
2) Within the installation directory, modify the makefile as needed for the fortran compiler you have; here we use gfortan. Run "make".
3) A file called "disort_driver_mat” should be created. This is the code that will call DISORT.
4) Move disort_driver_mat into the main directory (where sample_run.m is), or create an alias there.
5) Test the code by running sample_run.m from Matlab or Octave.




A note on precision: some variables in DISORT are single precision. Sensitivity studies indicate that best accuracy is achieved when all layer optical depths in typical model atmospheres are above 10^-5. For example, for a total optical depth of 0.5 and temperatures near 240 K, our studies indicate that including a layer optical depth of 10^-6 results in round-off errors in zenith downwelling radiance of ~0.2 mW/(m2 sr cm-1), whereas omitting these layers from the calculation causes an error of only ~0.004 mW/(m2 sr cm-1). Furthermore, omitting these layers saves computational time. For this reason, the code excludes upper atmospheric layers with optical depths below 10^-5. If such thin layers need to be included, it is possible to compile entirely in double-precision at a computational cost probably less than 20% (see DISORT documentation). As described by Istvan Laszlo (personal communication), “To run DISORT in double-precision you should use DISORTsp.f, which is currently only available in pre-version 3 distributions, like in DISORT2.0beta. To create the double-precision version it is, however, not sufficient to simply auto-double DISORTsp.f at compile time. You would first need to change all instances of R1MACH to D1MACH. Then you should compile all files with the auto-double option, except one file: RDI1MACH.f should be compiled without this option. Finally you would need to link all compiled files to create the executable.” We plan to include a double-precision option soon.


Notes:

1) The DISORT code is DISORTv2beta. Please be sure to acknowledge use of DISORT. 

2) Although provided here, the DISORT (Fortran) code is identical to the original code except for two changes. In disort.f, the maximum number of layers, MXCLY has been increased from 6 to 120. A header that was printed to the screen has been suppressed. Thus we recommend running run_disort.sh, provided with the DISORT codes, and comparing output to their test output. While not all of the DISORT codes are needed here, we have included them so that this code can be run to ensure DISORT is running properly on your system. See README.txt.

3) disort_driver_mat is based on disotest.f, provided with the DISORT codes.



REFERENCES for DISORT code: 

Stamnes, K., Tsay, S.-C., Wiscombe, W. J., and Jayaweera, K.: Numerically stable algorithm for discrete-ordinate-method radiative transfer in multiple scattering and emitting layered media, Appl. Optics, 27(12), 2502–2509, doi:10.1364/AO.27.002502, 1988. Stamnes, K., Tsay, S.-C., Wiscombe, W., and Laszlo, I.: DISORT, a general-purpose Fortran program for discrete-ordinate-method radiative transfer in scattering and emitting layered media: Documentation of methodology, Tech. rep., Dept. of Physics and Engineering Physics, Stevens Institute of Technology, Hoboken, NJ 07030, 2000. REFERENCE FOR THE DISORT DRIVER for matlab: P.M. Rowe, S. Neshyba, and V. P. Walden. A reference is pending, so please contact Penny Rowe (prowe@harbornet.com for proper referencing)



For more information, contact
Penny M. Rowe (prowe@harbornet.com)
March 4, 2015



