runDisort_mat
==============

This code runs DISORT from MATLAB or Octave to perform radiative transfer calculations for layered model atmospheres, including absorption and scattering by gases and clouds (aerosols are not currently included but could be added). Inputs include gaseous optical depths, atmospheric parameters, and cloud properties. Radiative transfer is performed using DISORT 2.0 Beta (from ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple_Scatt/). This is a work in progress, and may have bugs. It is shared without guarantees of any kind. I am not maintaining this code anymore as I have switched to Python (please see runDisort_py), but I welcome feedback and will try to respond.

Installation Instructions
1) To use this code, first copy these files to your directory. 
2) Within the installation directory, modify the makefile as needed for the fortran compiler you have; here we use gfortan. Run "make" to compile disort_driver_mat (a fortran executable). This is the code that will call disort.
3) Make an alias of disort_driver_mat or copy it into the main directory (where sample_run.m is) 
4) To try out the code (see next step), you will need these files: "pmom_water_T300_DW_S331.nc" and "pmom_ice_sphere_T266_Warren_S331.nc. They are available here:
https://people.nwra.com/rowe/public/singleScatteringParams/pmom_water_T300_DW_S331.nc and https://people.nwra.com/rowe/public/singleScatteringParams/pmom_ice_sphere_T266_Warren_S331.nc. There are other single scattering parameter sets there as well. Please be sure to reference the use of these.
5) In matlab, you can try out the code using "sample_run.m" in the sampleRun folder.  Make sure you add the path to the folder where "run_disort.m" is to your path and change the directories as needed in "sample_run.m."


Upgrades that would be nice: 
1) A new version of DISORT, DISORT3 is available. 
2) It would speed up the code to call DISORT directly from Matlab or Octave (as is done in the python version), rather than making an external call to the system and using files as intermediates.
3) This code requires gaseous optical depths as input. We use the line-by-line radiative transfer model (LBLRTM; http://rtweb.aer.com/lblrtm.html) for this purpose. A github directory for a codes for running LBLRTM from MATLAB or Octave, runLBLRTM_mat would be nice.
6) Currently on some systems the following may appear when running the code, "Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL." This is not an error but a consequence of DISORT's operation. It would be great to supress this warning.


Acknowledging use of this code (please let us know if any are missing):

1) If this code is used in work leading to a publication, please include an acknowledgement to: Rowe, P. M., Fergoda, M., and Neshyba, S. (2020). Temperature‐Dependent Optical Properties of Liquid Water From 240 to 298 K. Journal of Geophysical Research: Atmospheres, 125(17), e2020JD032624.

2) DISORT references: 

    Stamnes, K., Tsay, S.-C., Wiscombe, W. J., and Jayaweera, K.: Numerically stable algorithm for discrete-ordinate-method radiative transfer in multiple scattering and emitting layered media, Appl. Optics, 27(12), 2502Ð2509, doi:10.1364/AO.27.002502, 1988. 

    Stamnes, K., Tsay, S.-C., Wiscombe, W., and Laszlo, I.: DISORT, a general-purpose Fortran program for discrete-ordinate-method radiative transfer in scattering and emitting layered media: Documentation of methodology, Tech. rep., Dept. of Physics and Engineering Physics, Stevens Institute of Technology, Hoboken, NJ 07030, 2000. 

3) Temperature-dependent single-scattering parameters of liquid water (https://people.nwra.com/rowe/single_scatter.shtml):

    Rowe, P.M., S. Neshyba, and V. P. Walden, Radiative consequences of low-temperature infrared refractive indices for supercooled water clouds, Atmos. Chem. Phys., 13, 11925Ð11933, 2013. www.atmos-chem-phys.net/13/11925/2013/ doi:10.5194/acp-13-11925-2013 as well as Zasetsky et al., 2005 and Wagner et al., 2005 as given within.

4) Other single scattering parameters of liquid and ice (https://people.nwra.com/rowe/single_scatter.shtml): Please reference Rowe, P. M., Fergoda, M., and Neshyba, S. (2020). Temperature‐Dependent Optical Properties of Liquid Water From 240 to 298 K. Journal of Geophysical Research: Atmospheres, 125(17), e2020JD032624. as well as:
- Liquid water at 300 K: Downing, H. D. and Williams, D.: Optical constants of water in the Infrared, J. Geophys. Res., 80, 1656Ð1661, 1975, as well as 
- Ice spheres: Warren, S. G. and Brandt, R. E.: Optical constants of ice from the ultraviolet to the microwave: a revised compilation, J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744, 2008. 
- Ice rough solid column: Yang, P., Bi, L., Baum, B. A., Liou, K. N., Kattawar, G. W., Mishchenko, M. I., and Cole, B. (2013). Spectrally Consistent Scattering, Absorption, and Polarization Properties of Atmospheric Ice Crystals at Wavelengths from 0.2 to 100 μm. J. Atmos. Sci., 70, 330–347, doi:10.1175/JAS-D-12-039.1. 

6) For use of the solar irradiance spectra: Kurucz, R.L., Synthetic infrared spectra, in Infrared Solar Physics, IAU Symp. 154, edited by D.M. Rabin and J.T. Jefferies, Kluwer, Acad., Norwell, MA, 1992.

7) For use of cloud_2012012006_cld1.mat (used by sample_run.m), please reference: Cox, C., Rowe, P. M., Neshyba, S., & Walden, V. P. (2016). A synthetic data set of high-spectral resolution infrared spectra for the Arctic atmosphere. Earth System Science Data Discussions, 1Ð29. http://doi.org/10.5194/essd-2015-40, in review for Earth System Science Data. Please also see our database of atmospheric profiles and cloudy and clear sky up and downwelling radiances characteristic of the Arctic at the Arctic Observing Network (AON) Arctic data repository (at https://www.aoncadis.org/dataset/AAIRO_spectra.html; doi:10.5065/D61J97TT).




A note on precision: some variables in DISORT are single precision. Sensitivity studies indicate that best accuracy is achieved when all layer optical depths in typical model atmospheres are above 10^-5. For example, for a total optical depth of 0.5 and temperatures near 240 K, our studies indicate that including a layer optical depth of 10^-6 results in round-off errors in zenith downwelling radiance of ~0.2 mW/(m2 sr cm-1), whereas omitting this layer from the calculation causes an error of only ~0.004 mW/(m2 sr cm-1). Furthermore, omitting extremely thin layers saves computational time. For this reason, the code excludes upper atmospheric layers with optical depths below 10^-5 (this is done on a wavenumber-by-wavenumber basis). If such thin layers need to be included, it is possible to compile entirely in double-precision at a computational cost probably less than 20% (see DISORT documentation). As described by Istvan Laszlo (personal communication), "To run DISORT in double-precision you should use DISORTsp.f, which is currently only available in pre-version 3 distributions, like in DISORT2.0beta. To create the double-precision version it is, however, not sufficient to simply auto-double DISORTsp.f at compile time. You would first need to change all instances of R1MACH to D1MACH. Then you should compile all files with the auto-double option, except one file: RDI1MACH.f should be compiled without this option. Finally you would need to link all compiled files to create the executable."


Credits:

1) The DISORT code is DISORT2.0beta, from ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple_Scatt/. Please be sure to acknowledge use of DISORT. Although provided here, it is identical to the original code except for two changes. In disort.f, the maximum number of layers, MXCLY has been increased from 6 to 120. A header that was printed to the screen has been suppressed.
2) Development of this Matlab/Octave suite of codes for running DISORT for layered atmospheres including both clouds and gases is by Penny Rowe, Von Walden, and Steven Neshyba, with helpful insights provided by Christopher Cox. While a small project, work has proceeded over several years, and thus was funded from a variety of sources, including the National Aeronautics and Space Administration (NASA) Research Opportunities in Space and Earth Sciences program (contract NNX08AF79G), the National Science Foundation (NSF) Idaho Experimental Program to Stimulate Competitive Research (EPSCoR), and NSF award ARC-1108451. Rowe also acknowledges support from USACH-DICYT and Neshyba acknowledges support from the University of Puget Sound and the National Science Foundation under grant CHE Ð 1306366.
3) Matlab code is written by Steven Neshyba and Penny Rowe (except as noted within). 



Please email me with any questions or problems: prowe@harbornet.com

Thank you and hope it is useful!
Dr. Penny M. Rowe
