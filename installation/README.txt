This folder includes the standard DISORT version 2.0 available 
at Wiscombe's NASA ftp site: ftp://climate1.gsfc.nasa.gov/wiscombe/. 

The purpose is to make the program available in a single file, and to 
provide a script to automatically compile and run the program as well 
as an example Makefile.

No modifications were made to the program itself; it is the standard 
version 2.0 beta release code. If you have any issues/comments/questions 
please send me (Snorre Stamnes) an email at sstamnes@stevens.edu.

The following scripts expect the gfortran compiler (versions 4.x and 3.x)
to be installed. But if you have a different compiler it will most 
likely work (substitute gfortran for your compiler, e.g. ifort.)

This package should run on Linux, Mac OS X, *BSD, Solaris, etc.

Compile and run DISORT 2.0 BETA by executing the run_disort.sh script:
 ./run_disort.sh

To test your system execute the test_disort.sh script
(assuming no modifications were made to DISOTEST.f or any other files):
 ./test_disort.sh

Alternatively, for convenience a Makefile has been created. 
To compile using the Makefile type "gmake" 
("make" may be used if gmake is not available)
--------------------------------------
The standard DISORT Readme file:

README:  DISORT version 2.0 (beta release)

This new version of DISORT is being released as a beta version.  While we have
every confidence that it is working correctly, and while it passes all the
test problems, there may still be subtle problems that we would appreciate
hearing about (our e-mail addresses are in DISORT.doc).

The major changes from version 1.2 are:

* Nakajima intensity corrections, allowing use of many fewer streams (NSTR) to
obtain the same intensity accuracy

* a proper lower-boundary bidirectional reflectance option

A version 2.1, incorporating LAPACK replacements for LINPACK, may be issued
soon.
--------------------------------------
The standard DISORT RDIMACH Readme file:

{R,D,I}1MACH:  The routines we needed but hated
W. Wiscombe (wiscombe@gsfc.nasa.gov)
July 1998

   The machine-constant routines R1MACH, D1MACH, I1MACH caused more 
problems for me and for users of my programs than any others.  Their 
functions were simple, but people just had a hard time getting the 
versions distributed on netlib (and which I re-distributed) to work 
correctly.

   At this point in time, it no longer makes sense to distribute or 
use these routines.  Fortran-90 contains intrinsic functions which 
contain all the functionality of R1MACH and D1MACH, and almost all of 
I1MACH.  Eric Grosse of Bell Labs has been kind enough to provide versions 
of {R,D,I}1MACH which use these new F-90 intrinsic functions.  I have 
slightly edited his routines.  The package is called RDI1MACH.f and 
is self-documenting.

   Fortran-90 compilers have matured on most platforms and we 
highly recommend buying/using them (see http://www.fortran.com/fortran/).
There are even some free Fortran-90 subset compilers (http://www.lahey.com/).
Soon, Fortran-77 compilers may no longer be supported.  Fortran-90 
compilers can compile any Fortran-77 programs or routines, so once they 
became reliable it is inevitable that support for f77 compilers will wither.

   Since Fortran-90 is entirely backward compatible with Fortran-77, 
you need not use any other feature of Fortran-90; you can just use 
RDI1MACH.f inside an old Fortran-77 program.  Remember that the .f
extension only means the file is in fixed source form, not that it is
pure Fortran-77.  In fact, any .f file can contain f90 constructs.

  Those without access to Fortran-90 compilers can obtain the old 
versions of {R,D,I}1MACH at http://www.netlib.org/.
