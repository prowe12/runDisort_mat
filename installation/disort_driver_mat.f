c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      PROGRAM  disort_driver

c    Runs DISORT using a wrapper to get input and write output data.

c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c    BDREF:    Sets bidirectional reflectance of lower boundary

c    GETMOM:   Sets phase function Legendre coefficients

c    PRTFIN:   Prints fluxes and intensities and their ratios to
c              the correct values

c    CHEKDO:   Data block containing correct fluxes and intensities

c    RATIO :   Ratio of calculated to correct value with underflow
c              and overflow protection (kept in file DISORT.f)

c+---------------------------------------------------------------------+
c
c                 ** DISORT I/O specifications **
c
c     Increased MAXUMU, PMR 10 May, 2012

      INTEGER  MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, CLDLYR, NCLDLYR
      PARAMETER ( MAXCLY = 120, MAXMOM = 1000, MAXPHI = 1, MAXULV = 2,
     &            MAXUMU = 10 )
      CHARACTER  HEADER*127
      LOGICAL  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU, DEBUG
      INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     &         PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ), 
     &         PMOM_CLD(0:MAXMOM),
     &         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )

c+---------------------------------------------------------------------+

c                     ** Correct answers **

      INTEGER  MXPROB, MXCASE, MAXTAU, MAXMU, MAXAZ
      PARAMETER  ( MXPROB = 9, MXCASE = 8, MAXTAU = 5, MAXMU = 10,
     &             MAXAZ = 3 )
      REAL  TSTFIR( MAXTAU, MXCASE, MXPROB ),
     &      TSTFDN( MAXTAU, MXCASE, MXPROB ),
     &      TSTFUP( MAXTAU, MXCASE, MXPROB ),
     &      TSTDFD( MAXTAU, MXCASE, MXPROB ),
     &      TSTUU ( MAXTAU, MAXMU, MAXAZ, MXCASE, MXPROB )
      COMMON / DOCHEK / TSTFIR, TSTFDN, TSTFUP, TSTDFD, TSTUU

c+---------------------------------------------------------------------+

      INTEGER  MXTAU, MXMU, MXPHI
      PARAMETER     ( MXTAU = 5, MXMU = 10, MXPHI = 3 )
      CHARACTER*1   ABC(18)*1, TITLE*100, BLANKS*3
      CHARACTER*9   dumStr(45)*9
      LOGICAL       DOPROB( 13 )
      INTEGER       ICAS, IOD, ISS, IU, J, K, LC, LENTIT, LU, NPROB
      REAL          CMPFIR( MXTAU ), CMPFDN( MXTAU ), CMPFUP( MXTAU ),
     &              CMPDFD( MXTAU ), CMPUU ( MXTAU, MXMU, MXPHI ),
     &              PI
      CHARACTER*100 OUTFILE, INFILE

c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX
c     ..
      DATA  ABC / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
     &            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r' /,
     &      PRNT / .FALSE., 3*.FALSE., .FALSE. /, ACCUR / 0.0 /,
     &      BLANKS / '   ' /,  DOPROB / .TRUE., 12*.FALSE. /
c     &      badval / -999 /
c



c   ..Namelist stuff..
c     Added NUMU, UMU, NTAU, UTAU, PMR 10 May 2012
      NAMELIST /DISORTINPUT/ 
     &                 NLYR, DTAUC, SSALB, NMOM, NCLDLYR,
     &                 TEMPER, WVNMLO, WVNMHI, USRTAU, 
     &                 NUMU, UMU, NTAU, UTAU,
     &                 NSTR, USRANG, NPHI, PHI, IBCND,
     &                 FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO,
     &                 BTEMP,TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, 
     &                 OUTFILE, DEBUG
      NAMELIST /PMOMINPUT/ 
     &                 CLDLYR, PMOM_CLD

	 
c    .. Initialize variables to zeros
      DO J=1,MAXCLY
        SSALB(J) = 0.0
      ENDDO	  
	  DO J=0,MAXCLY
	    TEMPER(J) = 0.0
	  ENDDO
	  DO J=1,MAXCLY
	    DTAUC(J) = 0.0
	  ENDDO
	  DO J=0,MAXCLY
	    TEMPER(J) = 0.0
	  ENDDO
      DO I=0,MAXMOM
        DO J=1,MAXCLY
          PMOM(I,J) = 0.0
        ENDDO
      ENDDO

c   ..Set everything in namelist to default values..
c      NLYR   = badvalbadval
c      DTAUC(1)  = YMDHms
c      SSALB(1)  = badval
c      NMOM(1)   = badval
c      PMOM(0,1)   = badval
c      TEMPER(1) = badval
c      WVNMLO = badval
c      WVNMHI = badval
c      USRTAU, NTAU, UTAU, NSTR,
c      USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
c      UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
c      TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
c      HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
c      MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
c      ALBMED, TRNMED )


c   ..Get the input data via namelist..
      ICOUNT = COMMAND_ARGUMENT_COUNT()
      IF (ICOUNT .GT. 0) THEN
        CALL GETARG(1,INFILE)
      ELSE
        INFILE = 'disortinput.nml'
      ENDIF
c      WRITE (6,*) INFILE
      OPEN(7,FILE= INFILE)
	  
	  
c   ..Read the main namelist..	  
      READ(7,NML=DISORTINPUT)


c   ..Write out the namelist values
c      WRITE(6,NML=DISORTINPUT)
	  
     
c   ..Fill in the moments from the namelist
      DO ICLDLYR=1,NCLDLYR
	    READ(7,NML=PMOMINPUT)
c		WRITE(6,*) '---------'
c		WRITE(6,*) 'For cldlyr = ',CLDLYR
        DO I=0,NMOM
           PMOM(I, CLDLYR) = PMOM_CLD(I)
c           WRITE(6,*) PMOM(I, CLDLYR)
        ENDDO
      ENDDO


c     Close the namelist file
      CLOSE(7)
	
c   ..Some debugging
c      if(debug) then
c      	write(6,*) 'UTAU = ', UTAU
c      	write(6,*) 'sum(DTAUC) = ', sum(DTAUC(1:NLYR))
c      endif
			
c   ..Don't let number streams exceed number of legendre moments
      IF (NSTR .GT. NMOM) THEN
		NMOM = NSTR
c		WRITE (6,*) 'Reducing NSTR to NMOM-1 ...'
      ENDIF

c    .. *Do Not* Hardwire it to give values at TOA and surface
c    .. or viewing angles 
c    .. PMR May 12, 2012
c      NTAU = 2
c      UTAU(1) = sum(DTAUC(1:NLYR))
c	  UTAU(2) = 0
c	  NUMU    = 2
c	  UMU(1)  = -1
c	  UMU(2)  = 1
c      WRITE(6,*) NUMU
c	  WRITE(6,*) UMU


c     Specify the header, although it won't be printed out.
      HEADER = 'Run from namelist.'

c
c   ..Call disort..
      CALL  DISORT( 
     &                 NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )
c

c      print (6,*), TEMIS,TEMPER
c
c    ..Check if any inputs were not set, and bomb if not..
c      IF ( (NLYR.eq.badval) .or. (TEMPER(1).eq.badval) ) then
c        print (6,*), TEMIS
c      ENDIF


c   ..Output..
c   ..Note: UU( MAXUMU, MAXULV, MAXPHI ), MAXULV hardwired to 2
      OPEN(8,FILE=OUTFILE,STATUS='unknown',FORM='FORMATTED')
      WRITE (8,*) RFLDN(1), '% RFLDN: Diffuse down flux, sfc'
      WRITE (8,*) FLUP(2), '% FLUP: Diffuse upward flux, toa'
      DO I=1,NUMU
        WRITE (8,*) UU(I,1,1), '% UU: down radiance(z. ang), sfc'
	  ENDDO
      DO I=1,NUMU
        WRITE (8,*) UU(I,2,1), '% UU: up radiance(z. ang), toa'
	  ENDDO
      CLOSE(8)
      STOP
      END

c     OTHER POSSIBLE OPTIONS (not included now)
c      WRITE (8,*) RFLDN(2), '% RFLDN: Diffuse down flux, toa'
c      WRITE (8,*) FLUP(1), '% FLUP: Diffuse upward flux, sfc'
c      WRITE (8,*) UU(2,1,1), '% UU: Upwelling Radiance, sfc'
c      WRITE (8,*) UU(1,2,1), '% UU: Downwelling Radiance, toa'
c      WRITE (8,*) RFLDIR(1), '% RFLDIR: Direct down flux, sfc'
c      WRITE (8,*) RFLDIR(2), '% RFLDIR: Direct down flux, toa'



