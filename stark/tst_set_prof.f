	PROGRAM TST_SET_PROF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
	INTEGER*4 ND
	INTEGER*4, PARAMETER :: ND_MAX=5
	INTEGER*4, PARAMETER :: NFREQ=201	!Must be odd
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: LU_STK=20
!
	REAL*8 ED_IN(ND_MAX)
	REAL*8 TEMP_IN(ND_MAX)
	REAL*8 CHIL(ND_MAX)
	REAL*8 VDOP_IN(ND_MAX)
	REAL*8 ZERO_VEC(ND_MAX)
!
	REAL*8 PROF(NFREQ,ND_MAX)
	REAL*8 PRO_VEC(ND_MAX)
	REAL*8 NU(NFREQ),LAM(NFREQ)
	REAL*8 VEL_KMS(NFREQ),LOGLAM(NFREQ)
	REAL*8 NORM(ND_MAX)
	REAL*8 Z_IN,AMASS_IN
	REAL*8 C_KMS
	REAL*8 VTURB
	REAL*8 VTURB_FIX
	REAL*8 GAM_NAT
	REAL*8 GAM_COL
	REAL*8 DOP_PROF_LIMIT
	REAL*8 VOIGT_PROF_LIMIT
!
	INTEGER*4 PROF_LIST_LOCATION
	INTEGER*4 NL,NUP
	CHARACTER*12 ION_ID
	CHARACTER*12 PROF_TYPE
!
	REAL*8 LAM_VAC,SPEED_OF_LIGHT
	EXTERNAL LAM_VAC,SPEED_OF_LIGHT
!
	REAL*8 START_FREQ
	REAL*8 DEL_NU
	REAL*8 NU_ZERO
	REAL*8 LAMBDA
!
	REAL*8 T1
	INTEGER*4 I,J,ML,ML_ST,ML_CUR
	LOGICAL PLT
!
! Set default values
!
	Z_IN=1.0
	AMASS_IN=1.0
	ND=ND_MAX
	DO I=1,ND_MAX
	  ED_IN(I)=10.0D0**( 9+I )
	END DO
	TEMP_IN(1:ND_MAX)=2.0
	GAM_COL=0.0D0
	GAM_NAT=0.0D0
	PROF_TYPE='LIST'
!
	CHIL(:)=1.0D+10
	ZERO_VEC(:)=0.0D0
!
	CALL INIT_PROF_MODULE(ND_MAX,10,NFREQ)
	CALL RD_STRK_LIST(20)
!
	ION_ID='HeI'
	AMASS_IN=4.0
	Z_IN=1.0
!
1000	CONTINUE
!
	CALL GEN_IN(ION_ID,'Ion identification [e.g., HeI]')
	IF(ION_ID .EQ. "")STOP
	CALL GEN_IN(AMASS_IN,'Atomic mass for species')
	CALL GEN_IN(Z_IN,'Effective atomic charge')
	CALL GEN_IN(LAMBDA,'Exact wavelength of transition [-ve for air]')
	IF(LAMBDA .LT. 0)THEN
	  LAMBDA=ABS(LAMBDA); LAMBDA=LAM_VAC(LAMBDA)
	END IF
	CALL GEN_IN(PROF_TYPE,'Profile type: LIST, DOPPLER, VOIGT')
	WRITE(6,*)'PROF_TYPE=',PROF_TYPE
!
	ED_IN(1:ND)=LOG10(ED_IN(1:ND))
	CALL GEN_IN(ED_IN,ND,ND_MAX,'Log10(Electron density[cm^-3])')
	ED_IN(1:ND)=10**(ED_IN(1:ND))
	CALL GEN_IN(TEMP_IN,ND,ND_MAX,'T (10^4) K')
!
	CALL GEN_IN(VTURB,'VTURB (km/s)')
	VDOP_IN(1:ND)=VTURB
	VTURB_FIX=VTURB
	WRITE(6,*)'Number of profiles to be computed=',ND
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	NU_ZERO=0.01D0*C_KMS/LAMBDA
!
! NL and NUP are not used.
!
	NL=2
	NUP=3
!
! Compute the frequecny grid
!
	DEL_NU=1.5/C_KMS/SQRT(AMASS_IN)
	T1=1.1
	NU(NFREQ/2+1)=0.0D0
	I=NFREQ/2
	J=NFREQ/2+2
	DO ML=1,NFREQ/2
	  DEL_NU=DEL_NU*T1
	  NU(I)=NU(I+1)+DEL_NU ; I=I-1
	  NU(J)=NU(J-1)-DEL_NU ; J=J+1
	END DO
	LAM(:)=-0.01D0*C_KMS/NU_ZERO*NU(:)/(1.0D0+NU(:))
	NU(:)=NU_ZERO*(1.0D0+NU(:))
	VEL_KMS(:)=C_KMS*(NU_ZERO-NU(:))/NU_ZERO
!
	WRITE(6,*)'Call SET_PROF_LIMITS'
	CALL SET_PROF_LIMITS_V2(START_FREQ,VTURB_FIX,
	1                     CHIL,ED_IN,TEMP_IN,VDOP_IN,ND,
	1                     PROF_TYPE,PROF_LIST_LOCATION,
	1                     NU_ZERO,NL,NUP,
	1                     ION_ID,AMASS_IN,Z_IN,
	1                     GAM_NAT,GAM_COL,VTURB_FIX,
	1                     DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,L_FALSE)
	WRITE(6,*)'Ended call SET_PROF_LIMITS'
	WRITE(6,*)'PROF_TYPE=',PROF_TYPE
	WRITE(6,*)START_FREQ,PROF_LIST_LOCATION
	WRITE(6,*)NU_ZERO,NU(1),NU(NFREQ)
!
	CALL TUNE(1,'SET_PROF')
	PROF(:,:)=0.0D0
	ML_ST=1
	DO ML=1,NFREQ
	  ML_CUR=ML
	  CALL SET_PROF_V3(PRO_VEC,NU,ML_CUR,ML_ST,NFREQ,
	1               ED_IN,ZERO_VEC,ZERO_VEC,TEMP_IN,VDOP_IN,ND,
	1               PROF_TYPE,PROF_LIST_LOCATION,
	1               NU_ZERO,NL,NUP,AMASS_IN,Z_IN,
	1               GAM_NAT,GAM_COL,VTURB_FIX,L_FALSE,L_FALSE,LU_STK)
!	1               GAM_NAT,GAM_COL,VTURB_FIX,L_FALSE,L_TRUE,LU_STK)
	  PROF(ML,1:ND)=PRO_VEC(1:ND)
	END DO
	CALL TUNE(2,'SET_PROF')
	CALL TUNE(3,' ')
C
C Check profile normalization.
C
	DO I=1,ND
	  NORM(I)=0.0D0
	  DO ML=1,NFREQ-1
	    NORM(I)=NORM(I)+0.5D0*(NU(ML)-NU(ML+1))*
	1                               (PROF(ML,I)+PROF(ML+1,I))
   	  END DO
	  WRITE(6,'(A,I2,A,1P,E15.8)')' I=',I,'   T1=',1.0D+15*NORM(I)
	END DO
C
	PROF(:,:)=PROF(:,:)*1.0E+12
	DO I=1,ND
	  CALL DP_CURVE(NFREQ,VEL_KMS,PROF(1,I))
	END DO
	CALL GRAMON_PGPLOT('\gDV/c','Prof',' ',' ')
!
	PLT=.FALSE.
	CALL GEN_IN(PLT,'Plot log of profile?')
	IF(PLT)THEN
	  PROF(:,:)=LOG10( PROF(:,:)+1.0D-40 )-12
	  DO I=1,ND
	    CALL DP_CURVE(NFREQ,VEL_KMS,PROF(1,I))
	  END DO
	  CALL GRAMON_PGPLOT('\gDV/c','Prof',' ',' ')
	END IF
!
	PLT=.FALSE.
	CALL GEN_IN(PLT,'Plot log of profile in wavelength space?')
	IF(PLT)THEN
	  DO I=1,ND
	    CALL DP_CURVE(NFREQ,LAM,PROF(1,I))
	  END DO
	  CALL GRAMON_PGPLOT('Lam(\A)','Prof',' ',' ')
	END IF
!
	GOTO 1000
!
!	LOGLAM(4:NFREQ)=LOG10(LAM(4:NFREQ))
!	DO I=1,ND
!	  CALL DP_CURVE(NFREQ-3,LOGLAM(4),PROF(4,I))
!	END DO
!	CALL GRAMON_PGPLOT('Log \gl(\A)','Prof',' ',' ')
C
	END
