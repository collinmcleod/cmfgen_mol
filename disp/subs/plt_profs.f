!
! Subroutine designed to plot intrinsinc line profiles for a range of
! temperatures and densities.
!
! FULL_STRK_LIST and Stark files must be linked for the routine to work.
!
	SUBROUTINE PLT_PROFS
	USE SET_KIND_MODULE
	USE MOD_DISP
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created: 20-May-2015 : Based on $cmfdist/stark/tst_set_prof.f
!
	INTEGER ND
	INTEGER, PARAMETER :: ND_MAX=6		!Maximum number of profiles to be plotted.
	INTEGER, PARAMETER :: NFREQ=201		!Must be odd
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL(KIND=LDP) ED_IN(ND_MAX)
	REAL(KIND=LDP) TEMP_IN(ND_MAX)
	REAL(KIND=LDP) CHIL(ND_MAX)
	REAL(KIND=LDP) VTURB_IN(ND_MAX)
	REAL(KIND=LDP) ZERO_VEC(ND_MAX)
!
	REAL(KIND=LDP) PROF(NFREQ,ND_MAX)
	REAL(KIND=LDP) PRO_VEC(ND_MAX)
	REAL(KIND=LDP) NU(NFREQ),LAM(NFREQ)
	REAL(KIND=LDP) VEL_KMS(NFREQ),LOGLAM(NFREQ)
	REAL(KIND=LDP) NORM(ND_MAX)
	REAL(KIND=LDP) Z_IN
	REAL(KIND=LDP) AMASS_IN
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) VTURB
	REAL(KIND=LDP) VTURB_FIX
	REAL(KIND=LDP) GAM_NAT
	REAL(KIND=LDP) GAM_COL
	REAL(KIND=LDP) DOP_PROF_LIMIT
	REAL(KIND=LDP) VOIGT_PROF_LIMIT
	REAL(KIND=LDP) VEC_DOP_VMIN
	REAL(KIND=LDP) MAX_PROF_ED
	REAL(KIND=LDP) V_PROF_LIM
	REAL(KIND=LDP) LOC_C4
	REAL(KIND=LDP) LOC_ARAD
!
	INTEGER PROF_LIST_LOCATION
	INTEGER NL,NUP
	CHARACTER*12 LOC_ION_ID
	CHARACTER*12 PROF_TYPE
!
	REAL(KIND=LDP) LAM_VAC,LAMVACAIR,SPEED_OF_LIGHT
	EXTERNAL LAM_VAC,LAMVACAIR,SPEED_OF_LIGHT
!
	REAL(KIND=LDP) START_FREQ
	REAL(KIND=LDP) DEL_NU
	REAL(KIND=LDP) NU_ZERO
	REAL(KIND=LDP) LAMBDA
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TWO_PI
	REAL(KIND=LDP) WAVE
	INTEGER ID,ISPEC
	INTEGER I,J,ML,ML_ST,ML_CUR
	INTEGER LOOP_COUNT
	INTEGER LU_STK
	LOGICAL PLT
	LOGICAL FOUND_TRANS
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	CHARACTER(LEN=5) CHECK_DF
!
! Set default values
!
	TWO_PI=8.0D0*ATAN(1.0D0)
	Z_IN=1.0D0
	AMASS_IN=1.0D0
	ND=ND_MAX
	NL=1; NUP=2
	DO I=1,ND_MAX
	  ED_IN(I)=10.0D0**( 11+I )
	END DO
	TEMP_IN(1:ND_MAX)=2.0
	GAM_COL=0.0D0
	GAM_NAT=0.0D0
	PROF_TYPE='LIST_VGT'
	CALL GET_LU(LU_STK,'PLT_PROFS')
!
	CHIL(:)=1.0D+10
	ZERO_VEC(:)=0.0D0
!
	IF(FIRST_TIME)THEN
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'Make sure line profile data files are linked'
	  CALL GEN_IN(CHECK_DF,'Input any character to continue')
	  WRITE(6,*)DEF_PEN
	  CALL INIT_PROF_MODULE(ND_MAX,10,NFREQ)
	  CALL RD_STRK_LIST(LU_STK)
	  FIRST_TIME=.FALSE.
	END IF
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LOC_ION_ID='HeI'
	AMASS_IN=4.0
	Z_IN=1.0
!
1000	CONTINUE
!
	CALL GEN_IN(LOC_ION_ID,'Ion identification [e.g., HeI; EX(IT) or "" to quit]')
	IF(LOC_ION_ID .EQ. "" .OR. LOC_ION_ID(1:2) .EQ. 'EX')RETURN
	ID=0
	DO I=1,NUM_IONS
	  IF(LOC_ION_ID .EQ. ION_ID(I))THEN
	    ID=I
	    EXIT
	  END IF
	END DO
	IF(ID .EQ. 0)GOTO 1000
!
	ISPEC=SPECIES_LNK(ID)
	AMASS_IN=AT_MASS(ISPEC)
	Z_IN=ATM(ID)%ZxZV
!
	CALL GEN_IN(NL,'Lower level of transition [0 to for LAM promt]')
	LOOP_COUNT=1
	DO WHILE(NL .EQ. 0)
	  CALL GEN_IN(WAVE,'Approximate wavelength of line')
	  FOUND_TRANS=.FALSE.
	  DO WHILE(.NOT. FOUND_TRANS)
	    DO I=1,ATM(ID)%NXzV_F
	      DO J=I+1,ATM(ID)%NXzV_F
	        NU_ZERO=ATM(ID)%EDGEXzV_F(I)-ATM(ID)%EDGEXzV_F(J)
	        T1=0.01D0*C_KMS/NU_ZERO
	        IF(ABS(WAVE-T1) .LT. 5.0D0*LOOP_COUNT .AND. ATM(ID)%AXzV_F(I,J) .NE. 0)THEN
	          WRITE(6,'(2I5,2F10.2,3X,3A)')I,J,T1,LAMVACAIR(NU_ZERO),
	1                 TRIM(ATM(ID)%XzVLEVNAME_F(I)),'-',
	1                 TRIM(ATM(ID)%XzVLEVNAME_F(J))
	          FOUND_TRANS=.TRUE.
	        END IF
	      END DO
	    END DO
	    LOOP_COUNT=LOOP_COUNT+1
	  END DO
	  CALL GEN_IN(NL,'Lower level of transition')
	END DO
!
	CALL GEN_IN(NUP,'Upper level of transition')
	IF(  NL  .LT. 0 .OR. NL  .GT. ATM(ID)%NXzV_F .OR.
	1    NUP .LT. 0 .OR. NUP .GT. ATM(ID)%NXzV_F)THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A,I5)')' Error -- invalid NL or NUP: Max vlaue is',ATM(ID)%NXzV_F
	  WRITE(6,'(A)')DEF_PEN
	  GOTO 1000
	END IF
	NU_ZERO=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	IF(ATM(ID)%AXzV_F(NL,NUP) .EQ. 0.0D0)THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')' Error -- oscillator strength is zero'
	  WRITE(6,'(A)')DEF_PEN
	  GOTO 1000
	END IF
	PROF_TYPE='LIST_VGT'
	CALL GEN_IN(PROF_TYPE,'Profile type: LIST, LIST_VGT, DOPPLER, VOIGT')
!
	ED_IN(1:ND)=LOG10(ED_IN(1:ND))
	CALL GEN_IN(ED_IN,ND,ND_MAX,'Log10(Electron density[cm^-3])')
	ED_IN(1:ND)=10**(ED_IN(1:ND))
	CALL GEN_IN(TEMP_IN,ND,ND_MAX,'T (10^4) K')
!
	VTURB=0.0D0
	CALL GEN_IN(VTURB_FIX,'VTURB (km/s)')
	IF(VTURB_FIX .EQ. 0.0D0 .AND. LOC_ION_ID(1:1) .EQ. 'H')THEN
	   VTURB_FIX=1.0D0
	   WRITE(6,*)'We set VTURB to be 1 km/s as it cannot be zero for H or He'
	END IF
	VTURB_IN(1:ND)=VTURB_FIX
!
	T1=MINVAL(TEMP_IN(1:ND))
	VEC_DOP_VMIN=12.85*SQRT(T1/AMASS_IN+VTURB_FIX*VTURB_FIX)
	MAX_PROF_ED=1.0D+20
	V_PROF_LIM=5000.0D0/MIN(SQRT(AMASS_IN),5.0D0)
!
	LOC_ARAD=ATM(ID)%ARAD(NL)+ATM(ID)%ARAD(NUP)
	LOC_C4=ATM(ID)%GAM2(NL)+ATM(ID)%GAM2(NUP)
!
! Compute the frequency grid
!
	DEL_NU=1.5/C_KMS/SQRT(AMASS_IN)
	T1=1.1D0
	DO I=1,10
	  T1=(V_PROF_LIM*(T1-1))**(2.0D0/NFREQ)
	END DO
	WRITE(6,*)' '
	WRITE(6,*)'Frequency scaling factor for frequency grid is',T1
	WRITE(6,*)' '
!
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
	CALL SET_PROF_LIMITS_V4(START_FREQ,VEC_DOP_VMIN,
	1                     CHIL,ED_IN,TEMP_IN,VTURB_IN,ND,
	1                     PROF_TYPE,PROF_LIST_LOCATION,
	1                     NU_ZERO,NL,NUP,LOC_ION_ID,AMASS_IN,Z_IN,
	1                     LOC_ARAD,LOC_C4,TEMP_IN(1),AMASS_IN,VTURB_FIX,
	1                     DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,V_PROF_LIM,MAX_PROF_ED,L_FALSE)
	WRITE(6,*)'Ended call SET_PROF_LIMITS'
!
	WAVE=0.01D0*C_KMS/NU_ZERO
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(X,A,T8,2A)')LOC_ION_ID, 'PROF_TYPE=',PROF_TYPE
	WRITE(6,'(T8,4A,F9.2,A)')TRIM(ATM(ID)%XzVLEVNAME_F(NL)),' - ',
	1                     TRIM(ATM(ID)%XzVLEVNAME_F(NUP)),
	1                       ' at ',WAVE,'Ang'
	IF(PROF_TYPE .EQ. 'VOIGT')THEN
	  IF(LOC_C4 .NE. 0.0D0)THEN
	    T1=1.55D+04*1.0D+16*(ABS(LOC_C4)**0.667D0)
	    WRITE(6,'(A,ES8.2,11X,A,ES9.2,A,ES8.2)')' Gam(A)=',LOC_ARAD,'C4=',LOC_C4,'  GAM_COL(10^16)=',T1
	    T2=1.0D-15*C_KMS*LOC_ARAD/TWO_PI/NU_ZERO
	    T1=1.0D-15*C_KMS*T1/TWO_PI/NU_ZERO
	    WRITE(6,'(A,ES8.2,A,ES8.2,A)')'  dv(A)=',T2,'    dv(10^16)=',T1,'   (km/s)'
	  ELSE
	    T2=1.0D-15*C_KMS*LOC_ARAD/TWO_PI/NU_ZERO
	    WRITE(6,'(A,ES8.2,A,ES8.2,A)')' Gam(A)=',LOC_ARAD,'  dv(A)=',T2,'km/s'
	  END IF
	END IF
!
	CALL TUNE(1,'SET_PROF')
	PROF(:,:)=0.0D0
	ML_ST=1
	DO ML=1,NFREQ
	  ML_CUR=ML
	  CALL SET_PROF_V5(PRO_VEC,NU,ML_CUR,ML_ST,NFREQ,
	1               ED_IN,ZERO_VEC,ZERO_VEC,TEMP_IN,VTURB_IN,ND,
	1               PROF_TYPE,PROF_LIST_LOCATION,
	1               NU_ZERO,NL,NUP,AMASS_IN,Z_IN,
	1               LOC_ARAD,LOC_C4,
	1               TEMP_IN(1),VTURB_FIX,AMASS_IN,MAX_PROF_ED,
	1               L_FALSE,L_FALSE,LU_STK)
	  PROF(ML,1:ND)=PRO_VEC(1:ND)
	END DO
	CALL TUNE(2,'SET_PROF')
	CALL TUNE(3,' ')
!
! Check profile normalization.
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,*)'Checking profile normalization',RED_PEN
	DO I=1,ND
	  NORM(I)=0.0D0
	  DO ML=1,NFREQ-1
	    NORM(I)=NORM(I)+0.5D0*(NU(ML)-NU(ML+1))*
	1                               (PROF(ML,I)+PROF(ML+1,I))
   	  END DO
	  WRITE(6,'(A,I2,A,1P,E15.8)')'   I=',I,'    Int Phi(v) dv=',1.0D+15*NORM(I)
	END DO
	WRITE(6,*)DEF_PEN
!
	PROF(:,:)=PROF(:,:)*1.0E+12
	DO I=1,ND
	  CALL DP_CURVE(NFREQ,VEL_KMS,PROF(1,I))
	END DO
	CALL GRAMON_PGPLOT('V(kms)','10\u12\d \gF(\gn)',' ',' ')
!
	PLT=.FALSE.
	CALL GEN_IN(PLT,'Plot profile in wavelength space?')
	IF(PLT)THEN
	  DO I=1,ND
	    CALL DP_CURVE(NFREQ,LAM,PROF(1,I))
	  END DO
	  CALL GRAMON_PGPLOT('Lam(\A)','10\u12\d \gP(\n)',' ',' ')
	END IF
!
	GOTO 1000
!
	END
