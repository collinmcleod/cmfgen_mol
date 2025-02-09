!
! Routinine to compute the frequency grid for CMFGEN. Routine also
! allocates the vectors needed for the line data, sets the line data,
! and puts the line data into numerical order.
!
! The option FREQ_GRID_OPTION was introduced to allow an old frequency grid to be used.
! This could be useful when testing new options, and tesing the code against
! earlier versions. Its still possible that the grids may not be absolutely
! identical.
!
! If IOPT is 0, and anything else except 1, or 2 (only 1 is curently implmented),
! the default frequncy grid will used.
!
	SUBROUTINE SET_FREQUENCY_GRID(NU,FQW,LINES_THIS_FREQ,NU_EVAL_CONT,
	1               NCF,NCF_MAX,N_LINE_FREQ,
	1               OBS_FREQ,OBS,N_OBS,LUIN,IMPURITY_CODE)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE LINE_VEC_MOD
	IMPLICIT NONE
!
! Created 8-Jun-2004
!
	INTEGER NCF_MAX
	INTEGER NCF				!Total number of continuum points.
	INTEGER N_LINE_FREQ                     !Number of lines
	INTEGER N_OBS
	INTEGER LUIN
!
	REAL(KIND=LDP) NU(NCF_MAX)
	REAL(KIND=LDP) FQW(NCF_MAX)
	REAL(KIND=LDP) NU_EVAL_CONT(NCF_MAX)
	INTEGER LINES_THIS_FREQ(NCF_MAX)
!
	REAL(KIND=LDP) OBS_FREQ(NCF_MAX)
	REAL(KIND=LDP) OBS(NCF_MAX)
!
	LOGICAL IMPURITY_CODE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) LAMVACAIR
!
	REAL(KIND=LDP) NU_MAX_OBS
	REAL(KIND=LDP) NU_MIN_OBS
	REAL(KIND=LDP) T1,T2
!
	INTEGER I,J,K
	INTEGER ID
	INTEGER ML
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	LOGICAL FIRST
!
	CHARACTER(LEN=132) TEMP_CHAR
!
	LUER=ERROR_LU()
!
!
!
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species.
!
	ML=0			!Initialize line counter.
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNUP=2,ATM(ID)%NXzV_F
	      NUP= ATM(ID)%F_TO_S_XzV(MNUP)+ ATM(ID)%EQXzV-1
	      DO MNL=1,MNUP-1
	        NL= ATM(ID)%F_TO_S_XzV(MNL)+ ATM(ID)%EQXzV-1
	        IF(ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)ML=ML+1
	      END DO
	    END DO
	  END IF
	END DO
	N_LINE_FREQ=ML
!
! Now that we have the number of lines, we can allocate the needed memory.
! VEC_TRANS_NAME is allocated temporaruly so that we can output the full
! transitions name to TRANS_INFO.
!
	ALLOCATE (VEC_TRANS_NAME(N_LINE_FREQ),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_FREQ(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_OSCIL(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_EINA(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_DP_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_INDX(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_ID(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_NL(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_NUP(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_MNL_F(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_MNUP_F(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_INT_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_SPEC(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_CHAR_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_TRANS_TYPE(N_LINE_FREQ) ,STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( LINE_ST_INDX_IN_NU(N_LINE_FREQ) ,STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( LINE_END_INDX_IN_NU(N_LINE_FREQ) ,STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( LINE_LOC(N_LINE_FREQ) ,STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Unable to allocate memory in SET_FREQUENCY_GRID'
	  STOP
	END IF
!
! Now get the lines
!
	ML=0			!Initialize line counter.
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNUP=2,ATM(ID)%NXzV_F
	      NUP= ATM(ID)%F_TO_S_XzV(MNUP)+ ATM(ID)%EQXzV-1
	      DO MNL=1,MNUP-1
	        NL= ATM(ID)%F_TO_S_XzV(MNL)+ ATM(ID)%EQXzV-1
	        IF( ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)THEN
	          ML=ML+1
	          VEC_FREQ(ML)= ATM(ID)%EDGEXzV_F(MNL)- ATM(ID)%EDGEXzV_F(MNUP)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_ID(ML)=ID
	          VEC_NL(ML)=NL
	          VEC_NUP(ML)=NUP
	          VEC_MNL_F(ML)=MNL
	          VEC_MNUP_F(ML)=MNUP
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL,MNUP)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP,MNL)
	          VEC_TRANS_TYPE(ML)=ATM(ID)%XzV_TRANS_TYPE
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1           '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP))//'-'//
	1           TRIM( ATM(ID)%XzVLEVNAME_F(MNL))//')'
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
!
! 
!
! Get lines and arrange in numerically decreasing frequency. This will
! allow us to consider line overlap, and to include lines with continuum
! frequencies so that the can be handled automatically.
!
	N_LINE_FREQ=ML
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_ID,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
!
! Output all lines to TRANS_INFO. Usefule for diagnostic purposes.
!
	I=160	!Record length - allow for long names
	CALL GEN_ASCI_OPEN(LUIN,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LUIN,*)'!'
	  WRITE(LUIN,*)' Wavelengths are in air for Lambda > 2000 A'
	  WRITE(LUIN,*)'!'
	  WRITE(LUIN,*)
	1     '     I    NL_F  NUP_F        Nu',
	1     '       Lam(A)    /\V(km/s)    Transition'
	  WRITE(LUIN,
	1    '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,16X,A)')
	1         IONE,VEC_MNL_F(1),VEC_MNUP_F(1),
	1         VEC_FREQ(1),LAMVACAIR(VEC_FREQ(1)),
	1         TRIM(VEC_TRANS_NAME(VEC_INDX(1)))
	  DO ML=2,N_LINE_FREQ
	    T1=LAMVACAIR(VEC_FREQ(ML))
	    T2=2.998E+05_LDP*(VEC_FREQ(ML-1)-VEC_FREQ(ML))/VEC_FREQ(ML)
	    IF(T2 .GT. 2.998E+05_LDP)T2=2.998E+05_LDP
	    IF(T1 .LT. 1.0E+04_LDP)THEN
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    ELSE
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,1X,1P,E11.4,0P,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    END IF
	  END DO
	CLOSE(UNIT=LUIN)
	DEALLOCATE (VEC_TRANS_NAME)
!
! 
!
! GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
! The local species setting only takes precedence when it is set to NONE.
!
	IF(GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE')THEN
	  DO I=1,N_LINE_FREQ
	    VEC_TRANS_TYPE(I)=GLOBAL_LINE_SWITCH
	  END DO
	ELSE
	  DO I=1,N_LINE_FREQ
	    CALL SET_CASE_UP(VEC_TRANS_TYPE(I),IZERO,IZERO)
	  END DO
	END IF
!
! If desired, we can set transitions with:
!      wavelengths > FLUX_CAL_LAM_END (in A) to the SOBOLEV option.
!      wavelengths < FLUX_CAL_LAM_BEG (in A) to the SOBOLEV option.
!
! The region defined by FLUX_CAL_LAM_BEG < LAM < FLUX_CAL_LAM_END will be computed using
! transition types determined by the earlier species and global options.
!
! Option has 2 uses:
!
! 1. Allows use of SOBOLEV approximation in IR where details of radiative
!    transfer is unimportant. In this case FLUX_CAL_LAM_BEG should be set to zero.
! 2. Allows a full flux calculation to be done in a limited wavelength region
!    as defined by FLUX_CAL_LAM_END and FLUX_CAL_LAM_BEG.
!
	IF(SET_TRANS_TYPE_BY_LAM)THEN
	  IF(FLUX_CAL_LAM_END .LT. FLUX_CAL_LAM_BEG)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)'FLUX_CAL_LAM_END must be > FLUX_CAL_LAM_BEG'
	    STOP
	  END IF
	  IF( (.NOT. FLUX_CAL_ONLY) .AND. FLUX_CAL_LAM_BEG .NE. 0)THEN
	    WRITE(LUER,*)'WARNING in CMFGEN'
	    WRITE(LUER,*)'WARNING in CMFGEN'
	    WRITE(LUER,*)'FLUX_CAL_LAM_BEG is normally zero for non-FLUX'
	    WRITE(LUER,*)'calculations:'
	  END IF
	  GLOBAL_LINE_SWITCH='NONE'
	  T1=SPEED_OF_LIGHT()*1.0E-07_LDP
	  DO I=1,N_LINE_FREQ
	    IF(T1/VEC_FREQ(I) .GE. FLUX_CAL_LAM_END)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	    IF(T1/VEC_FREQ(I) .LE. FLUX_CAL_LAM_BEG)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	  END DO
	END IF
!
!
! Read in or calculate continuum frequency values.
!
	IF(RD_CONT_FREQ .OR. IMPURITY_CODE)THEN
	  CALL GEN_ASCI_OPEN(LUIN,'CFDAT','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening CFDAT in CMFGEN, IOS=',IOS
	      STOP
	    END IF
	    TEMP_CHAR=' '
	    DO WHILE(INDEX(TEMP_CHAR,'!Number of continuum frequencies')
	1                                              .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)TEMP_CHAR
	      IF(IOS .NE. 0)THEN
                WRITE(LUER,*)'Error reading in number of continuum ',
	1                        'frequencies in CMFGEN'
	        STOP
	      END IF
	    END DO
	    READ(TEMP_CHAR,*)NCF
	    IF(NCF .GT. NCF_MAX)THEN
	      WRITE(LUER,*)'Error - NCF > NCF_MAX in CMFGEN'
	      STOP
	    END IF
	    READ(LUIN,*)(NU(I),I=1,NCF)
	  CLOSE(UNIT=LUIN)
!
! Now need to set continuum frequencies.
!
	ELSE
	  FIRST=.TRUE.		!Check cross section at edge is non-zero.
	  NCF=0 		!Initialize number of continuum frequencies.
!
	  DO ID=1,NUM_IONS-1
	    CALL SET_EDGE_FREQ_V3(ID,OBS,NCF,NCF_MAX,
	1            ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1            ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV, ATM(ID)%N_XzV_PHOT)
	  END DO
!
	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS-1
	      CALL SET_X_FREQ_V2(OBS,NCF,NCF_MAX, MAX_CONT_FREQ,
	1             AT_NO(SPECIES_LNK(ID)),ATM(ID)%ZXzV,
	1             ATM(ID)%XzV_PRES, ATM(ID+1)%XzV_PRES)
	    END DO
	  END IF
!
! Now insert addition points into frequency array. WSCI is used as a
! work array - okay since of length NCF_MAX, and zeroed in QUADSE.
! OBS contains the bound-free edges - its contents are zero on
! subroutine exit. J is used as temporary variable for the number of
! frequencies transmitted to SET_CONT_FREQ. NCF is returned as the number
! of frequency points. FQW is used a an integer array for the sorting ---
! we know it has the correct length since it is the same size as NU.
! LUIN --- Used as temporary LU (opened and closed).
!
	  J=NCF
	  IF(FREQ_GRID_OPTION .EQ. 1)THEN
	    CALL SET_CONT_FREQ(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        J,NCF,NCF_MAX,LUIN)
	  ELSE IF(FREQ_GRID_OPTION .EQ. 2)THEN
	    CALL SET_CONT_FREQ_V3(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_XRAY_END,
	1                        J,NCF,NCF_MAX,LUIN)
	  ELSE
!
! Checks if BIG_FREQ_AMP was set to old definition.
!
	    IF(BIG_FREQ_AMP .LT. 0.0_LDP .OR. BIG_FREQ_AMP .GT. 1.0_LDP)THEN
	      BIG_FREQ_AMP=0.5_LDP
	    END IF
	    CALL SET_CONT_FREQ_V4(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_XRAY_END,
	1                        J,NCF,NCF_MAX,LUIN)
	  END IF
!
	END IF              !End set continuum if
!
! 
!
! We have found all lines. If we are doing a blanketing calculation for this
! line we insert them into the continuum frequency set, otherwise the
! line is not included.
!
	DO ML=1,NCF
	  FQW(ML)=NU(ML)	!FQW has temporary storage of continuum freq.
	END DO
	V_DOP=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
	CALL INS_LINE_V4(  NU,LINES_THIS_FREQ,I,NCF_MAX,
	1		  VEC_FREQ,VEC_TRANS_TYPE,
	1                 LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,
	1                 N_LINE_FREQ,FQW,NCF,
	1   		  V_DOP,FRAC_DOP,MAX_DOP,VINF,
	1                 dV_CMF_PROF,dV_CMF_WING,
	1                 ES_WING_EXT,R_CMF_WING_EXT  )
!
	K=NCF		!# of continuum frequencies: Need for DET_MAIN...
	NCF=I		!Revised
	CALL DET_MAIN_CONT_FREQ(NU,NCF,FQW,K,NU_EVAL_CONT,
	1         V_DOP,DELV_CONT,COMPUTE_ALL_CROSS)
!
	WRITE(LUER,'(A,T40,I6)')' Number of line frequencies is:',N_LINE_FREQ
	WRITE(LUER,'(A,T40,I6)')' Number of frequencies is:',NCF
	WRITE(LUER,*)' '
!
! Redefine frequency quadrature weights.
!
	IF(FREQ_GRID_OPTION.EQ. 1)THEN
	  CALL SMPTRP(NU,FQW,NCF)
	ELSE
	  CALL TRAPUNEQ(NU,FQW,NCF)
	END IF
	DO ML=1,NCF
	  FQW(ML)=FQW(ML)*1.0E+15_LDP
	END DO
!
! Revise number of lines so only those in frequency grid are included.
!
	I=N_LINE_FREQ
	DO ML=I,1,-1
	  N_LINE_FREQ=ML
	  IF(VEC_FREQ(ML) .GT. MIN_CONT_FREQ*(1.0_LDP+1.1_LDP*VINF/2.998E+05_LDP))EXIT
	END DO
	I=I-N_LINE_FREQ
	IF(I .NE. 0)THEN
	  WRITE(LUER,*)' Warning from SET_FREQUENCY_GRID'
	  WRITE(LUER,'(1X,I5,A,A)')I,' weak lines in ',
	1        'extreme IR will be ignored as outside continuum range.'
	  WRITE(LUER,*)'Min(Nu_CONT)=',NU(NCF)
	  WRITE(LUER,*)'Min(Nu_LINE)=',VEC_FREQ(I+N_LINE_FREQ)
	END IF
!
! Set observers frequencies. The slight fiddling in setting NU_MAX and NU_MIN
! is done so that the CMF frequencies encompass all observers frame
! frequencies. This allows computation of all observers fluxes allowing for
! velocity effects.
!
! We insert lines into the observers frame frequencies if they are being
! treated in blanketed mode, or if SOBN_FREQ_IN_OBS is set to TRUE.
!
	NU_MAX_OBS=NU(3)		!3 Ensures OBS_FREQ contained inside CMF band.
	T1=NU(NCF)*(1.0_LDP+2.0_LDP*VINF/2.998E+05_LDP)
	NU_MIN_OBS=MAX(NU(NCF-3),T1)
	CALL INS_LINE_OBS_V3(OBS_FREQ,N_OBS,NCF_MAX,
	1               VEC_FREQ,VEC_TRANS_TYPE,N_LINE_FREQ,SOB_FREQ_IN_OBS,
	1		NU_MAX_OBS,NU_MIN_OBS,VINF,
	1               dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP)
!
	RETURN
	END
