!
! Subroutine to compute the CMF line frequencies. The spacing of the line frequencies
! is determined by passed parameters which specify the resonance zone extent, spacing
! in Doppler widths etc.
!
! The CMF line frequencies are merged with the previously adopted continuum frequency set.
!  Unnecessary frequencies are eliminated. Edge frequencies are retained.
!
	SUBROUTINE INS_LINE_V4(
	1		FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1		NU_LINE,TRANS_TYPE,
	1               LINE_ST_INDX,LINE_END_INDX,N_LINES,
	1		NU_CONT,NCF,
	1		V_DOP,FRAC_DOP,MAX_DOP,VINF,dV_CMF_PROF,
	1               dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 19-Aug-2004 : Bug fix with definition of continuum edges. This will not have
!                         a large effect on models, but will change the final adopted
!                         frequency grid. We also set the test for edge frequencies to
!                         to 1 part in 10^{-8} which is adequate since 2*10^{-11} is
!                         used in SET_CONT_FREQ_V# to define the edges. Previously we
!                         adopted 10^{-6}. Practically this will have little effect.
!
! Altered 12-Dec-1997 : If no lines are to be inserted, FREQ is set to
!                         NU_CONT, and there is an immediate return.
!                         All edge frequencies [defined by
!                             ABS( NU(ML)/MU(M+1)-1 ) < 0.000001 ]
!                         are now retained in the list.
!                         Changes based/requested by PACO.
!
! Altered 15-Aug-1997 : Bug fix: Inserting to many pints in profile zone
!                         when some lines not in BLANK mode. LST_LN_INDX
!                         variable now used.
! Altered 05-Dec-1996 : NDOP made integer (bug fix)
! Altered 30-May-1996 : Warning about lines below minimum continuum frequency
!                         installed.
! Altered 15-May-1996 : ES_WING_EXT introduced (Now V4). May be zero if
!                         coherent scattering (in km/s).
!                       Some cleaning done to minimize frequencies which are
!                         unnecessarily close.
!                       R_CMF_WING_EXT now refers only to the wing caused by
!                         coherent scattering.
!
! Altered 28-Feb-1995 : dv_CMF_WING,CMF_WING_EXT inserted (_V2)
! Altered 02-Feb-1995 : Not correctly inserting all continuum points.
!
	INTEGER NCF,NFREQ_MAX,N_LINES
	INTEGER NFREQ				!Returned
!
! Vectors returned by subroutine:
!
! Line+continuum frequencies
!
	REAL(KIND=LDP) FREQ(NFREQ_MAX)			!Continuum frequencies
	INTEGER LINES_THIS_FREQ(NFREQ_MAX)	!Indicates that this frequency
						!  has line contributions,
!
	INTEGER LINE_ST_INDX(N_LINES)		!Start index for the line
						!  in the NEW frequency array.
	INTEGER LINE_END_INDX(N_LINES)	!End index for the line
						! in the NEW frequency array.
!
! Passed vectors.
!
	REAL(KIND=LDP) NU_CONT(NCF)		!Continuum frequencies
	REAL(KIND=LDP) NU_LINE(N_LINES)		!Line frequencies
	CHARACTER*(*) TRANS_TYPE(N_LINES)
!
! Passed constants:
	REAL(KIND=LDP) VINF		!Terminal velocity of wind.
	REAL(KIND=LDP) V_DOP		!Doppler velocity (km/s).
	REAL(KIND=LDP) FRAC_DOP		!Indicates dNU across line in Doppler widths.
	REAL(KIND=LDP) MAX_DOP		!Half the extent of intrinsic profile
				!  in Doppler widths,
	REAL(KIND=LDP) dV_CMF_PROF	!Indicate spacing in profile but outside
                                !  resonance zone (in km/s).
	REAL(KIND=LDP) dV_CMF_WING	!Indicate spacing in wings (i.e. outside
				!  intrinsic profile) (in km/s).
!
! R_CMF_WING_EXT indicates how far profile should extend beyond red edge
! of RESONANCE zone. This will normally be greater than 2VINF as electron
! scattering by cold electrons broadens the line to the red side. Value
! is set on the basis of cold electrons.
!
! ES_WING_EXT is useful when have non-coherent electron scattering.
! Used for both blue and red sides of the line profile.
!
	REAL(KIND=LDP) ES_WING_EXT
	REAL(KIND=LDP) R_CMF_WING_EXT
!
! Local variables.
!
	REAL(KIND=LDP) dNU_on_NU	!Actual spacing used across intrinsic line
				!  profile given by dNU =NU*dNU_on_NU
	INTEGER NDOP		!Number of frequencies across intrinsic
				!  profile.
	REAL(KIND=LDP) RES_EXTENT	!Maximum frequency in line is NU*RES_EXTENT
!
	REAL(KIND=LDP) BLUE_WING_EXT	!In km/s
	REAL(KIND=LDP) RED_WING_EXT
	REAL(KIND=LDP) APP_RES_EXT	!No units.
	REAL(KIND=LDP) APP_ESBW_EXT
	REAL(KIND=LDP) EDGE_SEP_FAC
	REAL(KIND=LDP) MIN_FREQ_RAT
!
	INTEGER INDX		!Current frequency index.
	INTEGER LN_INDX	!Current line whose frequencies we are
				!   installing.
	INTEGER LST_LN_INDX	!Last line whose frequencies we
				!   installed.
!
	INTEGER ML		!Continuum frequency index
	INTEGER I,J,K		!Miscellaneous loop variables.
	INTEGER LU_ER
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) DELF
	REAL(KIND=LDP) MIN_FREQ
	REAL(KIND=LDP) SWITCH_FREQ
	REAL(KIND=LDP) TEMP_FREQ
	REAL(KIND=LDP) T1
!
	LOGICAL EDGE_FREQ(NCF)
!
! External functions
!
	INTEGER ERROR_LU
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
	LU_ER=ERROR_LU()
!
! Check validity of some of the parameters.
!
	IF(R_CMF_WING_EXT .LT. 2.0_LDP .OR. R_CMF_WING_EXT .GT. 20.0_LDP)THEN
	  WRITE(LU_ER,*)'Invalid value for R_CMF_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'R_CMF_WING_EXT=',R_CMF_WING_EXT
	  STOP
	END IF
	IF(ES_WING_EXT .LT. 0.0_LDP .OR. ES_WING_EXT .GT. 20000.0_LDP)THEN
	  WRITE(LU_ER,*)'Invalid value for ES_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'ES_WING_EXT=',ES_WING_EXT
	  STOP
	END IF
	IF(MAX_DOP .LT. 2.0_LDP .OR. MAX_DOP .GT. 20.0_LDP)THEN
	  WRITE(LU_ER,*)'Invalid value for MAX_DOP in INS_LINE'
	  WRITE(LU_ER,*)'MAX_DOP=',MAX_DOP
	  STOP
	END IF
	IF(dV_CMF_PROF .LT. 1.0_LDP .OR. dV_CMF_PROF .GT. 1.0E+05_LDP)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_PROF in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_PROF=',dV_CMF_PROF
	  STOP
     	END IF
	IF(dV_CMF_WING .LT. 1.0_LDP .OR. dV_CMF_WING .GT. 1.0E+05_LDP)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_WING in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_WING=',dV_CMF_WING
	  STOP
	END IF
	IF(NU_CONT(NCF) .GT. NU_LINE(N_LINES))THEN
	  LN_INDX=N_LINES-1
	  DO WHILE(NU_CONT(NCF) .GT. NU_LINE(LN_INDX))
	    LN_INDX=LN_INDX-1
	  END DO
	  WRITE(LU_ER,*)'Warning from INS_LINE_V4'
	  WRITE(LU_ER,'(1X,I5,A,A)')N_LINES-LN_INDX,' weak lines in ',
	1        'extreme IR will be ignored as outside continuum range.'
	  WRITE(LU_ER,*)'Min(Nu_CONT)=',NU_CONT(NCF)
	  WRITE(LU_ER,*)'Min(Nu_LINE)=',NU_LINE(N_LINES)
	END IF
!
! 
!
! We assume that both lines and continuum are ordered from highest to
! lowest frequencies.
!
	dNU_on_NU=FRAC_DOP*V_DOP/C_KMS
	NDOP=NINT(MAX_DOP/FRAC_DOP)
	RES_EXTENT=(1.0_LDP+dNU_on_NU)**NDOP
!
! To avoid numerical instabilities in the iteration procedure when solving
! for the corrections we ensure that the frequencies bracketing a bound-free
! edge are EDGE_SEP_FAC*FRAC_DOP Doppler widths apart. We adjust the lower
! frequency to ensure this. MIN_FREQ_RAT is the minimum ratio allowed between
! successive frequencies.
!
	EDGE_SEP_FAC=0.1_LDP
	MIN_FREQ_RAT=1.0_LDP+EDGE_SEP_FAC*dNU_on_NU
!
! Define approximate edges of the resonance zone and the e.s blue wing.
! These limit getting frequencies unnecessarily close.
!
	BLUE_WING_EXT=ES_WING_EXT+(NDOP+2)*V_DOP*FRAC_DOP	!In km/s
	RED_WING_EXT=ES_WING_EXT+R_CMF_WING_EXT*VINF
!
	APP_RES_EXT=RES_EXTENT*(1.0_LDP+0.5_LDP*dNU_on_NU)	!No units
	APP_ESBW_EXT=1.0_LDP+(BLUE_WING_EXT+0.2_LDP*dV_CMF_WING)/C_KMS
!
! Determine continuum frequencies bracketing bound-free edges. We keep
! these in our final continuum list. To avoid numerical instabilities
! in the iteration procedure when solving for the corrections we ensure
! that the frequencies bracketing a bound-free edge are EDGE_SEP_FAC*FRAC_DOP
! Doppler widths apart. We adjust the lower frequency to ensure this.
!
	EDGE_FREQ(1:NCF)=.FALSE.
	I=2
	DO WHILE (I .LT. NCF)
	  IF( ABS(NU_CONT(I-1)/NU_CONT(I)-1.0_LDP) .LT. 1.0E-08_LDP)THEN
	    IF(NU_CONT(I)/MIN_FREQ_RAT .GT. MIN_FREQ_RAT*NU_CONT(I+1))THEN
	      EDGE_FREQ(I-1)=.TRUE.
	      EDGE_FREQ(I)=.TRUE.
	      NU_CONT(I)=NU_CONT(I)/MIN_FREQ_RAT
	      I=I+2
	    ELSE
	      EDGE_FREQ(I-1)=.TRUE.
	      K=1
	      DO WHILE ( NU_CONT(I)/MIN_FREQ_RAT .LT. NU_CONT(I+K))
	        K=K+1
	      END DO
	      EDGE_FREQ(I+K)=.TRUE.
	      I=I+K+1
	    END IF
	  ELSE
	    I=I+1
	  END IF
	END DO
!
! Find the first line that is to be included as a blanketed line.
!
	LN_INDX=1
	DO WHILE(LN_INDX .LE. N_LINES .AND.
	1          TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	  LN_INDX=LN_INDX+1
	END DO
!
! If no lines to be included, set FREQ directly to continuum frequencies
! and return.
!
	IF(LN_INDX .GT. N_LINES)THEN
	  FREQ(1:NCF)=NU_CONT(1:NCF)
	  NFREQ=NCF
	  LU_ER=ERROR_LU()
	  WRITE(LU_ER,*)'Warning --- no line inserted in INS_LINE_V4'
	  RETURN
	END IF
! 
!
! Compute the frequency grid. We use a combination of CONTINUUM and LINE
! frequencies. Continuum frequencies are inserted simultaneously with the
! line frequencies (rather than after the line grid is created) to avoid
! the need for extra temporary storage for LINE_THIS_FREQ, and to avoid
! having to alter LINE_ST_INDX etc.
!
! All edge frequencies are included.
!
	ML=1
	INDX=1
	FREQ(1)=NU_CONT(1)
	LINES_THIS_FREQ(1)=0
!
	DO WHILE (ML .LE. NCF)
          IF( NU_CONT(ML) .GE. FREQ(INDX) )THEN
	     ML=ML+1			!Ignore as past this freq.
!
! If continuum frequency is within 0.2 Doppler widths of last set frequency,
! there is no need to use it, unless it is a bound-free edge frequency.
!
          ELSE IF( NU_CONT(ML) .GE. FREQ(INDX)/(1.0_LDP+0.2_LDP*dNU_on_NU)
	1                       .AND. .NOT. EDGE_FREQ(ML) )THEN
	     ML=ML+1			!Use current set frequency.
!
	  ELSE IF( LN_INDX .GT. N_LINES .OR. NU_CONT(ML) .GT.
	1          NU_LINE(LN_INDX)*MAX(APP_RES_EXT,APP_ESBW_EXT) )THEN
!
! No nearby line, so continuum point becomes next frequency point.
!
	    INDX=INDX+1
	    IF(INDX .GT. NFREQ_MAX)GOTO 9999
	    FREQ(INDX)=NU_CONT(ML)
	    LINES_THIS_FREQ(INDX)=0
	    ML=ML+1
	  ELSE
!
! Add in additional frequencies for this line.
!
! First, we allow for electron scattering on the blue side of the line profile.
! Strictly only need if scattering is non-coherent. Insertion of frequencies
! only occurs until we hit the resonance zone.
!
! Insert points at beginning (high frequency side) of electron scattering
! wing if needed.
!
	    IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*APP_ESBW_EXT)THEN
	      T1=NU_LINE(LN_INDX)*(1+BLUE_WING_EXT/C_KMS)
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND.
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	        ML=ML+1
	      END DO
              IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=0
	      END IF
	    END IF
!
! Now insert points in the blue e.s. wing, until we reach the resonance
! zone.
!
!
	    T1=FREQ(INDX)/(1.0_LDP+1.1_LDP*dV_CMF_WING/C_KMS)
	    IF(T1 .GT. NU_LINE(LN_INDX)*RES_EXTENT)THEN
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND.
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	        ML=ML+1
	      END DO
              IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=0
	      END IF
!
! Only add in extra points if resonance zone is 1.2dV_CMF_WING away.
!
	      I=INT( (FREQ(INDX)-NU_LINE(LN_INDX)*RES_EXTENT)/
	1               FREQ(INDX)*C_KMS/dV_CMF_WING - 0.1_LDP)
	      DELF=(FREQ(INDX)-NU_LINE(LN_INDX)*RES_EXTENT)/(I+1)
	      DO J=1,I
	        T1=FREQ(INDX)-DELF
	        DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND.
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=NU_CONT(ML)
	            LINES_THIS_FREQ(INDX)=0
	          END IF
	          ML=ML+1
	        END DO
	        IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=T1
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	      END DO
	    END IF
!
! Now begin the resonance zone.
!
	    LINE_ST_INDX(LN_INDX)=INDX+1
	    DO I=-NDOP,NDOP
	      T1=NU_LINE(LN_INDX)*(1.0_LDP+dNU_on_NU)**(-I)
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND.
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=1
	        END IF
	        ML=ML+1
	      END DO
	      IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=1
	      END IF
	    END DO
	    LINE_END_INDX(LN_INDX)=INDX
!
! Ready for next line or continuum point.
!
	    LST_LN_INDX=LN_INDX
	    LN_INDX=LN_INDX+1
	    DO WHILE(LN_INDX .LE. N_LINES .AND.
	1                           TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	      LN_INDX=LN_INDX+1
	    END DO
!
! Check if overlapping line [Overlap of intrinsic profile only!]
!
100	    CONTINUE
	    IF(LN_INDX .LE. N_LINES)THEN
	      IF(FREQ(INDX) .LT. NU_LINE(LN_INDX)*RES_EXTENT)THEN
!
! Check how far line extends back
!
	        J=INDX
	        DO WHILE(FREQ(J) .LT. NU_LINE(LN_INDX)*RES_EXTENT)
	          LINES_THIS_FREQ(J)=LINES_THIS_FREQ(J)+1
	          J=J-1
	        END DO
	        LINE_ST_INDX(LN_INDX)=J+1
!
	        DO WHILE(FREQ(INDX) .GT. NU_LINE(LN_INDX)/RES_EXTENT)
	          T1=FREQ(INDX)/(1.0_LDP+dNU_on_NU)
	          DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	          IF(EDGE_FREQ(ML) .AND.
	1                NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	              INDX=INDX+1
	              IF(INDX .GT. NFREQ_MAX)GOTO 9999
	              FREQ(INDX)=NU_CONT(ML)
	              LINES_THIS_FREQ(INDX)=1
	            END IF
	            ML=ML+1
	          END DO
	          IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=T1
	            LINES_THIS_FREQ(INDX)=1
	          END IF
	        END DO
!
	        LINE_END_INDX(LN_INDX)=INDX
	        LST_LN_INDX=LN_INDX
	        LN_INDX=LN_INDX+1
	        DO WHILE(LN_INDX .LE. N_LINES .AND.
	1                           TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	          LN_INDX=LN_INDX+1
	        END DO
!
! Check another line
!
	        GOTO 100
	      END IF
	    END IF
!
! Now we need to extend line to allow for the outflow velocity. We only extend
! until the blue edge of the next line.
!
! We extend the frequencies with spacing dV_CMF_PROF (in km/s) from the
! resonance zone edge for another 2Vinf/kms. Beyond this freq (called
! SWITCH_FREQ) we use until dV_VMF_WING, which is the spacing in the e.s. wings.
!
	    MIN_FREQ=NU_LINE(LST_LN_INDX)/RES_EXTENT/
	1              (1.0_LDP+RED_WING_EXT/C_KMS)
!
! The V_CMF_PROF is added to allow for some bleeding.
!
	    SWITCH_FREQ=NU_LINE(LST_LN_INDX)/RES_EXTENT/
	1              (1.0_LDP+(2.0_LDP*VINF+dV_CMF_PROF)/C_KMS)
!
! We check that the minimum frequency does not extend beyond the
! resonance zone of the next line. As we will put a frequency at
! the beginning of the resonance zone, we back of a little bit.
!
	    T1=RES_EXTENT*(1.0_LDP+0.3_LDP*MIN(dV_CMF_PROF,dV_CMF_WING)/C_KMS)
	    IF(LN_INDX .GT. N_LINES)THEN
	    ELSE IF(MIN_FREQ .LT. NU_LINE(LN_INDX)*T1)THEN
	      MIN_FREQ=NU_LINE(LN_INDX)*T1
	    END IF
	    DO WHILE(FREQ(INDX) .GT. MIN_FREQ)
	      IF( FREQ(INDX) .GT. SWITCH_FREQ)THEN
	        TEMP_FREQ=FREQ(INDX)/(1.0_LDP+dV_CMF_PROF/C_KMS)
	      ELSE
	        TEMP_FREQ=FREQ(INDX)/(1.0_LDP+dV_CMF_WING/C_KMS)
	      END IF
!
! We need to check again, since we didn't know the frequency step size.
!
	      IF(TEMP_FREQ .GT. MIN_FREQ)THEN
	        DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. TEMP_FREQ)
	        IF(EDGE_FREQ(ML) .AND.
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=NU_CONT(ML)
	            LINES_THIS_FREQ(INDX)=0
	          END IF
	          ML=ML+1
	        END DO
                IF(TEMP_FREQ .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=TEMP_FREQ
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	      ELSE
	        GOTO 200	!Exit from this section
	      END IF
	    END DO
200	    CONTINUE
!
	  END IF
	END DO
!
! Set the number of frequencies.
!
	NFREQ=INDX
!
! Test for monotocity of frequencies, and determine MINIMUM frequency
! spacing in velocity space.
!
	T1=10000.0_LDP
	DO ML=1,NFREQ-1
	  T1=MIN( T1 , C_KMS*(FREQ(ML)-FREQ(ML+1))/FREQ(ML) )
	  IF(FREQ(ML) .LE. FREQ(ML+1))THEN
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_LINE_V4'
	    WRITE(LU_ER,*)' ML=',ML
	    DO I=MAX(1,ML-20),MIN(NCF,ML+20)
	      WRITE(LU_ER,*)I,FREQ(I)
	    END DO
	    STOP
	  END IF
	END DO
	WRITE(LU_ER,'(1X,A,1PE9.2,A)')
	1          'Minimum frequency spacing is:',T1,'km/s'
!
	DO I=2,NFREQ
	  WRITE(63,'(1X,I6,1P,2E12.4)')I,FREQ(I),
	1            3.0E+05_LDP*(FREQ(I-1)-FREQ(I))/FREQ(I-1)
	END DO
!
!	DO I=1,N_LINES
!	  J=LINE_ST_INDX(I)
!	  K=LINE_END_INDX(I)
!          WRITE(83,'(1P,3E16.6,0P,2I7,A)')NU_LINE(I),FREQ(J),FREQ(K),J,K,TRIM(TRANS_TYPE(I))
!        END DO
!
	RETURN
!
9999	CONTINUE
	LU_ER=ERROR_LU()
	WRITE(LU_ER,*)'Error --- insufficient frequencies to store'//
	1               ' both line and continuum frequencies'
	WRITE(LU_ER,*)'ML= ',ML
	WRITE(LU_ER,*)'LN_INDX= ',LN_INDX
	WRITE(LU_ER,*)'INDX= ',INDX
	WRITE(LU_ER,*)'NFREQ_MAX= ',NFREQ_MAX
	STOP
	END
