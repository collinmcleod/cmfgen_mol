!
! Routine to generate frequencies for radiative transfer program given
! a set of bound-free edges, and maximum and minimum frequencies.
! The frequency spacing is contolled by three passable parameters.
!
! Every bound-free edged is defined by two frequencies, except
! where two bound-free edges coincide to 1 part in 10^11. Genrally
! points are inserted so that simpsons rule can be used, although in
! some regions where the bound-free edges are close the trapazoidal
! rule will be used.
!
	SUBROUTINE SET_CONT_FREQ_V2(NEW_FREQ,EDGE,TYPE,INDX,
	1                        SMALL_RAT,BIG_AMP,DNU_MAX,
	1                        MAX_FREQ,MIN_FREQ,
	1                        dV_LEV,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        dV_CONT,dV_DOP,
	1                        N,NCF,NCF_MAX,LUOUT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 26-May-1996 : ERROR_LU installed.
! Cleaned
! Created 29-Mar-1990 (Based on GEN_FREQ).
!
	INTEGER N,NCF,NCF_MAX,LUOUT,INDX(NCF_MAX)
	REAL(KIND=LDP) EDGE(NCF_MAX)
	REAL(KIND=LDP) NEW_FREQ(NCF_MAX)
	CHARACTER*(*) TYPE(NCF)
!
	REAL(KIND=LDP) MAX_FREQ		!Maximum continuum frequency
	REAL(KIND=LDP) MIN_FREQ		!Minimum continuum frequency
!
	REAL(KIND=LDP) DNU_MAX		!Maximum frequency spacing near/above
				!  bound-free edge: i.e. dNU <  DNU_MAX
	REAL(KIND=LDP) BIG_AMP  	!Amplification. dNU increases by a factor
				! BIG_AMP as we move away from the b.f. edge.
				! for frequencies above SWITCH_FREQ.
   	REAL(KIND=LDP) SMALL_RAT	!Used to define frequency spacing for
				! frequencies less than SWITCH_FREQ.
				! dNU/NU=SMALL_RAT-1
!
! Parameters for installing etra frequencies near bound-free edges
! (low frequency side) to allow for level dissolution.
!
	REAL(KIND=LDP) dV_LEV			!Spacing near b-f edge.
	REAL(KIND=LDP) AMP_DIS			!Amplification factor for dNU as we
					!  move to smaller frequencies.
	REAL(KIND=LDP) MIN_FREQ_LEV_DIS		!Indicates that the extra frequencies
					!  should only be installed for
					!  frequencies above MIN_FREQ_LEV_DIS.
	REAL(KIND=LDP) dV_DOP			!Minimum spacing in Doppler line profiles
                                        !  (km/s). Used to set minimum frequecny
                                        !  spacing.
	REAL(KIND=LDP) dV_CONT                  !Maximum spacing in arbitrary section of
                                        !  continuum (km/s)
!
	INTEGER, PARAMETER :: RZERO=0.0_LDP
	INTEGER, PARAMETER :: RONE=1.0_LDP
	INTEGER, PARAMETER :: RTWO=2.0_LDP
!
! N+2 as we insert MAX_FREQ and MIN_FREQ into the array.
!
	CHARACTER*10   CHAR_WRK(N+2)
!
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) SPEED_OF_LIGHT
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) RN
	REAL(KIND=LDP) COINCIDENT_FRAC
	REAL(KIND=LDP) DEL_NU
	REAL(KIND=LDP) DEL_NU_TO_NEXT_EDGE
!
	INTEGER INDX_DIS
	INTEGER I,J,K
	LOGICAL NUMER,EQUAL
	REAL(KIND=LDP) FAC,EQUAL_FAC
	REAL(KIND=LDP) SWITCH_FREQ
!
	LUER=ERROR_LU()
	COINCIDENT_FRAC=0.2_LDP
!
! Sort frequencies into numerical order. NEW_FREQ is used as a
! work array.
!
	DO I=1,N
	   WRITE(168,*)I,EDGE(I),TYPE(I)
	END DO
	N=N+2
	EDGE(N-1)=MIN_FREQ;   TYPE(N-1)=' '
	EDGE(N)=MAX_FREQ;     TYPE(N)=' '
	NUMER=.TRUE.
	CALL INDEXX(N,EDGE,INDX,NUMER)
	CALL SORTDP(N,EDGE,INDX,NEW_FREQ)
	CALL SORTCHAR(N,TYPE,INDX,CHAR_WRK)
!
	DO I=1,N
	  WRITE(167,*)I,EDGE(I),TYPE(I)
	END DO
!
! Do a quick check that sort was done correctly.
! Not equal as MIN_FREQ, and MAX_FREQ have been inserted in FREQ array.
!
	IF(MIN_FREQ .NE. EDGE(1))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - MIN_FREQ too big'
	  WRITE(LUER,*)'MIN_FREQ =',MIN_FREQ
	  WRITE(LUER,*)'EDGE(1) =',EDGE(1)
	  STOP
	END IF
	IF(MAX_FREQ .NE. EDGE(N))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - MAX_FREQ too small'
	  WRITE(LUER,*)'MAX_FREQ =',MAX_FREQ
	  WRITE(LUER,*)'EDGE(N) =',EDGE(N)
	  STOP
	END IF
!
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
	EQUAL_FAC=COINCIDENT_FRAC*dV_DOP/C_KMS
	FAC = RONE + EQUAL_FAC
!
! SMAL_FAC is the ratio used to set the frequency spacing for
! small frequencies (frequencies less than 1.0 approximately).
! Two points are inserted so that simpsons rule can be used.
!
	IF(SMALL_RAT .LE. RONE .OR. SMALL_RAT .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - SMALL_RAT Outside ',
	1           'valid range ',SMALL_RAT
	  STOP
	END IF
!
! DNU_MAX is the maximum frequecy spacing for frequencies adjacent
! to a bound-free edge. For large v, we require constant spacing since
! error in integral is a function of {del v }.
!
	IF(DNU_MAX .LE. RZERO .OR. DNU_MAX .GT. RTWO)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - DNU_MAX Outside ',
	1           'valid range ',DNU_MAX
	  STOP
	END IF
!
! BIG_AMP is used to amplify the spacing for large frequencies. Close
! to the edge, the spacing is DNU_MAX, but for every
! points, the spacing increases by a factor of BIG_AMP
!
	IF(BIG_AMP .LT. RONE .OR. BIG_AMP .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - BIG_AMP ',
	1                 'Outside valid range ',BIG_AMP
	  STOP
	END IF
!
! Determin the freuqncy where we switch from using DNU_MAX to
! SMALL_RAT.
!
	SWITCH_FREQ=DNU_MAX/(SMALL_RAT-RONE)
!
! Now begin inserting extra points into frequency array.
!
	INDX_DIS=0
	K=1
	NEW_FREQ(1)=MIN_FREQ
	I=2
	DO WHILE(I .LE. N)
!
! Get next edge for which we will consider level dissolution.
!
          IF(I .GT. INDX_DIS)THEN
	    INDX_DIS=I
	    DO WHILE( INDEX(TYPE(INDX_DIS),'D') .EQ. 0 .AND. INDX_DIS .LT. N)
	      INDX_DIS=INDX_DIS+1
 	    END DO
	  END IF
!
! Determine whether next EDGE frequency is too close to selected frequency
! already.
!
	  DEL_NU=EDGE(I)/FAC-NEW_FREQ(K)
	  DEL_NU_TO_NEXT_EDGE=DEL_NU
	  IF(DEL_NU/NEW_FREQ(K) .LT. EQUAL_FAC)THEN
	     NEW_FREQ(K)=EDGE(I)*FAC
	     I=I+1
	     GOTO 1000
	  END IF
!
! dV_CONT is the minimum spacing (velocity space) at which we have to sample the
! photoioinization cross-sections.
!
	  T1=NEW_FREQ(K)*dV_CONT/C_KMS
	  IF(T1 .LT. DEL_NU)THEN
	     DEL_NU=T1
	  END IF
!
! For frequency close to a bound-free edge (but on the high side) we
! use a spacing DNU_MAX near the edge. As we move blueward,
! the spacing increaes by BIG_AMP on each succesive insertion.
!
	  IF(DEL_NU .GT. DNU_MAX .AND. NEW_FREQ(K) .GT. SWITCH_FREQ)THEN
	    T1=(NEW_FREQ(K)-EDGE(I-1))/DNU_MAX
	    RN=LOG( 1.0_LDP+T1*(BIG_AMP-1.0_LDP))/LOG(BIG_AMP) - 1.0_LDP
	    T1=DNU_MAX*(BIG_AMP**RN)
	    IF(T1 .LT. DEL_NU)DEL_NU=T1
	  END IF
!
! Small frequencies where we use a fixed velocity spacing. dV_CONT
! has essentially the same meaning, and may supercede SMALL_RAT.
!
	  T1=NEW_FREQ(K)*(SMALL_RAT-1.0_LDP)
	  IF(DEL_NU .GT. T1 .AND. NEW_FREQ(K) .LE. SWITCH_FREQ)THEN
	    DEL_NU=T1
	  END IF
!
! Now check for level dissolution. We only apply level dissolution to
! those levels which have a 'D' specified in the TYPE variable.
!
	  IF(EDGE(INDX_DIS) .GT. MIN_FREQ_LEV_DIS .AND. INDX_DIS .NE. N)THEN
	    T2=EDGE(INDX_DIS)*dV_LEV/C_KMS
	    T1=(EDGE(INDX_DIS)-NEW_FREQ(K))/T2
	    RN=LOG( 1.0_LDP+T1*(AMP_DIS-1.0_LDP))/LOG(AMP_DIS) - 1.0_LDP
	    T1=T2*(AMP_DIS**RN)
	    IF(T1 .LT. DEL_NU)DEL_NU=T1
	  END IF
!
! Insert frequency edge.
!
	  IF(DEL_NU .EQ. DEL_NU_TO_NEXT_EDGE)THEN
	    IF(I .EQ .N)THEN
	      K=K+1
	      NEW_FREQ(K)=EDGE(N)
	    ELSE
	      K=K+1
	      NEW_FREQ(K)=EDGE(I)/FAC
	      K=K+1
	      NEW_FREQ(K)=EDGE(I)*FAC
	    END IF
	    I=I+1
	  ELSE
	    K=K+1
	    NEW_FREQ(K)=NEW_FREQ(K-1)+DEL_NU
	  END IF
!
	  IF(K .GT. NCF_MAX)THEN
	    WRITE(LUER,*)'Error NCF too small in SET_CONT_FREQ'
	    WRITE(LUER,*)N,NCF_MAX,SMALL_RAT,BIG_AMP,DNU_MAX
	    OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	      DO J=1,NCF
	        WRITE(LUOUT,100)NEW_FREQ(J)
	      END DO
	    CLOSE(LUOUT)
	    STOP
	  END IF
1000	  CONTINUE
!
	END DO				!I loop
!
	NCF=K
!
! Now sort frequencies into numerical decreasing order.
!
	DO I=1,NCF/2
	  T1=NEW_FREQ(I)
	  NEW_FREQ(I)=NEW_FREQ(NCF-I+1)
	  NEW_FREQ(NCF-I+1)=T1
	END DO
!
	OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	  DO I=1,NCF
	    WRITE(LUOUT,100)NEW_FREQ(I)
	  END DO
100	  FORMAT(1X,F22.16)
	CLOSE(LUOUT)
!
! Check frequency array is monotonic.
!
	DO I=1,NCF-1
	  IF( NEW_FREQ(I) .LE. NEW_FREQ(I+1) )THEN
	    WRITE(LUER,*)'Error in SET_CONT_FREQ - frequency array is',
	1             ' not monotonic.'
	    STOP
	  END IF
	END DO
!
! Ensure FREQ array is zeroed (as probably will use OBSF).
!
	DO I=1,NCF
	  EDGE(I)=0.0_LDP
	END DO
	WRITE(6,*)'Returning'
!
	RETURN
	END
