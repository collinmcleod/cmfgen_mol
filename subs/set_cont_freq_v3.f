C
C Routine to generate frequencies for radiative transfer program given
C a set of bound-free edges, and maximum and minimum frequencies.
C The frequency spacing is contolled by three passable parameters.
C
C Every bound-free edged is defined by two frequencies, except
C where two bound-free edges coincide to 1 part in 10^11. Genrally
C points are inserted so that simpsons rule can be used, although in
C some regions where the bound-free edges are close the trapazoidal
C rule will be used.
C
	SUBROUTINE SET_CONT_FREQ_V3(NEW_FREQ,FREQ,INDX,
	1                        SMALL_RAT,BIG_AMP,DNU_MAX,MAX_FREQ,MIN_FREQ,
	1                        dV_LEV,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_END_XRAY,N,NCF,NCF_MAX,LUOUT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 26-May-1996 : ERROR_LU installed.
C Cleaned
C Created 29-Mar-1990 (Based on GEN_FREQ).
C
	INTEGER N,NCF,NCF_MAX,LUOUT,INDX(NCF_MAX)
	REAL(KIND=LDP) FREQ(NCF_MAX),NEW_FREQ(NCF_MAX)
C
	REAL(KIND=LDP) MAX_FREQ		!Maximum continuum frequency
	REAL(KIND=LDP) MIN_FREQ		!Minimum continuum frequency
	REAL(KIND=LDP) DNU_MAX		! Twice the maximum frequency spacing near/above
				!  bound-free edge: i.e. dNU < 0.5* DNU_MAX
	REAL(KIND=LDP) BIG_AMP  	!Amplification. dNU increases by a factor
				! BIG_AMP as we move away from the b.f. edge.
				! for frequencies above SWITCH_FREQ.
   	REAL(KIND=LDP) SMALL_RAT	!Used to define frequency spacing for
				! frequencies less than SWITCH_FREQ.
				! dNU/NU=SMALL_RAT-1
C
C Parameters for installing etra frequencies near bound-free edges
C (low frequency side) to allow for level dissolution.
C
	REAL(KIND=LDP) dV_LEV			!Spacing near b-f edge.
	REAL(KIND=LDP) AMP_DIS			!Amplification factor for dNU as we
					!  move to smaller frequencies.
	REAL(KIND=LDP) MIN_FREQ_LEV_DIS		!Indicates that the extra frequencies
					! should only be installed for
					! frequencies above MIN_FREQ_LEV_DIS.
!
	REAL(KIND=LDP) DELV_CONT
	REAL(KIND=LDP) DELV_XRAY
	REAL(KIND=LDP) NU_END_XRAY
C
        REAL(KIND=LDP), PARAMETER :: RZERO=0.0_LDP
        REAL(KIND=LDP), PARAMETER :: RHALF=0.5_LDP
        REAL(KIND=LDP), PARAMETER :: RONE=1.0_LDP
        REAL(KIND=LDP), PARAMETER :: RTWO=2.0_LDP
C
        REAL(KIND=LDP) SPEED_OF_LIGHT
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
C
	REAL(KIND=LDP) T1,dV_NEW,dV
	REAL(KIND=LDP) C_KMS
C
	INTEGER I,J,K,L,ML,NPTS,K_BEG
	LOGICAL NUMER,EQUAL
	REAL(KIND=LDP) FAC,EQUAL_FAC
	REAL(KIND=LDP) UP,LOW,DELF,DIFF,RAT,SWITCH_FREQ
C
	LUER=ERROR_LU()
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
C
C Sort frequencies into numerical order. New freq is used as a
C work array.
C
	N=N+2
	FREQ(N-1)=MIN_FREQ
	FREQ(N)=MAX_FREQ
	NUMER=.TRUE.
	CALL INDEXX(N,FREQ,INDX,NUMER)
	CALL SORTDP(N,FREQ,INDX,NEW_FREQ)
C
C Not equal as MIN_FREQ, and MAX_FREQ have been inserted in FREQ array.
C
	IF(MIN_FREQ .NE. FREQ(1))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - MIN_FREQ too big'
	  WRITE(LUER,*)'MIN_FREQ =',MIN_FREQ
	  WRITE(LUER,*)'FREQ(1) =',FREQ(1)
	  STOP
	END IF
	IF(MAX_FREQ .NE. FREQ(N))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - MAX_FREQ too small'
	  WRITE(LUER,*)'MAX_FREQ =',MAX_FREQ
	  WRITE(LUER,*)'FREQ(N) =',FREQ(N)
	  STOP
	END IF
C
	EQUAL_FAC=1.0E-11_LDP
	FAC = RONE + EQUAL_FAC
C
C SMAL_FAC is the ratio used to set the frequency spacing for
C small frequencies (frequencies less than 1.0 approximately).
C Two points are inserted so that simpsons rule can be used.
C
	IF(SMALL_RAT .LE. RONE .OR. SMALL_RAT .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - SMALL_RAT Outside ',
	1           'valid range ',SMALL_RAT
	  STOP
	END IF
C
C DNU_MAX is the maximum frequecy spacing for frequencies adjacent
C to a bound-free edge. For large v, we require constant spacing since
C error in integral is a function of {del v }.
C
	IF(DNU_MAX .LE. RZERO .OR. DNU_MAX .GT. RTWO)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - DNU_MAX Outside ',
	1           'valid range ',DNU_MAX
	  STOP
	END IF
C
C G is used to amlify the spacing for large frequencies. Close
C to the edge, the spacing is less than DNU_MAX, but for every
C two points, the spacing increases by a factor of G.
C
	IF(BIG_AMP .LT. RONE .OR. BIG_AMP .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - BIG_AMP ',
	1                 'Outside valid range ',BIG_AMP
	  STOP
	END IF
C
	SWITCH_FREQ=DNU_MAX/(SMALL_RAT-RONE)
C
C Now begin inserting extra points into frequency array.
C
	K=1
	NEW_FREQ(1)=MIN_FREQ
	DO I=1,N-1
	  IF( EQUAL(FREQ(I),FREQ(I+1),EQUAL_FAC) )THEN
	    FREQ(I+1)=FREQ(I)
	  ELSE
	    LOW=FREQ(I)
	    UP=FREQ(I+1)
	    IF(I .NE. 1)THEN
	      LOW=LOW*FAC
	      K=K+1
	      NEW_FREQ(K)=LOW
	    END IF
	    IF(I .NE. N-1)UP=UP/FAC
	    RAT=UP/LOW
	    DIFF=UP-LOW
C
C The insertion of frequencies is treated differently depending on whether the
C EDGE_FREQ (stored in FREQ) is greater than, or less than, SWITCH_FREQ.
C
C If EDGE_FREQ < SWITCH_FREQ  we insert fequencies so that FREQ(K+1)/FREQ(K)
C < SMALL_RAT.
C
C If EDGE_FREQ > SWITCH_FREQ  we insert fequencies so that FREQ(K+1)-FREQ(K)
C > DNU_MAX at the  bound-free edge. The spacing is ramped up by a factor
C G as we move to higher frequencies away from the bound-free edge.
C
C Note the spacing is actually a factor of 2 better than implied by DNU_MAX
C and SMALL_RAT becuase we insert points midway so that Simpson's rule can
C be used.
C
	    IF(LOW .LT. SWITCH_FREQ)THEN
	      IF(RAT .GT. SMALL_RAT)THEN
	        NPTS=LOG(RAT)/LOG(SMALL_RAT) + 1
	        DELF=(RAT)**(RONE/NPTS)
	        DO J=1,NPTS-1
	          K=K+1
	          NEW_FREQ(K)=NEW_FREQ(K-1)*DELF
	        END DO
	        K=K+1
	        NEW_FREQ(K)=UP
	      ELSE
	        K=K+1
	        NEW_FREQ(K)=UP
	      END IF
	    ELSE IF(DIFF .GT. DNU_MAX)THEN
	      NPTS=LOG( RONE - (UP-LOW)*(RONE-BIG_AMP)/DNU_MAX )/
	1                 LOG(BIG_AMP) + 1
	      DELF=(UP-LOW)*(RONE-BIG_AMP)/(RONE-(BIG_AMP**NPTS))
	      DO J=1,NPTS-1
	        K=K+1
	        NEW_FREQ(K)=NEW_FREQ(K-1)+DELF
	        DELF=DELF*BIG_AMP
	      END DO
	      K=K+1
	      NEW_FREQ(K)=UP
	    ELSE
	      K=K+1
	      NEW_FREQ(K)=UP
	    END IF			!< switch_freq if
	  END IF			!Equal if
C
C Now we will insert a few extra points before the bound-free edge.
C This will allow for blending of levels due to level dissolution near the
C bound-free edge.
C
C This is only done for frequencies above MIN_FREQ_LEV_DIS.
C
	  IF(NEW_FREQ(K) .GT. MIN_FREQ_LEV_DIS)THEN
	    K_BEG=K-1
500	    CONTINUE
	    dV=C_KMS*(NEW_FREQ(K)-NEW_FREQ(K_BEG))/NEW_FREQ(K)
	    IF(dV .GT. dV_LEV .AND. I .LT . N-1)THEN
	      T1= LOG( dV/dV_LEV*(AMP_DIS-RONE) + RONE) / LOG(AMP_DIS)
	      J=NINT(T1)-1
	      IF(K_BEG .NE. 1)THEN
	        IF( dV_LEV*AMP_DIS**(J+1) .LT. C_KMS*
	1               (NEW_FREQ(K_BEG)-NEW_FREQ(K_BEG-1))/NEW_FREQ(K))THEN
	          K_BEG=K_BEG-1
	          GOTO 500
	        END IF
	      END IF
	      dV_NEW=dV*(AMP_DIS-RONE)/(AMP_DIS**(J+1) -RONE)
	      UP=NEW_FREQ(K)
	      NEW_FREQ(K_BEG+J+1)=UP
	      T1=0.0_LDP
	      DO L=J,1,-1
	        T1=T1+dV_NEW*(AMP_DIS**(J-L))/C_KMS
	        NEW_FREQ(K_BEG+L)=UP*(RONE-T1)
	      END DO
	      K=K_BEG+J+1
	    END IF	
	  END IF
C
	  IF(K .GT. NCF_MAX)THEN
	    WRITE(LUER,*)'Error NCF too small in SET_CONT_FREQ_V3'
	    WRITE(LUER,*)N,NCF_MAX,SMALL_RAT,BIG_AMP,DNU_MAX
	    OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	      DO J=1,NCF
	        WRITE(LUOUT,100)NEW_FREQ(J)
	      END DO
	    CLOSE(LUOUT)
	    STOP
	  END IF
C
	END DO				!I loop
C
	NCF=K
C
C Now insert points midway to allow use of Simpson rule if desired. If the
C points are closer than 30km/s (i.e. T1*c), no points are inserted.
C
	DO I=1,NCF
	 FREQ(I)=NEW_FREQ(I)	
	END DO
	K=1
	T1=1.0E-04_LDP
	DO I=2,NCF
	  IF( EQUAL(FREQ(I),FREQ(I-1), T1) )THEN
	    K=K+1
	    NEW_FREQ(K)=FREQ(I)
	  ELSE
	    K=K+1
	    NEW_FREQ(K)=RHALF*(FREQ(I)+FREQ(I-1))
	    K=K+1
	    NEW_FREQ(K)=FREQ(I)
	  END IF
	END DO
	NCF=K
C
C Now sort frequencies into numerical decreasing order. Then check
C frequency array is monotonic.
C
	DO I=1,NCF/2
	  DELF=NEW_FREQ(I)
	  NEW_FREQ(I)=NEW_FREQ(NCF-I+1)
	  NEW_FREQ(NCF-I+1)=DELF
	END DO
!
! Add in extra frequencies to ensure adequate sampling of continuum.
! NU_EVAL is used as a temporary aray.
!
	FREQ(1:NCF)=NEW_FREQ(1:NCF)	
        K=1
        DO ML=2,NCF
          T1=C_KMS*(FREQ(ML-1)-FREQ(ML))/FREQ(ML)
          IF(T1 .GT. 1.25_LDP*DELV_XRAY .AND. FREQ(ML) .GT. NU_END_XRAY)THEN
            J=T1/DELV_XRAY/1.2_LDP
            DO L=1,J
              K=K+1
              IF(K .GT. NCF_MAX)EXIT
              NEW_FREQ(K)=FREQ(ML-1)-L*(FREQ(ML-1)-FREQ(ML))/(J+1.0_LDP)
            END DO
          ELSE IF(T1 .GT. 1.25_LDP*DELV_CONT .AND. FREQ(ML) .LE. NU_END_XRAY)THEN
            J=T1/DELV_CONT/1.2_LDP
            DO L=1,J
              K=K+1
              IF(K .GT. NCF_MAX)EXIT
              NEW_FREQ(K)=FREQ(ML-1)-L*(FREQ(ML-1)-FREQ(ML))/(J+1.0_LDP)
            END DO
          END IF
          K=K+1
          IF(K .GT. NCF_MAX)THEN
            WRITE(LUER,*)'Error NCF_MAX is too small in SET_CONT_FREQ_V3'
            WRITE(LUER,*)'NCF_MAX=',NCF_MAX
            STOP
          END IF
          NEW_FREQ(K)=FREQ(ML)
        END DO
        NCF=K
C
	OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	  DO I=1,NCF
	    WRITE(LUOUT,100)NEW_FREQ(I)
	  END DO
100	  FORMAT(1X,F22.16)
	CLOSE(LUOUT)
C
	DO I=1,NCF-1
	  IF( NEW_FREQ(I) .LE. NEW_FREQ(I+1) )THEN
	    WRITE(LUER,*)'Error in SET_CONT_FREQ_V3 - frequency array not',
	1             ' monotonic'
	    STOP
	  END IF
	END DO
C
C Ensure FREQ array is zeroed (as probably will use OBSF).
C
	DO I=1,NCF
	  FREQ(I)=0.0_LDP
	END DO
C
	RETURN
	END
