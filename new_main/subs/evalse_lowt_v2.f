!
! Subroutine to increment the statistical equilibrium equations for each
! depth point for recombinations. Routine designed only for the case
! when HNST is very high, and has been set to zero. That is, it is designed
! to handle species such as OVII when T is very low.
!
	SUBROUTINE EVALSE_LOWT_V2(RJ,NU,FQW,COMPUTE_BA,NT,ND,FIRST_FREQ)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 04-Apr-2011 - Changed check on whether minval is zero (added XzVLTE_MIN).
! Created 26-Nov-2010
!
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) FQW
!
	LOGICAL COMPUTE_BA
	LOGICAL FIRST_FREQ
!	
	INTEGER ID              !Ionic species identifier.
	INTEGER ION_LEV	        !Levl ID of target in DI (i.e. the ION).
	INTEGER ION_EQ
	INTEGER ION_V
	INTEGER NIV
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	INTEGER I
	INTEGER J
	INTEGER K
	INTEGER IPR
!
! REV_HNST referes to the LTE population  of the level defined with respect
! to the actual destination (target) level.
!
	REAL(KIND=LDP) REV_HNST
	REAL(KIND=LDP) LOG_B_RAT
	REAL(KIND=LDP) DI_FAC
	REAL(KIND=LDP) ED_FAC
	REAL(KIND=LDP) T_FAC
	REAL(KIND=LDP) JREC
	REAL(KIND=LDP) dRR			!Increment to radiative recombinaton rate.
!
	REAL(KIND=LDP) SUM_SE
	REAL(KIND=LDP) SUM_VK_R
	REAL(KIND=LDP) T1,T2,T3,T4
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: XzVLTE_MIN(:)
!
	IF(FIRST_FREQ .AND. .NOT. ALLOCATED(XzVLTE_MIN))ALLOCATE(XzVLTE_MIN(NUM_IONS))
	IF(FIRST_FREQ)THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	       XzVLTE_MIN(ID)=MINVAL(ATM(ID)%XzVLTE)
	    ELSE
	       XzVLTE_MIN(ID)=0.0_LDP
	    END IF
	  END DO
	END IF
!
! REV_HNST= HNST * B(ION_LEV)/B(1) where b is the deparure coefficient.
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ID,IPR,K,J,
!$OMP1          NIV,ION_LEV,ION_EQ,ION_V,T1,T2,T3,T4,
!$OMP1          JREC,SUM_SE,LOG_B_RAT,dRR,ED_FAC,DI_FAC,T_FAC)
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES .AND. XzVLTE_MIN(ID) .EQ. 0.0_LDP)THEN
	    NIV=SE(ID)%N_IV
            DO IPR=1,ATM(ID)%N_XzV_PHOT
	      ION_LEV=ATM(ID)%XzV_ION_LEV_ID(IPR)
	      ION_EQ=SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV)
	      ION_V=ION_EQ
	      IF(ION_LEV .EQ. 0)EXIT
!
	      T1=TWOHCSQ*(NU**3)
	      DO K=1,ND
	        JREC=(T1+RJ(K))*FQW/NU
	        SUM_SE=0.0_LDP
	        LOG_B_RAT=(ATM(ID+1)%LOG_XzVLTE(1,K)-ATM(ID+1)%LOG_XzVLTE(ION_LEV,K)) +
	1                    LOG(ATM(ID+1)%XzV(ION_LEV,K)/ATM(ID+1)%XzV(1,K))
	        DO J=1,ATM(ID)%NXzV
	          IF(ATM(ID)%XzVLTE(J,K) .EQ. 0.0_LDP .AND. ATM(ID)%WSXzV(J,K,IPR) .GT. 0.0_LDP)THEN
	            dRR=JREC*ATM(ID)%WSXzV(J,K,IPR)*EXP( ATM(ID)%LOG_XzVLTE(J,K)+LOG_B_RAT-HDKT*NU/T(K) )
	            SE(ID)%STEQ(J,K)=SE(ID)%STEQ(J,K)+dRR
	            SUM_SE=SUM_SE+dRR
	          END IF
	        END DO
	        SE(ID)%STEQ(ION_EQ,K)=SE(ID)%STEQ(ION_EQ,K)-SUM_SE
	      END DO
!
	      IF(COMPUTE_BA)THEN
	        T1=TWOHCSQ*(NU**3)
	        DO K=1,ND
	          JREC=(T1+RJ(K))*FQW/NU
	          LOG_B_RAT=(ATM(ID+1)%LOG_XzVLTE(1,K)-ATM(ID+1)%LOG_XzVLTE(ION_LEV,K)) +
	1                       LOG(ATM(ID+1)%XzV(ION_LEV,K)/ATM(ID+1)%XzV(1,K))
	          DO J=1,ATM(ID)%NXzV
	            IF(ATM(ID)%XzVLTE(J,K) .EQ. 0.0_LDP .AND. ATM(ID)%WSXzV(J,K,IPR) .GT. 0.0_LDP)THEN
	              T3=JREC*EXP( ATM(ID)%LOG_XzVLTE(J,K)+LOG_B_RAT-HDKT*NU/T(K) )
	              T4=T3*ATM(ID)%dWSXzVdT(J,K,IPR)
	              T3=T3*ATM(ID)%WSXzV(J,K,IPR)
!
	              DI_FAC=T3/ATM(ID+1)%XzV(ION_LEV,K)
	              ED_FAC=T3/ED(K)
	              T_FAC=T4 + T3*( ATM(ID)%dlnXzVLTE_dlnT(J,K) +
	1               (ATM(ID+1)%dlnXzVLTE_dlnT(1,K)-ATM(ID+1)%dlnXzVLTE_dlnt(ION_LEV,K)) + HDKT*NU/T(K) )/T(K)
!
	              SE(ID)%BA_PAR(J,ION_V,K) =SE(ID)%BA_PAR(J,ION_V,K)   + DI_FAC
	              SE(ID)%BA_PAR(J,NIV-1,K) =SE(ID)%BA_PAR(J,NIV-1,K)   + ED_FAC
	              SE(ID)%BA_PAR(J,NIV,K)   =SE(ID)%BA_PAR(J,NIV,K)     + T_FAC
!
! Include ionizations/recombinations implicitly in the rate equation
! of the target ion (eg He++(gs) for He+ ion/recoms ).
!
	              SE(ID)%BA_PAR(ION_EQ,ION_V,K) =SE(ID)%BA_PAR(ION_EQ,ION_V,K)  - DI_FAC
	              SE(ID)%BA_PAR(ION_EQ,NIV-1,K) =SE(ID)%BA_PAR(ION_EQ,NIV-1,K)  - ED_FAC
	              SE(ID)%BA_PAR(ION_EQ,NIV,K)   =SE(ID)%BA_PAR(ION_EQ,NIV,K)    - T_FAC
	            END IF
	          END DO
	        END DO
	      END IF
	    END DO
	  END IF
	END DO
!
	RETURN
	END
