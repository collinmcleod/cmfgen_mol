C
C General routine to
C
C   1:  compute the LTE populations of the levels in the FULL atom given the
C       POPULATIONS in the super level atoms,
C   2:  compute the populations of the levels in the full atom given the
C       POPULATIONS in the supel level atoms,
C   3:  compute the populations of the levels in the super-level atom,
C   4:  compute the occupation probailities, and
C   5:  routine the ion population used to compute the LTE population in the
C       FULL model atom.
C
C for C2.
C
C Routine is written for any 2 successive ionization stages --- not just
C C2 and CIII.
C
C Notation:
C
C         We use _F to denote populations and variables for the FULL atom,
C            with all terms and levels treated separately.
C	  We use _S to denote populations and variables for the SMALL model
C            atom, with many terms and levels treated as one (i.e using
C            SUPER levels).
C
	SUBROUTINE SUP_TO_FULL_V3(
	1   C2_F,C2LTE_F,W_C2_F,EDGEC2_F,
	1   GC2_F,F_TO_S_MAP_C2,INT_SEQ_C2,NC2_F,DIC2_F,GIONC2_F,
	1   C2_S,C2LTE_S,NC2_S,DIC2_S,ZC2,C2_PRES,
	1   EDGECIII_F,GCIII_F,F_TO_S_MAP_CIII,NCIII_F,
	1   CIII_PRES,T,ED,ND)
	IMPLICIT NONE
C
C Altered 27-may-1996 : Dynamic arrays installed for GION,EDGE_S, SUM and CNT.
C
C Altered 02-Jan-1996. INT_SEQ_C2 inserted in call.
C                      This variable will us to estimate non constant
C                        departure coeficients in a super level. This
C                        has been done specifically to improve the treatment
C                        of hign n lines in HeII and HeI.
C
C Altered 24-Oct-1995. GIONC2_S deleted from call.
C
	INTEGER*4 ND
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature 10^4K
	REAL*8 DIC2_S(ND)		!Ion density (Super levels)
	REAL*8 DIC2_F(ND)		!Ion density (Full model atom)
	REAL*8 ZC2			!Ion charge
C
	INTEGER*4 NC2_F
	REAL*8 C2_F(NC2_F,ND)
	REAL*8 C2LTE_F(NC2_F,ND)
	REAL*8 W_C2_F(NC2_F,ND)
	REAL*8 EDGEC2_F(NC2_F)
	REAL*8 GC2_F(NC2_F)
	INTEGER*4 F_TO_S_MAP_C2(NC2_F)
	INTEGER*4 INT_SEQ_C2(NC2_F)
	REAL*8 GIONC2_F
C
	INTEGER*4 NC2_S
	REAL*8 C2_S(NC2_S,ND)
	REAL*8 C2LTE_S(NC2_S,ND)
C
	LOGICAL C2_PRES,CIII_PRES
C
	INTEGER*4 NCIII_F
	REAL*8 EDGECIII_F(NCIII_F)
	REAL*8 GCIII_F(NCIII_F)
	INTEGER*4 F_TO_S_MAP_CIII(NCIII_F)
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER*4 ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER*4 I,K,L,M,LU_ER
	REAL*8 B,B1,B2,T1
	REAL*8 X,Y,RGU
	INTEGER*4 INT_SL
C
C Dynamic arrays.
C
	REAL*8 GION(ND)
	REAL*8 EDGE_S(NC2_S)
	REAL*8 SUM(NC2_S)
	INTEGER*4 CNT(NC2_S)
C
	IF(.NOT. C2_PRES)RETURN
C
C Compute the partition function of the CIII g.s.. NB. If CIII is not present
C the statistical weight of the ion ground state must just be GIONC2_F.
C
	IF(CIII_PRES)THEN
	  DO K=1,ND
	    GION(K)=0.0D0
	    DO I=1,NCIII_F
	      IF(F_TO_S_MAP_CIII(I) .EQ. 1)THEN
	        GION(K)=GION(K)+GCIII_F(I)*EXP( HDKT*(EDGECIII_F(I)-
	1                                         EDGECIII_F(1))/T(K) )
	      END IF
	    END DO
	  END DO
	ELSE
	  DO K=1,ND
	    GION(K)=GIONC2_F
	  END DO
	END IF
C
C Compute the ion density used to compute LTE populations in the full atom.
C This is essentially the ground-state population.
C
	DO K=1,ND
	  DIC2_F(K)=DIC2_S(K)*GIONC2_F/GION(K)
	END DO
C
C Compute the occupation probabilities.
C
	CALL OCCUPATION_PROB(W_C2_F,EDGEC2_F,ZC2,NC2_F,ND)
C
C Since no the effective statistical weight, can now compute the LTE
C populations of the levels in the full atom.
C
	T1=2.07078D-22
	RGU=LOG(T1)
	DO K=1,ND
	  X=HDKT/T(K)
	  Y=ED(K)*DIC2_S(K)*( T(K)**(-1.5) )/GION(K)
	  DO I=1,NC2_F
	    C2LTE_F(I,K)=W_C2_F(I,K)*GC2_F(I)*Y*EXP(EDGEC2_F(I)*X+RGU)
	  END DO
	END DO
C
C Compute the LTE pops in the atom with super levels, after initializing them.
C
	DO K=1,ND
	  DO I=1,NC2_S
	    C2LTE_S(I,K)=0.0D0
	  END DO
	END DO
C
	DO K=1,ND
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    C2LTE_S(L,K)=C2LTE_S(L,K)+C2LTE_F(I,K)
	  END DO
	END DO
C
C Can now compute populations in full atom.
C
	DO K=1,ND
C
	  DO L=1,NC2_S
	    CNT(L)=0
	    EDGE_S(L)=0.0D0
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    EDGE_S(L)=EDGE_S(L)+EDGEC2_F(I)*C2LTE_F(I,K)
	    CNT(L)=CNT(L)+1
	  END DO
	  DO L=1,NC2_S
	    EDGE_S(L)=EDGE_S(L)/C2LTE_S(L,K)
	  END DO
C
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    IF(CNT(L) .EQ. 1 .OR. INT_SEQ_C2(I) .EQ. 0)THEN
	      C2_F(I,K)=C2LTE_F(I,K)*(C2_S(L,K)/C2LTE_S(L,K))
	    ELSE
C
C Find the closest level of the same interpolating sequence. We first
C attempt to bracket (in energy) the enegry of the level whose population
C is to be determined.
C
	      INT_SL=0
	      IF(EDGEC2_F(I) .LE. EDGE_S(L))THEN
	        M=I+1
	        DO WHILE(M .LE. NC2_F .AND. INT_SL .EQ. 0)
	          IF(INT_SEQ_C2(M) .EQ. INT_SEQ_C2(I) .AND.
	1             F_TO_S_MAP_C2(I) .NE. F_TO_S_MAP_C2(M))THEN
	            INT_SL=M
	          END IF
	          M=M+1
	        END DO
	      END IF
	      IF(INT_SL .EQ. 0)THEN
	        M=I-1
	        DO WHILE(M .GE. 1 .AND. INT_SL .EQ. 0)
	          IF(INT_SEQ_C2(M) .EQ. INT_SEQ_C2(I) .AND.
	1             F_TO_S_MAP_C2(I) .NE. F_TO_S_MAP_C2(M))THEN
	            INT_SL=M
	          END IF
	          M=M-1
	        END DO
	      END IF
	      IF(INT_SL .EQ. 0)THEN
	        LU_ER=ERROR_LU()
	        WRITE(LU_ER,*)'Error in SUP_TO_FULL_V3'
	        WRITE(LU_ER,*)'No interpolating sequence found'
	        WRITE(LU_ER,*)'N_S=',NC2_S
	        WRITE(LU_ER,*)'N_F=',NC2_F
	        WRITE(LU_ER,*)'Full level is:',I
	        STOP
	      END IF
	      INT_SL=F_TO_S_MAP_C2(INT_SL)
C
	      T1=LOG(EDGE_S(L)/EDGEC2_F(I)) /
	1                  LOG(EDGE_S(L)/EDGE_S(INT_SL))
	      B1=C2_S(INT_SL,K)/C2LTE_S(INT_SL,K)
	      B2=C2_S(L,K)/C2LTE_S(L,K)
              B=T1*B1 + (1.0D0-T1)*B2
C
C Constrain the interpolation. Hopefully this is not necessary.
C
	      IF(B1 .LE. 1. .AND. B2 .LE. 1 .AND. B .GT. 1)B=1.0
	      IF(B1 .GE. 1. .AND. B2 .GE. 1 .AND. B .LT. 1)B=1.0
	      IF(B .LT. 0)B=B1
	      C2_F(I,K)=C2LTE_F(I,K)*B
	    END IF
	  END DO
C
C The revised departure coefficents now needs adjusting so that the total
C population of the C2 levels in the full atom matches the corresponding
C super level. This correction will generally be small.
C
	  DO L=1,NC2_S
	    SUM(L)=0.0D0
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    SUM(L)=SUM(L)+C2_F(I,K)
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    C2_F(I,K)=C2_F(I,K)*(C2_S(L,K)/SUM(L))
	  END DO
C
	END DO
C
	RETURN
	END
