!
! Subroutine to determine correspondence between SPECIES in the CHARGE exchange
! reactions, and the corresponding program variables. It is advised 
! (and program checks) that levels involved in charge exchange reactions should
! be distinct super-levels. 
!
! Program assumes (and checks) that each charge reaction involves full terms.
!
! Charge exchange reactions must have the form (and ordering)
!
!     Y(n+) + X([m-1]+)  <--> Y([n-1]+) + X(M+)
!
	SUBROUTINE SET_CHG_LEV_ID_V4(ND,LUOUT)
	USE MOD_CMFGEN
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
	INTEGER ND,LUOUT
!
! Altered  11-Apr-2002 : Changed to allow automatic splitting of charge exchange rates
!                          among a single LS state.
! Altered  09-Oct-1999 : Error reporting inproved.
! Created  20-Aug-1998 : BASED on V1. Very different calls and a change in
!                          philosiphy of how super levels are managed.
!
! Local variables
!
	INTEGER, PARAMETER :: IZERO=0
!
	INTEGER ID
	INTEGER, ALLOCATABLE :: ID_POINTER(:,:)
	INTEGER, ALLOCATABLE :: LEV_CNT(:,:)
	INTEGER, ALLOCATABLE :: CHG_ID(:,:)
	REAL, ALLOCATABLE :: G_SUM(:,:)
!
	INTEGER NOUT
	INTEGER NIN
	INTEGER J,K,L
	INTEGER IPOS
	INTEGER I_S,I_F
	INTEGER LST
	INTEGER II
	INTEGER IP
	INTEGER IOS
	INTEGER N_CHG_OMITTED
	CHARACTER(LEN=30) LOC_NAME,NAME1,NAME2
	CHARACTER(LEN=132) STRING
	LOGICAL LEVEL_SET
!
	IF(.NOT. DO_CHG_EXCH)RETURN
!
!	WRITE(LUER,*)'Entered SET_CHG_LEV_ID_V4'
	ALLOCATE (LEV_CNT(N_CHG_RD,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (G_SUM(N_CHG_RD,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ID_POINTER(N_CHG_RD,4),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in SET_CHG_LEV_ID_V4'
	  WRITE(LUER,*)'Unable to allocate memory'
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
!
! Perform iniilizations.
!
        LEV_CNT(:,:)=0
        G_SUM(:,:)=0
	ID_POINTER(:,:)=0
!
! Count actual number of reactions. This is to allow for split J levels.
!
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO J=1,N_CHG_RD
	      DO K=1,4
	        IF(SPEC_ID_CHG_RD(J,K) .EQ. ION_ID(ID))THEN
	          ID_POINTER(J,K)=ID
!
! We first check if species 1 or 4 corresponds to the final ioization state.
!
	          IF( (K .EQ. 2 .OR. K .EQ. 3) .AND. 
	1           ATM(ID)%EQXzV+ATM(ID)%NXzV .EQ. EQ_SPECIES(SPECIES_LNK(ID)) )THEN
	            IF(K .EQ. 3)THEN
	              LEV_CNT(J,1)=1
	              G_SUM(J,1)=ATM(ID)%GIONXzV_F
	            ELSE
	              LEV_CNT(J,4)=1
	              G_SUM(J,4)=ATM(ID)%GIONXzV_F
	            END IF
	          END IF
!
! We now check against regular levels.
! 
	          I_S=0
	          DO I_F=1,ATM(ID)%NXzV_F
	            LOC_NAME=ATM(ID)%XzVLEVNAME_F(I_F)
	            IF(INDEX(LEV_NAME_CHG_RD(J,K),'[') .EQ. 0)THEN
	              IPOS=INDEX(LOC_NAME,'[')
	              IF(IPOS .NE. 0)LOC_NAME=LOC_NAME(1:IPOS-1)
	            END IF
	            IF(LOC_NAME .EQ. LEV_NAME_CHG_RD(J,K) .OR.
	1               LOC_NAME .EQ. ALT_LEV_NAME_CHG_RD(J,K))THEN
	              IF(I_S .EQ. 0)THEN
                        I_S=ATM(ID)%F_TO_S_XZv(I_F)
	                LEV_CNT(J,K)=1
	              ELSE IF(I_S .EQ.  ATM(ID)%F_TO_S_XZv(I_F))THEN
!
! Do nothing as same Super Level.
!
	              ELSE
!
! Another SL
!
	                LEV_CNT(J,K)=LEV_CNT(J,K)+1
	              END IF
	              G_SUM(J,K)=G_SUM(J,K)+ATM(ID)%GXzV_F(I_F)
	            END IF
	          END DO
!
	        END IF			!Species match
	      END DO			!Over variable in charge exch. reaction
	    END DO			!Over charge reaction
	  END IF                        !Species present
	END DO				!Over species
!
! Now determine total number ocharge reactions.
!
	N_CHG=0
	N_CHG_OMITTED=0
	DO J=1,N_CHG_RD
          K=LEV_CNT(J,1)*LEV_CNT(J,2)*LEV_CNT(J,3)*LEV_CNT(J,4)
	  IF(K .EQ. 0)THEN
	    N_CHG_OMITTED=N_CHG_OMITTED+1
	    CHG_INCLUDED_RD(J)=.FALSE.
          ELSE
	    N_CHG=N_CHG + K 
	  END IF
	END DO
	WRITE(LUER,*)'Number of charge transitions read in is',N_CHG_RD
	WRITE(LUER,*)'Number of charge transitions omitted is ',N_CHG_OMITTED
	WRITE(LUER,*)'Number of revixsed charge transitions, including split state, is',N_CHG
!
	ALLOCATE (TYPE_CHG(N_CHG),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TLO_CHG(N_CHG),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (THI_CHG(N_CHG),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (COEF_CHG(N_CHG,N_COEF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SPEC_ID_CHG(N_CHG,4),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE (CHG_REACTION_AVAILABLE(N_CHG),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ID_ION_CHG(N_CHG,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LEV_IN_POPS_CHG(N_CHG,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LEV_IN_ION_CHG(N_CHG,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (G_CHG(N_CHG,4),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (Z_CHG(N_CHG,4),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (AI_AR_CHG(N_CHG,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (dlnAI_AR_CHG_dlnT(N_CHG,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (COOL_CHG(N_CHG,ND),STAT=IOS)
!
        IF(IOS .EQ. 0)ALLOCATE (CHG_ID(N_CHG,ND),STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in SET_CHG_LEV_ID_V4'
	  WRITE(LUER,*)'Unable to allocate main charge exchange memory'
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
	WRITE(LUER,*)'Done memory allocations for charge exchange reactions'
!
! Perform inilizations.
!
	ID_ION_CHG(:,:)=0
	LEV_IN_POPS_CHG(:,:)=0
	LEV_IN_ION_CHG(:,:)=0
!
! Real variables
!
	Z_CHG(:,:)=0.0D0
        AI_AR_CHG(:,:)=0.0D0
        dlnAI_AR_CHG_dlnT(:,:)=0.0D0
        COOL_CHG(:,:)=0.0D0
	INITIALIZE_ARRAYS=.TRUE.
!
! Determine ionization stage and levels for species involved in charge 
! exchange reactions. Now determine whether the present species is in 
! the CHARGE exchange reaction list.
!
	L=0
	DO J=1,N_CHG_RD
	  WRITE(117,*)'Operating on charge exchange reaction J=',J
!
! NOUT and NIN are use to loop over ALL possible charge exchnage reactions.
! For LS coupling, we assume that the reaction rates are independent of
! the initial states, and proportional to the statistical weight of the final
! state.
!
	  NOUT=1
	  NIN=LEV_CNT(J,1)*LEV_CNT(J,2)*LEV_CNT(J,3)*LEV_CNT(J,4)
	  WRITE(117,*)NIN,LEV_CNT(J,1),LEV_CNT(J,2),LEV_CNT(J,3),LEV_CNT(J,4)
	  IF(NIN .NE. 0)THEN
	    LST=L
	    DO K=1,4
	      L=LST
	      NIN=NIN/LEV_CNT(J,K)
!
! Level 1 or 4 may be the final ionization stage, which we first check.
!
	      ID=1
	      IF(K .EQ. 1)ID=ID_POINTER(J,3)
	      IF(K .EQ. 4)ID=ID_POINTER(J,2)
	      I_S=(K-1)*(K-4)         
	      IF( I_S .EQ. 0 .AND. ATM(ID)%EQXzV+ATM(ID)%NXzV .EQ. EQ_SPECIES(SPECIES_LNK(ID)) )THEN
	        DO IP=1,NOUT
	          L=L+1
	          Z_CHG(L,K)=ATM(ID)%ZXzV
	          ID_ION_CHG(L,K)=ID+1
	          LEV_IN_POPS_CHG(L,K)=ATM(ID)%EQXzV+ATM(ID)%NXzV
	          LEV_IN_ION_CHG(L,K)=1
	          CHG_ID(L,K)=J
	          G_CHG(L,K)=ATM(ID)%GIONXzV_F
	          WRITE(117,*)J,K,L,Z_CHG(L,K),ID+1
	        END DO
	      ELSE
	        ID=ID_POINTER(J,K)
	        DO IP=1,NOUT
	          I_S=0
	          DO I_F=1,ATM(ID)%NXzV_F
	            LOC_NAME=ATM(ID)%XzVLEVNAME_F(I_F)
	            IF(INDEX(LEV_NAME_CHG_RD(J,K),'[') .EQ. 0)THEN
	              IPOS=INDEX(LOC_NAME,'[')
	              IF(IPOS .NE. 0)LOC_NAME=LOC_NAME(1:IPOS-1)
	            END IF
	            IF(LOC_NAME .EQ. LEV_NAME_CHG_RD(J,K) .OR.
	1               LOC_NAME .EQ. ALT_LEV_NAME_CHG_RD(J,K))THEN
	              IF(I_S .EQ. 0 .OR. I_S .NE. ATM(ID)%F_TO_S_XzV(I_F))THEN
	                DO II=1,NIN
	                  L=L+1
	                  ID_ION_CHG(L,K)=ID     
	                  I_S=ATM(ID)%F_TO_S_XzV(I_F)
	                  LEV_IN_ION_CHG(L,K)=I_S
	                  LEV_IN_POPS_CHG(L,K)=ATM(ID)%EQXzV+I_S-1
	                  Z_CHG(L,K)=ATM(ID)%ZXzV-1.0D0
	                  CHG_ID(L,K)=J
	                  G_CHG(L,K)=ATM(ID)%GXzV_F(I_F)
	                  WRITE(117,*)'B',J,K,L,Z_CHG(L,K),ID+1
	                END DO
	              ELSE
	                DO II=1,NIN
	                  G_CHG(L,K)=G_CHG(L,K)+ATM(ID)%GXzV_F(I_F)
	                END DO
	                WRITE(117,*)'C',J,K,L,Z_CHG(L,K),ID+1
	              END IF
	            END IF
	          END DO
	        END DO
	      END IF
	      NOUT=NOUT*LEV_CNT(J,K)
!
	    END DO                   !Loop over K
!
! Save the reaction rates for each of the new reactions.
!
	    DO K=LST+1,L
	      TYPE_CHG(K)=TYPE_CHG_RD(J)
	      TLO_CHG(K)=TLO_CHG_RD(J)
	      THI_CHG(K)=THI_CHG_RD(J)
	      SPEC_ID_CHG(K,1:4)=SPEC_ID_CHG_RD(J,1:4)
	      COEF_CHG(K,1:N_COEF_MAX)=COEF_CHG_RD(J,1:N_COEF_MAX)
	      COEF_CHG(K,1)=COEF_CHG(K,1)*G_CHG(K,3)*G_CHG(K,4)/G_SUM(J,3)/G_SUM(J,4)
	    END DO
          END IF                    !Is reaction available
!	    
	END DO                      !Loop over J (reactions read in)
!
! Check that the SL's involved in a charge exchange reaction correspond to a
! single LS state.
!
	CALL GEN_ASCI_OPEN(LUOUT,'CHG_EXCH_CHK','UNKNOWN',' ','WRITE',IZERO,IOS)
	WRITE(117,'(A)')
	WRITE(117,'(A)')' Summary of charge exchange info'
	WRITE(117,'(A)')
	DO J=1,N_CHG
	  STRING=' '
	  DO K=1,4
	    ID=ID_ION_CHG(J,K)
	    IF(ATM(ID)%XzV_PRES)THEN
	      I_S=LEV_IN_ION_CHG(J,K)
	      NAME1=' '
	      DO I_F=1,ATM(ID)%NXzV_F
	        IF(ATM(ID)%F_TO_S_XzV(I_F) .EQ. I_S)THEN
                  IF(NAME1 .EQ. ' ')THEN
	            NAME1=ATM(ID)%XzVLEVNAME_F(I_F)
	            L=LEN_TRIM(STRING)
	            IF(K .NE. 1)L=MAX(L,(K-1)*25+5*MOD(K+1,2))
	            STRING(L+1:)=TRIM(SPEC_ID_CHG(J,K))//'{'//TRIM(NAME1)//'}'
	            IPOS=INDEX(NAME1,'[')
	            IF(IPOS .NE. 0)NAME1=NAME1(1:IPOS-1)
	          ELSE
	            NAME2=ATM(ID)%XzVLEVNAME_F(I_F)
	            IPOS=INDEX(NAME2,'[')
	            IF(IPOS .NE. 0)NAME2=NAME2(1:IPOS-1)
	            IF(NAME2 .NE. NAME1)THEN
	               WRITE(LUER,*)'Error in SET_CHG_LEV_ID'
	               WRITE(LUER,*)'Each charge exchange reaction must go to a SL corresponding to a single LS state'
	               WRITE(LUER,*)J,K,ID,I_S
	               WRITE(LUER,*)'NAME1=',NAME1
	               WRITE(LUER,*)'NAME2=',NAME2
	               STOP
	            END IF
	          END IF
	        END IF
	      END DO
	      IF(NAME1 .EQ. ' ')THEN
	        WRITE(LUER,*)'Error in SET_CHG_LEV_ID - no name match'
	        WRITE(LUER,*)J,K,ID,I_S
	        STOP
	      END IF
	    ELSE
	      L=LEN_TRIM(STRING)
	      IF(K .NE. 1)L=MAX(L,(K-1)*25+5*MOD(K+1,2))
	      STRING(L+1:)=TRIM(SPEC_ID_CHG(J,K))//'{ion}'
	    END IF
	  END DO
	  WRITE(LUOUT,'(1X,I3,3X,A,3X,ES12.4)')J,STRING(1:100),COEF_CHG(J,1)
	  WRITE(117,'(4I5,F5.2,I4,F7.2,A)')
	1       (TYPE_CHG(J),ID_ION_CHG(J,K),LEV_IN_ION_CHG(J,K),LEV_IN_POPS_CHG(J,K),Z_CHG(J,K),
	1       CHG_ID(J,K),G_CHG(J,K),TRIM(SPEC_ID_CHG(J,K)), K=1,4)
	END DO
!
	IF(N_CHG_OMITTED .NE. 0)THEN
          WRITE(LUOUT,'(A)')
	  WRITE(LUOUT,'(A)')' The following reactions were not included.'
	  WRITE(LUOUT,'(A)')' This may have been because a species was not included in this model.'
	  WRITE(LUOUT,'(A)')' Alternatively check species and level names for consistency with model atoms.'
	  WRITE(LUOUT,'(A)')
	  DO J=1,N_CHG_RD
	    IF(.NOT. CHG_INCLUDED_RD(J))THEN
	       STRING=' '
	       DO K=1,4
	         STRING(1+(K-1)*30:)=TRIM(SPEC_ID_CHG_RD(J,K))//'{'//TRIM(LEV_NAME_CHG_RD(J,K))//'}'
	       END DO
	       WRITE(LUOUT,'(1X,A)')TRIM(STRING)
	    END IF
	  END DO
	END IF
	CLOSE(LUOUT)
!
	DEALLOCATE (LEV_CNT)
	DEALLOCATE (CHG_ID)
	DEALLOCATE (G_SUM)
	DEALLOCATE (ID_POINTER)
!
	RETURN
	END
