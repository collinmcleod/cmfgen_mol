!
! Subroutine to compute the contribution to the opacity AND emissivity
! by FREE-FREE and BOUND-FREE processes for a general ion. The
! contribution is added directly to the opacity CHI and emissivity ETA.
!
	SUBROUTINE COMP_FREE_FREE_V2(CHI_FF,ETA_FF,VCHI_FF,VETA_FF,CONT_FREQ,FREQ,
	1                 INIT,USE_EHB,DO_VAR,ND,NT)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 04-May-2022 : CHange HMI to H0
! Created 11-Jul-2019 : Based on COMP_OPAC and GENOPAETA_V10
!
	INTEGER ND
	INTEGER NT
	REAL*8 FREQ
	REAL*8 CONT_FREQ
	REAL*8 CHI_FF(ND)			!Opacity
	REAL*8 ETA_FF(ND)			!Emissivity
	REAL*8 VCHI_FF(NT,ND)			!Opacity
	REAL*8 VETA_FF(NT,ND)			!Emissivity
	LOGICAL INIT
	LOGICAL USE_EHB
	LOGICAL DO_VAR
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Vectors to save computational effort.
!
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 GFF_VAL(ND)		!g(ff) as a function of depth
!
	REAL*8, ALLOCATABLE, SAVE :: POP_SUM(:,:)
	REAL*8, ALLOCATABLE, SAVE :: GFF_STORE(:,:)
!
! Local constants.
!
	INTEGER ID
	INTEGER EQION
	INTEGER I,K,L
!
	REAL*8 T1
	REAL*8 NU
	REAL*8 ALPHA,TCHI1,TETA1
	REAL*8 EMIS
	REAL*8 HNUONKT
!
	IF(INIT)THEN
	  IF(.NOT. ALLOCATED(POP_SUM))ALLOCATE(POP_SUM(ND,NUM_IONS))
	  IF(.NOT. ALLOCATED(GFF_STORE))ALLOCATE(GFF_STORE(ND,20))
	  GFF_STORE=0.0D0
	  POP_SUM=0.0D0
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      POP_SUM(1:ND,ID)=SUM(ATM(ID+1)%XzV,1)
	      I=NINT(ATM(ID)%ZXzV)
	      IF(I .GT. 20)THEN
	        WRITE(6,*)'Insufficient storage for GFF_STORE in COMP_FREE_FREE_V2'
	        STOP
	      ELSE IF(I .LE. 0)THEN
	      ELSE IF(GFF_STORE(1,I) .EQ. 0.0D0)THEN
	        CALL GFF_VEC(GFF_STORE(1,I),NU,T,ATM(ID)%ZXzV,ND)
	      END IF
	    END IF
	  END DO
	END IF
!
	NU=FREQ
        T1=-HDKT*NU
        DO I=1,ND
          EMHNUKT(I)=EXP(T1/T(I))
          CHI_FF(I)=0.0D0
          ETA_FF(I)=0.0D0
        END DO
	VCHI_FF=0.0D0
	VETA_FF=0.0D0
!
! Add in free-free contribution. Because SN can be dominated by elements other
! than H and He, we now sum over all levels. To make sure that we only do this
! once, we only include the FREE-FREE contribution for the ion when PHOT_ID is one.
!
!!!$OMP  PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(+:CHI_FF,ETA_FF)
!!!$OMP1 PRIVATE(ID,I,GFF_VAL,TCHI1,TETA1,ALPHA)
	DO ID=1,NUM_IONS
	  IF(.NOT. ATM(ID)%XzV_PRES)THEN
	  ELSE IF(ION_ID(ID) .EQ. 'H0')THEN
	    I=7                                 !Used to read in data on first entry (will not be used here.)
	    CALL DO_H0_FF(ETA_FF,CHI_FF,POP_SUM(1,ID),ED,T,EMHNUKT,NU,I,ND)
	  ELSE IF(NINT(ATM(ID)%ZXzV) .GT. 0.0D0)THEN
!
! Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
!
	    GFF_VAL=GFF_STORE(:,NINT(ATM(ID)%ZXzV))
	    CALL FF_RES_GAUNT(GFF_VAL,NU,T,ID,ATM(ID)%GIONXzV_F,ATM(ID)%ZXzV,ND)
!
! We use POP_SUM as a temporary vector containing the sum of all level populations in
! the ion at each depth.
!
	    TCHI1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV/(FREQ*FREQ*FREQ)
	    TETA1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV*TWOHCSQ
	    DO I=1,ND
	      ALPHA=ED(I)*POP_SUM(I,ID)*GFF_VAL(I)/SQRT(T(I))
	      CHI_FF(I)=CHI_FF(I)+TCHI1*ALPHA*(1.0D0-EMHNUKT(I))
	      ETA_FF(I)=ETA_FF(I)+TETA1*ALPHA*EMHNUKT(I)
	    END DO
	  END IF
	END DO
!!!$OMP END PARALLEL DO
!
	IF(USE_EHB .AND. DO_VAR)THEN
	  DO ID=1,NUM_IONS
	    IF(.NOT. ATM(ID)%XzV_PRES)THEN
	    ELSE
	      GFF_VAL=GFF_STORE(:,NINT(ATM(ID)%ZXzV))
	      CALL FF_RES_GAUNT(GFF_VAL,NU,T,ID,ATM(ID)%GIONXzV_F,ATM(ID)%ZXzV,ND)
	      EQION=ATM(ID+1)%EQXzV
	      TCHI1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV/(FREQ*FREQ*FREQ)
	      TETA1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV*TWOHCSQ
	      DO I=1,ND
	        ALPHA=POP_SUM(I,ID)*GFF_VAL(I)/SQRT(T(I))
	        VCHI_FF(NT-1,I)=VCHI_FF(NT-1,I)+TCHI1*ALPHA*(1.0D0-EMHNUKT(I))
	        VETA_FF(NT-1,I)=VETA_FF(NT-1,I)+TETA1*ALPHA*EMHNUKT(I)
!
	        ALPHA=TCHI1*ED(I)*GFF_VAL(I)*(1.0D0-EMHNUKT(I))/SQRT(T(I))
	        EMIS=TETA1*ED(I)*GFF_VAL(I)*EMHNUKT(I)/SQRT(T(I))
	        DO L=1,ATM(ID+1)%NXzV
	           K=EQION+L-1
	           VCHI_FF(K,I)=VCHI_FF(K,I)+ALPHA
	           VETA_FF(K,I)=VETA_FF(K,I)+EMIS
	        END DO
!
	        HNUONKT=HDKT*NU/T(I)
	        ALPHA=ED(I)*POP_SUM(I,ID)*GFF_VAL(I)/SQRT(T(I))/T(I)
		VCHI_FF(NT,I)=VCHI_FF(NT,I)-TCHI1*ALPHA*(0.5D0+(HNUONKT-0.5D0)*EMHNUKT(I))
		VETA_FF(NT,I)=VETA_FF(NT,I)+TETA1*ALPHA*(HNUONKT-0.5D0)*EMHNUKT(I)
	      END DO
	    END IF
	  END DO
	END IF
!
	RETURN
	END
