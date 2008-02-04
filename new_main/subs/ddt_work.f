!
! Routine to compute change in gas enthalpy.
! Routine should be compatible with eval_temp_ddt_v2.
!
	SUBROUTINE DDT_WORK(WORK,POPS,T_FROM_J,OLD_T_RET,TIME_SEQ_NO,ND,NT)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 13-Dec-2005 : Based on EVAL_ADIABATIC_V3
!
	INTEGER NT
	INTEGER ND
!
! Output:
!
	REAL*8 WORK(ND)
	REAL*8 T_FROM_J(ND)
	REAL*8 OLD_T_RET(ND)
!
! Input:
!
	REAL*8 POPS(NT,ND)
	INTEGER TIME_SEQ_NO
!
! Local vectors.
!
	REAL*8 TCUR(ND)
	REAL*8 EK_VEC(ND)
	REAL*8 EI_VEC(ND)
	REAL*8 P_VEC(ND)
	REAL*8 GAMMA(ND)
	REAL*8 INT_EN(ND)
!
	REAL*8 OLD_POPS(NT,ND)
	REAL*8 OLD_R(ND)
	REAL*8 OLD_T(ND)
	REAL*8 OLD_ED(ND)
	REAL*8 OLD_GAMMA(ND)
	REAL*8 OLD_POP_ATOM(ND)
	REAL*8 OLD_INT_EN(ND)
!
	REAL*8 ION_EN(NT)
	REAL*8 TOT_ENERGY(NT)
	REAL*8 AVE_ENERGY(NT)
!
	REAL*8, ALLOCATABLE, SAVE :: T_STORE(:,:)
	REAL*8, ALLOCATABLE, SAVE :: INT_EN_STORE(:,:)
	LOGICAL MATCH
!
! Local variables.
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	INTEGER ERROR_LU
	REAL*8 BOLTZMANN_CONSTANT,FUN_PI
	EXTERNAL BOLTZMANN_CONSTANT,FUN_PI,ERROR_LU
!
	REAL*8 SCALE
	REAL*8 T1,T2,PI
	REAL*8 DELTA_T_SECS
	INTEGER I,J,K,L
	INTEGER LUER
	INTEGER ISPEC
	INTEGER ID
	INTEGER LU
	INTEGER LUIN
	LOGICAL WRITE_CHK
!
	LOGICAL, SAVE :: FIRST=.TRUE.
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	CALL GET_LU(LU)
	IF(FIRST)THEN
	  OPEN(UNIT=LU,FILE='DDT_WORK_CHK',STATUS='UNKNOWN')
	ELSE
	  OPEN(UNIT=LU,FILE='DDT_WORK_CHK',STATUS='OLD',POSITION='APPEND')
	END IF
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')'Beginning new DDT_WORK cycle'
	WRITE(LU,'(A)')' '
!
	IF(T_FROM_J(1) .EQ. 0.0D0)THEN
	  TCUR(1:ND)=T(1:ND)
	  T_FROM_J(1:ND)=T(1:ND)
	ELSE
	  TCUR(1:ND)=T_FROM_J(1:ND)
	END IF
!
	IF(.NOT. ALLOCATED(T_STORE))THEN
	  ALLOCATE (T_STORE(ND,20))
	  ALLOCATE (INT_EN_STORE(ND,20))
	  T_STORE=0.0D0
	  INT_EN_STORE=0.0D0
	END IF
!
! Define the average energy of each super level. At present this is
! depth independent, which should be adequate for most models.
! This average energy is used to scale the line cooling rates in
! the radiative equilibrium equation so that is more consistent
! with the electron cooling rate. The need for this scaling
! arises when levels within a super level have a 'relatively large'
! energy separation, and the dominat rates are scattering.
!
        AVE_ENERGY(:)=0.0D0
	DO ID=1,NUM_IONS-1
	  CALL AVE_LEVEL_ENERGY(AVE_ENERGY, ATM(ID)%EDGEXzV_F,
	1         ATM(ID)%GXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%EQXzV,
	1         ATM(ID)%NXzV,   ATM(ID)%NXzV_F, NT, ATM(ID)%XzV_PRES)
	END DO
!
! Compute the total excitation energy of each level.
!
	TOT_ENERGY(1:NT)=0.0D0
	ION_EN(1:NT)=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  T1=0.0D0
	  T2=0.0D0
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    T2=T2+AVE_ENERGY(ATM(ID)%EQXzV)
	    DO I=1,ATM(ID)%NXzV
	      J=ATM(ID)%EQXzV+I-1
	      TOT_ENERGY(J)=(AVE_ENERGY(ATM(ID)%EQXzV)-AVE_ENERGY(J))+T1
	      ION_EN(J)=T2
	    END DO
	    J=ATM(ID)%EQXzV
	    T1=T1+AVE_ENERGY(J)			!Adding on ionization energy
	  END DO
	  ID=SPECIES_END_ID(ISPEC)-1
	  IF(ID .GT. 0)THEN
	    J=ATM(ID)%EQXzV
	    TOT_ENERGY(J+ATM(ID)%NXzV)=T1
	    ION_EN(J+ATM(ID)%NXzV)=T2
	  END IF
	END DO
!
! Get the populations at the previous time step. These are put onto the same
! V grid as the curent model. The returned populations are NOT corected
! for advection, radioactive decays, and are not normalized.
!
	LUIN=7
	CALL GET_POPS_AT_PREV_TIME_STEP_V4(OLD_POPS,OLD_R,
	1      L_FALSE,L_FALSE,L_FALSE,TIME_SEQ_NO,ND,NT,LUIN)
	OLD_ED(:)=OLD_POPS(NT-1,:)
	OLD_T(:)=OLD_POPS(NT,:)
	OLD_POP_ATOM=0.0D0
	DO I=1,ND
	  DO J=1,NT-2
	    OLD_POP_ATOM(I)=OLD_POP_ATOM(I)+OLD_POPS(J,I)
	  END DO
	END DO
!
	WRITE(LU,'(A,6ES14.4)')'OLD_T',OLD_T(ND-5:ND)
	WRITE(LU,'(A,6ES14.4)')'    T',T(ND-5:ND)
	WRITE(LU,'(A,6ES14.4)')'CUR_T',TCUR(ND-5:ND)
!
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    T1=0.0D0
	    T2=0.0D0
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        T1=T1+POPS(K,J)
	        T2=T2+OLD_POPS(K,J)
	      END DO
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      I=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      T1=T1+POPS(I,J)
	      T2=T2+OLD_POPS(I,J)
	    END IF
	  END DO
	END DO
!
! Compute time step. The factor of 10^5 arises because R is in units of 10^10 cm, and
! V is in units of km/s.
!
	DELTA_T_SECS=1.0D+05*(R(ND)-OLD_R(ND))/V(ND)
!
! Compute the mean energy per atom. At first it is units of 10^15Hz.
!
	INT_EN(:)=0.0D0
	OLD_INT_EN(:)=0.0D0
	DO I=1,ND
	  DO J=1,NT-2
	     INT_EN(I)=INT_EN(I)+POPS(J,I)*TOT_ENERGY(J)
	     OLD_INT_EN(I)=OLD_INT_EN(I)+OLD_POPS(J,I)*TOT_ENERGY(J)
	  END DO
	END DO
	INT_EN=HDKT*INT_EN/POP_ATOM
	OLD_INT_EN=HDKT*OLD_INT_EN/OLD_POP_ATOM
!
! If new INT_EN calculation, store for interpolation purposes.
!
	DO I=1,ND
	  MATCH=.FALSE.
	  DO J=1,20
	    IF(T(I) .EQ. T_STORE(I,J))THEN
	      MATCH=.TRUE.
	      EXIT
	    END IF
	  END DO
	  IF(.NOT. MATCH)THEN
	    DO J=1,20
	      IF(T_STORE(I,J) .EQ. 0.0D0)THEN
	        T_STORE(I,J)=T(I)
	        INT_EN_STORE(I,J)=INT_EN(I)
	        EXIT
	      ELSE IF(T(I) .LT. T_STORE(I,J))THEN
	             T_STORE(I,20:J+1:-1)=     T_STORE(I,19:J:-1)
	        INT_EN_STORE(I,20:J+1:-1)=INT_EN_STORE(I,19:J:-1)
	        T_STORE(I,J)=T(I)
	        INT_EN_STORE(I,J)=INT_EN(I)
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
!
	DO I=1,ND
	  IF(T_STORE(I,2) .EQ. 0)THEN
	    IF(ABS(OLD_T(I)/T(I)-1.0D0) > 0.02)THEN
	      T1=(TCUR(I)-OLD_T(I))/(T(I)-OLD_T(I))
	      INT_EN(I)=INT_EN(I)*T1+OLD_INT_EN(I)*(1.0D0-T1)
	    END IF
	  ELSE IF(TCUR(I) .LT. T_STORE(I,1))THEN
	    T1=(TCUR(I)-T_STORE(I,1))/(T_STORE(I,2)-T_STORE(I,1))
	    INT_EN(I)=INT_EN_STORE(I,2)*T1+INT_EN_STORE(I,1)*(1.0D0-T1)
	  ELSE
	    DO J=2,20
	      IF(T_STORE(I,J) .EQ. 0.0D0)THEN
	        T1=(TCUR(I)-T_STORE(I,J-2))/(T_STORE(I,J-1)-T_STORE(I,J-2))
	        INT_EN(I)=INT_EN_STORE(I,J-1)*T1+INT_EN_STORE(I,J-2)*(1.0D0-T1)
	        EXIT
	      ELSE IF(TCUR(I) .LE. T_STORE(I,J))THEN
	        T1=(TCUR(I)-T_STORE(I,J-1))/(T_STORE(I,J)-T_STORE(I,J-1))
	        INT_EN(I)=INT_EN_STORE(I,J)*T1+INT_EN_STORE(I,J-1)*(1.0D0-T1)
                EXIT
	      END IF   
 	    END DO
	  END IF
	END DO
!
! We now compute constants for each of the 4 terms. These make
! it simpler and cleaner for the evaluation of the linearization.
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the adiabatic cooling rate by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. We need to introduce another factor of 10^4 since T is 
! in units of 10^4K.
!
	PI=FUN_PI()
	SCALE=1.0D+14*BOLTZMANN_CONSTANT()/4.0D0/PI
	DO I=1,ND
	  EK_VEC(I)=1.5D0*SCALE*POP_ATOM(I)/DELTA_T_SECS
	  EI_VEC(I)=SCALE*POP_ATOM(I)/DELTA_T_SECS
          P_VEC(I)=-SCALE*(POP_ATOM(I)+ED(I))*TCUR(I)/DELTA_T_SECS
	END DO
!
	DO I=1,ND
	  GAMMA(I)=ED(I)/POP_ATOM(I)
	  OLD_GAMMA(I)=OLD_ED(I)/OLD_POP_ATOM(I)
	END DO
!
	DO I=1,ND
 	    WORK(I)=EK_VEC(I)*( (1.0D0+GAMMA(I))*TCUR(I)- (1.0D0+OLD_GAMMA(I))*OLD_T(I) ) +
	1           EI_VEC(I)*(INT_EN(I)-OLD_INT_EN(I))      +
	1           P_VEC(I)*LOG(POP_ATOM(I)/OLD_POP_ATOM(I)) 
	END DO
	OLD_T_RET(1:ND)=OLD_T(1:ND)
!
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')'Term comparisons'
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,11(4X,A))')'    I','  Cur. T','   Old T','  Ek(cur)','  Ek(old)','  IE(cur)','  IE(old)',
	1               'Rek(cur)','Rek(old)','RIE(cur)','RIE(old)','   Pterm'
	T1=16.0D0*1.0D+16*5.67D-05/2.998D+10
	DO I=1,ND
	    T2=1.0D+04*BOLTZMANN_CONSTANT()*POP_ATOM(I)
	    WRITE(45,'(I4,15ES12.4)')I,TCUR(I),OLD_T(I),
	1           T1*TCUR(I)**4,T1*OLD_T(I)**4,
	1           1.5D0*T2*(1.0D0+GAMMA(I))*TCUR(I), 1.5D0*T2*(1.0D0+OLD_GAMMA(I))*OLD_T(I),
	1           T2*INT_EN(I), T2*OLD_INT_EN(I),
	1           EK_VEC(I)*(1.0D0+GAMMA(I))*TCUR(I), EK_VEC(I)*(1.0D0+OLD_GAMMA(I))*OLD_T(I),
	1           EI_VEC(I)*INT_EN(I), EI_VEC(I)*OLD_INT_EN(I),
	1           P_VEC(I)*LOG(POP_ATOM(I)/OLD_POP_ATOM(I))
	END DO
	WRITE(LU,'(A,11(4X,A))')'    I','  Cur. T','   Old T','  Ek(cur)','  Ek(old)','  IE(cur)','  IE(old)',
	1               'Rek(cur)','Rek(old)','RIE(cur)','RIE(old)','   Pterm'
	CLOSE(UNIT=LU)
!
	FIRST=.FALSE.
	RETURN
	END
