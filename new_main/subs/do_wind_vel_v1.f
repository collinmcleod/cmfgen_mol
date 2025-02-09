	MODULE DO_WIND_VEL_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP), ALLOCATABLE :: GRAD_STORE(:)
	REAL(KIND=LDP), ALLOCATABLE :: GRAD_RD(:)
	REAL(KIND=LDP), ALLOCATABLE :: GELEC_RD(:)
	REAL(KIND=LDP), ALLOCATABLE :: VdVdR_RD(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: R_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: GRAD_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: GELEC_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: VdVdR_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC_OLD(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: RHAT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_NEW(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_NEW(:)
	REAL(KIND=LDP), ALLOCATABLE :: GELEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: GRAD(:)
	REAL(KIND=LDP), ALLOCATABLE :: GRAD_MIN(:)
	REAL(KIND=LDP), ALLOCATABLE :: dFdR(:)
	REAL(KIND=LDP), ALLOCATABLE :: F_CORR(:)
	REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)
!
	LOGICAL FIRST/.TRUE./
	INTEGER LST_V_UPDATE/0/
	INTEGER NUM_V_ITS/0/
	SAVE
!
	END MODULE DO_WIND_VEL_MOD
!
CONTAINS
	SUBROUTINE DO_WIND_VEL_V1(R_RD,V_RD,SIGMA_RD,CLUMP_FAC_RD,MSTAR,
	1                DONE_V_REVISION,MAIN_COUNTER,ND_RD)
	USE SET_KIND_MODULE
	USE DO_WIND_VEL_MOD
	USE UPDATE_KEYWORD_INTERFACE
	IMPLICIT NONE
!
	INTEGER ND_RD
	INTEGER MAIN_COUNTER
	REAL(KIND=LDP) MSTAR
	REAL(KIND=LDP) R_RD(ND_RD)
	REAL(KIND=LDP) V_RD(ND_RD)
	REAL(KIND=LDP) SIGMA_RD(ND_RD)
	REAL(KIND=LDP) CLUMP_FAC_RD(ND_RD)
	LOGICAL DONE_V_REVISION
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
        INTEGER, PARAMETER :: LUIN=7
        INTEGER, PARAMETER :: LUSCR=8
	INTEGER ND
	INTEGER ND_OLD
	INTEGER NINS       !Number of points to be inserted into grid (>_ 0)
!
	REAL(KIND=LDP) GCRIT
	REAL(KIND=LDP) RSTAR
	REAL(KIND=LDP) MSUN
	REAL(KIND=LDP) GRAV_CONST
	REAL(KIND=LDP) VCRIT
	REAL(KIND=LDP) VSOUND
	REAL(KIND=LDP) RSOUND
	REAL(KIND=LDP) RX		!Start radius for integration
	REAL(KIND=LDP) VX 		!Start velocity for integration
!
	REAL(KIND=LDP) LAMBERT_WM_FUN
	EXTERNAL LAMBERT_WM_FUN
!
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) GINT
	REAL(KIND=LDP) Z
	REAL(KIND=LDP) LNZ
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) ERROR, ALLOWED_CHANGE_IN_GRAD
	INTEGER I,J,K,IT
	INTEGER IOS
	INTEGER NUM_ITS
	INTEGER MAX_ITS
	INTEGER IT_FREQ
	INTEGER FIRST_IT
!
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) FILENAME
	LOGICAL F_FUN_V
	LOGICAL AVERAGE_V
	LOGICAL FILE_PRES
!
	DONE_V_REVISION=.FALSE.
	GRAV_CONST=6.67259E-08_LDP
	MSUN=1.989E+33_LDP
!
! Get parameters controlling the corrections to the velcoity law.
!
	CALL GEN_ASCI_OPEN(LUIN,'DO_V_PARAMS','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
!	  WRITE(6,*)'Error opening DO_V_PARAMS in WIND_HYD, IOS=',IOS
	  RETURN
	END IF
!
	ALPHA=0.5_LDP; ALLOWED_CHANGE_IN_GRAD=0.02_LDP; IT_FREQ=10; NINS=0
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
	CALL RD_STORE_INT(NUM_ITS,'N_GRAD_ITS',L_TRUE,'Number of iterations adjusting g(rad)')
	CALL RD_STORE_INT(ND,'DST',L_TRUE,'Depth to integrate from')
	CALL RD_STORE_INT(FIRST_IT,'FST_IT',L_TRUE,'First iteration')
	CALL RD_STORE_INT(MAX_ITS,'MAX_IT',L_TRUE,'Maximum number of iterations')
	CALL RD_STORE_INT(IT_FREQ,'IT_FREQ',L_FALSE,'Frequency of revisions of V')
	CALL RD_STORE_DBLE(ALLOWED_CHANGE_IN_GRAD,'dGRAD',L_FALSE,'Maximum fractional change in GRAD allowed')
	CALL RD_STORE_DBLE(VSOUND,'VSOUND',L_TRUE,'Sound speed in km/s')
	CALL RD_STORE_DBLE(ALPHA,'ALPHA',L_FALSE,'Exponent relating grad with vdvdr')
	CALL RD_STORE_INT(NINS,'NINS',L_FALSE,'Number of points/interval to be inserted into grid')
!	
	F_FUN_V=.FALSE.
	AVERAGE_V=.FALSE.
	CALL CLEAN_RD_STORE()
	CLOSE(UNIT=LUIN)
	CLOSE(UNIT=LUSCR)
!
	IF(MAIN_COUNTER .LT. FIRST_IT .OR. NUM_V_ITS .GE. MAX_ITS)THEN
	  RETURN
	ELSE IF(MAIN_COUNTER .EQ. FIRST_IT)THEN
	ELSE IF(MAIN_COUNTER .LT. LST_V_UPDATE+IT_FREQ)THEN
	  RETURN
	END IF
!
! Allocate all vectors.
!
	ND_OLD=ND_RD+NINS*(ND_RD-1)
	IF(FIRST)THEN
	  ALLOCATE(GELEC_RD(ND_RD))
	  ALLOCATE(GRAD_STORE(ND_RD)); GRAD_STORE=0.0_LDP
	  ALLOCATE(GRAD_RD(ND_RD))
	  ALLOCATE(VdVDR_RD(ND_RD))
!
	  ALLOCATE(R_OLD(ND_OLD))
	  ALLOCATE(V_OLD(ND_OLD))
	  ALLOCATE(SIGMA_OLD(ND_OLD))
	  ALLOCATE(GELEC_OLD(ND_OLD))
	  ALLOCATE(GRAD_OLD(ND_OLD))
	  ALLOCATE(VdVdR_OLD(ND_OLD))
	  ALLOCATE(CLUMP_FAC_OLD(ND_OLD))
!
	  ALLOCATE(RHAT(ND_OLD))
	  ALLOCATE(GELEC(ND_OLD))
	  ALLOCATE(GRAD(ND_OLD))
	  ALLOCATE(GRAD_MIN(ND_OLD))
	  ALLOCATE(CLUMP_FAC(ND_OLD))
	  ALLOCATE(dfdR(ND_OLD))
	  ALLOCATE(F_CORR(ND_OLD))
!
	  ALLOCATE(V_NEW(ND_OLD))
	  ALLOCATE(SIGMA_NEW(ND_OLD))
	  FIRST=.FALSE.
!
	END IF
!
! Get the hydrodynamic data.
!
! Note: V is in km/s, R in units of 10^10 cm.
!
! Note the VdV/dR is the same in program units as in cgs units.
! Thus dVdR_OLD will be in PORGAM units!
!
	OPEN(UNIT=10,FILE='HYDRO',STATUS='OLD',ACTION='READ')
	READ(10,'(A)')STRING
	DO I=1,ND_RD
	  READ(10,*)T1,T1,T1,VdVdR_RD(I),T1,T1,GRAD_RD(I),GELEC_RD(I)
	END DO
	CLOSE(UNIT=10)
	RSTAR=R_RD(ND_RD)
!
	ERROR=0.0_LDP
	DO I=1,ND_RD
	  T1=2.0_LDP*ABS(GRAD_RD(I)-GRAD_STORE(I))/( ABS(GRAD_RD(I))+ABS(GRAD_STORE(I)) )
	  ERROR=MAX(ERROR,T1)
	END DO
	IF(ERROR .GT. ALLOWED_CHANGE_IN_GRAD)THEN
	  WRITE(6,'(A)')'Change in GRAD too large to allow an revision of the wind velocity'
	  GRAD_STORE=GRAD_RD
	  RETURN
	END IF
	WRITE(6,'(A,F9.3)')' Maximum % change in GRAD since last iteration is',100.0*ERROR
!
! Interpolate onto the fine grid.
!
	WRITE(6,*)'Defining the new R grid'
	J=0
	DO I=1,ND_RD-1
	  J=J+1; R_OLD(J)=R_RD(I)
	  DO K=1,NINS
	    J=J+1; R_OLD(J)=R_OLD(J-1)+(R_RD(I+1)-R_RD(I))/(NINS+1)
	  END DO
	END DO
	J=J+1; R_OLD(J)=R_RD(ND_RD)
	WRITE(6,*)'Starting the interpolations',J,ND_OLD
	CALL MON_INTERP(V_OLD,        ND_OLD,IONE,R_OLD,ND_OLD,V_RD,       ND_RD,R_RD,ND_RD)
	CALL MON_INTERP(GRAD_OLD,     ND_OLD,IONE,R_OLD,ND_OLD,GRAD_RD,    ND_RD,R_RD,ND_RD)
	CALL MON_INTERP(GELEC_OLD,    ND_OLD,IONE,R_OLD,ND_OLD,GELEC_RD,   ND_RD,R_RD,ND_RD)
	CALL MON_INTERP(CLUMP_FAC_OLD,ND_OLD,IONE,R_OLD,ND_OLD,CLUMP_FAC_RD,ND_RD,R_RD,ND_RD)
	WRITE(6,*)'Done the interpolations'
!
!	DO I=1,ND_OLD
!	  WRITE(45,'(I3,6ES14.4)')I,R_OLD(I),V_OLD(I),GRAD_OLD(I),GELEC_OLD(I),CLUMP_FAC_OLD(I)
!	END DO
!
	ALLOCATE (COEF(ND_OLD,4))
        CALL MON_INT_FUNS_V2(COEF,V_OLD,R_OLD,ND_OLD)
          DO I=1,ND_OLD
	    WRITE(94,*)I,V_OLD(I),COEF(I,3),V_OLD(I)*COEF(I,3),
	1              V_OLD(I)*(V_OLD(I)-V_OLD(I+1))/(R_OLD(I)-R_OLD(I+1))
	    T1=COEF(I,3)
            SIGMA_OLD(I)=T1*R_OLD(I)/V_OLD(I)-1.0_LDP
	    VdVdR_OLD(I)=V_OLD(I)*T1
	  END DO
!
!
! Set parameters of model.
!
!	VSOUND=20.0D0; CALL GEN_IN(VSOUND,'Sound speed in km/s')
!	NUM_ITS=3; CALL GEN_IN(NUM_ITS,'Number of iterations adjusting g(rad)')
!	F_FUN_V=.TRUE.; CALL GEN_IN(F_FUN_V,'Is the clumping factor a function of V?')
!	AVERAGE_V=.TRUE.; CALL GEN_IN(AVERAGE_V,'Average Vnew with Vold?')
!
! Find the locaton of the sonic point.
!
! Currently read depth in.
!
	DO I=1,ND_OLD
	  IF(V_OLD(I) .LT. VSOUND)THEN
	    J=I-1
	    T1=(VSOUND-V_OLD(J+1))/(V_OLD(J)-V_OLD(J+1))
	    RSOUND=R_OLD(J)*T1 + (1.0_LDP-T1)*R_OLD(J+1)
	    EXIT
	  END IF
	END DO
!	IF(ND < 0)THEN
!	  ND=-ND; ND=(ND-1)*NINS+ND
	  ND=(ND-1)*NINS+ND
!	ELSE
!	  ND=J
!	END IF
!
!	CALL DERIVCHI(dFdR,CLUMP_FAC,R_OLD,ND_OLD,'LINMON')
!	F_CORR=dFdR*VSOUND*VSOUND/CLUMP_FAC
!	DO I=1,ND_OLD
!	  WRITE(25,'(I5,10ES14.4)')I,R_OLD(I),V_OLD(I),CLUMP_FAC_OLD(I),dFdR(I),F_CORR(I),
!	1        VSOUND*VSOUND/R_OLD(I)
!	END DO
!
	WRITE(6,'(A)')' '
	WRITE(STRING,'(I4)')ND_RD
	WRITE(6,'(A,T40,A)')' Number of depth points read is:',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(I4)')ND_OLD
	WRITE(6,'(A,T40,A)')' Number of depth points set is:',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(I4)')ND
	WRITE(6,'(A,T40,A)')' Sarting depth for integration is: ',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(ES12.6)')RSTAR
	WRITE(6,'(A,T40,A)')'     Stellar  radius is: ',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(ES12.6)')RSOUND
	WRITE(6,'(A,T40,A)')'     Critical radius is: ',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(ES12.6)')R_OLD(ND)
	WRITE(6,'(A,T40,A)')'   Connection radius is: ',TRIM(ADJUSTL(STRING))
	WRITE(STRING,'(ES12.6)')V_OLD(ND)
	WRITE(6,'(A,T40,A)')' Connection velocity is: ',TRIM(ADJUSTL(STRING))
!
! Change R to normalized coordinates.
!
	RHAT=R_OLD/RSTAR
	RSOUND=RSOUND/RSTAR
!
	VCRIT=1.0E-10_LDP*SQRT(GRAV_CONST*(MSTAR*MSUN)/RSTAR)
	WRITE(STRING,'(F10.2)')VCRIT
	WRITE(6,'(A,T40,A)')' VCRIT is:',TRIM(ADJUSTL(STRING))
	VCRIT=VCRIT/VSOUND
!
! We limit the minimum value GRAD to GRAD_MIN to avoid a deaccelerating flow.
!
        DO I=1,ND_OLD
          GRAD_MIN(I)=GRAD(I)
          IF(V_OLD(I) .GT. VSOUND)GRAD_MIN(I)=1.02E-20_LDP*GRAV_CONST*(MSTAR*MSUN)/R_OLD(I)**2
	END DO
!
	GELEC_OLD=GELEC_OLD*RSTAR/VSOUND/VSOUND
	GRAD_OLD=GRAD_OLD*RSTAR/VSOUND/VSOUND
	GRAD_MIN=GRAD_MIN*RSTAR/VSOUND/VSOUND
	F_CORR=0.0_LDP
	GELEC=GELEC_OLD
	GRAD=GRAD_OLD
!
	RX=RHAT(ND)
	VX=V_OLD(ND)/VSOUND
	OPEN(UNIT=16,FILE='LAMBERT_CHK',STATUS='UNKNOWN')
	WRITE(6,*)' '
	WRITE(6,*)'Beginning Lambert integration. '
	WRITE(6,*)'Diagnostic information written to LAMBERT_CHK'
	WRITE(6,*)' '
	WRITE(16,'(2X,A,10X,A,10X,A,7X,A,11X,A,13X,A,5(10X,A))')
	1            'Depth','RHAT','GINT','EXP_ARG','LNZ','Z','V_OLD','V_NEW','GMIN','GRAD','FCOR'
	WRITE(16,'(A)')' '
	DO IT=1,NUM_ITS
	  GINT=0.0_LDP
	  DO I=ND-1,1,-1
	    WRITE(50,*)I,GRAD(I)*RHAT(I)/VCRIT/VCRIT,GRAD(I+1)*RHAT(I+1)/VCRIT/VCRIT
	    GINT=GINT+0.5_LDP*(RHAT(I)-RHAT(I+1))*(MAX(GRAD_MIN(I),GRAD(I)+F_CORR(I))+
	1                                       MAX(GRAD_MIN(I+1),GRAD(I+1)+F_CORR(I+1)))
	    Z=-(RX/RHAT(I))**4
	    T1=2*LOG(VX)-VX*VX-2*VCRIT*VCRIT*(1.0_LDP/RHAT(I)-1.0_LDP/RX)-2.0_LDP*GINT
	    LNZ=LOG(-Z)+T1
	    Z=Z*EXP(T1)
	    T2=LAMBERT_WM_FUN(Z,LNZ)
	    V_NEW(I)=SQRT(-T2)*VSOUND
	    WRITE(16,'(I7,10ES14.4)')I,RHAT(I),GINT,T1,LNZ,Z,V_OLD(I),V_NEW(I),GRAD_MIN(I),GRAD(I),F_CORR(I)
	  END DO
	  V_NEW(ND:ND_OLD)=V_OLD(ND:ND_OLD)
	  K=ND/5
	  IF(IT .EQ. 1)THEN
	    WRITE(6,'(A,I3,6ES14.4)')' Iteration: ',0,(V_OLD(I),I=1,ND,K)
	  END IF
	  WRITE(6,'(A,I3,6ES14.4)')' Iteration: ',IT,(V_NEW(I),I=1,ND,K)
!
	  IF(AVERAGE_V)THEN
	    V_NEW=0.5_LDP*(V_NEW+V_OLD)
	  END IF
!
	  IF(F_FUN_V)THEN
	    CALL MON_INTERP(F_CORR,ND_OLD,IONE,V_NEW,ND_OLD,CLUMP_FAC_OLD,ND_OLD,V_OLD,ND_OLD)
	    CLUMP_FAC=0.5_LDP*(CLUMP_FAC_OLD+F_CORR)
	    CALL DERIVCHI(dFdR,CLUMP_FAC,R_OLD,ND_OLD,'LINMON')
	    F_CORR=dFdR*VSOUND*VSOUND/CLUMP_FAC
	    F_CORR=F_CORR*RSTAR/VSOUND/VSOUND
	  END IF
!
! When we iterate, we can correct for the fact that G(rad) depends on the velocity
! gradient.
!
          CALL MON_INT_FUNS_V2(COEF,V_NEW,R_OLD,ND_OLD)
          DO I=1,ND_OLD
            SIGMA_NEW(I)=COEF(I,3)
	    WRITE(33,'(I5,2ES14.4)')I,SIGMA_NEW(I),VdVdR_OLD(I)
	    T1=V_NEW(I)*SIGMA_NEW(I)
	    T2=VdVdR_OLD(I)
	    GRAD(I)=GELEC(I)+(GRAD_OLD(I)-GELEC(I))*(T1/T2)**ALPHA
	    WRITE(77,*)T2,T1,GRAD_OLD(I),GRAD(I)
            SIGMA_NEW(I)=R_OLD(I)*SIGMA_NEW(I)/V_NEW(I)-1.0_LDP
	  END DO
!
	END DO
!
	INQUIRE(FILE='RVSIG_COL',EXIST=FILE_PRES)
	IF(FILE_PRES)THEN
	  WRITE(STRING,*)MAIN_COUNTER
	  STRING='cp RVSIG_COL RVSIG_COL_'//ADJUSTL(STRING)
	  CALL SYSTEM(STRING)
	END IF
!
	WRITE(6,*)'New RV grid being written to RVSIG_COL'
	WRITE(6,*)' '
        OPEN(UNIT=10,FILE='RVSIG_COL',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(10,'(A)')'!'
        WRITE(10,'(A,7X,A,9X,10X,A,11X,A,3X,A)')'!','R','V(km/s)','Sigma','Depth'
        WRITE(10,'(A)')'!'
        WRITE(10,'(A)')' '
        WRITE(10,'(I4,20X,A)')ND_RD,'!Number of depth points`'
        WRITE(10,'(A)')' '
	J=0
        DO I=1,ND_OLD,NINS+1
          WRITE(10,'(F18.8,ES20.10,F17.7,4X,I4)')R_OLD(I),V_NEW(I),SIGMA_NEW(I),(I-1)/(NINS+1)+1
	  J=J+1
	  R_RD(J)=R_OLD(I)
	  V_RD(J)=V_NEW(I)
	  SIGMA_RD(J)=SIGMA_NEW(I)
        END DO
        CLOSE(UNIT=10)
	R_RD(ND_RD/2)=R_RD(ND_RD/2)*(1.0_LDP+1.0E-12_LDP)
!
	DONE_V_REVISION=.TRUE.
	DEALLOCATE (COEF)
	LST_V_UPDATE=MAIN_COUNTER
	NUM_V_ITS=NUM_V_ITS+1
!	CALL UPDATE_KEYWORD(V_RD(1),'[VINF]','VADAT',L_TRUE,L_TRUE,LUIN)
!
	RETURN
	END
