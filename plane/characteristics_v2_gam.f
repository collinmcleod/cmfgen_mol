!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE CHARACTERISTICS_V2_GAM(R_GRID,P,V_GRID,VDOP_VEC,VDOP_FRAC,ND,NC,NP)
	USE SET_KIND_MODULE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This subroutine uses the Runge-Kutta method (solves sets of ordinary
! differential equations) to determine the characteristic rays needed
! to make the comoving frame radiative transfer equation a perfect
! differential. The "rk4" routine used below is from Numerical Recipes
! (92), chapter 16.
!
! Written    3/95 DLM
! Altered    8/95 DLM Inserted into transfer program and changed to calculate
!                     on the given radial grid.  Renamed characteristics from
!                     runge.
! Altered  2/21/96 DLM Updated to F90 standard
! Altered  4/8/97  DLM Put path length and angle variables in module.  Added
!                     use MOD_SPACE_GRID and allocation of variables
! Altered  5/22/97 DLM Added calculation of b_p and b_m for relativistic
!                     transfer terms
! Altered 27/12/04 DJH Cleaned and new velocity law routine inserted.
! Altered 26/11/08 DJH Changed to version 2.
!                  Altered to allow creation of finer grid along each ray.
!
!--------------------------------------------------------------------
!
      USE MOD_SPACE_GRID_V2_GAM
      IMPLICIT NONE
!
! Grid size variables: passed to routine
!
      INTEGER :: ND,NC,NP
!
! Grid variables
!
      REAL(KIND=LDP), DIMENSION(ND) :: R_GRID
      REAL(KIND=LDP), DIMENSION(NP) :: P
      REAL(KIND=LDP), DIMENSION(ND) :: V_GRID
      REAL(KIND=LDP), DIMENSION(ND) :: VDOP_VEC(ND)
      REAL(KIND=LDP)  VDOP_FRAC
!
! 'ns' is number of steps to take between each grid point when
! integrating wrt r.  'mu_step' is mu step to use when integrating wrt s.
! It is used to determine ds.
!
      INTEGER, PARAMETER :: NS=20
      INTEGER, PARAMETER :: IONE=1
!
      REAL(KIND=LDP), PARAMETER ::  MU_STEP=0.005_LDP
!
! Number of equations in Runge-Kutta solution
!
      INTEGER, PARAMETER :: EQU=2
!
! Derivative functions for r and s
!
      EXTERNAL DERIVS,DERIVR
!
! Loop variables
!
      INTEGER :: I,J,IP,ID,JD,IS
      INTEGER :: IDM1
      INTEGER :: NP_LIMIT
!
! Determine direction of characteristic ray (+1,-1)
!
      REAL(KIND=LDP) :: DIRECTION
!
! Radius variables
!
      REAL(KIND=LDP) dZ,dV,dR
      REAL(KIND=LDP) :: RMIN,RMAX,DS,SS,RR,VEL,A_OLD,DR_OLD,DS_1,DR_1
      REAL(KIND=LDP) :: MU_E,MU_H
!
! Velocity variables
!
      REAL(KIND=LDP) :: BETAMAX,BETA,DBETADR,GAMMA
!
! Runge-Kutta variables
!
      REAL(KIND=LDP), DIMENSION(EQU) :: A,DADS,DADR
!
	REAL(KIND=LDP), ALLOCATABLE :: Z(:)
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: S(:)
	REAL(KIND=LDP), ALLOCATABLE :: MU(:)
	REAL(KIND=LDP), ALLOCATABLE :: B(:)
!
	INTEGER NINS
	INTEGER NRAY
	INTEGER NRAY_MAX
	INTEGER NRAY_SM
!
	INTEGER ERROR_LU
	INTEGER LUER
	INTEGER IOS
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
!
	WRITE(LUER,*)'   Inside CHARACTERISTICS_V2_GAM'
	DO I=1,ND
	  WRITE(101,'(I3,2X,3ES14.4)')I,R_GRID(I),V_GRID(I),VDOP_VEC(I)
	END DO
	WRITE(101,*)'P '
	DO I=1,NP
	  WRITE(101,'(I3,1X,ES14.4)')I,P(I)
	END DO
!
	OPEN(UNIT=88,FILE='ray_check_for_grays.dat',STATUS='UNKNOWN',ACTION='WRITE',
	1	IOSTAT=IOS)
!
! KDW
      WRITE(88,'(3(A3,1X,I3,3X))')'NP:',NP,'ND:',ND,'NC:',NC
!
!
! If we are defining a new grid, we need to deallocate variables defined
! along each ray.
!
	IF(ALLOCATED(RAY(1)%S_P))THEN
	  DO IP=1,NP
	    DEALLOCATE (RAY(IP)%S_P, RAY(IP)%MU_P, RAY(IP)%B_P, RAY(IP)%R_RAY)
	    DEALLOCATE (RAY(IP)%S_M, RAY(IP)%MU_M, RAY(IP)%B_M)
	  END DO
	END IF
!
! Allocate temporary storage used for each ray. This memory is deallocated at the
! end of the subroutine.
!
	NRAY_MAX=100*ND
	ALLOCATE (Z(NRAY_MAX),STAT=IOS)
!	DO I=1,NP
!	  ALLOCATE (RAY(I)%Z_SM(ND))
!	  ALLOCATE (RAY(I)%Z_LG(NRAY_MAX))
!	END DO
	IF(IOS .EQ. 0)ALLOCATE (R(NRAY_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(NRAY_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (S(NRAY_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (MU(NRAY_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (B(NRAY_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Unable to allocate R, V, S etc in CHARACTERISTICS_V2'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Allocate path length and direction variables. P stands for +S (+MU)
! direction,  M for -ve S (mu) direction.
!
!
! Allocate variables for frequency independent parts of advection and
! abberation terms in relativistic transfer equation.  They will be used
! in solve_rel_formal.
!
!  b_*=gamma*((1-mu_*^2)*beta/r+gamma^2*mu_**(mu_*+beta)*dbetadr)
!
!      where _* is _p or _m for ray in plus or minus direction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Integrate in the positive direction during the first pass and in the
! negative direction on the second pass.  This is controlled by the
! variable 'direction'.  Set to +1 now, and will set to -1 at the end
! of the ip do loop
!
      DIRECTION=1.0_LDP
      RAY_POINTS_INSERTED=.FALSE.       !Will be reset, if necessary
!
!--------------------------------------------------------------------
!----------Loop over p-rays------------------------------------------
!--------------------------------------------------------------------
!
! Integrate with s as the independent variable till first depth point
! is reached
!
! Must first calculate along central ray, mu=1, dmu/dr=0.
! Because of the relativistic terms, the path length varies.
!
 100  CONTINUE
!
!      IF(.NOT. ALLOCATED(RAY(1)%NINS))ALLOCATE (RAY(1)%NINS(ND))
!      RAY(1)%NINS=0
!
      DO I=1,ND
        Z(I)=R_GRID(I)
!	RAY(1)%Z_SM(I)=R_GRID(I)
      END DO
      R(1)=R_GRID(1)
!      RAY(1)%Z_LG(1)=RAY(1)%Z_SM(1)
      ID=1
      DO I=1,ND-1
        dV=ABS(V_GRID(I)-V_GRID(I+1))
	IF(dV .GT. 1.25_LDP*VDOP_FRAC*VDOP_VEC(I))THEN
          RAY_POINTS_INSERTED=.TRUE.
	  WRITE(LUER,'(A,1X,I1)') "INSERTING RAY POINTS for ray:",1
	  NINS=dV/(VDOP_FRAC*VDOP_VEC(I))
!	  RAY(1)%NINS(I)=NINS
          dR=(R_GRID(I)-R_GRID(I+1))/(NINS+1)
          IF(ID+NINS .GT. NRAY_MAX)THEN
	    WRITE(LUER,*)'Error in CHRACTERISTICS_V2_GAM'
	    WRITE(LUER,*)'NRAY_MAX (',NRAY_MAX,') is too small'
	    STOP
	   END IF
          DO J=1,NINS
            ID=ID+1
            R(ID)=R(ID-1)-dR
!	    RAY(1)%Z_LG(ID)=R(ID-1)-dR
          END DO
        END IF
        ID=ID+1
        IF(ID .GT. NRAY_MAX)THEN
	  WRITE(LUER,*)'Error in CHRACTERISTICS_V2'
	  WRITE(LUER,*)'NRAY_MAX (',NRAY_MAX,') is too small'
	  STOP
	 END IF
	 R(ID)=R_GRID(I+1)
!	 RAY(1)%Z_LG(ID)=R_GRID(I+1)
      END DO
      NRAY=ID
      CALL MON_INTERP(V,NRAY,IONE,R,NRAY,V_GRID,ND,R_GRID,ND)
!
      IP=1
      IDM1=NRAY-1
      ID=NRAY
      RAY(IP)%NZ=NRAY
!
! KDW
      IF(IP .GT. NC)THEN
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'NP-IP+1:',NP-IP+1
      ELSE
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A3,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'ND:',ND
      END IF
!
!
      RR=R(ND)
      CALL VELOCITY_LAW(RR,IDM1,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
      S(ID)=DIRECTION*R(ID)
      MU(ID)=DIRECTION
      B(ID)=GAMMA*((1.0_LDP-MU(ID)**2)*BETA/R(ID)+GAMMA**2*MU(ID)*(MU(ID)+BETA)*DBETADR)
!
      A(1)=0.0_LDP
      A(2)=DIRECTION
      SS=DIRECTION*R(ID)
      DO JD=IDM1,1,-1
        RR=R(JD+1)
        DR=(R(JD)-R(JD+1))/10.0_LDP
        DO IS=1,10
          CALL VELOCITY_LAW(RR,JD,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
          CALL DERIVR(RR,A,DADR,BETA,DBETADR,GAMMA)
          SS=SS+DADR(1)*DR
          RR=RR+DR
        END DO
!
! Save s and mu for this grid point
!
        CALL VELOCITY_LAW(RR,JD,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
        S(JD)=SS
        MU(JD)=DIRECTION
        B(JD)=GAMMA*((1.0_LDP-MU(JD)**2)*BETA/R(JD)+
     *          GAMMA**2*MU(JD)*(MU(JD)+BETA)*DBETADR)
      END DO
!
      IF(DIRECTION .GT. 0)THEN
	ALLOCATE (RAY(IP)%S_P(NRAY),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%MU_P(NRAY),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%B_P(NRAY),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%R_RAY(NRAY),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in characteristics.f'
	  WRITE(LUER,*)'Unable to allocate memory for P with IP=',1
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
        DO I=1,NRAY
          RAY(IP)%R_RAY(I)=R(I)
          RAY(IP)%S_P(I)=S(I)
          RAY(IP)%MU_P(I)=MU(I)
          RAY(IP)%B_P(I)=B(I)
        END DO
      ELSE
	ALLOCATE (RAY(IP)%S_M(NRAY),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%MU_M(NRAY),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%B_M(NRAY),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in characteristics.f'
	  WRITE(LUER,*)'Unable to allocate memory for P with IP=',1
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
        DO I=1,NRAY
          RAY(IP)%S_M(I)=S(I)
          RAY(IP)%MU_M(I)=MU(I)
          RAY(IP)%B_M(I)=B(I)
        END DO
      END IF
!
! Calculate the rest of the rays
!
      NP_LIMIT=NP-1
      IF(ND+NC .GT. NP)NP_LIMIT=NP
      DO IP=2,NP_LIMIT
!
! Define a the fine grid such that dV < VDOP_FRAC*VTURB along ray.
! For simplicity we ignor aberration terms in defining the new grid.
!
        NRAY_SM=ND
        IF(IP .GT. NC)NRAY_SM=ND-(IP-NC-1)
!	IF(.NOT. ALLOCATED(RAY(IP)%NINS))ALLOCATE (RAY(IP)%NINS(NRAY_SM))
!	RAY(IP)%NINS=0
        DO I=1,NRAY_SM
          Z(I)=SQRT( (R_GRID(I)-P(IP))*(R_GRID(I)+P(IP)) )
!          RAY(IP)%Z_SM(I)=SQRT((R_GRID(I)-P(IP))*(R_GRID(I)+P(IP)))
	END DO
        R(1)=R_GRID(1)
!	RAY(IP)%Z_LG(1)=RAY(IP)%Z_SM(1)
        ID=1
        DO I=1,NRAY_SM-1
          dV=Z(I)*V_GRID(I)/R_GRID(I)-Z(I+1)*V_GRID(I+1)/R_GRID(I+1)
          IF(dV .GT. 1.25_LDP*VDOP_FRAC*VDOP_VEC(I))THEN
            RAY_POINTS_INSERTED=.TRUE.
	    WRITE(LUER,'(A,1X,I3)') "INSERTING RAY POINTS for ray:",IP
            NINS=dV/(VDOP_FRAC*VDOP_VEC(I))
!	    RAY(IP)%NINS(I)=NINS
            dZ=(Z(I)-Z(I+1))/(NINS+1)
            IF(ID .GT. NRAY_MAX)THEN
	      WRITE(LUER,*)'Error in CHRACTERISTICS_V2'
	      WRITE(LUER,*)'NRAY_MAX (',NRAY_MAX,') is too small'
	      STOP
	    END IF
            DO J=1,NINS
              ID=ID+1
              R(ID)=SQRT( P(IP)*P(IP)+(Z(I)-J*dZ)**2 )
!	      RAY(IP)%Z_LG(ID)=RAY(IP)%Z_LG(ID-1)-dZ
            END DO
          END IF
          ID=ID+1
          IF(ID .GT. NRAY_MAX)THEN
	    WRITE(LUER,*)'Error in CHRACTERISTICS_V2'
	    WRITE(LUER,*)'NRAY_MAX (',NRAY_MAX,') is too small'
	    STOP
	  END IF
          R(ID)=R_GRID(I+1)
!	  RAY(IP)%Z_LG(ID)=RAY(IP)%Z_SM(I+1)
        END DO
	NRAY=ID
        CALL MON_INTERP(V,NRAY,IONE,R,NRAY,V_GRID,ND,R_GRID,ND)
!
! If core ray then integrate WRT s till r(nd).  If non-core ray
! then integrate WRT s till r(ni)
!
	IF(IP .LE. NC)THEN
	  ID=NRAY
          RAY(IP)%NZ=NRAY
!
! KDW
      IF(IP .GT. NC)THEN
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'NP-IP+1:',NP-IP+1
      ELSE
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A3,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'ND:',ND
      END IF
!
!
	  DS_1=DIRECTION*R(ID)/FLOAT(NS)
	ELSE
          ID=NRAY-1
          RR=R(ID+1)
          CALL VELOCITY_LAW(RR,ID,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
          RAY(IP)%NZ=NRAY
!
! KDW
      IF(IP .GT. NC)THEN
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'NP-IP+1:',NP-IP+1
      ELSE
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A3,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'ND:',ND
      END IF
!
!
          S(ID+1)=0.0_LDP
          MU(ID+1)=-BETA
          B(ID+1)=GAMMA*((1.0_LDP-MU(ID+1)**2)*BETA/R(ID+1)+
     *       GAMMA**2*MU(ID+1)*(MU(ID+1)+BETA)*DBETADR)
          DS_1=DIRECTION*(R(ID)-R(ID+1))/FLOAT(NS)
	END IF
!
! Determine ds by using dmu*ds/dmu at mu=-beta
!
        DADS(2)=(1.0_LDP-BETA*BETA)**1.5_LDP/P(IP)
        DS=DIRECTION*MU_STEP/DADS(2)
        IF(ABS(DS_1).LT.ABS(DS))DS=DS_1
!
        SS=0.0_LDP
!
! Boundary conditions on r and dr/ds
!
        A(1)=P(IP)
        DADS(1)=0.0_LDP
!
! Boundary conditions on mu and dmu/ds
!
        CALL VELOCITY_LAW(A(1),MIN(ID,NRAY-1),R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
        A(2)=-BETA
        DADS(2)=(1.0_LDP-BETA*BETA)**1.5_LDP/P(IP)
!
! Storage of change in r [a(1)_new-a(2)_new)]
!
        DR_OLD=0.0_LDP
        A_OLD=A(1)
!
! Integrate with s as independent variable.  Must do several
! steps away from s=0 in order to find dmu/dr.
!
        DO WHILE((R(ID)-A(1)-DR_OLD).GT.(DR_OLD*1.0E-3_LDP))
          CALL RUNGE_KUTTA(A,DADS,EQU,SS,DS,BETA,DBETADR,GAMMA,DERIVS)
          DR_OLD=A(1)-A_OLD
          A_OLD=A(1)
          CALL VELOCITY_LAW(A(1),MIN(ID,NRAY-1),R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
          CALL DERIVS(SS,A,DADS,BETA,DBETADR,GAMMA)
          DS=DIRECTION*MU_STEP/DADS(2)
          IF(ABS(DS_1).LT.ABS(DS))DS=DS_1
        END DO
!
! Put initial conditions into vectors for passage to subroutine
!
        RR=A(1)
        A(1)=SS
        DADR(1)=1.0_LDP/DADS(1)
        DADR(2)=DADS(2)/DADS(1)
!
! Choose increment dr to finish integration to r(id).
!
        DR=R(ID)-RR
!
! Integrate with r as the independent variable
!
        CALL RUNGE_KUTTA(A,DADR,EQU,RR,DR,BETA,DBETADR,GAMMA,DERIVR)
        CALL VELOCITY_LAW(RR,MIN(ID,NRAY-1),R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
        CALL DERIVR(RR,A,DADR,BETA,DBETADR,GAMMA)
        S(ID)=A(1)
        MU(ID)=A(2)
        B(ID)=GAMMA*((1.0_LDP-MU(ID)**2)*BETA/R(ID)+
     *       GAMMA**2*MU(ID)*(MU(ID)+BETA)*DBETADR)
!
! Integrate along the p-ray to determine s and mu at each grid point r.
! Thus integrate with r as the independent variable.  For a ray there will
! be nz(ip)-1 such increments
!
        DR_OLD=0.0_LDP
        A_OLD=R(ID)
!
        DO JD=ID-1,1,-1
!
! Integrate in increments between r(id+1) and r(id)
!
          DR_1=(R(JD)-R(JD+1))/FLOAT(NS)
          DR=ABS(MU_STEP/DADR(2))
          IF(DR_1.LT.DR)DR=DR_1
!
          DO WHILE((R(JD)-RR-DR_OLD).GT.(DR_OLD*1.0E-3_LDP))
            CALL RUNGE_KUTTA(A,DADR,EQU,RR,DR,BETA,DBETADR,GAMMA,DERIVR)
            DR_OLD=RR-A_OLD
            A_OLD=RR
            CALL VELOCITY_LAW(RR,JD,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
            CALL DERIVR(RR,A,DADR,BETA,DBETADR,GAMMA)
            DR=ABS(MU_STEP/DADR(2)/4.0_LDP)
            IF(DR_1.LT.DR)DR=DR_1
          END DO
!
! Choose increment dr to finish integration to r(id).
!
          DR=R(JD)-RR
          CALL RUNGE_KUTTA(A,DADR,EQU,RR,DR,BETA,DBETADR,GAMMA,DERIVR)
          RR=R(JD)
          CALL VELOCITY_LAW(RR,JD,R,V,NRAY,VEL,BETA,DBETADR,GAMMA)
          CALL DERIVR(RR,A,DADR,BETA,DBETADR,GAMMA)
!
! Save s and mu for this grid point
!
          S(JD)=A(1)
          MU(JD)=A(2)
          B(JD)=GAMMA*((1.0_LDP-MU(JD)**2)*BETA/R(JD)+
     *          GAMMA**2*MU(JD)*(MU(JD)+BETA)*DBETADR)
!
        END DO
!
	IF(DIRECTION .GT. 0)THEN
	  ALLOCATE (RAY(IP)%S_P(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%R_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%MU_P(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%B_P(NRAY),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in characteristics.f'
	    WRITE(LUER,*)'Unable to allocate memory for P with IP=',IP
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
          DO I=1,NRAY
	    RAY(IP)%S_P(I)=S(I)
	    RAY(IP)%MU_P(I)=MU(I)
	    RAY(IP)%B_P(I)=B(I)
	    RAY(IP)%R_RAY(I)=R(I)
	  END DO
	ELSE
	  ALLOCATE (RAY(IP)%S_M(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%MU_M(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(IP)%B_M(NRAY),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in characteristics.f'
	    WRITE(LUER,*)'Unable to allocate memory for M with IP=',IP
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
          DO I=1,NRAY
	    RAY(IP)%S_M(I)=S(I)
	    RAY(IP)%MU_M(I)=MU(I)
	    RAY(IP)%B_M(I)=B(I)
	  END DO
        END IF
      END DO		!Loop over NP
!
      IF(NP .EQ. NC+ND)THEN
	RAY(NP)%NZ=1
        NRAY=1
!
! KDW
      IF(IP .GT. NC)THEN
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'NP-IP+1:',NP-IP+1
      ELSE
        WRITE(88,'(A3,1X,I3,3X,A5,1X,I3,4X,A3,1X,I3)')'IP:',IP,'NRAY:',
	1	RAY(IP)%NZ,'ND:',ND
      END IF
!
!
        IF(DIRECTION .GT. 0)THEN
	  ALLOCATE (RAY(NP)%S_P(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(NP)%R_RAY(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(NP)%MU_P(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(NP)%B_P(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in characteristics.f'
	    WRITE(LUER,*)'Unable to allocate memory for P with IP=',NP
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
         RAY(NP)%S_P(1)=0.0_LDP
         RAY(NP)%MU_P(1)=-BETA
         RAY(NP)%B_P(1)=GAMMA*((1.0_LDP-RAY(NP)%MU_P(1)**2)*BETA/R(1)+
     *       GAMMA**2*RAY(NP)%MU_P(1)*(RAY(NP)%MU_P(1)+BETA)*DBETADR)
         RAY(NP)%R_RAY(1)=R(1)
       ELSE
	  ALLOCATE (RAY(NP)%S_M(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(NP)%MU_M(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RAY(NP)%B_M(NRAY),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in characteristics.f'
	    WRITE(LUER,*)'Unable to allocate memory for M with IP=',NP
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
          RAY(NP)%S_M(1)=0.0_LDP
          RAY(NP)%MU_M(1)=-BETA
          RAY(NP)%B_M(1)=GAMMA*((1.0_LDP-RAY(NP)%MU_M(1)**2)*BETA/R(1)+
     *       GAMMA**2*RAY(NP)%MU_M(1)*(RAY(NP)%MU_M(1)+BETA)*DBETADR)
        END IF
      END IF
!
!--------------------------------------------------------------------
!----------End loop over core p-rays---------------------------------
!--------------------------------------------------------------------
!
! Do loop a second time for s<0
!
      IF(DIRECTION .GT. 0.0_LDP)THEN
        DIRECTION=-1.0_LDP
        GOTO 100
      END IF
!
! Free up memory.
!
      DEALLOCATE (R,V,Z,S,MU,B)
!
      CLOSE(UNIT=88)
      OPEN(UNIT=7,FILE='nray_info.dat',STATUS='UNKNOWN',ACTION='WRITE')
      WRITE(7,'(A4,2X,A4)')'IP','NRAY'
      DO I=1,NP
	WRITE(7,'(I4,2X,I4)')I,RAY(I)%NZ
      END DO
      CLOSE(UNIT=7)
      RETURN
      END
