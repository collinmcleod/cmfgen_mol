!
! Program to modify RVSIG_COL. Various options are available.
! Ideal for revising grid etc.
!
	PROGRAM REVISE_RVSIG
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NMAX=200
	INTEGER, PARAMETER :: IONE=1
!
	REAL*8 OLD_R(NMAX)
	REAL*8 OLD_V(NMAX)
	REAL*8 OLD_SIGMA(NMAX)
!
	REAL*8 R(NMAX)
	REAL*8 V(NMAX)
	REAL*8 SIGMA(NMAX)
!
	REAL*8 TMP_R(NMAX)
	REAL*8 X1(NMAX)
	REAL*8 X2(NMAX)	
!
	REAL*8, ALLOCATABLE :: COEF(:,:)
!
	INTEGER ND_OLD
	INTEGER ND
!
	REAL*8 RX
	REAL*8 T1,T2
	REAL*8 BETA
	REAL*8 VINF
	REAL*8 FAC
	REAL*8 V_MAX
	REAL*8 V_MIN 
!
	REAL*8 R_TRANS
	REAL*8 V_TRANS
	REAL*8 dVdR_TRANS
	REAL*8 SCALE_HEIGHT
	REAL*8 RO
	REAL*8 dVdR
	REAL*8 MDOT
	REAL*8 OLD_MDOT
	REAL*8 TOP,BOT 
	REAL*8 dTOPdR,dBOTdR 
	INTEGER TRANS_I
	INTEGER VEL_TYPE
!
	INTEGER NN
	INTEGER NX
	INTEGER N_ADD
	INTEGER NX_IN
	INTEGER NX_OUT
!
	INTEGER I,J,K
	INTEGER I_ST,I_END
!
	CHARACTER(LEN=10) OPTION
	CHARACTER(LEN=80) OLD_RVSIG_FILE
	CHARACTER(LEN=80) NEW_RVSIG_FILE
	CHARACTER(LEN=80) STRING

	OLD_RVSIG_FILE='RVSIG_COL_OLD'
	CALL GEN_IN(OLD_RVSIG_FILE,'File containing old R, V and sigma values')
	OPEN(UNIT=10,FILE=OLD_RVSIG_FILE,STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE (INDEX(STRING,'!Number of depth points') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)ND_OLD
	  DO I=1,ND_OLD
	    READ(10,*)OLD_R(I),OLD_V(I),OLD_SIGMA(I)
	  END DO
	CLOSE(UNIT=10)
!
	OPTION='NEW_ND'
	WRITE(6,'(A)')
	WRITE(6,'(A)')'Current options are:'
	WRITE(6,'(A)')'       NEW_ND: revise number of depth points'
	WRITE(6,'(A)')'         ADDR: add extra grid points'
	WRITE(6,'(A)')'         EXTR: extend grid to larger radii'
	WRITE(6,'(A)')'         MDOT: change mass-loss rate'
	WRITE(6,'(A)')'         NEWG: revise grid'
	WRITE(6,'(A)')

	CALL GEN_IN(OPTION,'Enter option for revised RVSIG file')
	IF(OPTION .EQ. 'NEW_ND')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows a new R grid to be output'
	  WRITE(6,'(A)')'Grid spacing is similar to input grid.'
	  WRITE(6,'(A)')' '
	  ND=70
	  CALL GEN_IN(ND,'Number of depth points')
	  DO I=1,ND_OLD
	    X1(I)=I
	  END DO
	  T1=DFLOAT(ND_OLD-1)/DFLOAT(ND-1)
	  DO I=1,ND
	    X2(I)=(I-1)*T1+1
	  END DO
	  X2(1)=X1(1)
	  X2(ND)=X1(ND_OLD)
	  
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
!
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0D0*COEF(J,2)+3.0*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'ADDR')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows extra grid points to be added'
	  WRITE(6,'(A)')'The grid points are specfied by the user as requested'
	  WRITE(6,'(A)')'Please insert the grid points largest to smallest'
	  N_ADD=1
	  CALL GEN_IN(N_ADD,'Number of grid points to be added')
	  DO I=1,N_ADD
	    CALL GEN_IN(X1(I),'New grid point:')
	  END DO
	  DO I=1,N_ADD-1
	    IF(X1(I) .LE. X1(I+1))THEN
	      WRITE(6,*)'Error: new grid points are invalid -- equal or wrong order'
	      STOP
	    END IF
	  END DO
	  IF(X1(1) .GE. OLD_R(1))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside defined grid'
	    WRITE(6,*)'OLD_R(1)=',OLD_R(1)
	    WRITE(6,*)'R_INS=',X1(1)
	    STOP
	  END IF
	  IF(X1(N_ADD) .LE. OLD_R(ND_OLD))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside defined grid'
	    WRITE(6,*)'OLD_R(ND_OLD)=',OLD_R(ND_OLD)
	    WRITE(6,*)'R_INS=',X1(N_ADD)
	    STOP
	  END IF
!
	  ND=N_ADD+ND_OLD
	  J=1
	  I=2
	  K=1
	  R(J)=OLD_R(1)
	  DO WHILE (K .LE. N_ADD)
	    IF(X1(K) .LE. R(J) .AND. X1(K) .GT. OLD_R(I))THEN
	      J=J+1
	      R(J)=X1(K)
	      K=K+1
	    ELSE
	      J=J+1
	      R(J)=OLD_R(I)
	      I=I+1
	    END IF
	  END DO
	  DO K=J+1,ND
	    R(K)=OLD_R(I)
	    I=I+1
	  END DO
!
! Check ordering etc.
!
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(6,*)'Error: R values not monotonic'
	    WRITE(6,*)'I=',I
	    WRITE(6,*)'R(I)=',R(I)
	    WRITE(6,*)'R(I+1)=',R(I+1)
	    STOP
	  END IF
	END DO
!
! Now compute the revised V and SIGMA.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0D0*COEF(J,2)+3.0*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'NEWG')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows ithe grid to be redeifined'
!
	  CALL GEN_IN(V_MAX,'Maximum velocity for grid refinement')
	  CALL GEN_IN(V_MIN,'Initial velocity for grid refinement')
	  IF(V_MAX .GE. OLD_V(1) .OR. V_MAX .LE. OLD_V(ND_OLD))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	    STOP
	  END IF
	  IF(V_MIN .GE. OLD_V(1) .OR. V_MIN .LE. OLD_V(ND_OLD))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	    STOP
	  END IF
	  CALL GEN_IN(N_ADD,'Number of points in interval (exclusive)')
!
! Find interval
!
	  DO I=1,ND_OLD-1
	    IF(V_MAX .GT. OLD_V(I+1))THEN
	      I_ST=I
	      EXIT
	    END IF
	  END DO
	  IF(V_MAX-OLD_V(I+1) .LT. OLD_V(I)-V_MAX)I_ST=I+1
!
	  DO I=1,ND_OLD-1
	    IF(V_MIN .GT. OLD_V(I+1))THEN
	      I_END=I
	      EXIT
	    END IF
	  END DO
	  IF(V_MIN-OLD_V(I+1) .LT. OLD_V(I)-V_MIN)I_END=I+1
	  WRITE(6,*)' I_ST=',I_ST,OLD_V(I_ST)
	  WRITE(6,*)'I_END=',I_END,OLD_V(I_END)
!
	  ND=N_ADD+ND_OLD-(I_END-I_ST-1)
	  V(1:I_ST)=OLD_V(1:I_ST)
	  T1=EXP(DLOG(V_MAX/V_MIN)/(N_ADD+1))
	  DO I=I_ST+1,I_ST+N_ADD
	   V(I)=V(I-1)/T1
	  END DO
	  V(I_ST+N_ADD+1:ND)=OLD_V(I_END:ND_OLD)
!
	  CALL MON_INTERP(R,ND,IONE,V,ND,OLD_R,ND_OLD,OLD_V,ND_OLD)
!
! Now compute the revised SIGMA. V has already been computed.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
!	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0D0*COEF(J,2)+3.0*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
	  
	ELSE IF(OPTION .EQ. 'EXTR')THEN
	  FAC=2.0D0
	  CALL GEN_IN(FAC,'Factor to extend RMAX by')
          NX_OUT=2
	  VINF=1000.0D0
	  CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	  BETA=1.0D0
	  CALL GEN_IN(BETA,'Beta for velocity law')
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Outputing old grid near outer boudary for NX estimate.'
	  WRITE(6,*)'NB: The grid ratio will differ for NX+1 values.'
	  WRITE(6,*)' '
	  DO I=1,7
	    WRITE(6,*)I,OLD_R(I)/OLD_R(I+1)
	  END DO
	  NX_OUT=2
	  CALL GEN_IN(NX_OUT,'Number of points used at outer boundary to refine grid')
!
	  RX=OLD_R(1)*(1.0D0-(OLD_V(1)/VINF)**(1.0D0/BETA))
	  WRITE(6,*)'RX=',RX
!
! Set up a rough grid so we can define a better grid, equally spaced
! in log density.
!
	  TMP_R(1)=FAC*OLD_R(1)
	  NX=20 
	  T1=EXP(LOG(FAC)/NX)
	  TMP_R(1)=FAC*OLD_R(1)
	  DO I=2,NX
	    TMP_R(I)=TMP_R(I-1)/T1
	  END DO
	  TMP_R(NX+1)=OLD_R(1)
	  DO I=1,NX+1
	    V(I)=VINF*(1.0D0-RX/TMP_R(I))**BETA
	    X1(I)=TMP_R(I)*TMP_R(I)*V(I)
	  END DO
	  NX=NX+1
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Outputing old V.r^2 near outer boudary for GRID estimate.'
	  WRITE(6,*)' '
	  DO I=1,7
	    WRITE(6,*)I,(OLD_V(I)/OLD_V(I+1))*(OLD_R(I)/OLD_R(I+1))**2
	  END DO
	  FAC=SQRT(OLD_R(5)*OLD_R(5)*OLD_V(5)/OLD_R(7)/OLD_R(7)/OLD_V(7))
	  CALL GEN_IN(FAC,'Grid spacing ratio')
	  T1=LOG(X1(1)/X1(NX))
	  NN=NINT(T1/LOG(FAC))
	  T1=EXP(T1/NN)
	  WRITE(6,*)'Number of points to be add is',NN
	  WRITE(6,*)'Grid spacing factor is ',T1
!
	  X2(1)=X1(1)
	  DO I=2,NN
	    X2(I)=X2(I-1)/T1
	  END DO
	  CALL MON_INTERP(R,NN,IONE,X2,NN,TMP_R,NX,X1,NX)
!
	  DO I=NN,2,-1
	    R(I+NX_OUT)=R(I)
	  END DO
	  IF(NX_OUT .EQ. 1)THEN
	    R(2)=R(1)-0.2*(R(1)-R(3))
	  ELSE IF(NX_OUT .EQ. 2)THEN
	    R(2)=R(1)-0.1*(R(1)-R(4))
	    R(3)=R(1)-0.4*(R(1)-R(4))
	  ELSE IF(NX_OUT .EQ. 3)THEN
	    R(2)=R(1)-0.05*(R(1)-R(5))
	    R(3)=R(1)-0.15*(R(1)-R(5))
	    R(4)=R(1)-0.40*(R(1)-R(5))
	  ELSE IF(NX_OUT .EQ. 3)THEN
	    R(2)=R(1)-0.015*(R(1)-R(6))
	    R(3)=R(1)-0.05*(R(1)-R(6))
	    R(4)=R(1)-0.15*(R(1)-R(6))
	    R(5)=R(1)-0.40*(R(1)-R(6))
	  END IF
!
	  DO I=1,NN+NX_OUT
	    V(I)=VINF*(1.0D0-RX/R(I))**BETA
	    SIGMA(I)=BETA*V(I)*RX/R(I)/R(I)/(1.0D0-RX/R(I))-1.0D0
	  END DO
!
! We remove the fine grid in the old model.
!
	  J=NN+NX_OUT+1
	  R(J)=OLD_R(1)
	  V(J)=OLD_V(1)
	  SIGMA(J)=OLD_SIGMA(1)
!
	  DO I=NX_OUT+2,ND_OLD
	    J=NN+I
	    R(J)=OLD_R(I)
	    V(J)=OLD_V(I)
	    SIGMA(J)=OLD_SIGMA(I)
	  END DO
	  ND=NN+ND_OLD
	ELSE IF(OPTION .EQ. 'MDOT')THEN
!
	  ND=ND_OLD
	  CALL GEN_IN(OLD_MDOT,'Old mass-loss rate in Msun/yr')
	  CALL GEN_IN(MDOT,'New mass-loss rate in Msun/yr')
	  CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	  CALL GEN_IN(BETA,'Beta for velocity law')
	  CALL GEN_IN(V_TRANS,'Connection velocity in km/s')
!
	  WRITE(6,'(A)')
	  WRITE(6,'(A)')'Type 1: W(r).V(r) = 2V(t) + (Vinf-2V(t))*(1-r(t)/r))**BETA'
	  WRITE(6,'(A)')'Type 2: W(r).V(r) = Vinf*(1-rx/r)**BETA'
	  WRITE(6,'(A)')'        with W(r) = 1.0D0+exp( (r(t)-r)/h )'
	  WRITE(6,'(A)')
!
	  CALL GEN_IN(VEL_TYPE,'Velocity law to be used: 1 or 2')
!
! Find conection velocity and index.
!
	  DO I=1,ND_OLD
	    IF(V_TRANS .LE. OLD_V(I) .AND. V_TRANS .GE. OLD_V(I+1))THEN
	      TRANS_I=I
	      EXIT
	    END IF
	  END DO
	  IF( OLD_V(TRANS_I)-V_TRANS .GT. V_TRANS-OLD_V(TRANS_I+1))TRANS_I=TRANS_I+1
	  V_TRANS=OLD_V(TRANS_I)
	  R(1:ND)=OLD_R(1:ND_OLD)
!
! In the hydrostatic region, the velocity is simply scaled by the change in
! mass-loss rate. This preserves the density. Only valid if wind does not have
! a significant optical depth.
!
	  DO I=TRANS_I,ND
	    V(I)=MDOT*OLD_V(I)/OLD_MDOT
	    SIGMA(I)=OLD_SIGMA(I)
	  END DO
!
! Now do the new wind law, keeping the same radius grid.
!
	  R_TRANS=R(TRANS_I)
	  V_TRANS=MDOT*V_TRANS/OLD_MDOT
	  dVdR_TRANS=(SIGMA(TRANS_I)+1.0D0)*V_TRANS/R_TRANS
!
	  IF(VEL_TYPE .EQ. 1)THEN
	    RO = R_TRANS * (1.0D0 - (2.0D0*V_TRANS/VINF)**(1.0D0/BETA) )
	    T1= R_TRANS * dVdR_TRANS / V_TRANS
	    SCALE_HEIGHT =  0.5D0*R_TRANS / (T1 - BETA*RO/(R_TRANS-RO) )
! 
	    WRITE(6,*)'  Transition radius is',R_TRANS
	    WRITE(6,*)'Transition velocity is',V_TRANS
	    WRITE(6,*)'                 R0 is',RO
	    WRITE(6,*)'       Scale height is',SCALE_HEIGHT
!
	    DO I=1,TRANS_I-1
              T1=RO/R(I)
              T2=1.0D0-T1
              TOP = VINF* (T2**BETA)
              BOT = 1.0D0 + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
              V(I) = TOP/BOT
                                                                                
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
                                                                                
              dTOPdR = VINF * BETA * T1 / R(I) * T2**(BETA - 1.0D0)
              dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT )  / SCALE_HEIGHT
              dVdR = dTOPdR / BOT  + V(I)*dBOTdR/BOT
              SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    END DO
	  ELSE
	    SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	    WRITE(6,*)'  Transition radius is',R_TRANS
	    WRITE(6,*)'Transition velocity is',V_TRANS
	    WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	    DO I=1,TRANS_I-1
	      T1=R_TRANS/R(I)
	      T2=1.0D0-T1
	      TOP = 2.0D0*V_TRANS + (VINF-2.0D0*V_TRANS) * T2**BETA
	      BOT = 1.0D0 + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
	      V(I) = TOP/BOT
                                                                                
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
                                                                                
	      dTOPdR = (VINF - 2.0D0*V_TRANS) * BETA * T1 / R(I) * T2**(BETA - 1.0D0)
	      dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
	      dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
              SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    END DO
	  END IF
	ELSE
	  WRITE(6,'(A)')'Option not recognixed'
	  STOP
	END IF
!
	NEW_RVSIG_FILE='RVSIG_COL_NEW'
	CALL GEN_IN(NEW_RVSIG_FILE,'File for new R, V and sigma values')
	OPEN(UNIT=10,FILE=NEW_RVSIG_FILE,STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A,7X,A,9X,10X,A,11X,A,3X,A)')'!','R','V(km/s)','Sigma','Depth'
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')' '
	  WRITE(10,'(I4,20X,A)')ND,'!Number of depth points`'
	  WRITE(10,'(A)')' '
	  DO I=1,ND
	    WRITE(10,'(F18.8,ES17.7,F17.7,4X,I4)')R(I),V(I),SIGMA(I),I
	  END DO
	CLOSE(UNIT=10)
!
	STOP
	END
