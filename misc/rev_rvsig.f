!
! Program to modify RVSIG_COL. Various options are available.
! Ideal for revising grid etc.
!
	PROGRAM REVISE_RVSIG
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 24-Oct-2013: Added VEL_TYPE=3 so that can do a velocity law with 2 components.
! Altered 14-Mar-2011: Improved header output to RVSIG_COL.
! Altered 15-Oct-2010: Fixed bug with SIGMA computation for the extra
!                        points added with the EXTR option.
! Altered 27-Aug-2007: Revised file read so as all ! (1st character) 
! comments are ignored.
!
	INTEGER, PARAMETER :: NMAX=5000
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
	REAL*8 RTMP(NMAX)
	REAL*8 OLD_TAU(NMAX)
	REAL*8 TAU_SAV(NMAX)
	REAL*8 TAU(NMAX)
	REAL*8 CHI_ROSS(NMAX)
	REAL*8 CLUMP_FAC(NMAX)
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
	REAL*8 NEW_RSTAR
	REAL*8 T1,T2,T3
	REAL*8 BETA
	REAL*8 RP2
	REAL*8 VEXT
	REAL*8 BETA2
	REAL*8 VINF
	REAL*8 FAC
	REAL*8 V_MAX
	REAL*8 V_MIN 
!
	REAL*4 XVAL,YVAL
	REAL*8 R_TRANS
	REAL*8 V_TRANS
	REAL*8 dVdR_TRANS
	REAL*8 SCALE_HEIGHT
	REAL*8 RO
	REAL*8 dVdR
	REAL*8 MDOT
	REAL*8 OLD_MDOT
	REAL*8 LSTAR
	REAL*8 OLD_LSTAR
	REAL*8 TOP,BOT 
	REAL*8 dTOPdR,dBOTdR 
	REAL*8 ALPHA
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
	INTEGER N_HEAD
	INTEGER IOS
	INTEGER PGCURS
	INTEGER CURSERR
!
	LOGICAL ROUND_ERROR
	LOGICAL RD_MEANOPAC
	LOGICAL REPLOT
!
        CHARACTER*30 UC
        EXTERNAL UC
!
	CHARACTER(LEN=1) CURSVAL
	CHARACTER(LEN=10) OPTION
	CHARACTER(LEN=80) OLD_RVSIG_FILE
	CHARACTER(LEN=80) NEW_RVSIG_FILE
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) OLD_HEADER(30)
!
	OLD_RVSIG_FILE='RVSIG_COL_OLD'
	CALL GEN_IN(OLD_RVSIG_FILE,'File containing old R, V and sigma values')
	OPEN(UNIT=10,FILE=OLD_RVSIG_FILE,STATUS='OLD',ACTION='READ')
	  STRING=' '
	  N_HEAD=0
	  DO WHILE (INDEX(STRING,'!Number of depth points') .EQ. 0)
	    READ(10,'(A)')STRING
	    N_HEAD=N_HEAD+1
	    N_HEAD=MIN(30,N_HEAD)
	    OLD_HEADER(N_HEAD)=STRING
	  END DO
	  N_HEAD=N_HEAD-1
	  READ(STRING,*)ND_OLD
	  STRING=' '
	  DO WHILE (STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)OLD_R(1),OLD_V(1),OLD_SIGMA(1)
	  DO I=2,ND_OLD
	    READ(10,*)OLD_R(I),OLD_V(I),OLD_SIGMA(I)
	  END DO
	CLOSE(UNIT=10)
!
	OPTION='NEW_ND'
	WRITE(6,'(A)')
	WRITE(6,'(A)')'Current options are:'
	WRITE(6,'(A)')'       NEW_ND: revise number of depth points by simple scaling'
	WRITE(6,'(A)')'         ADDR: explicitly add extra grid points'
	WRITE(6,'(A)')'         EXTR: extend grid to larger radii'
	WRITE(6,'(A)')'         FGOB: insert N points at outer boundary (to make finer grid)'
	WRITE(6,'(A)')'         MDOT: change the mass-loss rate or velocity law'
	WRITE(6,'(A)')'         NEWG: revise grid between two velocities or depth indicies'
	WRITE(6,'(A)')'         SCLR: scale radius of star to new value'
	WRITE(6,'(A)')'         SCLV: scale velocity law to new value'
	WRITE(6,'(A)')'         PLOT: plot V and SIGMA from old RVSIG file'
	WRITE(6,'(A)')

	CALL GEN_IN(OPTION,'Enter option for revised RVSIG file')
	OPTION=UC(TRIM(OPTION))
!
! Read in optical depth scale. Needed for SPP and TAU options. Also used
! for checking purposes (if available) for some other options.
!
! We use TAU_SAV for dTAU, and is used to increase the precision of the TAU
! scale (as insufficient digits may be print out).
!
	ROUND_ERROR=.FALSE.
	OPEN(UNIT=20,FILE='MEANOPAC',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(20,'(A)')STRING
	    DO I=1,ND_OLD
	      READ(20,*)RTMP(I),J,OLD_TAU(I),TAU_SAV(I),T1,CHI_ROSS(I)
	      J=MAX(I,2)
	      T1=OLD_R(J-1)-OLD_R(J)
	      IF( ABS(RTMP(I)-R(I))/T1 .GT. 2.0D-03 .AND. .NOT. ROUND_ERROR)THEN
	        WRITE(6,*)' '
	        WRITE(6,*)'Possible eror with MEANOPAC -- inconsistent R grid'
	        WRITE(6,*)'Error could simply be a lack of sig. digits in MEANOPAC'
	        WRITE(6,*)' RMO(I)=',RTMP(I)
	        WRITE(6,*)'   R(I)=',OLD_R(I)
	        WRITE(6,*)' R(I+1)=',OLD_R(I+1)
	        ROUND_ERROR=.TRUE.
	        CALL GEN_IN(ROUND_ERROR,'Continue as only rounding error?')
	        IF(.NOT. ROUND_ERROR)STOP
	      END IF
	    END DO
	    DO I=8,1,-1
              OLD_TAU(I)=OLD_TAU(I+1)-TAU_SAV(I)
            END DO
	    RD_MEANOPAC=.TRUE.
	    WRITE(6,*)'Successfully read MEANOPAC'
	  ELSE
	    RD_MEANOPAC=.FALSE.
	    WRITE(6,*)'MEANOPAC has not been read'
	    IF(OPTION .EQ. 'SPP' .OR. OPTION .EQ. 'TAU')THEN
	      WRITE(6,*)'MENAOPAC  is equired for the TAU option'
	      STOP
	    END IF
	    CLOSE(UNIT=20)
	  END IF
	  IF(ROUND_ERROR .AND. RD_MEANOPAC)THEN
	     RTMP(1:ND_OLD)=OLD_R(1:ND_OLD)
	  END IF 
!
	IF(OPTION .EQ. 'SPP')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option takes a plane-parallel model and oututs a spherical model'
	  WRITE(6,'(A)')'The MEANOPAC from the plane-parallel model is required'
	  WRITE(6,'(A)')'The VADAT file is also required'
	END IF
!
	IF(OPTION .EQ. 'TAU')THEN
	  WRITE(6,*)'Option still under development'
	  I=7
          CALL RD_SING_VEC_RVTJ(CLUMP_FAC,ND_OLD,'Clumping Factor','RVTJ',I,IOS)
	  ND=ND_OLD
	  CALL GEN_IN(ND,'Number of depth points')
	  T1=MAXVAL(OLD_SIGMA(1:ND_OLD))
	  WRITE(6,*)'Maximum value of SIGMA is',T1
	  DO I=1,ND_OLD
!	    X2(I)=CLUMP_FAC(I)*CHI_ROSS(I)*(0.5D0+OLD_SIGMA(I)/10.0D0)   !/(1.0D0+SQRT(OLD_V(I))) !*(1.0D0+OLD_SIGMA(I)/T1)
	    X2(I)=CLUMP_FAC(I)*CHI_ROSS(I)*(1.0D0+(1.0D0+OLD_SIGMA(I)))
!	    X2(I)=CLUMP_FAC(I)*CHI_ROSS(I)
	  END DO
	  X1(1)=X2(1)*OLD_R(1)
	  DO I=2,ND_OLD
	    X1(I)=X1(I-1)+ 0.5D0*(OLD_R(I-1)-OLD_R(I))*(X2(I-1)+X2(I))
	  END DO
!
!	  DO I=1,ND_OLD
!	    T1=MIN(OLD_V(I),300.0D0)/0.1D0
!	    X1(I)=X1(I)/MAX(1.0D0,SQRT(T1))
!	  END DO
!
	  T1=EXP(LOG(X1(ND_OLD)/X1(1))/(ND-4))
	  WRITE(6,*)'TAU step size ratio is',T1
	  X2(1)=X1(1)
	  X2(2)=X1(1)+0.02*X1(1)*(T1-1.0D0)
	  X2(3)=X1(1)*T1
	  DO I=4,ND-3
	    X2(I)=X2(I-1)*T1
	  END DO
	  X2(ND-2)=X2(ND-3)*SQRT(T1)
	  X2(ND-1)=X2(ND-2)*(T1**0.25D0)
	  X2(ND)=X1(ND_OLD)
	  WRITE(6,*)'Defined X2 grid'
! 
	  DO I=1,ND_OLD
	    WRITE(6,*)I,OLD_R(I),X2(I)
	  END DO
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
	  WRITE(6,*)'Defined new R grid by interpolation'
	  DO I=1,ND
	    WRITE(6,*)I,R(I),X2(I)
	  END DO
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
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0D0*COEF(J,2)+3.0*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
	  WRITE(6,*)'Determined R and V on the new grid'
	  CALL MON_INTERP(TAU,ND,IONE,R,ND,OLD_TAU,ND_OLD,OLD_R,ND_OLD)
	  WRITE(6,*)'Determined TAU on the new grid'
	  I=1
	  WRITE(30,'(3X,A,6X,A,T28,A,T40,A,T53,A,T66,A,T78,A,T93,A,T101,A)')
	1           'I','R','V','SIGMA','TAU','dTAU','dTRAT','X2','V(I)/V(I+1)'
	  WRITE(30,'(1X,I3,F16.9,7ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I),TAU(I+1)-TAU(I),T2,X2(I),V(I)/V(I+1)
	  FLUSH(UNIT=30)
	  T2=0.00
	  DO I=2,ND-1
	    T1=TAU(I+1)-TAU(I)
	    T2=T1/(TAU(I)-TAU(I-1))
	    WRITE(30,'(1X,I3,F16.9,7ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I),T1,T2,X2(I),V(I)/V(I+1)
	    FLUSH(UNIT=30)
	  END DO
	  I=ND
	  WRITE(30,'(1X,I3,F16.9,6ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I),0.0D0,0.0D0,X2(I)
	  FLUSH(UNIT=30)
!
	ELSE IF(OPTION .EQ. 'NEW_ND')THEN
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
	ELSE IF(OPTION .EQ. 'SCLV')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option simply scales V'
	  WRITE(6,'(A)')'Sigma will not be changed by this operation.'
	  WRITE(6,'(A)')' '
	  T1=1.1D0
	  ND=ND_OLD
	  CALL GEN_IN(T1,'Factor to sclae velocity')
	  V(1:ND)=T1*OLD_V(1:ND)
	  R(1:ND)=OLD_R(1:ND)
	  SIGMA(1:ND)=OLD_SIGMA(1:ND)
!
	ELSE IF(OPTION .EQ. 'FG')THEN
          I_ST=1; I_END=ND
          CALL GEN_IN(I_ST,'Start index for fine grid')
          CALL GEN_IN(I_END,'End index for fine grid')
          WRITE(6,*)'Number of points in requeted interval is',I_END-I_ST-1
          NX=I_END-I_ST
          CALL GEN_IN(NX,'New number of grid points for this interval')
          ND=I_ST+NX+(ND_OLD-I_END)+1
!
            DO I=1,I_ST
              X2(I)=I
            END DO
            T1=(I_END-I_ST)/(NX+1.0D0)
            DO I=1,NX
              X2(I_ST+I)=I_ST+I*T1
            END DO
            DO I=I_END,ND_OLD
              X2(I_ST+NX+I+1-I_END)=I
            END DO
           DO I=1,ND
             X1(I)=I
           END DO
           CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
	   WRITE(6,*)'Done interpolation'
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
	  WRITE(6,'(A)')'This option allows the grid to be redefined over a specified interval'
	  WRITE(6,'(A)')'Originally defined for small ranges'
	  WRITE(6,'(A)')'Data points in the requested interval are equally spaced on log r^2. V'
	  WRITE(6,'(A)')' '
!
	  CALL GEN_IN(V_MAX,'Maximum velocity for grid refinement (-ve for depth index)')
	  CALL GEN_IN(V_MIN,'Initial velocity for grid refinement (-ve for depth index)')
	  IF(V_MAX*V_MIN .LT. 0.0D0)THEN
	    WRITE(6,*)'Error: either velocities of depth index but not mixture'
	    STOP
	  END IF
!
	  IF(V_MAX .LT. 0.0D0)THEN
	    V_MAX=ABS(V_MAX); V_MIN=ABS(V_MIN)
	    I_ST=NINT(MIN(V_MIN,V_MAX))
	    I_END=NINT(MAX(V_MIN,V_MAX))
	  ELSE
	    IF(V_MAX .GE. OLD_V(1) .OR. V_MAX .LE. OLD_V(ND_OLD))THEN
	      WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	      STOP
	    END IF
	    IF(V_MIN .GE. OLD_V(1) .OR. V_MIN .LE. OLD_V(ND_OLD))THEN
	      WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	      STOP
	    END IF
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
	    IF(V_MIN-OLD_V(I_END+1) .LT. OLD_V(I_END)-V_MIN)I_END=I_END+1
	  END IF
!
	  WRITE(6,*)' I_ST=',I_ST,OLD_V(I_ST)
	  WRITE(6,*)'I_END=',I_END,OLD_V(I_END)
	  V_MAX=OLD_V(I_ST); V_MIN=OLD_V(I_END)
	  WRITE(6,*)'Current number of points in interval is',I_END-I_ST-1
	  CALL GEN_IN(N_ADD,'Number of points in interval (exclusive)')
!
	  X1(1:ND_OLD)=OLD_V(1:ND_OLD)*OLD_R(1:ND_OLD)*OLD_R(1:ND_OLD)
	  ND=N_ADD+ND_OLD-(I_END-I_ST-1)
	  X2(1:I_ST)=X1(1:I_ST)
	  T1=EXP(DLOG(V_MAX*OLD_R(I_ST)*OLD_R(I_ST)/V_MIN/OLD_R(I_END)/OLD_R(I_END))/(N_ADD+1))
	  DO I=I_ST+1,I_ST+N_ADD
	   X2(I)=X2(I-1)/T1
	  END DO
	  X2(I_ST+N_ADD+1:ND)=X1(I_END:ND_OLD)
!
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
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
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0D0*COEF(J,2)+3.0*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
	  
	ELSE IF(OPTION .EQ. 'FGOB')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows the grid to be redefined at the outer boundary'
	  WRITE(6,'(A)')'Adding N points decreases the spacing at the outer boundary by 3^N'
!
	  N_ADD=1
	  CALL GEN_IN(N_ADD,'Number of grid points to be added at outer boundary')
	  ND=ND_OLD+N_ADD
	  R(N_ADD+2:ND)=OLD_R(2:ND_OLD)
	  R(1)=OLD_R(1)
	  DO I=N_ADD,1,-1
	    R(I+1)=R(1)-0.3333D0*(R(1)-R(I+2))
	    WRITE(6,*)I,R(1),R(I+1)
	  END DO
!
! Now compute the revised V & SIGMA. The new points lie in the first interval.
!
	  V(1)=OLD_V(1); V(N_ADD+2:ND)=OLD_V(2:ND_OLD)
	  SIGMA(1)=OLD_SIGMA(1); SIGMA(N_ADD+2:ND)=OLD_SIGMA(2:ND_OLD)
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  DO I=2,N_ADD+1
	    T1=R(I)-OLD_R(1)
	    V(I)=COEF(1,4)+T1*(COEF(1,3)+T1*(COEF(1,2)+T1*COEF(1,1)))
	    SIGMA(I)=COEF(1,3)+T1*(2.0D0*COEF(1,2)+3.0*T1*COEF(1,1))
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
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
	    SIGMA(I)=BETA*RX/R(I)/(1.0D0-RX/R(I))-1.0D0
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
	ELSE IF(OPTION .EQ. 'SCLR')THEN
!
	  ND=ND_OLD
	  NEW_RSTAR=OLD_R(ND_OLD)
	  V_TRANS=4.0D0
	  CALL GEN_IN(NEW_RSTAR,'New radius')
	  CALL GEN_IN(V_TRANS,'Connection velocity in km/s')
	  CALL GEN_IN(OLD_MDOT,'Old mass-loss rate in Msun/yr')
	  MDOT=OLD_MDOT
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'The defaut mass assume Mdot/R^1.5 is to be preserved'
	  WRITE(6,*)DEF_PEN
	  MDOT=OLD_MDOT*(NEW_RSTAR/OLD_R(ND_OLD))**1.5D0
	  CALL GEN_IN(MDOT,'New mass-loss rate in Msun/yr')
!
	  CALL GEN_IN(OLD_LSTAR,'Old stellar luminosity in Lsun')
	  T1=OLD_LSTAR*(NEW_RSTAR/OLD_R(ND_OLD))**2
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(A,ES14.4)')'New luminosity if Teff is to be preserved',T1
	  WRITE(6,*)DEF_PEN
!
	  WRITE(6,'(A)')
	  WRITE(6,'(A)')'Type 1: W(r).V(r) = Vinf*(1-rx/r)**BETA'
	  WRITE(6,'(A)')'        with W(r) = 1.0D0+exp( (r(t)-r)/h )'
	  WRITE(6,'(A)')'Type 2: W(r).V(r) = 2V(t) + (Vinf-2V(t))*(1-r(t)/r))**BETA'
	  WRITE(6,'(A)')'        with W(r) = 1.0D0+exp( (r(t)-r)/h )'
	  WRITE(6,'(A)')
!
	  VEL_TYPE=1
	  CALL GEN_IN(VEL_TYPE,'Velocity law to be used: 1 or 2')
	  VINF=1000.0D0
	  CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	  BETA=1.0D0
	  CALL GEN_IN(BETA,'Beta for velocity law')
!
! Find conection velocity and index.
!
	  TRANS_I=0
	  DO I=1,ND_OLD
	    IF(V_TRANS .LE. OLD_V(I) .AND. V_TRANS .GE. OLD_V(I+1))THEN
	      TRANS_I=I
	      EXIT
	    END IF
	  END DO
	  IF(TRANS_I .EQ. 0)THEN
	    WRITE(6,'(/,1X,A)')'Error V_TRANS is outside range'
	    WRITE(6,'(1X,3(A,ES15.8,3X))')'V_TRANS=',V_TRANS,'OLD_V(1)=',OLD_V(1),'V(ND_OLD)=',OLD_V(ND_OLD)
	  END IF

	  IF( OLD_V(TRANS_I)-V_TRANS .GT. V_TRANS-OLD_V(TRANS_I+1))TRANS_I=TRANS_I+1
	  V_TRANS=OLD_V(TRANS_I)
	  R(1:ND)=OLD_R(1:ND_OLD)+(NEW_RSTAR-OLD_R(ND_OLD))
	  R(ND)=NEW_RSTAR
!
! In the hydrostatic region, the velocity is simply scaled so as to keep the
! density constant.
!
	  T1=R(ND)/OLD_R(ND_OLD)
	  DO I=TRANS_I,ND
	    V(I)=MDOT*OLD_V(I)/T1/T1/OLD_MDOT
	    SIGMA(I)=(OLD_SIGMA(I)+1.0D0)*R(I)/OLD_R(I)-1.0D0
	  END DO
	  V_TRANS=V(TRANS_I)
!
! Now do the new wind law, keeping the same radius grid.
!
	  R_TRANS=R(TRANS_I)
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
!
	    SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	    WRITE(6,*)'  Transition radius is',R_TRANS
	    WRITE(6,*)'Transition velocity is',V_TRANS
	    WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	    ALPHA=2.0D0
!	    IF(VEL_TYPE .EQ. 4)ALPHA=3.0D0
	    BETA2=BETA
	    CALL GEN_IN(BETA2,'Beta2 for velocity law')
	    DO I=1,TRANS_I-1
	      T1=R_TRANS/R(I)
	      T2=1.0D0-T1
	      T3=BETA+(BETA2-BETA)*T2
	      TOP = (VINF-ALPHA*V_TRANS) * T2**T3
	      BOT = 1.0D0 + (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT )

!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.

	      dTOPdR = (VINF - ALPHA*V_TRANS) * BETA * T1 / R(I) * T2**(T3-1.0D0) +
	1                  T1*TOP*(BETA2-BETA)*(1.0D0+LOG(T2))/R(I)
	      dBOTdR=  (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
!
	      TOP = ALPHA*V_TRANS + TOP
	      dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
	      V(I) = TOP/BOT
              SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    END DO
	  END IF
!
	ELSE IF(OPTION .EQ. 'MDOT')THEN
!
	  ND=ND_OLD
	  OLD_MDOT=0.0D0; VINF=0.0D0; BETA=1.0D0; V_TRANS=10.0D0
	  CALL GEN_IN(OLD_MDOT,'Old mass-loss rate in Msun/yr')
	  MDOT=OLD_MDOT
	  CALL GEN_IN(MDOT,'New mass-loss rate in Msun/yr')
!
	  WRITE(6,*)' '
	  WRITE(6,*)' If increasing Mdot, simply enter a number slightly smaller than the sound speed.'
	  WRITE(6,*)' If decreasing Mdot, you may need to enter a smaller number,'
	  WRITE(6,*)' since as you decrease Mdot, the extent of the photosphere increase'
!
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Type 1: W(r).V(r) with'
	  WRITE(6,'(A)')'            W(r) = VINF*(1-rt/r)**BETA'
	  WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')'Type 2: W(r).V(r) with'
	  WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*(1-rt/r)**BETA'
	  WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Type 3: W(r).V(r) with'
          WRITE(6,'(A)')'              X  = 1-rt/r'
	  WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*X**[BETA+(BETA2-BETA)*X]'
	  WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')'Type 4: W(r).V(r) with'
          WRITE(6,'(A)')'              X  = 1-rt/r'
	  WRITE(6,'(A)')'            W(r) = 3*VTRANS + (VINF-2*VTRANS)*X**[BETA+(BETA2-BETA)*X]'
	  WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Type 5: (W(r)+E(r))/V(r) with'
	  WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*(1-rt/r)**BETA'
	  WRITE(6,'(A)')'            E(r)= VEXT*(1-RP2/r)**BETA2 '
	  WRITE(6,'(A)')'        and V(r) = 1.0D0+2*exp( (rt-r)/h )'
	  WRITE(6,'(A)')DEF_PEN
!
	  VEL_TYPE=2
	  CALL GEN_IN(VEL_TYPE,'Velocity law to be used: 1, 2, 3, 4 or 5')
	  CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	  CALL GEN_IN(BETA,'Beta for velocity law')
	  CALL GEN_IN(V_TRANS,'Connection velocity in km/s')
	  V_TRANS=V_TRANS*OLD_MDOT/MDOT
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
	  ELSE IF(VEL_TYPE .EQ. 2)THEN
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
	  ELSE IF(VEL_TYPE .EQ. 3 .OR. VEL_TYPE .EQ. 4)THEN
!
	    SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	    WRITE(6,*)'  Transition radius is',R_TRANS
	    WRITE(6,*)'Transition velocity is',V_TRANS
	    WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	    ALPHA=2.0D0
	    IF(VEL_TYPE .EQ. 4)ALPHA=3.0D0
	    CALL GEN_IN(BETA2,'Beta2 for velocity law')
	    DO I=1,TRANS_I-1
	      T1=R_TRANS/R(I)
	      T2=1.0D0-T1
	      T3=BETA+(BETA2-BETA)*T2
	      TOP = (VINF-ALPHA*V_TRANS) * T2**T3
	      BOT = 1.0D0 + (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT )

!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.

	      dTOPdR = (VINF - ALPHA*V_TRANS) * BETA * T1 / R(I) * T2**(T3-1.0D0) +
	1                  T1*TOP*(BETA2-BETA)*(1.0D0+LOG(T2))/R(I)
	      dBOTdR=  (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
!
	      TOP = ALPHA*V_TRANS + TOP
	      dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
	      V(I) = TOP/BOT
              SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    END DO
	  ELSE IF(VEL_TYPE .EQ. 5)THEN
	    BETA2=BETA; VEXT=0.1D0*VINF; RP2=2.0D0
	    CALL GEN_IN(VEXT,'Additonal V component (add to VINF)- VEXT')
	    CALL GEN_IN(BETA2,'Beta2 for velocity law')
	    CALL GEN_IN(RP2,'RP2 for velocity law (in terms of R(ND)')
	    RP2=RP2*R(ND)
	    SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	    WRITE(6,*)'  '
	    WRITE(6,*)'  Transition radius is',R_TRANS
	    WRITE(6,*)'Transition velocity is',V_TRANS
	    WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	    WRITE(6,*)'  '
	    DO I=1,TRANS_I-1
	      T1=R_TRANS/R(I)
	      T2=1.0D0-T1
	      TOP = 2.0D0*V_TRANS + (VINF-2.0D0*V_TRANS) * T2**BETA
	      BOT = 1.0D0 + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
              IF(RP2/R(I) .LT. 1.0D0)THEN
                 TOP=TOP+ VEXT*(1.0D0-RP2/R(I))**BETA2
                 dTOPdR=(RP2/R(I)/R(I))*BETA2*VEXT*(1.0D0-RP2/R(I))**(BETA2-1)
              ELSE
                 dTOPdR=0.0D0
	      END IF
	      V(I) = TOP/BOT
                                                                                
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
                                                                                
	      dTOPdR = dTOPdR + (VINF - 2.0D0*V_TRANS) * BETA * T1 / R(I) * T2**(BETA - 1.0D0)
	      dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
	      dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
              SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    END DO
	  END IF
!
	ELSE IF(OPTION .EQ. 'CUR')THEN
!
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,*)'Use ''r'' to replace a data point'
	  WRITE(6,*)'Use ''a'' to add a data point'
	  WRITE(6,*)'Use ''d'' to delete a data point'
	  WRITE(6,*)'Use ''e'' to exit'
	  WRITE(6,'(A)')DEF_PEN
!
	  CALL DP_CURVE(ND_OLD,OLD_R,OLD_V)
 	  CALL GRAMON_PGPLOT('R/R\d*\u','V(km/s)',' ',' ')
	  WRITE(6,'(A)')'Cursor now available to modify data points'
	  ND=ND_OLD
	  R(1:ND)=OLD_R(1:ND); V(1:ND)=OLD_V(1:ND)
!
	  DO WHILE(1 .EQ. 1)
	    DO WHILE(1 .EQ. 1)
	      CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	      WRITE(6,*)XVAL,YVAL
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
	      IF(CURSVAL .EQ. 'r' .OR. CURSVAL .EQ. 'R')THEN
	        T1=XVAL
	        TMP_R(1:ND)=ABS(R(1:ND)-XVAL)
	        I=MINLOC(TMP_R(1:ND),IONE)
	        V(I)=YVAL
	        WRITE(6,*)'Replaced V for R=',R(I)
	      ELSE IF(CURSVAL .EQ. 'd' .OR. CURSVAL .EQ. 'D')THEN
	        T1=XVAL
	        TMP_R(1:ND)=ABS(R(1:ND)-T1)
	        I=MINLOC(TMP_R(1:ND),IONE)
	        R(I:ND-1)=R(I+1:ND)
	        V(I:ND-1)=V(I+1:ND)
	        ND=ND-1
	        WRITE(6,*)'Deleted R=',R(I),' from grid'
	      ELSE IF(CURSVAL .EQ. 'a' .OR. CURSVAL .EQ. 'A')THEN
	        DO I=1,ND-1
	          IF( (R(I)-XVAL)*(R(I+1)-XVAL) .LT. 0 )THEN
	            DO J=ND,I+1,-1
                      R(J+1)=R(J)
                      V(J+1)=V(J)
	            END DO
	            R(I+1)=XVAL; V(I+1)=YVAL
	            ND=ND+1
	            WRITE(6,*)'Add R=',R(I),' to grid'
	            EXIT
	          END IF
	        END DO
	      ELSE
	        WRITE(6,*)RED_PEN
	        WRITE(6,*)'Error - use r(eplace), a(dd), d(elete), e'
	        WRITE(6,*)DEF_PEN
	      END IF
	    END DO
	    CALL GEN_IN(REPLOT,'Replot to see revision to V')
	    IF(REPLOT)THEN
	      CALL DP_CURVE(ND_OLD,OLD_R,OLD_V)
	      CALL DP_CURVE(ND,R,V)
 	      CALL GRAMON_PGPLOT('R/R\d*\u','V(km/s)',' ',' ')
	    ELSE
	      EXIT
	    END IF
	  END DO
!
! Compute SIGMA
!
	  ALLOCATE (COEF(ND,4))
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    SIGMA(I)=COEF(I,3)
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'PLOT')THEN
	  ND=ND_OLD 
	  DO I=1,ND
	    R(I)=OLD_R(I)
	    V(I)=OLD_V(I)
	    SIGMA(I)=OLD_SIGMA(I)
	  END DO
	ELSE
	  WRITE(6,'(A)')'Option not recognixed'
	  STOP
	END IF
!
	IF(OPTION .NE. 'PLOT')THEN
	  NEW_RVSIG_FILE='RVSIG_COL_NEW'
	  CALL GEN_IN(NEW_RVSIG_FILE,'File for new R, V and sigma values')
	  OPEN(UNIT=10,FILE=NEW_RVSIG_FILE,STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(10,'(A)')'!'
	    IF(OPTION .EQ. 'MDOT')THEN
	      WRITE(10,'(A,ES14.6)')'! Old mass-loss rate in Msun/yr=',OLD_MDOT
	      WRITE(10,'(A,ES14.6)')'! New mass-loss rate in Msun/yr=',MDOT 
	      WRITE(10,'(A,ES14.6)')'! Velocity at infinity in km/s =',VINF
	      WRITE(10,'(A,ES14.6)')'! Beta for velocity law        =',BETA
	      WRITE(10,'(A,I3)'    )'! Velocity law (type)          =',VEL_TYPE
	      WRITE(10,'(A,ES14.6)')'! Transition velocity is       =',V_TRANS
	      IF(VEL_TYPE .EQ. 5)THEN
	        WRITE(10,'(A,ES14.6)')'! Beta for 2nd vel. component  =',BETA2
	        WRITE(10,'(A,ES14.6)')'! VEXT for 2nd vel. component  =',VEXT
	        WRITE(10,'(A,ES14.6)')'! RP2  for 2nd vel. component  =',RP2/R(ND)
	        ELSE IF(VEL_TYPE .EQ. 3 .OR. VEL_TYPE .EQ. 3)THEN
	      WRITE(10,'(A,ES14.6)')'! Beta in inner wind           =',BETA2
	        WRITE(10,'(A,ES14.6)')'! ALPHA                        =',ALPHA
	        END IF
	      WRITE(10,'(A,ES14.6)')'! R(1)/R(ND)                   =',R(1)/R(ND)
	    ELSE IF(OPTION .EQ. 'SCLR')THEN
	      WRITE(10,'(A,ES14.6)')'! Old stellar radius           =',OLD_R(ND_OLD)
	      WRITE(10,'(A,ES14.6)')'! New stellar rdaius           =',R(ND_OLD)
	      WRITE(10,'(A,ES14.6)')'! Velocity at infinity in km/s =',VINF
	      WRITE(10,'(A,ES14.6)')'! Beta for velocity law        =',BETA
	      WRITE(10,'(A,I3)'    )'! Velocity law (type)          =',VEL_TYPE
	      WRITE(10,'(A,ES14.6)')'! Transition velocity is       =',V_TRANS
	      WRITE(10,'(A,ES14.6)')'! R(1)/R(ND)                   =',R(1)/R(ND)
	    END IF
!
! Use a !! to inidicate from old file.
!
	    DO I=1,N_HEAD
	      WRITE(10,'(A,A)')'!',TRIM(OLD_HEADER(I))
	    END DO
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
	END IF
!
	T1=R(ND)
	R(1:ND)=R(1:ND)/T1
	T1=OLD_R(ND_OLD)
	OLD_R(1:ND_OLD)=OLD_R(1:ND_OLD)/T1
	WRITE(6,*)'Plotting V versus R/R*'
	CALL DP_CURVE(ND_OLD,OLD_R,OLD_V)
	CALL DP_CURVE(ND,R,V)
	CALL GRAMON_PGPLOT('R/R\d*\u','V(km/s)',' ',' ')
	WRITE(6,*)'Plotting Sigma versus R/R*'
	CALL DP_CURVE(ND_OLD,OLD_R,OLD_SIGMA)
	CALL DP_CURVE(ND,R,SIGMA)
	CALL GRAMON_PGPLOT('R/R\d*\u','SIGMA',' ',' ')
!
	STOP
	END
