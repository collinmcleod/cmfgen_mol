	SUBROUTINE BRANCH_RAT(BETA,XV,YV,XAXIS,YAXIS,XSPEC,NMAX,NV,ND)
	USE SET_KIND_MODULE
	USE MOD_DISP
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE MOD_LEV_DIS_BLK
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER NMAX
	INTEGER ND
	INTEGER NV
	REAL(KIND=LDP) BETA(NMAX,NMAX)
	REAL(KIND=LDP) XV(NV),YV(NV)
	CHARACTER(LEN=*) XSPEC
!
	CHARACTER(LEN=*) XAXIS
	CHARACTER(LEN=*) YAXIS
        CHARACTER(LEN=120) DEFAULT
!
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER I,J,K,IB
	INTEGER NL,NUP
	INTEGER ID
!
	LOGICAL FLAG
	LOGICAL RADIAL_TAU
	LOGICAL LINE_STRENGTH
	LOGICAL VALID_VALUE
	LOGICAL FOUND
!
        CHARACTER*30 UC
        EXTERNAL UC
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	DOUBLE PRECISION CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	REAL(KIND=LDP) ASUM(2000)
	REAL(KIND=LDP) GLDGU
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) TAU_CONSTANT
	REAL(KIND=LDP) ANG_TO_HZ
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER, PARAMETER :: NB=8
	REAL(KIND=LDP) RANGE_START(NB),RANGE_END(NB)
	DATA RANGE_START/0.99D0, 0.50, 0.25D0, 0.125D0, 0.0625D0,0.03125D0, 0.01, 0.001D0/ !.0.1D0,0.01D0,0.001D0,0.0001D0,0.00001D0/
!
	RANGE_END(1)=1.0D0
	RANGE_END(2:NB)=RANGE_START(1:NB-1)
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07      !10^8/10^15
!
	I=ND/2
	DEFAULT=WR_STRING(I)
	VALID_VALUE=.FALSE.
	DO WHILE(.NOT. VALID_VALUE)
	  CALL USR_OPTION(I,'DEPTH',DEFAULT,'Depth for plotting TAUL')
	  IF(I .GE. 1 .AND. I .LE. ND)THEN
	     WRITE(T_OUT,'(A,I4,A,ES12.4)')'     Radius at depth',I,' is',R(I)
	     WRITE(T_OUT,'(A,I4,A,ES12.4)')'   Velocity at depth',I,' is',V(I)
	     WRITE(T_OUT,'(A,I4,A,ES12.4)')'Temperature at depth',I,' is',T(I)
	     VALID_VALUE=.TRUE.
	  END IF
	END DO
	CALL USR_OPTION(RADIAL_TAU,'RD_TAU','TRUE',
	1      'Use radial (alt. is TANGENTIAL) direction to evaluate the Sobolev optical depth')
!
	TAU_CONSTANT=OPLIN*R(I)*2.998E-10/V(I)
	IF(RADIAL_TAU)TAU_CONSTANT=TAU_CONSTANT/(1.0D0+SIGMA(I))
!
! Compute line opacity and emissivity.
!
	WRITE(6,*)'XSPEC=',XSPEC
	DO IB=1,NB
	  J=0
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') )THEN
!
	      ASUM=0.0D0
	      DO NUP=2,ATM(ID)%NXzV_F
	        DO NL=1,NUP-1
	          T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	          FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	          GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	          T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	          IF(T2 .NE. 0.0D0)THEN
	            T2=T2*TAU_CONSTANT/FREQ
	            BETA(NUP,NL)=(1.0D0-EXP(-T2))/T2
	            ASUM(NUP)=ASUM(NUP)+ATM(ID)%AXzV_F(NUP,NL)*BETA(NUP,NL)
!	            WRITE(37,'(2I4,4ES12.4)')NL,NUP,ANG_TO_HZ/FREQ,BETA(NUP,NL),T2,ASUM(NUP)
	          ELSE
	            BETA(NUP,NL)=0.0D0
	          END IF
	        END DO
	      END DO
!
	      FLAG=.FALSE.
	      DO NL=1,ATM(ID)%NXzV_F
	        DO NUP=NL+1,ATM(ID)%NXzV_F
	         T1=1000
	         IF(ASUM(NUP) .GT. 0.0D0)T1=ATM(ID)%AXzV_F(NUP,NL)*BETA(NUP,NL)/ASUM(NUP)
	          IF(T1 .GE. RANGE_START(IB) .AND. T1 .LE. RANGE_END(IB))THEN
	            J=J+1
	            T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	            FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	            GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	            IF(FLAG)THEN
                      XV(J)=LOG10( ANG_TO_HZ/FREQ )
	            ELSE
                      XV(J)=ANG_TO_HZ/FREQ
	            END IF
	            T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	            WRITE(25,*)NL,NUP,T1,T2
	            WRITE(25,*)J,XV(J),T2
	            IF(T2 .NE. 0)THEN
	              YV(J)=LOG10( ABS(T2)*TAU_CONSTANT/FREQ )
	            ELSE
	              J=J-1
	            END IF
	            WRITE(25,*)J,XV(J),T2
	            FOUND=.TRUE.
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
!	!  IF(.NOT. FOUND)THEN
!	!    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
!	  END IF
	  IF(J  .NE. 0)CALL DP_CURVE(J,XV,YV)
	  IF(J .NE. 0)WRITE(T_OUT,'(I6,A,F7.5,A,F7.5)')J,' lines plotted with a braching ratio between ',RANGE_START(IB),
	1         ' and ',RANGE_END(IB)
	END DO
!
	IF(FLAG)THEN
	  XAXIS='Log(\gl(\A))'
	ELSE
	  XAXIS='\gl(\gV)'
	END IF
	YAXIS='Log(\gt)'
	WRITE(T_OUT,'(A)')
	WRITE(T_OUT,'(A)')'Use C option in PLOT_SPEC to draw vertical lines (V)'
	WRITE(T_OUT,'(A)')
!
	RETURN
	END
