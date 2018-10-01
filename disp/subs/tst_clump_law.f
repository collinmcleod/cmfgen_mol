!
! Routine to be used with DISPGEN to test different clumping laws.
!
	SUBROUTINE TST_CLUMP_LAW(CLUMP_FAC,R,VEL,TAU,ND)
	USE MOD_COLOR_PEN_DEF
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Created 01-Oct-2018
!
	INTEGER ND
	REAL*8 CLUMP_FAC(ND)
	REAL*8 R(ND)
	REAL*8 VEL(ND)
	REAL*8 TAU(ND)
!
	INTEGER, PARAMETER :: NPAR_MAX=5
	INTEGER, PARAMETER :: LUER=6
	REAL*8 CLUMP_PAR(NPAR_MAX)
	REAL*8 T1
	INTEGER NPAR
	INTEGER K
	CHARACTER(LEN=10) CLUMP_LAW
!
        CHARACTER(LEN=30) UC
        EXTERNAL UC
!
	WRITE(LUER,'(A)')BLUE_PEN
	WRITE(LUER,'(A)')'EXPO:   F=C1 + (1-C1-C3)EXP(-V/C2) + C3.EXP(-V/C4)'
	WRITE(LUER,'(A)')'MEXP:   F=(C1 + (1-C1-C3)EXP(-V/C2))EXP(-R/R(1))'
	WRITE(LUER,'(A)')'REXP:   F=C1 + (1-C1)EXP(-V/C2)    + (1-C1)EXP((V-VINF)/C3)'
	WRITE(LUER,'(A)')'POW:    F=C1 + (1-C1)(-V/VINF)**C2'
	WRITE(LUER,'(A)')DEF_PEN
!
	CALL USR_OPTION(CLUMP_LAW,'LAW','EXPO','Clumping law: EXPO, CREXP, POW')
	CALL USR_OPTION(NPAR,'NPAR','2','Number of clumping parameters [ <6 ]')
	CLUMP_LAW=UC(CLUMP_LAW)
!
	IF(CLUMP_LAW(1:4) .EQ. 'EXPO')THEN
!
! CLUMP_PAR(1) is the clumping factor at infinity.
! CLUMP_PAR(2) is a velocity, and determines how fast the clumping factor
! approach CLUMP_PAR(1).
!
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  IF(NPAR .EQ. 3 .OR. NPAR .EQ. 4)THEN
	    CALL USR_OPTION(CLUMP_PAR(3),'CLP3','0.1','Second clumping factor')
	    CALL USR_OPTION(CLUMP_PAR(4),'CLP4','1.0','Second velocity scale')
	  END IF
	  IF(CLUMP_PAR(3) .EQ. 0.0D0)CLUMP_PAR(4)=1.0D0
	  DO K=1,ND
	    CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0D0-CLUMP_PAR(1)-CLUMP_PAR(3))*
	1                   EXP(-VEL(K)/CLUMP_PAR(2))+
	1                   CLUMP_PAR(3)*EXP(-VEL(K)/CLUMP_PAR(4))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:4) .EQ. 'MEXP')THEN
!
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity is CLP1/(1+CLP3)')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  CALL USR_OPTION(CLUMP_PAR(3),'CLP3','1.0','Second clumping factor')
	  CALL USR_OPTION(CLUMP_PAR(4),'CLP4','0.5','Radius scale (< 1.0')
	  CALL USR_OPTION(CLUMP_PAR(5),'CLP5','0.5','Exponent ')
	  DO K=1,ND
	    T1=1.0D0-(MAX(1.0D0,CLUMP_PAR(4)*R(1)/R(K)))**CLUMP_PAR(5)
	    CLUMP_FAC(K)=(CLUMP_PAR(1)+(1.0D0-CLUMP_PAR(1))*EXP(-VEL(K)/CLUMP_PAR(2)))/
	1                     (1.0D0+CLUMP_PAR(3)*EXP(T1))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:4) .EQ. 'REXP')THEN
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  CALL USR_OPTION(CLUMP_PAR(3),'CLP3','0.1','Second velocity factor')
	  DO K=1,ND
	    CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0D0-CLUMP_PAR(1))*
	1         EXP(-VEL(K)/CLUMP_PAR(2))+
	1        (1.0D0-CLUMP_PAR(1))*EXP( (VEL(K)-VEL(1))/CLUMP_PAR(3))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:3) .EQ. 'POW')THEN
	  IF(NPAR .NE. 2)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)' WRONG VALUE NPAR=',NPAR
	    WRITE(LUER,*)'NPAR assumed to be 2'
	  END IF
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','1.0D0','Power law exponent (+ve)')
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0D0-(1.0D0-CLUMP_PAR(1))*(VEL(K)/VEL(1))**CLUMP_PAR(2)
	  END DO
	ELSE
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Invalid law for computing clumping factor'
	  WRITE(LUER,*)'Set it to 1'
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0D0
	  END DO
	END IF
!
	RETURN
	END
