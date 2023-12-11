!
! Program to test the solution of the large martices generated by
! CMFGEN. The routine reads in the file BA_ASCI_N_DX (X=depth) and generates
! the solution.
!
	PROGRAM TST_BA_MAT_SOL
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
!
! Altered 24-Jun-2009: Depth of matrix now input.
!
	IMPLICIT NONE
	INTEGER ID
	INTEGER NT,N
	INTEGER ILOW,IUP
	INTEGER IOS
	LOGICAL FILE_OPEN
	CHARACTER(LEN=200) STRING
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Total number of variables') .NE. 0)THEN
	        READ(STRING,*)NT
	        WRITE(6,'(A,I4)')' Number of variables in the model is:',NT
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  NT=0
	  WRITE(6,*)' Unable to open/read MODEL file to get # of variables'
	  CALL GEN_IN(NT,'Number of elements (==NT in MODEL)')
	END IF
	IF(NT .EQ. 0)THEN
	  WRITE(6,*)'Invalid number of elements'
	  WRITE(6,*)'Check NT in file MODEL'
	  STOP
	END IF
!
	ID=1
	CALL GEN_IN(ID,'Depth indicator: i.e., X in BA_ASCI_N_DX:')
!
	ILOW=1; IUP=NT
	CALL GEN_IN(ILOW,'Lower index of sub matrix')
	CALL GEN_IN(IUP,'Upper index of sub matrix')
	N=IUP-ILOW+1
!
	CALL SOLVE_BA_MAT_V2(NT,N,ILOW,IUP,ID)
!
	STOP
	END
! 
!
	SUBROUTINE SOLVE_BA_MAT_V2(NT,N,ILOW,IUP,ID)
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER N,NT
	INTEGER ILOW,IUP
	INTEGER ID
	INTEGER IEQ
!
	REAL(KIND=LDP) POPS_RD(NT)
	REAL(KIND=LDP) STEQ_RD(NT)
	REAL(KIND=LDP) CMAT_RD(NT,NT)
!
! Allocate and declare dynamic arrays and vectors.
!
	REAL(KIND=LDP) POPS(N)
	REAL(KIND=LDP) PLTE(N)
	REAL(KIND=LDP) STEQ(N)
	REAL(KIND=LDP) CMAT(N,N)
	REAL(KIND=LDP) STAT_WT(N)
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) POP_SOLS(N,5)
!
	REAL(KIND=LDP) SAV_CMAT(N,N)
	REAL(KIND=LDP) SAV_STEQ(N)
	REAL(KIND=LDP) RHS(N)
	REAL(KIND=LDP) LARGEST_VAL(N)
	REAL(KIND=LDP) OLD_SOL(N)
!
	REAL(KIND=LDP) NEW_POPS(N)
	REAL(KIND=LDP) NEW_LTE(N)
	REAL(KIND=LDP) ED_NEW,DI_NEW,T_NEW
!
	REAL(KIND=LDP) ROW_SF(N)
	REAL(KIND=LDP) COL_SF(N)
!
	REAL(KIND=LDP) ROW_CND
	REAL(KIND=LDP) COL_CND
	REAL(KIND=LDP) MAX_VAL
	REAL(KIND=LDP) T1,T2
!
	REAL(KIND=LDP) HDKT
	REAL(KIND=LDP) GION,RGU,X,Y
	REAL(KIND=LDP) IONIZATION_ENERGY
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER IFAIL
	INTEGER IPIVOT(N)
!
	CHARACTER*200 STRING,FILENAME
	CHARACTER*1, PARAMETER :: NO_TRANS='N'
	INTEGER I,J,K,L
	INTEGER, PARAMETER :: NSNG=1
	LOGICAL USE_DC
	LOGICAL FILE_EXISTS
!
	HDKT=4.7994145_LDP
!
	WRITE(FILENAME,*)ID
	FILENAME=ADJUSTL(FILENAME)
	FILENAME='BA_ASCI_N_D'//TRIM(FILENAME)
	INQUIRE(FILE=FILENAME,EXIST=FILE_EXISTS)
	IF(.NOT. FILE_EXISTS)THEN
	  WRITE(FILENAME,'(I3.3)')ID
	  FILENAME='BA_ASCI_N_D'//TRIM(FILENAME)
	END IF
	OPEN(UNIT=10, FILE=FILENAME,STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'POP') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(10,'(A)')STRING
	  END DO
	  DO I=1,N
	    READ(STRING(10:),*)POPS_RD(I)
	    READ(10,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read POPS'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'STEQ') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(10,'(A)')STRING
	  END DO
	  DO I=1,NT
	    READ(STRING(10:),*)STEQ_RD(I)
	    READ(10,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read STEQ'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'C_MAT') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
!
	  STRING=' '
	  DO K=1,NT,5
	    DO WHILE(STRING .EQ. ' ')
	      READ(10,'(A)')STRING
	    END DO
	    IF( MOD(K-1,50) .EQ. 0)WRITE(6,*)'Reading element=',K,' in CMAT'
	    DO I=1,NT
	      READ(STRING(10:),*)(CMAT_RD(I,J),J=K,MIN(K+4,NT))
	      READ(10,'(A)',END=100)STRING
	    END DO
100	    CONTINUE
	  END DO
	  WRITE(6,*)'Successfully read CMAT'
!
	CLOSE(UNIT=10)
!
	CMAT=CMAT_RD(ILOW:IUP,ILOW:IUP)
	STEQ=STEQ_RD(ILOW:IUP)
	POPS=POPS_RD(ILOW:IUP)
	WRITE(6,*)'POPS(N)=',POPS_RD(N)
!
	SAV_CMAT=CMAT
	SAV_STEQ=STEQ
!
! For pure H models, we swicth to using depparture coefficients to
! see if this changes the accuracy of the solutions. We only request
! this option if N is small.
!
	USE_DC=.FALSE.
	IF(N .LT. 100)THEN
	  CALL GEN_IN(USE_DC,'Use departue coefficients? (only pure H models)')
	END IF
!
	IF(USE_DC)THEN
	  IONIZATION_ENERGY=109678.7640_LDP
	  DO I=1,10
	    EDGE(I)=IONIZATION_ENERGY*SPEED_OF_LIGHT()*1.0E-15_LDP/I/I
	    STAT_WT(I)=2.0_LDP*I*I
	  END DO
	  GION=1.0_LDP
!
	  RGU=LOG(2.07078E-22_LDP)
          X=HDKT/POPS(N)
	  Y=POPS(N-1)*POPS(N-2)*( POPS(N)**(-1.5_LDP) )/GION
          WRITE(6,*)RGU,X,Y
          DO I=1,N-3
            PLTE(I)=STAT_WT(I)*Y*EXP(EDGE(I)*X+RGU)
	    WRITE(6,*)EDGE(I),PLTE(I),POPS(I)
	  END DO
!
	  DO I=1,N
	    DO J=1,10
	      CMAT(I,N-2)=CMAT(I,N-2)+CMAT(I,J)
	      CMAT(I,N-1)=CMAT(I,N-1)+CMAT(I,J)
	      CMAT(I,N)=CMAT(I,N)-(1.5_LDP+HDKT*EDGE(J)/POPS(N))*CMAT(I,J)
	    END DO
	  END DO
	  SAV_CMAT=CMAT
	  CALL WR2D_MA(CMAT,N,N,'C_MAT_D61',96)
	END IF
!
! Perform the LU decomposition using DGETRF. We first equilibrize the matrix
! using DGEEQU sot the the maximum row and column values are approximately
! unity.
!
	WRITE(6,*)'Normalizing the matrix to get more accurate solution'
	CALL DGEEQU(N,N,CMAT,N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	WRITE(6,*)'Normlalization successful'
        DO J=1,N
          STEQ(J)=STEQ(J)*ROW_SF(J)
          DO I=1,N
            CMAT(I,J)=CMAT(I,J)*ROW_SF(I)*COL_SF(J)
          END DO
        END DO
!
	WRITE(6,*)'Performing LU decomposition'
        CALL DGETRF(N,N,CMAT,N,IPIVOT,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error in solution using DGETRF --A1'
          WRITE(6,*)'IFAIL=',IFAIL
	  GOTO 500
	END IF
!
! Now perform the solution.
!
        WRITE(6,*)'Do the back substituton to get solution'
	CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,STEQ,N,IFAIL)
        DO J=1,N
          STEQ(J)=STEQ(J)*COL_SF(J)
        END DO
!
! Check to see whether solution is accurate, by computing RHS with
! the newly determined solutions.
!
	RHS=0.0_LDP
	LARGEST_VAL=1.0E-42_LDP
	DO J=1,N
	  DO I=1,N
	    T1=SAV_CMAT(I,J)*STEQ(J)
	    LARGEST_VAL(I)=MAX(LARGEST_VAL(I),ABS(T1))
	    RHS(I)=RHS(I)+T1
	    IF(I .EQ. NT)WRITE(17,'(2I5,4ES14.4)')I,J,STEQ(J),SAV_CMAT(I,J),SAV_CMAT(I,J)/POPS(J),T1
	  END DO
	END DO
!
	T1=0.0_LDP
	DO IEQ=335,391
	DO J=335,394     !1,N
	  IF(SAV_STEQ(IEQ) .NE. 0)THEN
	     WRITE(125,'(I5,3ES30.16)')J,STEQ(J),SAV_CMAT(IEQ,J),SAV_CMAT(IEQ,J)*STEQ(J)/SAV_STEQ(IEQ)
!	     WRITE(125,'(I5,3ES30.16)')J,STEQ(J),SAV_CMAT(IEQ,J)/SAV_STEQ(IEQ),SAV_CMAT(IEQ,J)*STEQ(J)/SAV_STEQ(IEQ)
	     T1=T1+SAV_CMAT(IEQ,J)*STEQ(J)/SAV_STEQ(IEQ)
	   ELSE
	     WRITE(125,'(I5,3ES30.16)')J,STEQ(J)
	  END IF
	END DO
	WRITE(125,*)'Sum is for equation',IEQ,'is',T1
	END DO
!
	DO I=1,N
	  WRITE(15,'(I5,3ES16.8)')I,STEQ(I),SAV_CMAT(1,I),STEQ(I)*SAV_CMAT(1,I)
	  WRITE(16,'(I5,3ES16.8)')I,STEQ(I),SAV_CMAT(N,I),STEQ(I)*SAV_CMAT(N,I)
	END DO
!
	IF(USE_DC)THEN
	  DI_NEW=POPS(N-2)*(1.0_LDP-STEQ(N-2))
	  ED_NEW=POPS(N-1)*(1.0_LDP-STEQ(N-1))
	  T_NEW=POPS(N)*(1.0_LDP-STEQ(N))
	  RGU=LOG(2.07078E-22_LDP)
          X=HDKT/T_NEW
	  Y=ED_NEW*DI_NEW*( T_NEW**(-1.5_LDP) )/GION
          DO I=1,N-3
            NEW_LTE(I)=STAT_WT(I)*Y*EXP(EDGE(I)*X+RGU)
	  END DO
	  DO I=1,N-3
	    NEW_POPS(I)=(POPS(I)/PLTE(I))*(1.0_LDP-STEQ(I))*NEW_LTE(I)
	  END DO
	  NEW_POPS(N-2)=DI_NEW
	  NEW_POPS(N-1)=ED_NEW
	  NEW_POPS(N)=T_NEW
	  DO I=1,N
	    WRITE(6,'(1X,I3,4ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I),
	1                          1.0D0-NEW_POPS(I)/POPS(I)
	  END DO
	
	ELSE
!
! Output solutions and checks for digestion.
!
	  IF(N .LT. 20)THEN
	    DO I=1,N
	      WRITE(6,'(1X,I3,3ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I)
	    END DO
	  END IF
	  WRITE(12,'(1X,A,5(4X,A))')
	1          'Depth','    Sol   ','Old RHS   ','Eval. RHS ','Lgest Term','Error     '
	  DO I=1,N
	    WRITE(12,'(1X,I5,5ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I),
	1             LARGEST_VAL(I),ABS(SAV_STEQ(I)-RHS(I))/LARGEST_VAL(I)
	  END DO
	END IF
	WRITE(6,*)'Check fort.12 for accuracy check'
!
	DO J=1,N
	  POP_SOLS(J,1)=POPS(J)*(1.0_LDP-STEQ(J))
	END DO
!
	I=164; T2=0.0_LDP
	DO J=1,N
	  T1=SAV_CMAT(I,J)*STEQ(J)
	  T2=T2+T1
	  WRITE(175,'(I6,6ES15.4E4)')J,STEQ(J),SAV_CMAT(I,J),T1,T2,POPS(J),POP_SOLS(J,1)
	END DO
	WRITE(175,*)I,'RHS(I)=',RHS(I)
	FLUSH(UNIT=175)
!
! The following computes the residual to the equations, and then solves
! for the values needed to make the residuals zero. The solutions
! are the  updated.
!
        DO J=1,N
          RHS(J)=(SAV_STEQ(J)-RHS(J))*ROW_SF(J)
	END DO
	CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,RHS,N,IFAIL)
        DO J=1,N
          STEQ(J)=STEQ(J)+RHS(J)*COL_SF(J)
        END DO
!
	RHS=0.0_LDP
	LARGEST_VAL=1.0E-42_LDP
	DO J=1,N
	  DO I=1,N
	    T1=SAV_CMAT(I,J)*STEQ(J)
	    LARGEST_VAL(I)=MAX(LARGEST_VAL(I),ABS(T1))
	    RHS(I)=RHS(I)+T1
	  END DO
	END DO
	WRITE(12,*)' '
	WRITE(12,*)'Improved solution after 1 iteration on the residuals'
	WRITE(12,*)' '
	WRITE(12,'(1X,A,5(4X,A))')
	1          'Depth','    Sol   ','Old RHS   ','Eval. RHS ','Lgest Term','Error     '
	DO I=1,N
	  WRITE(12,'(1X,I5,5ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I),
	1             LARGEST_VAL(I),ABS(SAV_STEQ(I)-RHS(I))/LARGEST_VAL(I)
	END DO
!
	DO J=1,N
	  POP_SOLS(J,2)=POPS(J)*(1.0_LDP-STEQ(J))
	END DO
!
! Try a different scaling of CMAT to see if this makes any difference to the
! soluton obtained.
!
	CMAT=SAV_CMAT
	STEQ=SAV_STEQ
	DO I=1,N
	  OLD_SOL(I)=MAX(1.0_LDP,ABS(STEQ(I)))
!	  OLD_SOL(I)=2.0D0
	END DO
!
	DO J=1,N
	  DO I=1,N
	    CMAT(I,J)=CMAT(I,J)*OLD_SOL(J)
	  END DO
	END DO
	CALL DGEEQU(N,N,CMAT,N,ROW_SF,COL_SF,ROW_CND,COL_CND,MAX_VAL,IFAIL)
        DO J=1,N
          STEQ(J)=STEQ(J)*ROW_SF(J)
          DO I=1,N
            CMAT(I,J)=CMAT(I,J)*ROW_SF(I)*COL_SF(J)
          END DO
        END DO
	WRITE(6,*)'Performing LU decomposition'
        CALL DGETRF(N,N,CMAT,N,IPIVOT,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error ins solution using DGETRF'
	END IF
!
! Now perform the solution.
!
        WRITE(6,*)'Do the back substituton to get solution'
	CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,STEQ,N,IFAIL)
        DO J=1,N
          STEQ(J)=STEQ(J)*COL_SF(J)*OLD_SOL(J)
        END DO
!
	DO J=1,N
	  POP_SOLS(J,3)=POPS(J)*(1.0_LDP-STEQ(J))
	END DO
!
! Check to see whether solution is accurate, by computing RHS with
! the newly determined solutions.
!
	RHS=0.0_LDP
	LARGEST_VAL=1.0E-42_LDP
	DO J=1,N
	  DO I=1,N
	    T1=SAV_CMAT(I,J)*STEQ(J)
	    LARGEST_VAL(I)=MAX(LARGEST_VAL(I),ABS(T1))
	    RHS(I)=RHS(I)+T1
	  END DO
	END DO
	WRITE(12,*)' '
	WRITE(12,*)'Solution using a different scaling of CMAT'
	WRITE(12,*)' '
	WRITE(12,'(1X,A,5(4X,A))')
	1          'Depth','    Sol   ','Old RHS   ','Eval. RHS ','Lgest Term','Error     '
	DO I=1,N
	  WRITE(12,'(1X,I5,5ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I),
	1             LARGEST_VAL(I),ABS(SAV_STEQ(I)-RHS(I))/LARGEST_VAL(I)
	END DO
!
!	DO I=1,NT
!	  WRITE(30,'(1X,I5,3ES14.4)')I,STEQ(I),STEQ(I)*CMAT_RD(NT,I)
!	END DO
!
500	CONTINUE
!
! Solve directly for the populatons.
!
	CMAT=SAV_CMAT
	STEQ=SAV_STEQ
!
	DO J=1,N
	  DO I=1,N
	    CMAT(I,J)=CMAT(I,J)/POPS(J)
	  END DO
	END DO
	IF(N .EQ. NT)THEN
	   CMAT(:,N-1)=0.0_LDP
	   CMAT(:,N)=0.0_LDP
	   CMAT(N-1,N-1)=1.0_LDP;
	   CMAT(N,N)=1.0_LDP;
	END IF
	WRITE(36,*)CMAT(N-1,1:N)
	WRITE(36,*)CMAT(N,1:N)
	CALL DGEEQU(N,N,CMAT,N,ROW_SF,COL_SF,ROW_CND,COL_CND,MAX_VAL,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error in solution using DGEEQU'
	END IF
        DO J=1,N
          STEQ(J)=STEQ(J)*ROW_SF(J)
          DO I=1,N
            CMAT(I,J)=CMAT(I,J)*ROW_SF(I)*COL_SF(J)
          END DO
        END DO
	WRITE(6,*)'Performing LU decomposition'
        CALL DGETRF(N,N,CMAT,N,IPIVOT,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error in solution using DGETRF'
	END IF
!
! Now perform the solution.
!
        WRITE(6,*)'Do the back substituton to get solution'
	CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,STEQ,N,IFAIL)
        DO J=1,N
          STEQ(J)=STEQ(J)*COL_SF(J)
        END DO
!
	DO J=1,N
	  POP_SOLS(J,4)=POPS(J)-STEQ(J)
	END DO
!
	DO J=1,N
	  WRITE(35,'(I5,5ES16.6)')J,POPS(J),(POP_SOLS(J,I),I=1,4)
	END DO
!
	CMAT=SAV_CMAT
	IF(N .EQ. NT)THEN
	   CMAT(:,N-1)=0.0_LDP
	   CMAT(:,N)=0.0_LDP
	   CMAT(N-1,N-1)=1.0_LDP;
	   CMAT(N,N)=1.0_LDP;
	END IF

	COL_SF=0.0_LDP
	DO I=1,N
	  DO J=1,N
	    COL_SF(I)=COL_SF(I)+CMAT(I,J)
	    WRITE(38,*)I,J,COL_SF(I)
	  END DO
	END DO
	DO J=1,N
	  WRITE(37,'(I5,5ES16.6)')J,POPS(J),SAV_STEQ(J),COL_SF(J),COL_SF(J)-SAV_STEQ(J)
	END DO
!
	STOP
	END
