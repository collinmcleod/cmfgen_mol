C
C Subroutine to determine those frequencies at which the continuum opacity 
C will be evaluated. At other frequenecies it can be held fixed.
C
	SUBROUTINE DET_MAIN_CONT_FREQ(NU,NCF,NU_CONT,NCF_CONT,
	1                NU_EVAL,DOPV,DELV_CONT,COMPUTE_ALL_CROSS)
	IMPLICIT NONE
C
C NCF represents the toal number of frequencies at which we solve the radiative
C transfer equation. The frequencies are stored in the vector NU which
C has been chosen to sample all lines and continuum edges adequately.
C
	INTEGER*4 NCF
	REAL*8 NU(NCF)
C
C For each frequncy NU(I), the continuum cross-section will be evaluated at
C freqency NU_EVAL(I).
C
	REAL*8 NU_EVAL(NCF)
C
C NCF_CONT represents the number of continuum frequencies before line 
C insertion. NU_CONT contains these frequencies, and has been chosen to 
C sample continuum cross-sections (with allowance for level dissolution and 
C important bound-free edges) adequately. 
C
	INTEGER*4 NCF_CONT
	REAL*8 NU_CONT(NCF_CONT)
C
	REAL*8 DOPV		!Doppler spacing across lines
C
C DELV_CONT is the maximum sparation between points at which the continuum
C   opacity is evaluated.
C
	REAL*8 DELV_CONT
C
	LOGICAL COMPUTE_ALL_CROSS
C
	REAL*8 SPEED_OF_LIGHT
	INTEGER*4 ERROR_LU
	EXTERNAL SPEED_OF_LIGHT,ERROR_LU
C
C Local variables:
C
	INTEGER*4 LU_OUT,L,K,ML,ML_ST,ML_END,LST_COMP
	REAL*8 T1,T2,T3,T4,DOP_RAT,VRAT,C_KMS
C
	LU_OUT=ERROR_LU()
	IF(COMPUTE_ALL_CROSS .OR. DELV_CONT .EQ. 0)THEN
	  NU_EVAL(:)=NU(:)
	  WRITE(LU_OUT,'(A)')'The continuum will be evaluated at all',
	1                        ' frequencies.'
	  RETURN
	END IF
C
	C_KMS=SPEED_OF_LIGHT()/1.0D+05
C
C We set NU_EVAL to zero. Subsequently a NU_EVAL value of zero for any index
C indicates that NU_EVAL still needs to be set.
C
	NU_EVAL(:)=0.0D0
C
C We first ensure that the continuum opacity is evaluated at all important
C bound_free edges etc. These are located in NU_CONT. We use DOPV to do this
C since all continuum edges should be in NU unless they were in DOPV km/s
C of an inserted line frequency.
C
	K=1
	DOP_RAT=DOPV/C_KMS
	DO ML=1,NCF_CONT
	  DO WHILE( (NU(K)-NU_CONT(ML))/NU_CONT(ML) .GT. DOP_RAT)
	    K=K+1
	  END DO
	  L=K
	  DO WHILE( ABS( (NU_CONT(ML)-NU(L))/NU_CONT(ML)) .LE. DOP_RAT)
	    NU_EVAL(L)=NU(L) 
	    L=L+1
	  END DO
	END DO
	NU_EVAL(1)=NU(1)
	NU_EVAL(NCF)=NU(NCF)
C
C Now we set intermediate frequencies, such that the continuum opacities
C and emissivities etc are evaluated at least every DELV_CONT km/s/
C
	ML_ST=1
	LST_COMP=1
	VRAT=DELV_CONT/C_KMS
	DO WHILE(ML_ST .LT. NCF)
	  DO WHILE(NU_EVAL(ML_ST+1) .NE. 0.0D0)
	    ML_ST=ML_ST+1
	    IF(ML_ST .EQ. NCF)GOTO 100
	  END DO
	  ML_END=ML_ST+1
	  DO WHILE(NU_EVAL(ML_END) .EQ. 0.0D0)
	    ML_END=ML_END+1
	  END DO
	  DO L=ML_ST,ML_END-1
	    T1=(NU(LST_COMP)-NU(L))/NU(L) 
	    T2=(NU(L)-NU(ML_END))/NU(L)
	    IF(T1+T2 .GT. 1.25D0*VRAT .AND. T1 .GT. 0.75D0*VRAT)THEN
	      NU_EVAL(L)=NU(L)
	      LST_COMP=L
	    ELSE IF(T1 .LT. T2)THEN
	      NU_EVAL(L)=NU(LST_COMP)
	    ELSE
	      NU_EVAL(L)=NU(ML_END)
	    END IF
	  END DO
	  ML_ST=ML_END
	END DO
C
100	CONTINUE
	K=0
	DO ML=1,NCF
	  IF(NU(ML) .EQ. NU_EVAL(ML))K=K+1
	END DO
C
	WRITE(LU_OUT,'(A,I5,A)')' The continuum will be evaluated at ',
	1                        K,' frequencies'
	WRITE(LU_OUT,'(A)')' '
!
	OPEN(UNIT=17,FILE='CONT_FREQ',STATUS='UNKNOWN')
	   WRITE(17,*)NCF,'        !NCF'
	   T1=NU(1)
	   DO ML=2,NCF
	     IF(NU_EVAL(ML) .EQ. NU(ML))THEN
	       WRITE(17,'(X,I6,3X,2ES14.5,F10.2)')ML,NU(ML),0.01D0*C_KMS/NU(ML),C_KMS*(NU(ML)-T1)/NU(ML)
	       T1=NU(ML)
	     END IF
	   END DO
!
	   WRITE(17,'(A)')' '
	   WRITE(17,'(A)')' Check on continuum evaluations relative to bound-free edges'
	   WRITE(17,'(A)')' '
	   L=1
	   T1=0.0D0; T2=0.0D0; T3=0.0D0; T4=0.0D0;
	   DO ML=2,NCF_CONT-1
	     DO WHILE(NU(L) .GT. NU_CONT(ML))
	       L=L+1
	     END DO
	     IF(L .GT. 2)T1=C_KMS*(NU_EVAL(L-2)-NU_CONT(ML))/NU_CONT(ML)
	     IF(L .GT. 1)T2=C_KMS*(NU_EVAL(L-1)-NU_CONT(ML))/NU_CONT(ML)
	     T3=C_KMS*(NU_EVAL(L)-NU_CONT(ML))/NU_CONT(ML)
	     IF(L .LT. NCF)T4=C_KMS*(NU_EVAL(L+1)-NU_CONT(ML))/NU_CONT(ML)
	     WRITE(17,'(X,I6,ES15.7,4F12.2)')ML,NU_CONT(ML),T1,T2,T3,T4
	   END DO
	   
	CLOSE(UNIT=17)
C
	RETURN
	END
