!
! Simple program to plot data from the SCRTEMP file. This can be used to check
!   Convergence of a proram.
!   Convergence as a function of depth etc.
!
	PROGRAM PLT_SCR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
	INTEGER*4 ND,NT,NIT
C
	REAL*8, ALLOCATABLE :: POPS(:,:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
C
	REAL*8, ALLOCATABLE :: X(:)			!NIT
	REAL*8, ALLOCATABLE :: Y(:)			!NIT
	REAL*8, ALLOCATABLE :: Z(:)			!NIT
!
	INTEGER*4, ALLOCATABLE :: I_BIG(:)		!NT
	REAL*8, ALLOCATABLE :: Z_BIG(:)			!NT
C
	INTEGER*4, PARAMETER :: T_OUT=6
C
	INTEGER*4 NPLTS
	INTEGER*4 IREC
	INTEGER*4 IVAR
	INTEGER*4 J
	INTEGER*4 ID
	INTEGER*4 IT
	INTEGER*4 NY
	INTEGER*4 NITSF
	INTEGER*4 LST_NG
	INTEGER*4 RITE_N_TIMES
	INTEGER*4 LUSCR
	INTEGER*4 NUM_TIMES
	INTEGER*4 K
	INTEGER*4 IOS
C
	LOGICAL NEWMOD
	LOGICAL DO_ABS
	CHARACTER*10 PLT_OPT
	CHARACTER*80 YLABEL
	CHARACTER*132 STRING
C
	REAL*8 T1,T2,T3
C
	LUSCR=26
	NUM_TIMES=1
	NEWMOD=.TRUE.
C
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       POINT1.DAT'
	WRITE(T_OUT,*)'                                       SCRTEMP.DAT'
	WRITE(T_OUT,*)'                                       MODEL.DAT'
	WRITE(T_OUT,*)' '
C
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
C
100	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	END IF 
	CLOSE(UNIT=12)
C
	OPEN(UNIT=12,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)READ(12,*,IOSTAT=IOS)K,NIT
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Possible error reading POINT1'
	    CALL GEN_IN(NIT,'Number of iterations')
	  END IF
	CLOSE(UNIT=12)
C
	ALLOCATE (POPS(NT,ND,NIT))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
C
	J=MAX(NIT,ND); J=MAX(NT,J)
	ALLOCATE (X(J))
	ALLOCATE (Y(J))
	ALLOCATE (Z(J))
!
	WRITE(6,*)'   Number of depth points is:',ND
	WRITE(6,*)'Number of variables/depth is:',NT
	WRITE(6,*)'     Number of iterations is:',NIT
!
	ALLOCATE (I_BIG(NT))
	ALLOCATE (Z_BIG(NT))
C
	DO IREC=1,NIT
	    CALL SCR_READ(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,
	1              NT,ND,LUSCR,NEWMOD)
	END DO
C
200	CONTINUE
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'ND  :: ',ND
	WRITE(T_OUT,*)'NT  :: ',NT
	WRITE(T_OUT,*)'NIT :: ',NIT
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'F   :: Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)'
	WRITE(T_OUT,*)'R   :: [Y(K+2)-Y(K+1)]/[Y(K+1)-Y(K)]'
	WRITE(T_OUT,*)'D   :: Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)'
	WRITE(T_OUT,*)'Y   :: Z(K)=Y(K)'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'MR  :: Z(K)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)'IR  :: Z(ID)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'E   :: EXIT'
	WRITE(T_OUT,*)' '
	PLT_OPT='R'
	CALL GEN_IN(PLT_OPT,'Plot option: R(atio), F(rac) or D(elta), Y, IR, MR or E(xit)')
	CALL SET_CASE_UP(PLT_OPT,0,0)
!
	IF(PLT_OPT(1:2) .EQ. 'E ' .OR. PLT_OPT(1:2) .EQ. 'EX')STOP
	IF(PLT_OPT(1:5) .EQ. 'MED_R')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(ABS(T2). GT. 1.0D-50*ABS(T1))THEN
	        Z_BIG(IVAR)=T1/T2
	      ELSE
	        Z_BIG(IVAR)=20000
	      END IF
	    END DO
!
! Now find median value. First order.
!
	    CALL INDEXX(NT,Z_BIG,I_BIG,.TRUE.)
	    J=(NT+1)/2
	    Y(ID)=100.0D0*(Z_BIG(I_BIG(J))-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(r\dmed\u-1)'
!	    IF(ID .EQ. 30)THEN
	      WRITE(74,*)'ID=',ID
	      DO J=1,NT
	        WRITE(74,*)J,I_BIG(J),Z_BIG(I_BIG(J))
	      END DO
!	    END IF
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'MR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO_ABS=.TRUE.; CALL GEN_IN(DO_ABS,'Use absolute value')
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        T1=T1/T2
	        IF(DO_ABS)T1=ABS(T1)
	        Y(ID)=Y(ID)+ T1
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR
	      END IF
	    END DO
	    Y(ID)=100.0D0*(Y(ID)/K-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'HR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT
	    CALL GEN_IN(IT,'Iteration #')
	  END DO 
	  ID=0
	  DO WHILE(ID .LT. 1 .OR. ID .GT. ND)
	    ID=ND/2; CALL GEN_IN(ID,'Depth index')
	  END DO
	  Y(ID)=0.0D0
!
	  K=MIN(201,NT)
	  T3=200.0D0/(K-1)
	  DO J=1,K
	    X(J)=-50.0D0+(J-1)*T3
	  END DO
	  Y(1:K)=0.0D0
!
	  DO IVAR=1,NT
	    T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	    T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	    IF(T2 .NE. 0)THEN
	      T1=(100*(T1/T2-1.0D0)+51)/T3
	      IF(T1 .LT. 1)T1=1
	      IF(T1 .GT. K)T1=K
	      Y(NINT(T1))=Y(NINT(T1))+1
	    ELSE
	      Y(1)=Y(1)+1
	    END IF
	   Ylabel='N(r)'
	  END DO
	  CALL DP_CURVE(K,X,Y)
	  CALL GRAMON_PGPLOT('100(r-1)',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'IR')THEN
	  ID=0
	  DO WHILE(ID .LT. 1 .AND. ID .GT. ND)
	    ID=ND; CALL GEN_IN(ID,'Depth index')
	  END DO
	  DO IT=3,NIT
	    Y(IT)=0.0D0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        Y(IT)=Y(IT)+ T1/T2
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR,'for iteration',IT
	      END IF
	    END DO
	    Y(IT)=100.0D0*(Y(IT)/NT-1.0D0)
	    X(IT)=FLOAT(IT)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  IT=NIT-2
	  CALL DP_CURVE(IT,X(3),Y(3))
	  CALL GRAMON_PGPLOT('Iteration',Ylabel,' ',' ')
	  GOTO 200
	END IF
C
	DO WHILE(0 .EQ. 0)	
	  NPLTS=0
	  IVAR=1
	  DO WHILE(IVAR .NE. 0)
500	    IVAR=0
	    WRITE(STRING,'(I5,A)')NT,'](0 to plot)'
	    DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:) ; END DO
	    STRING='Variable to be plotted ['//STRING
	    CALL GEN_IN(IVAR,STRING)
	    IF(IVAR .LT. 0 .OR. IVAR .GT. NT)GO TO 500
	    IF(IVAR .EQ. 0)GOTO 1000
C
600	    CONTINUE
	    WRITE(STRING,'(I5,A)')ND,'](0 to plot)'
	    DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:); END DO
	    STRING='Depth of variable to be plotted ['//STRING
	    CALL GEN_IN(ID,STRING)
	    IF(ID .LE. 0 .OR. ID .GT. ND)GO TO 600
C
	    Y(1:NIT)=POPS(IVAR,ID,1:NIT)
C
	    IF(PLT_OPT(1:1) .EQ. 'F')THEN
	      DO K=1,NIT-1
	        Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)
	        X(K)=FLOAT(K)
	      END DO
	      NY=NIT-1
	      T1=MAXVAL(ABS(Z(1:NY)))
	      IF(T1 .LT. 1.0E-02)THEN
	        Z(1:NY)=Z(1:NY)*1.0D+03
	        YLABEL='\gDY/Y(%)\d \ux10\u3\d'
	        WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
	      ELSE
	        YLABEL='\gDY/Y(%)'
	      END IF
!	      
	    ELSE IF(PLT_OPT(1:1) .EQ. 'R')THEN
	      DO K=1,NIT-2
	        T1=Y(K+2)-Y(K+1)
	        T2=Y(K+1)-Y(K)
	        IF(T2 .NE. 0)THEN
	           Z(K)=T1/T2
	        ELSE
	           Z(K)=10.0
	        END IF
	        X(K)=FLOAT(K)+2
	      END DO
	      NY=NIT-2
	      YLABEL='\gDY(K+1)/\gDY(K)'
	    ELSE IF(PLT_OPT(1:1) .EQ. 'D')THEN
	      DO K=1,NIT
	        Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)
	        X(K)=FLOAT(K)
	      END DO
	      NY=NIT-1
	      T1=MAXVAL(ABS(Z(1:NY)))
	      IF(T1 .LT. 1.0E-02)THEN
	        Z(1:NY)=Z(1:NY)*1.0D+03
	        YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]\d \ux10\u3\d'
	        WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
	      ELSE
	        YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]'
	      END IF
	    ELSE IF(PLT_OPT(1:1) .EQ. 'Y')THEN
	      DO K=1,NIT
	        Z(K)=Y(K)
	        X(K)=FLOAT(K)
	      END DO
	      NY=NIT
	      YLABEL='Y(K) [%]'
	    ELSE
	      WRITE(T_OUT,*)' Option not recognized: try again'
	      GOTO 1000
	    END IF
	    CALL DP_CURVE(NY,X,Z)
	    NPLTS=NPLTS+1
	  END DO
1000	  CONTINUE
	  IF(NPLTS .NE. 0)THEN
	    CALL GRAMON_PGPLOT('Iteration number K',Ylabel,' ',' ')
	  ELSE
	    GOTO 200
	  END IF
	END DO
C
	END
C
C 
C
C Altered  12-Jan-1991 - By using call to DIR_ACC_PARS this version is now
C                        compatible with both CRAY and VAX fortran.
C Altered  3-Apr-1989 - LST_NG installed. LU now transmitted in call.
C
	SUBROUTINE SCR_READ(R,V,SIGMA,POPS,IREC,NITSF,
	1                        NUM_TIMES,LST_NG,NT,ND,LU,NEWMOD)
	IMPLICIT NONE
C
	LOGICAL NEWMOD
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the nuber of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the numer of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER*4 ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,L,LUER,ERROR_LU
	INTEGER*4 IOS,IST,IEND
	EXTERNAL ERROR_LU
C
	NEWMOD=.FALSE.
	LUER=ERROR_LU()
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C These are computer dependent, hence call to DIR_ACC_PARS. NB.
C REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
C passed to the OPEN statement.

	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length is too small to output R,V and'//
	1              ' sigma in SCR_RITE'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
C
C		' OLD MODEL '
C
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',
	1       RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    CLOSE(UNIT=LU)
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    RETURN
	  END IF
C
C Note that NITSF= # of successful iterations so far.
C
	  READ(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)READ(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading R,V, SIGMA vectors in READ_SCRTEMP'
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    RETURN
	  END IF
	  RECS_FOR_RV=2
C
C Read in the population data.
C
500	CONTINUE		!Try to read an earlier record.
C
C IREC ignores the number of records that it takes to write each time.
C Hence in POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
	ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	DO L=1,NUMRECS
	  IST=(L-1)*N_PER_REC+1
	  IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	  READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error on Scratch Read'
	    IREC=IREC-1
C	    IF(IREC .GE. 0)GOTO 500	!Get another record if one is available.
	    NEWMOD=.TRUE.
	    NITSF=0
	    IREC=0
	    LST_NG=-1000
	    CLOSE(UNIT=LU)
	    RETURN
	  END IF
	END DO
C
C Successful Read !
C
	CLOSE(UNIT=LU)
	RETURN
C
	END
C
C 
C
C Routine  to save population data. There is no limit on the
C size of the POPS array.
C
	SUBROUTINE SCR_RITE(R,V,SIGMA,POPS,
	1                      IREC,NITSF,NUM_TIMES,LST_NG,
	1                      NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
C
	LOGICAL WRITFAIL
	INTEGER*4 IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
C
C Local variables.
C
C REC_SIZE     is the (maximum) record length in bytes.
C REC_LEN      is the record length in computer units.
C UNIT_SIZE    is the nuber of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the numer of bytes used to represent the number.
C NUMRECS      is the # of records required to output POPS.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 REC_SIZE,REC_LEN
	INTEGER*4 NUMRECS
	INTEGER*4 N_PER_REC
	INTEGER*4 ARRAYSIZE  	!Size of POPS array.
C	
	INTEGER*4 ST_REC_M1,RECS_FOR_RV
	INTEGER*4 I,K,L,LUER,ERROR_LU
	INTEGER*4 IOS,IST,IEND
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C These are computer dependent, hence call to DIR_ACC_PARS. NB.
C REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
C passed to the OPEN statement.
C
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length is too small to output R,V and'//
	1              ' sigma in SCR_RITE'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND			
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
C
C*************************************************************************
C
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP in WRITE_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      WRITFAIL=.TRUE.
	      CLOSE(UNIT=LU)
	      RETURN
	    ELSE
	      IF(IOS .EQ. 0)IREC=0		!Since new file.
	    END IF
	  END IF
C
	IF(IREC .EQ. 0)THEN		!Newfile or newmodel
	  WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	    WRITFAIL=.TRUE.
	    RETURN
	  END IF
	END IF
	RECS_FOR_RV=2
C
C WRITE in the population data.
C
	IREC=IREC+1		!Next record output
	DO K=1,NUM_TIMES	!Write out POPS NUM_TIMES for safety.
C
C IREC ignores the number of records that it takes to write each time.
C Hence in POINT it will correspond to the iteration number.
C ST_REC_M1 + 1 is the first output record.
C
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)THEN
	      WRITFAIL=.TRUE.
	      WRITE(LUER,*)'Error writing SCRTEMP in SCR_RITE'
	      WRITE(LUER,*)'IOS=',IOS
	      CLOSE(UNIT=LU)
	      RETURN
	    END IF
	  END DO
	END DO
C
C Successful write.
C
	WRITFAIL=.FALSE.
1000	CLOSE(UNIT=LU)
C
C Write pointer to data files. 
C
	OPEN(UNIT=LU,FILE='POINT1',STATUS='UNKNOWN')
	  WRITE(LU,'(X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	CLOSE(UNIT=LU)
	OPEN(UNIT=LU,FILE='POINT2',STATUS='UNKNOWN')
	  WRITE(LU,'(X,4(I6,4X))')
	1          IREC,NITSF,NUM_TIMES,LST_NG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A)')
	1          'IREC','NITSF','#_TIMES','LST_NG'
	CLOSE(UNIT=LU)
C
	RETURN
	END
