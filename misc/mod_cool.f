!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc. 
!
	PROGRAM MOD_COOL
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: MAX_RECS=1000
!
	CHARACTER*132 TMP_STR
	CHARACTER*132 STRING
	CHARACTER*132 STR_VEC(MAX_RECS)
	REAL*8 VALS(MAX_RECS,10)
	REAL*8 TA(MAX_RECS)
	INTEGER INDX(MAX_RECS)
!
	REAL*8, ALLOCATABLE :: BOUND(:)
	REAL*8, ALLOCATABLE :: SUM(:)
!
	INTEGER*4 ND
	INTEGER*4 NV
	INTEGER*4 I,J,K,ID
	INTEGER N_INIT_RECS
	INTEGER NRECS
	REAL*8 T1
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	WRITE(6,'(A)',ADVANCE='NO')'INPUT ND: '
	READ(5,*)ND
!
	ALLOCATE (SUM(ND))
	ALLOCATE (BOUND(ND))
!
	OPEN(UNIT=20,FILE='GENCOOL',STATUS='OLD',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='WRITE')
!
	DO I=1,1+(ND-1)/10
	   ID=1+(I-1)*10
!
	   READ(20,'(A)')STRING
	   WRITE(21,'(A)')TRIM(STRING)
	   DO J=1,3
	     READ(20,'(A)')TMP_STR
	     READ(20,'(A)')STRING
	     IF(J .EQ. 1)WRITE(21,'(A,T11,10I12)')'Depth',(K,K=ID,MIN(ID+9,ND))
	     WRITE(21,'(A)')TMP_STR(4:9)//'   '//TRIM(STRING)
	     READ(20,'(A)')STRING
	   END DO
!
	   T1=0.0D0
	   DO WHILE(1 .EQ. 1)
	      READ(20,'(A)')STRING
	      IF( INDEX(STRING,'Coll') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T10,A)')TMP_STR(1:K)//'COL ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Free-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T10,A)')TMP_STR(1:K)//'FF ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'K-shell') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T10,A)')TMP_STR(1:K)//'XKS ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Rate') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T10,A)')'Net C.R.',TRIM(STRING)
	     ELSE IF( INDEX(STRING,'Net') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T10,A)')'% C.R.',TRIM(STRING)
	        IF(I .NE. 1+(ND-1)/10)THEN
	          READ(20,'(A)')STRING
	          WRITE(21,'(A)')STRING
	        END IF
	        EXIT
	      ELSE IF( INDEX(STRING,'Bound-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        SUM=0.0D0
	        DO WHILE(1 .EQ. 1)
	          READ(20,'(A)')STRING
	          IF(STRING .EQ. ' ')THEN
	            TMP_STR=ADJUSTL(TMP_STR)
	            K=INDEX(TMP_STR,' ')
	            WRITE(21,'(A)')STRING
                    WRITE(21,'(A,T11,10ES12.4)')TMP_STR(1:K)//'BF ',(SUM(K),K=ID,MIN(ID+9,ND))
	            EXIT
	          ELSE
	            READ(STRING,*)(BOUND(K),K=ID,MIN(ID+9,ND))
	            SUM=SUM+BOUND
	          END IF
	        END DO
	      END IF
	   END DO
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'Cooling data has been written to GENCOOL_SUM'
	WRITE(6,'(A)')'Will now sort data to display most important terms'
	WRITE(6,'(A)')' '
!
	OPEN(UNIT=20,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SORT',STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(6,'(A)',ADVANCE='NO')'Input maximum number of recors to be output to sorted file: '
	READ(5,*)NRECS
!
	N_INIT_RECS=5
	VALS=0.0D0
	DO I=1,1+(ND-1)/10
	  ID=1+(I-1)*10
	  DO K=1,5
	    READ(20,'(A)')STR_VEC(K)
	  END DO
!
	  K=N_INIT_RECS
	  DO WHILE(1 .EQ. 1)
	    K=K+1
100	    READ(20,'(A)')STR_VEC(K)
	    IF(STR_VEC(K) .EQ. ' ')GOTO 100
	    READ(STR_VEC(K)(12:),*)(VALS(K,J),J=1,MIN(10,ND-(I-1)*10))
	    IF( INDEX(STR_VEC(K),'%') .NE. 0)THEN
	      IF(I .NE. 1+(ND-1)/10 )READ(20,'(A)')STRING			!Getting record with ^L
	      EXIT
	    END IF
	  END DO
!
	  NV=K-N_INIT_RECS-2
	  TA(:)=ABS(VALS(:,1))
	  CALL INDEXX(NV,TA(6),INDX,L_FALSE)
	  NV=NV+N_INIT_RECS+2
!
	  DO K=1,N_INIT_RECS
	    WRITE(21,'(A)')TRIM(STR_VEC(K))
	  END DO
	  WRITE(21,'(A)')' '
	  WRITE(21,'(A)')TRIM(STR_VEC(NV))
	  WRITE(21,'(A)')TRIM(STR_VEC(NV-1))
	  WRITE(21,'(A)')' '
	  DO K=1,MIN(NRECS,NV-N_INIT_RECS-2)
	    WRITE(21,'(A)')TRIM(STR_VEC(INDX(K)+N_INIT_RECS))
	  END DO
!
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'Sorted cooling data has been written to GENCOOL_SORT'
	WRITE(6,'(A)')' '
!
	STOP
	END
