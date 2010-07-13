!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc. 
!
	PROGRAM MOD_PRRR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered: 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered: 08-Feb-2008: Extra terms (such as V term) sheck and output.
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
	REAL*8 PHOT(10,2000)
	REAL*8 RECOM(10)
!
	INTEGER*4 ND
	INTEGER*4 NV
	INTEGER*4 I,J,K
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	INTEGER IST,IEND,NLEV
	REAL*8 T1
	LOGICAL FILE_OPEN
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	OPEN(UNIT=20,FILE='FeIPRRR',STATUS='OLD',ACTION='READ')
	OPEN(UNIT=21,FILE='FeIPRRR_SUM',STATUS='UNKNOWN',ACTION='WRITE')
!
	DO I=1,1+(ND-1)/10
	   IST=1; IEND=MIN(10,ND-(I-1)*10)
!
	   DO J=1,13
	     READ(20,'(A)')TMP_STR
	     WRITE(21,'(A)')TRIM(TMP_STR)
	   END DO
	   READ(20,'(A)')TMP_STR	!Photoionization header.
!
	   NLEV=0
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)')STRING
	     IF(STRING .NE. ' ')THEN
	       NLEV=NLEV+1
	       WRITE(6,*)STRING(1:10)
	       READ(STRING,*)(PHOT(K,NLEV),K=IST,IEND)
	     ELSE
	       EXIT
	     END IF
	   END DO
!
	   DO J=1,3
	     READ(20,'(A)')STRING
	     WRITE(21,'(A)')TRIM(STRING)
	   END DO
	   WRITE(21,'(A,I4,A,I4,A)')'   Total photoionization rate (d=',(I-1)*10+1,' to',MIN(I*10,ND),'):'
	   WRITE(21,'(X,10ES12.4)')(SUM(PHOT(K,:)),K=IST,IEND)
	   WRITE(21,'(A)')
!
	   READ(20,'(A)')STRING
	   WRITE(21,'(A,A)')'   Net ',TRIM(STRING(4:))
!
	   DO J=1,NLEV
	     READ(20,*)(RECOM(K),K=IST,IEND)
	     WRITE(21,'(X,10ES12.4)')(RECOM(K)-PHOT(K,J),K=IST,IEND)
	   END DO
!
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)',END=200)STRING
	     WRITE(21,'(A)')TRIM(STRING)
	     IF(STRING(1:1) .EQ. '1')EXIT
	   END DO
	   FLUSH(21)
	END DO
200	CONTINUE
!
	STOP
	END
