C
C Program reads in departire coefficient estimates from a file and uses
C linear interpolation in the log plane to lay these estimates
C on the "new" electron grid. The electron density is
C assumed to be known on entering the program. The ion density is also
C retruned. Program will generate an error if NDOLD > 150, of if NOLD > 50.
C
C NB. The variables NST and NFIN are used so that this routine can also read
C     in the HEI arrays. N specifies the first dimension of the population
C     array. 
C           e.g. for Hydrogen read, set NST=1, and NFIN=N.
C           e.g. for HeI(Triplets) read, set NST=NSING+1, and NFIN=N.
C
	SUBROUTINE REGRID_B_ON_NE(DHEN,ED,ION,EDLOG,N,NST,NFIN,ND,
	1                         LU,FILNAME)
	IMPLICIT NONE
C
C Altered 14-Apr-2002 : LUER not being set before error output.
C Altered 07-Jul-1997 - Subroutine reads in and utilizes CLUMP_FAC, if it
C                         is available.
C Altered 26-May-1996 - Access to SCRATCH block removed.
C                       Dynamic memory allocation installed.
C                         Generical calls for LOG and EXP.
C Created 11-Feb-1987 - Based on 8-Feb-88 version of REGRID_B_ON_NE. ION density
C                       now retruned.
C
	INTEGER*4 N,NST,NFIN,ND,LU
	REAL*8 DHEN(2-NST:N-NST+1,ND),ED(ND),ION(ND),EDLOG(ND)
	CHARACTER*(*) FILNAME
C
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	REAL*8, ALLOCATABLE :: HYD(:,:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: OLDION(:)
C
C NLEV is the number of levels whose departure coefficients will be set.
C
	INTEGER*4 NLEV
C
	INTEGER*4 I,J,NX,NXST,NX_END,NZ,NOLD,NDOLD,IOS
	REAL*8 CONST,TSTAROLD,RPOLD
	REAL*8 ION_FRAC,VEL,OLD_CLUMP_FAC
	LOGICAL CLMP_PRES
	CHARACTER*80 STRING
C
C Read in values from previous model.
C
	OPEN(UNIT=LU,STATUS='OLD',FILE=FILNAME,ACTION='READ',IOSTAT=IOS)
        IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRID_B_ON_NE'
          WRITE(LUER,*)'IOS=',IOS
          STOP
        END IF
C
C Check whether the clumping factor has also been written to file. The
C mere presence of a string containing '!Format date' indicates that it has.
C
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(LU,'(A)')STRING
	END DO
	IF( INDEX(STRING,'!Format date') .NE. 0)THEN
	  CLMP_PRES=.TRUE.
	ELSE
	  CLMP_PRES=.FALSE.
	  REWIND(LU)
	END IF

	READ(LU,*)RPOLD,TSTAROLD,NOLD,NDOLD
C
C Allocate necessary memory.
C
	I=MAX(NOLD,NDOLD,N,ND)
	ALLOCATE (HYD(NOLD,NDOLD))
	ALLOCATE (OLDED(I))
	ALLOCATE (TA(I))
	ALLOCATE (TB(I))
	ALLOCATE (OLDION(I))
C
C NZ defines the number of atomic levels which can be found by
C direct interpolation
C
	NLEV=NFIN-NST+1
	NZ=NLEV
	IF(NLEV .GT. NOLD)NZ=NOLD
C
C TA is used for both R and T but note that they
C are not used and hence may be overwritten.
C
	DO I=1,NDOLD
	  IF(CLMP_PRES)THEN
	    READ(LU,*)TA(I),OLDION(I),OLDED(I),TA(I),
	1               ION_FRAC,VEL,OLD_CLUMP_FAC
	    OLDED(I)=LOG(OLDED(I)*OLD_CLUMP_FAC)
	  ELSE
	    READ(LU,*)TA(I),OLDION(I),OLDED(I),TA(I)
	    OLDED(I)=LOG(OLDED(I))
	  END IF
	  READ(LU,*)(TA(J),J=1,NOLD)
	  DO J=1,NZ
	    HYD(J,I)=TA(J)
	  END DO
	  OLDION(I)=LOG(OLDION(I))
	END DO
C
C Check to see if departure coefficients, or one minus departure
C coefficients. Assumes that at depth, the highest level should
C have b=1.
C
	IF( ABS((HYD(NZ,NDOLD)-1.0D0)) .LT. 0.1)THEN
	  CONST=0.0D0
	ELSE IF( ABS(HYD(NZ,NDOLD)) .LT. 0.1) THEN
	  CONST=1.0D0
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Possible error on REGRID_B_ON_NE - b(n) not equal'
	  WRITE(LUER,*)'to 1.0/pm0.1 or 0.0/pm0.1 at depth'
	END IF
C
C Interpolations are performed in the log plane.
C
	DO I=1,ND
	  EDLOG(I)=LOG(ED(I))
	END DO
C
C Detrmine whether new mesh extends beyond oldmesh.
C NX is used to define the region over which we may use a linear
C interpolation. Beyond the old radius mesh we set the departure coefficients
C constant or equal to 1.
C
	NXST=1
	DO WHILE (EDLOG(NXST) .LT. OLDED(1))
	  NXST=NXST+1
	END DO
	NX_END=ND
	DO WHILE (EDLOG(NX_END) .GT. OLDED(NDOLD))
	  NX_END=NX_END-1
	END DO
	NX=NX_END-NXST+1
C
C Interpolate the populations assuming that the departure
C coefficients remain constant. Constant is one if
C read in b's-1.0
C
	DO J=1,NZ
	  DO I=1,NDOLD
	    TB(I)=LOG(HYD(J,I)+CONST)
	  END DO
	  CALL LINPOP(EDLOG(NXST),TA(NXST),NX,OLDED,TB,NDOLD)
	  DO I=NXST,NX_END
	    DHEN(J,I)=EXP(TA(I))
	  END DO
	  IF(NXST .NE. 1)THEN
	    DO I=1,NXST-1
	      DHEN(J,I)=EXP(TB(3)+(EDLOG(I)-OLDED(3))
	1        /(OLDED(1)-OLDED(3))*(TB(1)-TB(3)))
	    END DO
!
! Assume d.c. is constant in outer region. This option is more
! consistent with the assumption of constant T.
!
	    DO I=1,NXST-1
	      DHEN(J,I)=EXP(TB(1))
	    END DO
	  END IF
	  IF(NX_END .LT. ND)THEN
	    DO I=ND,NX_END+1,-1
	      DHEN(J,I)=DHEN(J,NX_END)
	    END DO
	  END IF
	END DO
C
C Compute departure coefficints for N>NZ . New levels are set to have the same
C departure coefficient as highest known level.
C
	IF(NLEV .GT. NZ)THEN
	  DO I=1,ND
	    DO J=NZ+1,NLEV
	      DHEN(J,I)=DHEN(NZ,I)
	    END DO
	  END DO
	END IF
C
C Regrid ion density.
C
	CALL LINPOP(EDLOG(NXST),ION(NXST),NX,OLDED,OLDION,NDOLD)
	IF(NXST .NE. 1)THEN
	  DO I=1,NXST-1
	    ION(I)=OLDION(3)+(EDLOG(I)-OLDED(3))
	1          *(OLDION(1)-OLDION(3)) / (OLDED(1)-OLDED(3))
	  END DO
	END IF
	DO I=1,NX_END
	  ION(I)=EXP(ION(I))
	END DO
C
C Assume ion density is fixed relative to ion density.
C
	IF(NX_END .LT. ND)THEN
	  DO I=ND,NX_END+1,-1
	    ION(I)=ED(I)*ION(NX_END)/ED(NX_END)
	  END DO
	END IF
C
C Free memory.
C
	DEALLOCATE (HYD)
	DEALLOCATE (OLDED)
	DEALLOCATE (TA)
	DEALLOCATE (TB)
	DEALLOCATE (OLDION)
C
	RETURN
	END
