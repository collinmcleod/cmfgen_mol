!
! Read in collsion strengths from a file. The values are tabulated as a 
! function of temperature.
!
! OMEGA_SET is used to set the collison strength for those transitions in which
!    (1) No collison strength was contained in the file and
!    (2) has a zero oscillator strength.
!
! OMEGA_SCALE is used by SUBCOL_MULTI to scale all collsion strengths by a
! constant value.
!
	SUBROUTINE GEN_OMEGA_RD_V2(OMEGA_TABLE,T_TABLE,
	1                      ID_LOW,ID_UP,ID_INDX,
	1                      OMEGA_SCALE,OMEGA_SET,
	1                      LEVNAME,STAT_WT,NLEV,FILE_NAME,
	1                      NUM_TRANS,NUM_TVALS,
	1                      MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE)
	IMPLICIT NONE
!
! Altered 02-Mar-1998 : Ientifying the connection betwen table and model 
!                         levels improved. Routine can now handle split and
!                         combined LS stars in both the model and input data.
! Altered 08-Feb-1998 : STRING length invreased from 132 to 500 for Fe2.
! Altered 19-Dec-1996 : Bugs fixed, and changed splitting approximation.
! Altered 07-Dec-1996 : Can read in multiplet Omega values when split
!                         levls are present. Spliting is approximate.
! Altered 29-May-1996 : Bug fix. Incorrect handling of collisional data for
!                         levels not included in the model atom.
! Altered 24-May-1996 : GEN_SEQ_OPEN replaced by F90 OPEN statement.
! Altered 01-Mar-1996 : Warning about invalid levels commented out.
!                       Warning was occuring just because of additional levels
!                       in data file.
! Created 13-Sep-1995 
!
	INTEGER*4 NLEV
	CHARACTER*(*) LEVNAME(NLEV),FILE_NAME
	REAL*8 STAT_WT(NLEV)
	REAL*8 OMEGA_SCALE,OMEGA_SET
!
! ID_LOW:  ID of lower level of transition
! ID_UP:   ID of upper level of transition
! ID_INDX: Start location of tabulated values in table (OMEGA_TABLE)
!
! T_TABLE: Tempearure values (10^4 K) at which collsion strengths tabulated.
!
	INTEGER*4 MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE
	INTEGER*4 ID_LOW(MAX_TRANS)
	INTEGER*4 ID_UP(MAX_TRANS)
	INTEGER*4 ID_INDX(MAX_TRANS)    
	REAL*8 T_TABLE(MAX_TVALS)
	REAL*8 OMEGA_TABLE(MAX_TAB_SIZE)
!
	INTEGER*4 NUM_TRANS,NUM_TVALS
!
! Local variables.
!
	REAL*8 COL_VEC(MAX_TVALS)
	REAL*8 GL_SUM,GU_SUM,NORM_FAC
	INTEGER*4 NL_VEC(10)
	INTEGER*4 NUP_VEC(10)
	INTEGER*4 NL_CNT
	INTEGER*4 NUP_CNT
	INTEGER*4 LOC_INDX(NLEV,NLEV)
!
	INTEGER*4 I		!Used for lower level.
	INTEGER*4 K
	INTEGER*4 LOOP
	INTEGER*4 TRANS_CNT	!Used for reading in the dtata.
	INTEGER*4 NL,NUP
	INTEGER*4 MAX_LEV_SEP
	INTEGER*4 L,LUIN,IOS
	INTEGER*4 TRANS_INDX
	INTEGER*4 TAB_INDX
	INTEGER*4 LST_SPLIT_LEVEL
	INTEGER*4 LST_NL
	INTEGER*4 LST_NUP
	INTEGER*4 LUER,ERROR_LU,ICHRLEN
	CHARACTER*80 STRIP_LD_BLANK
	EXTERNAL ERROR_LU,ICHRLEN
!
	CHARACTER*132 LOW_LEV,UP_LEV
	CHARACTER*500 STRING
	CHARACTER*60 LOCNAME(NLEV)
!
	LUER=ERROR_LU()
!
	DO NL=1,NLEV
	  LOCNAME(NL)=STRIP_LD_BLANK(LEVNAME(NL))
	END DO
	LST_SPLIT_LEVEL=0
	DO NL=1,NLEV
	  IF(INDEX(LOCNAME(NL),'[') .NE. 0)LST_SPLIT_LEVEL=NL
	END DO
	MAX_LEV_SEP=1
	DO NL=1,LST_SPLIT_LEVEL
	  K=INDEX(LOCNAME(NL),'[')
	  IF(K .NE. 0)THEN
	    DO NUP=NL+1,NLEV
	      IF(LOCNAME(NL)(1:K) .EQ. LOCNAME(NUP)(1:K))THEN
	        MAX_LEV_SEP=MAX(MAX_LEV_SEP,NUP-NL)
	      END IF
	    END DO
	  END IF
	END DO
	LOC_INDX(:,:)=0
	OMEGA_TABLE(:)=0.0D0
!
! Open file with collisional data.
!
	LUIN=7
	OPEN(LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',FILE_NAME
	  WRITE(LUER,*)'IOS=',IOS	
	  STOP
	END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of transitions')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Numer of transitions string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
!
	READ(STRING,*)NUM_TRANS
	IF(NUM_TRANS .GT. MAX_TRANS)THEN
	  WRITE(LUER,*)'Error reading collisonal data from '//FILE_NAME
	  WRITE(LUER,*)'MAX_TRANS too small'
	  STOP
	END IF
	STRING=' '
!
	DO WHILE( INDEX(STRING,'!Number of T values OMEGA tabulated at')  
	1               .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Numer of T values string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)NUM_TVALS
	IF(NUM_TVALS .GT. MAX_TVALS)THEN
	  WRITE(LUER,*)'Error reading collisonal data from '//FILE_NAME
	  WRITE(LUER,*)'NU_TVALS too small'
	  STOP
	END IF
!         
! OMEGA_SCALE provides a means of scaling the OMEGA that have not been read in.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Scaling factor for OMEGA') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Scaling factor for Omega string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)OMEGA_SCALE
!         
! OMEGA_SET provides a means to set OMEGA for transition for which
! no atomic data is availaable, and for which the oscilator strength
! is zero.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Value for OMEGA if f=0') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Value for OMEGA if f=0 string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)OMEGA_SET
!
! We do this check here so that OMEGA_SCALE and OMEGA_SET are defined.
!
	IF(NUM_TRANS .EQ. 0)THEN
	  CLOSE(LUIN)
	  RETURN			!i.e use approximate formulae only.
	END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'Transition\T')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'Transition\T string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	L=INDEX(STRING,'ion\T')
	READ(STRING(L+5:),*)(T_TABLE(I),I=1,NUM_TVALS )
!
	STRING=' '
	DO WHILE(ICHRLEN(STRING) .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO             
!
! TRANS_CNT is the index used to read in all collisional data. 
! TRANS_INDX is used to specify the number of collisional transitions that
! are used by the present model. It will differ from TRANS_CNT because
! collisional data may be available for levels not included in the model
! atom.
!
! NB: TRANS_INDX =< TRANS_CNT.
!
	TAB_INDX=1
	TRANS_INDX=0
	LST_NL=1
	LST_NUP=1
	DO TRANS_CNT=1,NUM_TRANS
!
! Read and determine correspondance of lower level to that of model atom.
!
! Get lower level name.
!
	  I=1
	  DO WHILE (STRING(I:I) .EQ. ' ')
	    I=I+1
	  END DO
	  L=INDEX(STRING,'-')          
	  LOW_LEV=STRING(I:L-1)
!
! Determine corespondance of lower level. NB: The MOD function allows us
! to loop over all levels beggining from LST_NL. Thus the SEARCH order in
! INDEX space is
!              LST_NL, LST_NL+1, .... NLEV, 1, 2, ... , LST_NL-1
!
	  K=INDEX(LOW_LEV,'[')
	  NL_CNT=0
	  GL_SUM=0.0D0
	  IF(K .NE. 0)THEN
!
! Only a single unique match is possible.
!
	    DO LOOP=1,NLEV
	      NL=MOD(LOOP+LST_NL-2,NLEV)+1
	      IF( LOW_LEV .EQ. LOCNAME(NL) .OR.
	1            LOW_LEV(1:K-1).EQ. LOCNAME(NL) )THEN
	        NL_CNT=NL_CNT+1
	        NL_VEC(NL_CNT)=NL
	        GL_SUM=GL_SUM+STAT_WT(NL)
	        GOTO 100
	      END IF
	    END DO
	    GOTO 500			!No match
	  ELSE
!
! As the input level has no [], multiple matches are possible if model
! atom levels are split.
!
	    DO LOOP=1,NLEV
	      NL=MOD(LOOP+LST_NL-2,NLEV)+1
	      K=INDEX(LOCNAME(NL),'[')
	      IF(K .EQ. 0)THEN
	        IF(LOW_LEV .EQ. LOCNAME(NL))THEN
	          NL_CNT=NL_CNT+1
                  NL_VEC(NL_CNT)=NL
	          GL_SUM=GL_SUM+STAT_WT(NL)
	          GOTO 100
	        END IF
	      ELSE
	        IF(LOW_LEV .EQ. LOCNAME(NL)(1:K-1))THEN
	          DO I=MAX(1,NL-MAX_LEV_SEP),MIN(NLEV,NL+MAX_LEV_SEP)
	            IF(LOW_LEV .EQ. LOCNAME(I)(1:K-1))THEN
	              NL_CNT=NL_CNT+1
                      NL_VEC(NL_CNT)=I
	              GL_SUM=GL_SUM+STAT_WT(I)
	            END IF
	          END DO
	          GOTO 100
	        END IF
	      END IF
	    END DO
	    GOTO 500		!No match
	  END IF
100	  CONTINUE
!
! Read and determine correspondance of upper level to that of model atom.
!
! Get upper level name.
!
	  STRING(1:)=STRING(L+1:)
	  I=1
	  DO WHILE (STRING(I:I) .EQ. ' ')
	    I=I+1                 
	  END DO
	  STRING(1:)=STRING(I:)
	  L=INDEX(STRING,'  ')
	  UP_LEV=STRING(1:L-1)
!
! Determine upper level corespondance.
!
	  K=INDEX(UP_LEV,'[')
	  NUP_CNT=0
	  GU_SUM=0.0D0
	  IF(NL_VEC(1) .NE. LST_NL)LST_NUP=NL_VEC(1)
	  IF(UP_LEV .EQ. 'I')THEN
!
! Collisional ionization terms.
!
	    GU_SUM=1.0D0
	    GOTO 200
	  ELSE IF(K .NE. 0)THEN
	    DO LOOP=1,NLEV
	      NUP=MOD(LOOP+LST_NUP-2,NLEV)+1
	      IF(UP_LEV .EQ. LOCNAME(NUP) .OR.
	1            UP_LEV(1:K-1) .EQ. LOCNAME(NUP))THEN
	        NUP_CNT=NUP_CNT+1
	        NUP_VEC(NUP_CNT)=NUP
	        GU_SUM=GU_SUM+STAT_WT(NUP)
	        GOTO 200
	      END IF
	    END DO
	    GOTO 500				!No match
	  ELSE
	    DO LOOP=1,NLEV
	      NUP=MOD(LOOP+LST_NUP-2,NLEV)+1
	      K=INDEX(LOCNAME(NUP),'[')
	      IF(K .EQ. 0)THEN
	        IF(UP_LEV .EQ. LOCNAME(NUP))THEN
	          NUP_CNT=NUP_CNT+1
                  NUP_VEC(NUP_CNT)=NUP
	          GU_SUM=GU_SUM+STAT_WT(NUP)
	          GOTO 200
	        END IF
	      ELSE
!
! Several matches possible
!
	        IF(UP_LEV .EQ. LOCNAME(NUP)(1:K-1))THEN
	          DO I=MAX(1,NUP-MAX_LEV_SEP),MIN(NLEV,NUP+MAX_LEV_SEP)
	            IF(UP_LEV .EQ. LOCNAME(I)(1:K-1))THEN
	              NUP_CNT=NUP_CNT+1
                      NUP_VEC(NUP_CNT)=I
	              GU_SUM=GU_SUM+STAT_WT(I)
	            END IF
	          END DO
	          GOTO 200
	        END IF
	      END IF
	    END DO
	    GOTO 500		!No match
	  END IF
!
200	  NORM_FAC=1.0D0/GL_SUM/GU_SUM
	  LST_NL=NL_VEC(1)
	  LST_NUP=NUP_VEC(1)
	  STRING(1:)=STRING(L:)
!
! The multiplet collision strengths are converted to individual collision 
! strengths using:
!
! Omega(i,J)= g(i)g(j) OMEGA(MULITPLET) / G(L)/G(U)
!
! This formula automatically recovers the correct result that
!
! OMEGA = SUM(i) SUM(J)  Omega(i,j)
!
	  READ(STRING,*)(COL_VEC(I),I=1,NUM_TVALS)
	  COL_VEC(1:NUM_TVALS)=COL_VEC(1:NUM_TVALS)*NORM_FAC
	  IF(UP_LEV .EQ. 'I')THEN
	    DO NL=1,NL_CNT
	      TRANS_INDX=TRANS_INDX+1
	      ID_LOW(TRANS_INDX)=NL_VEC(NL)
	      ID_UP(TRANS_INDX)=NL_VEC(NL)
	      ID_INDX(TRANS_INDX)=TAB_INDX
	      OMEGA_TABLE(TAB_INDX:TAB_INDX+NUM_TVALS-1)=
	1         OMEGA_TABLE(TAB_INDX:TAB_INDX+NUM_TVALS-1)+
	1            COL_VEC(1:NUM_TVALS)*STAT_WT(NL_VEC(NL))
	      TAB_INDX=TAB_INDX+NUM_TVALS
	    END DO
	  ELSE
	    DO NL=1,NL_CNT
	      DO NUP=1,NUP_CNT
	        IF(LOC_INDX(NL,NUP) .EQ. 0)THEN
	          TRANS_INDX=TRANS_INDX+1
	          ID_LOW(TRANS_INDX)=NL_VEC(NL)
	          ID_UP(TRANS_INDX)=NUP_VEC(NUP)
	          ID_INDX(TRANS_INDX)=TAB_INDX
	          I=TAB_INDX
	          TAB_INDX=TAB_INDX+NUM_TVALS
	        ELSE
	          I=LOC_INDX(NL,NUP)
	        END IF
	        OMEGA_TABLE(I:I+NUM_TVALS-1)=OMEGA_TABLE(I:I+NUM_TVALS-1)+
	1            COL_VEC(1:NUM_TVALS)*STAT_WT(NUP_VEC(NUP))*
	1                              STAT_WT(NL_VEC(NL))
	      END DO
	    END DO
	  END IF
!
500	  CONTINUE
!
! Read in the next collisional data set. If the record is blank, we
! read in the NEXT record, as BLANK records are not included in the NUM_TRANS
! count.
!
	  STRING=' '
	  DO WHILE(ICHRLEN(STRING) .EQ. 0 .AND. TRANS_CNT .LT. NUM_TRANS)
	    READ(LUIN,'(A)')STRING
	  END DO
!
	END DO
!
! NUM_TRANS is now set equal to the actual number of transitions read.
!
	NUM_TRANS=TRANS_INDX
!
	CLOSE(LUIN)
	RETURN
	END

	FUNCTION STRIP_LD_BLANK(STRING)
	IMPLICIT NONE
C
	INTEGER*4 I
	CHARACTER*(*) STRING
	CHARACTER*80 STRIP_LD_BLANK
	I=1
	DO WHILE(I .LT. LEN(STRING) .AND. STRING(I:I) .EQ. ' ')
	 I=I+1
	END DO
	STRIP_LD_BLANK=STRING(I:)
	RETURN
	END
