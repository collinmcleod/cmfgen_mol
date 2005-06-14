!
! General routine to guess the departure coefficients for  a new species or
! ionization stage. The following files are required:
!
!     EDDFACTOR
!     RVTJ
!     XzV_F_OSCDAT
!
! These can come from a similar model but the frequency range must be sufficient.
! Program estimate th e ground state departure coefficient, and use the same
! excitation coefficient for all other levels.
!
	PROGRAM GUESS_DC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER NCF
	INTEGER ND
	REAL*8, POINTER :: RJ(:,:)
	REAL*8, POINTER :: NU(:)
	CHARACTER*10 DATA_TYPE
	CHARACTER*40 FILE_DATE
	CHARACTER*80 FILENAME
!
	INTEGER NUM_FILES
	INTEGER ID
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8, ALLOCATABLE :: PHOT_SUM(:)
	REAL*8, ALLOCATABLE :: RECOM_SUM(:)
	REAL*8, ALLOCATABLE :: GS_DC(:)
	REAL*8, ALLOCATABLE :: DC(:,:)
	REAL*8, ALLOCATABLE :: T_EXC(:)
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	INTEGER ND_ATM,NC_ATM,NP_ATM
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)
!
	INTEGER, PARAMETER :: N_MAX=2000
	CHARACTER*30 NAME(N_MAX)
	REAL*8 FEDGE(N_MAX)
	REAL*8 ENERGY(N_MAX)
	REAL*8 G(N_MAX)
	REAL*8 ION_EN
	REAL*8 ZION
	CHARACTER*30 EN_DATE
	INTEGER NLEV
!
	CHARACTER*6 METHOD,TYPE_ATM
	CHARACTER*10 NAME_CONVENTION
!
        INTEGER ACCESS_F
        INTEGER, PARAMETER :: EDD_CONT_REC=3
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
        INTEGER REC_SIZE
        INTEGER UNIT_SIZE
        INTEGER WORD_SIZE
        INTEGER N_PER_REC

! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,ISAV
	INTEGER ST_REC
	INTEGER REC_LENGTH
	REAL*8 SCALE_FAC
	REAL*8 TEMP
	REAL*8 T1,T2,T3
	REAL*8 LAMC
	REAL*8 T_ELEC
	LOGICAL AIR_LAM
	LOGICAL USE_V
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For file I/O
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: LU_HEAD=12
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER GET_INDX_DP
!
	CHARACTER*80 RVTJ_FILE_NAME
	CHARACTER*5 SPECIES
!
	REAL*8 SPEED_OF_LIGHT
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,UC
!
! 
! Set constants.
!
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
!
	METHOD='LOGMON'
	TYPE_ATM=' '
!
        CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
!  Read in EDDFACTOR file. This is used to compute the photoionization and
!         recombination rates.
!
	FILENAME='EDDFACTOR'
5	CALL GEN_IN(FILENAME,'First data file')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening/reading INFO file: check format'
	  WRITE(T_OUT,*)'Also check error file or fort.2'
	  GOTO 5
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error opening ',TRIM(FILENAME)
	     WRITE(T_OUT,*)'IOS=',IOS
	     GOTO 5
	  END IF
	  READ(LU_IN,REC=3)ST_REC,NCF,ND
	  ND=ND; NCF=NCF
	  ALLOCATE (RJ(ND,NCF))
	  ALLOCATE (NU(NCF))
	  DO ML=1,NCF
	    READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(RJ(I,ML),I=1,ND),NU(ML)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading all frequencies'
	      NCF=ML-1
	      EXIT
	    END IF
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in ',TRIM(FILENAME),' file as MODEL A (default)'
	WRITE(T_OUT,*)'Number of depth points is',ND
	WRITE(T_OUT,*)'Number of frequencies is ',NCF
	WRITE(T_OUT,*)' '
!
!
!
! *************************************************************************
!
! Read in basic model [i.e. R, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	RVTJ_FILE_NAME='RVTJ'
	CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
	OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    GOTO 10
	  END IF
	CLOSE(LU_IN)
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND_ATM,NC_ATM,NP_ATM,RVTJ_FILE_NAME,LU_IN)
	ALLOCATE (R(ND_ATM))
	ALLOCATE (V(ND_ATM))
	ALLOCATE (SIGMA(ND_ATM))
	ALLOCATE (T(ND_ATM))
	ALLOCATE (ED(ND_ATM))
	ALLOCATE (ROSS_MEAN(ND_ATM))
	ALLOCATE (FLUX_MEAN(ND_ATM))
	ALLOCATE (POP_ATOM(ND_ATM))
	ALLOCATE (MASS_DENSITY(ND_ATM))
	ALLOCATE (POPION(ND_ATM))
	ALLOCATE (CLUMP_FAC(ND_ATM))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND_ATM,LU_IN)
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND_ATM))
	 ALLOCATE (TB(ND_ATM))
	 ALLOCATE (TC(ND_ATM))
	 ALLOCATE (TAU_ROSS(ND_ATM))
	 ALLOCATE (TAU_ES(ND_ATM))
	 IF(ROSS_MEAN(ND_ATM) .NE. 0)THEN
	   CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
	 ELSE
	  TAU_ROSS(1:ND)=0.0D0
	 END IF
	 TA(1:ND_ATM)=6.65D-15*ED(1:ND_ATM)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
!
! 
!
! Loop to do many different species/ionization stages.
!
	SPECIES=' '
5000	CONTINUE
	WRITE(T_OUT,*)' '
	CALL GEN_IN(SPECIES,'Species (e.g., OIV, OVI) to guess d.c.''s for (or exit)')
	IF(UC(SPECIES(1:2)) .EQ. 'EX')STOP
	FILENAME=TRIM(SPECIES)//'_F_OSCDAT'
!
! Open Oscillator file.
!
	IOS=100
	DO WHILE(IOS .NE. 0)
	  CALL GEN_ASCI_OPEN(LU_HEAD,'HEAD_INFO','UNKNOWN',' ','WRITE',IZERO,IOS)
	  CALL GEN_IN(FILENAME,'Oscillator file')
	  CALL RD_ENERGY(NAME,G,ENERGY,FEDGE,NLEV,N_MAX,
	1       ION_EN,ZION,EN_DATE,FILENAME,LU_IN,LU_HEAD,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error occurred reading Oscillator file: try again'
	  END IF
	CLOSE(LU_HEAD)
	END DO
	WRITE(T_OUT,*)'Successfully read oscillator file'
	WRITE(T_OUT,*)' '
!
	ML=1
	DO WHILE(FEDGE(1) .LT. NU(ML))
	  ML=ML+1
	END DO
	J=ML-1
!
! Compute the recombination and photoionization rate to the ground state.
! For simplicity we assume the ground state population is set by a balance
! between photoionizations and recombinations. We ignore the frequency
! dependence of the photoionization cross-section, and not that its numerical
! value does not matter.
!
	IF(.NOT. ALLOCATED(PHOT_SUM))THEN
	  ALLOCATE (PHOT_SUM(1:ND))
	  ALLOCATE (RECOM_SUM(1:ND))
	END IF
	PHOT_SUM(1:ND)=0.0D0
	RECOM_SUM(1:ND)=0.0D0
	DO ML=1,J-1
	  DO I=1,ND
	    PHOT_SUM(I)=PHOT_SUM(I)+0.5D0*(NU(ML)-NU(ML+1))*(RJ(I,ML)+RJ(I,ML+1))
	    T1=(TWOHCSQ*(NU(ML)**3)+RJ(I,ML))*EXP(-HDKT*NU(ML)/T(I))
	    T2=(TWOHCSQ*(NU(ML+1)**3)+RJ(I,ML+1))*EXP(-HDKT*NU(ML+1)/T(I))
	    RECOM_SUM(I)=RECOM_SUM(I)+0.5D0*(NU(ML)-NU(ML+1))*(T1+T2)
	  END DO
	END DO
!
! Compute the ground-state departure coefficient, and the excitation temperature
! of the ground state.
!
	IF(.NOT. ALLOCATED(GS_DC))THEN
	  ALLOCATE (GS_DC(1:ND))
	  ALLOCATE (T_EXC(1:ND))
	END IF
	GS_DC(1:ND)=RECOM_SUM(1:ND)/PHOT_SUM(1:ND)
	T1=5.0D0
	DO I=1,ND
	  DO J=1,10
	    T1=LOG( GS_DC(I)*(T1/T(I))**1.5 )/HDKT/FEDGE(1)+1.0/T(I)
	    T1=1.0/T1
	  END DO
	  T_EXC(I)=T1
	  IF(I .EQ. 1 .OR. I .EQ. ND)THEN
	    WRITE(T_OUT,'(I4,2X,A,F8.3,3X,A,F8.3)')I,'T=',T(I),'T_EXC=',T_EXC(I)
	  END IF
	END DO	
	WRITE(T_OUT,*)GS_DC(1),FEDGE(1)
!
! Compute departure coefficients of all levels. We assume that they have the same
! departure coefficient as the ground state.
!
	IF(ALLOCATED(DC))DEALLOCATE(DC)
	ALLOCATE (DC(NLEV,ND))
	DO I=1,ND
	   DO J=1,NLEV
	     DC(J,I)=((T(I)/T_EXC(I))**1.5 )*EXP(HDKT*FEDGE(J)*(1.0D0/T_EXC(I)-1.0D0/T(I)))
	   END DO
	END DO
!
	FILENAME=TRIM(SPECIES)//'_IN'
	CALL GEN_IN(FILENAME,'Output file for DCs --- old file will be overwritten')
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'UNKNOWN',' ','WRITE',IZERO,IOS)
	  WRITE(LU_OUT,'(/,X,A,T40,A)')'07-Jul-1997','!Format date'
	  WRITE(LU_OUT,2120)R(ND),RLUM,NLEV,ND
	  DO I=1,ND
	    WRITE(LU_OUT,2122)R(I),1.0D-100,ED(I),T(I),0.0,V(I),CLUMP_FAC(I)
	    WRITE(LU_OUT,'(1X,1P,5E17.7)')(DC(J,I),J=1,NLEV)
	  END DO
	CLOSE(LU_OUT)
2120	FORMAT(/,X,F11.6,5X,1PE12.6,5X,0P,I4,5X,I4)
2122	FORMAT(/,X,1P,E16.8,6E17.8)
!
! Loop back to do addition species and ionization stages.
!
	SPECIES='EXIT'
	GOTO 5000
	END
