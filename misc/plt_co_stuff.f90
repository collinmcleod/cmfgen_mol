      PROGRAM PLT_CO_STUFF
      USE SET_KIND_MODULE
      USE GEN_IN_INTERFACE
      USE MOD_COLOR_PEN_DEF
      IMPLICIT NONE
!
! Simple plotting routine designed to visualize the CO cooling as a function of temperature
! Created Apr-20-2022
!
      INTEGER I,J,K,L,ND,NMOD,DEPTH,IOS
      INTEGER NRXNS,NRXNS_DISTINCT,CHOSEN_RXN
      INTEGER LU_CO,MIN_ITERS,MAX_ITERS,DEPTH_POINT
      LOGICAL IS_OPEN
      REAL(KIND=LDP), ALLOCATABLE :: T(:,:),CO_COOLING(:,:),CO_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: COMI_POPS(:,:),COM2_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CI_POPS(:,:),OI_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2_POPS(:,:),O2_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CO2MI_POPS(:,:),CO2M2_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2MI_POPS(:,:),C2M2_POPS(:,:) 
      REAL(KIND=LDP), ALLOCATABLE :: O2MI_POPS(:,:),O2M2_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: ED_POPS(:,:),COMP_POPS(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: COMI_D(:,:),COMI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: COM2_D(:,:),COM2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CI_D(:,:),CI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2_D(:,:),C2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: OI_D(:,:),OI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: O2_D(:,:),O2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2MI_D(:,:),C2MI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2M2_D(:,:),C2M2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: O2MI_D(:,:),O2MI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: O2M2_D(:,:),O2M2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CO2MI_D(:,:),CO2MI_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CO2M2_D(:,:),CO2M2_D_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: COMI_RR(:,:),COMI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CI_RR(:,:),CI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2_RR(:,:),C2_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: OI_RR(:,:),OI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: O2_RR(:,:),O2_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: C2MI_RR(:,:),C2MI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: O2MI_RR(:,:),O2MI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CO2MI_RR(:,:),CO2MI_PR(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RXN_RATES(:,:,:),TOT_COOLING(:,:),NET_COOLING(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RXN_RATES1(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: X(:),Y(:),Z(:)
      REAL(KIND=LDP), ALLOCATABLE :: X1(:),Y1(:)
      REAL(KIND=LDP), ALLOCATABLE :: CO_SCR_DATA(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RAD_HEAT(:,:),RHO(:)
      INTEGER, ALLOCATABLE :: INT_ARRAY(:)
      REAL(KIND=LDP) :: T1,T2
      CHARACTER(50) :: FILENAME,MODEL_LIST_FILE
      CHARACTER(50) :: FILE1,FILE2
      CHARACTER(100), ALLOCATABLE :: DIR_LIST(:),RXNS_LIST(:)
      CHARACTER(50) :: QUERY
      CHARACTER(200) :: STRING,TSTR1,TSTR2,TSTR3
      CHARACTER(30) :: PLT_OPT,XLABEL,YLABEL
      CHARACTER(30) :: CO_OPT,SCR_OPT,RCT1,RCT2
      CHARACTER(30) UC
      LOGICAL :: BOOL1,BOOL2,BOOL3
      INTEGER :: VI,V2,V3,J1,J2,J3
      INTEGER, ALLOCATABLE :: X_ARGSORT(:)
      
      TYPE ENERGY_LEVEL
      CHARACTER(33) :: LEVEL_NAME
      REAL(KIND=LDP) :: LEVEL_G
      REAL(KIND=LDP) :: ENERGY_CM1
      REAL(KIND=LDP) :: FREQ_1015HZ
      REAL(KIND=LDP) :: ENERGY_BELOW_ION_EV
      REAL(KIND=LDP) :: LAMBDA_ANG
      INTEGER :: LEVEL_ID
      REAL(KIND=LDP) :: ARAD
      END TYPE

      TYPE TRANSITION
      CHARACTER(28) :: LOWER_LEVEL_NAME
      CHARACTER(32) :: UPPER_LEVEL_NAME
      REAL(KIND=LDP) :: OSCILLATOR_STR
      REAL(KIND=LDP) :: EINSTEIN_A
      REAL(KIND=LDP) :: LAMBDA_ANG
      INTEGER :: LOWER_LEVEL_ID,UPPER_LEVEL_ID,TRANS_ID
      END TYPE
      
      TYPE(ENERGY_LEVEL), DIMENSION(:), ALLOCATABLE :: COMI_LEVEL_LIST
      TYPE(TRANSITION), DIMENSION(:), ALLOCATABLE :: COMI_TRANS_LIST
      TYPE(ENERGY_LEVEL) :: TEMPLEV1, TEMPLEV2
      TYPE(TRANSITION) :: TEMPTRANS1 ,TEMPTRANS2
      INTEGER :: NLEVELS, NTRANS
      REAL(KIND=LDP) :: COMI_ION_EV

      REAL(KIND=LDP), ALLOCATABLE :: THIN_COOL(:),LTE_POPS(:)
      REAL(KIND=LDP) :: PART_F,LOC_POP
      REAL(KIND=LDP), PARAMETER :: CM1TOEV=0.00012398470543514755D0
      REAL(KIND=LDP), PARAMETER :: PLANCKEV=4.135667696D-15
      REAL(KIND=LDP), PARAMETER :: LIGHTCM=29979245800.0D0
      REAL(KIND=LDP), PARAMETER :: KBEV=8.617333262D-5
      REAL(KIND=LDP), PARAMETER :: FACTOR1=1.602176634D-4 !useful conversion factor when calculating line cooling
                            ! Converts from eV*cm/A -> erg
      REAL(KIND=LDP) :: PI_NUM,HCON4PI,EMIT_CONST,HC
      
      EXTERNAL UC
!
      PI_NUM = 4.0_LDP*ATAN(1.0_LDP)
      HCON4PI = (PLANCKEV*LIGHTCM)/(4.0_LDP*PI_NUM)
      HC = PLANCKEV*LIGHTCM
      EMIT_CONST=HC*FACTOR1 !To reduce total number of calculations
!
      WRITE(6,'(A)') 'A simple program for plotting data related to CO cooling from'
! '
      WRITE(6,*)
      WRITE(6,'(A)') 'Reading CO level and transition data: make sure the data file is linked in&
      & the current directory as COMI_F_OSCDAT'
!'
      OPEN(UNIT=12,FILE='COMI_F_OSCDAT',ACTION='READ')

      DO
         READ(12,'(A)') TSTR1
         IF (INDEX(TSTR1,'!Date') .NE. 0) EXIT
      END DO
      
      READ(12,*) NLEVELS,TSTR1
      READ(12,*) COMI_ION_EV,TSTR1
      READ(12,'(A)') TSTR1
      READ(12,*) NTRANS,TSTR1

      ALLOCATE(COMI_LEVEL_LIST(NLEVELS))
      ALLOCATE(COMI_TRANS_LIST(NTRANS))
      
      READ(12,'(A)') TSTR1
      DO I=1,NLEVELS
         READ(12,'(A33,F5.1,F15.4,F12.5,F9.3,ES13.3,I8,ES12.3,ES12.3,&
         &ES12.3)') COMI_LEVEL_LIST(I)%LEVEL_NAME,&
         &COMI_LEVEL_LIST(I)%LEVEL_G,COMI_LEVEL_LIST(I)%ENERGY_CM1,&
         &COMI_LEVEL_LIST(I)%FREQ_1015HZ,COMI_LEVEL_LIST(I)%ENERGY_BELOW_ION_EV,&
         &COMI_LEVEL_LIST(I)%LAMBDA_ANG,COMI_LEVEL_LIST(I)%LEVEL_ID,&
         &COMI_LEVEL_LIST(I)%ARAD,T1,T2
!'
      END DO

      DO I=1,10
         READ(12,'(A)') TSTR1
      END DO

      WRITE(6,*) 'TSTR1:',TSTR1
      CALL FLUSH(6)

      DO I=1,NTRANS
         READ(12,'(A28,A1,A32,ES10.4,ES12.4,ES12.3,I8,A1,I4,I9)') COMI_TRANS_LIST(I)%LOWER_LEVEL_NAME,&
              TSTR1,COMI_TRANS_LIST(I)%UPPER_LEVEL_NAME,COMI_TRANS_LIST(I)%OSCILLATOR_STR,&
              COMI_TRANS_LIST(I)%EINSTEIN_A,COMI_TRANS_LIST(I)%LAMBDA_ANG,COMI_TRANS_LIST(I)%LOWER_LEVEL_ID,&
              TSTR2,COMI_TRANS_LIST(I)%UPPER_LEVEL_ID,COMI_TRANS_LIST(I)%TRANS_ID
      END DO

      CLOSE(12)

      QUERY='M'
      WRITE(6,*) ' '
      WRITE(6,*) 'S: Plot cooling data for a single model'
      WRITE(6,*) 'M: Plot cooling data for multiple models'
      CALL GEN_IN(QUERY,'SCR or MULT?')
      IF (UC(QUERY) .EQ. 'M') THEN
         QUERY='N'
         CALL GEN_IN(QUERY,'Set model list?')
         IF (INDEX(QUERY,'Y') .NE. 0 .OR. INDEX(QUERY,'y') .NE. 0) THEN
!
            NMOD=1
            CALL GEN_IN(NMOD,'Number of models')
            ND=65
            CALL GEN_IN(ND,'Number of depth points in each model')
            ALLOCATE(DIR_LIST(NMOD))
!
            TSTR1='DIR'
            DO I=1,NMOD
               WRITE(TSTR2,'(A,I3.3)') 'Location of model ',I
               CALL GEN_IN(TSTR1,TSTR2)
!
! Format the directory strings as './directory/'
!
               IF (INDEX(TSTR1,'./') .NE. 1) THEN
                  TSTR1 = './'//TSTR1
               END IF
               J=LEN_TRIM(TSTR1)
               IF (INDEX(TSTR1(J:J),'/') .eq. 0) THEN
                  TSTR1=TRIM(ADJUSTL(TSTR1))//'/'
               END IF
               DIR_LIST(I)=TSTR1
            END DO
!
            MODEL_LIST_FILE='MODELS_FOR_PLT_CO'
            OPEN(1,FILE=TRIM(ADJUSTL(MODEL_LIST_FILE)),STATUS='REPLACE',IOSTAT=IOS,ACTION='WRITE')
            IF (IOS .NE. 0) THEN
               WRITE(6,*) 'Error opening file ',MODEL_LIST_FILE
               STOP
            END IF
            WRITE(1,'(A)') 'N_MODELS:'
            WRITE(1,*) NMOD
            WRITE(1,'(A)') 'ND:'
            WRITE(1,*) ND
            WRITE(1,'(A)') 'DIR_LIST:'
            DO I=1,NMOD
               WRITE(1,'(A)') TRIM(ADJUSTL(DIR_LIST(I)))
            END DO
            CLOSE(1)
!
         ELSE
            MODEL_LIST_FILE='MODELS_FOR_PLT_CO'
            OPEN(1,FILE=TRIM(ADJUSTL(MODEL_LIST_FILE)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
            IF (IOS .NE. 0) THEN
               WRITE(6,*) 'Error opening file ',MODEL_LIST_FILE
               STOP
            END IF
!     
            READ(1,'(A)') TSTR1
            READ(1,*) NMOD
            ALLOCATE(DIR_LIST(NMOD))
            READ(1,'(A)') TSTR1
            READ(1,*) ND
            READ(1,'(A)') TSTR1
            DO I=1,NMOD
               READ(1,'(A)') DIR_LIST(I)
            END DO
            CLOSE(1)
         END IF
!
         ALLOCATE(T(NMOD,ND))
         T(:,:) = 0.0D0
         ALLOCATE(CO_COOLING(NMOD,ND))
         CO_COOLING(:,:)=0.0D0
         ALLOCATE(CO_POPS(NMOD,ND))
         CO_POPS(:,:)=0.0D0
         ALLOCATE(COMI_POPS(NMOD,ND))
         COMI_POPS(:,:)=0.0D0
         ALLOCATE(COM2_POPS(NMOD,ND))
         COM2_POPS(:,:)=0.0D0
         ALLOCATE(CI_POPS(NMOD,ND))
         CI_POPS(:,:)=0.0D0
         ALLOCATE(C2_POPS(NMOD,ND))
         C2_POPS(:,:)=0.0D0
         ALLOCATE(OI_POPS(NMOD,ND))
         OI_POPS(:,:)=0.0D0
         ALLOCATE(O2_POPS(NMOD,ND))
         O2_POPS(:,:)=0.0D0
         ALLOCATE(C2MI_POPS(NMOD,ND))
         C2MI_POPS(:,:)=0.0D0
         ALLOCATE(C2M2_POPS(NMOD,ND))
         C2M2_POPS(:,:)=0.0D0
         ALLOCATE(O2MI_POPS(NMOD,ND))
         O2MI_POPS(:,:)=0.0D0
         ALLOCATE(O2M2_POPS(NMOD,ND))
         O2M2_POPS(:,:)=0.0D0
         ALLOCATE(CO2MI_POPS(NMOD,ND))
         CO2MI_POPS(:,:)=0.0D0
         ALLOCATE(CO2M2_POPS(NMOD,ND))
         CO2M2_POPS(:,:)=0.0D0
         ALLOCATE(RXN_RATES(NMOD,3,ND))
         RXN_RATES(:,:,:)=0.0D0
         ALLOCATE(RXN_RATES1(NMOD,ND))
         RXN_RATES1(:,:)=0.0D0
         ALLOCATE(TOT_COOLING(NMOD,ND))
         TOT_COOLING(:,:)=0.0D0
         ALLOCATE(NET_COOLING(NMOD,ND))
         NET_COOLING(:,:)=0.0D0
         ALLOCATE(THIN_COOL(NMOD))
         ALLOCATE(LTE_POPS(NLEVELS))
         LTE_POPS(:)=0.0D0
         THIN_COOL(:)=0.0D0
         ALLOCATE(X(2*ND))
         ALLOCATE(Y(2*ND))
         ALLOCATE(X1(2*ND))
         ALLOCATE(Y1(2*ND))
         ALLOCATE(Z(2*ND))
         ALLOCATE(X_ARGSORT(2*ND))
         ALLOCATE(ED_POPS(NMOD,ND))
         ALLOCATE(RAD_HEAT(NMOD,ND))
         ED_POPS(:,:)=0.0D0
         RAD_HEAT(:,:)=0.0D0
!
         ALLOCATE(COMI_D(NMOD,ND))
         ALLOCATE(COMI_D_PR(NMOD,ND))
         COMI_D(:,:)=0.0D0
         COMI_D_PR(:,:)=0.0D0
         ALLOCATE(COM2_D(NMOD,ND))
         ALLOCATE(COM2_D_PR(NMOD,ND))
         COM2_D(:,:)=0.0D0
         COM2_D_PR(:,:)=0.0D0
         ALLOCATE(CI_D(NMOD,ND))
         ALLOCATE(CI_D_PR(NMOD,ND))
         CI_D(:,:)=0.0D0
         CI_D_PR(:,:)=0.0D0
         ALLOCATE(C2_D(NMOD,ND))
         ALLOCATE(C2_D_PR(NMOD,ND))
         C2_D(:,:)=0.0D0
         C2_D_PR(:,:)=0.0D0
         ALLOCATE(OI_D(NMOD,ND))
         ALLOCATE(OI_D_PR(NMOD,ND))
         OI_D(:,:)=0.0D0
         OI_D_PR(:,:)=0.0D0
         ALLOCATE(O2_D(NMOD,ND))
         ALLOCATE(O2_D_PR(NMOD,ND))
         O2_D(:,:)=0.0D0
         O2_D_PR(:,:)=0.0D0
         ALLOCATE(C2MI_D(NMOD,ND))
         ALLOCATE(C2MI_D_PR(NMOD,ND))
         C2MI_D(:,:)=0.0D0
         C2MI_D_PR(:,:)=0.0D0
         ALLOCATE(C2M2_D(NMOD,ND))
         ALLOCATE(C2M2_D_PR(NMOD,ND))
         C2M2_D(:,:)=0.0D0
         C2M2_D_PR(:,:)=0.0D0
         ALLOCATE(O2MI_D(NMOD,ND))
         ALLOCATE(O2MI_D_PR(NMOD,ND))
         O2MI_D(:,:)=0.0D0
         O2MI_D_PR(:,:)=0.0D0
         ALLOCATE(O2M2_D(NMOD,ND))
         ALLOCATE(O2M2_D_PR(NMOD,ND))
         O2M2_D(:,:)=0.0D0
         O2M2_D_PR(:,:)=0.0D0
         ALLOCATE(CO2MI_D(NMOD,ND))
         ALLOCATE(CO2MI_D_PR(NMOD,ND))
         CO2MI_D(:,:)=0.0D0
         CO2MI_D_PR(:,:)=0.0D0
         ALLOCATE(CO2M2_D(NMOD,ND))
         ALLOCATE(CO2M2_D_PR(NMOD,ND))
         CO2M2_D(:,:)=0.0D0
         CO2M2_D_PR(:,:)=0.0D0
!
         ALLOCATE(COMI_RR(NMOD,ND))
         ALLOCATE(COMI_PR(NMOD,ND))
         COMI_RR(:,:)=0.0D0
         COMI_PR(:,:)=0.0D0
         ALLOCATE(CI_RR(NMOD,ND))
         ALLOCATE(CI_PR(NMOD,ND))
         CI_RR(:,:)=0.0D0
         CI_PR(:,:)=0.0D0
         ALLOCATE(C2_RR(NMOD,ND))
         ALLOCATE(C2_PR(NMOD,ND))
         C2_RR(:,:)=0.0D0
         C2_PR(:,:)=0.0D0
         ALLOCATE(OI_RR(NMOD,ND))
         ALLOCATE(OI_PR(NMOD,ND))
         OI_RR(:,:)=0.0D0
         OI_PR(:,:)=0.0D0
         ALLOCATE(O2_RR(NMOD,ND))
         ALLOCATE(O2_PR(NMOD,ND))
         O2_RR(:,:)=0.0D0
         O2_PR(:,:)=0.0D0
         ALLOCATE(C2MI_RR(NMOD,ND))
         ALLOCATE(C2MI_PR(NMOD,ND))
         C2MI_RR(:,:)=0.0D0
         C2MI_PR(:,:)=0.0D0
         ALLOCATE(O2MI_RR(NMOD,ND))
         ALLOCATE(O2MI_PR(NMOD,ND))
         O2MI_RR(:,:)=0.0D0
         O2MI_PR(:,:)=0.0D0
         ALLOCATE(CO2MI_RR(NMOD,ND))
         ALLOCATE(CO2MI_PR(NMOD,ND))
         CO2MI_RR(:,:)=0.0D0
         CO2MI_PR(:,:)=0.0D0
         ALLOCATE(RHO(ND))
         RHO(:)=0.0D0
         ALLOCATE(INT_ARRAY(ND))
         INT_ARRAY(:)=0

!
!     Read in temperature and ED data for all models
!
         DO I=1,NMOD
            FILE1 = TRIM(ADJUSTL(DIR_LIST(I)))//'RVTJ'
!     WRITE(6,*) 'FILE1:',FILE1
            OPEN(2,FILE=TRIM(ADJUSTL(FILE1)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
            IF (IOS .NE. 0) THEN
               WRITE(6,*) 'Error opening file ',FILE1
               STOP
            END IF
!
            BOOL1=.FALSE.
            BOOL2=.FALSE.
            IF(I .EQ. 1) BOOL3=.FALSE.
            DO 
               READ(2,'(A)') TSTR1
               IF (INDEX(TSTR1,'Temperature (10^4K)') .NE. 0) THEN
                  READ(2,*) (T(I,J),J=1,ND)
                  BOOL1=.TRUE.
               END IF
               IF (INDEX(TSTR1,'Electron density') .NE. 0) THEN
                  READ(2,*) (ED_POPS(I,J),J=1,ND)
                  BOOL2=.TRUE.
               END IF
               IF (I .EQ. 1) THEN
                  IF (INDEX(TSTR1,'Mass Density') .NE. 0) THEN
                     READ(2,*) (RHO(J),J=1,ND)
                     BOOL3=.TRUE.
!                     WRITE(*,*) RHO
                  END IF
               END IF
               IF (BOOL1 .AND. BOOL2 .AND. BOOL3) EXIT
            END DO
            CLOSE(2)
         END DO
!
!     Read in cooling data for all models
!
         DO I=1,NMOD
            FILE1 = TRIM(ADJUSTL(DIR_LIST(I)))//'GENCOOL'
            OPEN(2,FILE=TRIM(ADJUSTL(FILE1)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
            IF (IOS .NE. 0) THEN
               WRITE(6,*) 'Error opening file ',FILE1
               STOP
            END IF
            J=1
            DO
               IF (J .GT. ND) EXIT
               DO
                  READ(2,'(A)') TSTR1
                  IF (INDEX(TSTR1,'Radiative decay') .NE. 0) EXIT
               END DO
               IF (ND-J .GE. 10) THEN
                  READ(2,'(1X,1P,10E12.4)') RAD_HEAT(I,J:J+9)
               ELSE
                  READ(2,'(1X,1P,10E12.4)') RAD_HEAT(I,J:ND)
               END IF
               DO
                  READ(2,'(A)') TSTR1
                  IF (INDEX(TSTR1,'Carbon Monoxide') .NE. 0) EXIT
               END DO
!     WRITE(6,*) TSTR1
               IF (ND-J .GE. 10) THEN
                  READ(2,'(1X,1P,10E12.4)') CO_COOLING(I,J:J+9)
               ELSE
                  READ(2,'(1X,1P,10E12.4)') CO_COOLING(I,J:ND)
               END IF
               DO
                  READ(2,'(A)') TSTR1
                  IF (INDEX(TSTR1,'Net Cooling') .NE. 0) EXIT
               END DO
!     WRITE(6,*) TSTR1
               IF (ND-J .GE. 10) THEN
                  READ(2,'(1X,1P,10E12.4)') NET_COOLING(I,J:J+9)
                  J=J+10
               ELSE
                  READ(2,'(1X,1P,10E12.4)') NET_COOLING(I,J:ND)
                  J=J+10
               END IF
            END DO
            CLOSE(2)
         END DO
         TOT_COOLING(:,:)=NET_COOLING(:,:)+RAD_HEAT(:,:)
!
!     Read in populations
!
         DO I=1,NMOD
            FILE1 = TRIM(ADJUSTL(DIR_LIST(I)))//'MOL_RXNS_BY_SPECIES'
            OPEN(2,FILE=TRIM(ADJUSTL(FILE1)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
            DO J=1,ND
               K=1
               DO WHILE(K .LE. 12)

                  READ(2,'(A)') TSTR1
!
                  IF (INDEX(TSTR1,'CI(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A7,ES15.5,A)') TSTR2,CI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,CI_PR(I,J),TSTR3,CI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'C2(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A7,ES15.5,A)') TSTR2,C2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,C2_PR(I,J),TSTR3,C2_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'OI(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A7,ES15.5,A)') TSTR2,OI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,OI_PR(I,J),TSTR3,OI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'O2(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A7,ES15.5,A)') TSTR2,O2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,O2_PR(I,J),TSTR3,O2_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'COMI(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,COMI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,COMI_PR(I,J),TSTR3,COMI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'COM2(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,COM2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'C2MI(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,C2MI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,C2MI_PR(I,J),TSTR3,C2MI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'C2M2(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,C2M2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'O2MI(POP=') .EQ. 1) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,O2MI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,O2MI_PR(I,J),TSTR3,O2MI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'O2M2(POP=') .EQ. 1) THEN
                     READ(TSTR1,'(A9,ES15.5,A)') TSTR2,O2M2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'CO2MI(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A10,ES15.5,A)') TSTR2,CO2MI_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A5,ES15.5,T50,A5,ES15.5)') TSTR2,CO2MI_PR(I,J),TSTR3,CO2MI_RR(I,J)
                     K=K+1
!
                  ELSE IF (INDEX(TSTR1,'CO2M2(POP=') .NE. 0) THEN
                     READ(TSTR1,'(A10,ES15.5,A)') TSTR2,CO2M2_POPS(I,J),TSTR3
                     READ(2,'(A)') TSTR2
                     DO L=1,5
                        READ(2,'(A)') TSTR2
                     END DO
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     READ(2,'(A)') TSTR2
                     K=K+1
!
                  END IF
               END DO
            END DO
            CLOSE(2)
         END DO
!     
!     Make the list of reactions
!
         FILE1=TRIM(ADJUSTL(DIR_LIST(1)))//'MOLECULAR_REACTIONS_SORTED'
         OPEN(2,FILE=TRIM(ADJUSTL(FILE1)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
         NRXNS=0
         NRXNS_DISTINCT=0
         DO
            READ(2,'(A)') TSTR1
            IF (INDEX(TSTR1,'T (10^4 K)') .NE. 0) EXIT
         END DO
         READ(2,'(A)') TSTR1
         DO 
            READ(2,'(A)') TSTR2
            IF (INDEX(TSTR2,':') .EQ. 0) EXIT
            NRXNS=NRXNS+1
            IF (TSTR2(1:27) .NE. TSTR1(1:27)) NRXNS_DISTINCT=NRXNS_DISTINCT+1
            TSTR1 = TSTR2
         END DO
         ALLOCATE(RXNS_LIST(NRXNS_DISTINCT))
         DO
            READ(2,'(A)') TSTR1
            IF (INDEX(TSTR1,'T (10^4 K)') .NE. 0) EXIT
         END DO
         READ(2,'(A)') TSTR1
         J=1
         DO
            READ(2,'(A)') TSTR2
            IF (INDEX(TSTR2,':') .EQ. 0) EXIT
            IF (INDEX(TSTR2,TSTR1(1:27)) .EQ. 0) THEN
               RXNS_LIST(J) = TSTR2(1:27)
               J=J+1
            END IF
            IF (J .GT. NRXNS_DISTINCT) EXIT
            TSTR1=TSTR2
         END DO
!     WRITE(6,*) 'NRXNS: ',NRXNS
!     WRITE(6,*) 'NRXNS_DISTINCT: ',NRXNS_DISTINCT
!     DO I=1,NRXNS_DISTINCT
!     WRITE(6,*) TRIM(ADJUSTL(RXNS_LIST(I)))
!     END DO
!
! Set defaults
!
         DEPTH=1
         XLABEL='Temperature (10\u4\d K)'
         DO I=1,NMOD
            X(I) = T(I,DEPTH)
         END DO
         X1(1:NMOD)=X(1:NMOD)
         ! Now we sort X; the sorting indices are stored in X_ARGSORT;
         ! X1 is destroyed
         DO I=1,NMOD
            T1=MAXVAL(X1(1:NMOD))+1.0D0
            L=MINLOC(X1(1:NMOD),1)
            X(I)=X1(L)
            X_ARGSORT(I)=L
            X1(L)=T1
         END DO
         CHOSEN_RXN=1
!
         DO
            WRITE(*,*) ' '
!
            PLT_OPT='P'
!
            WRITE(*,'(A,I5)') 'Depth ',DEPTH
            CALL GEN_IN(PLT_OPT,'Plot Option')
!     
            IF (INDEX(UC(PLT_OPT),'D_') .EQ. 1) THEN
               TSTR1=PLT_OPT(3:)
               READ(TSTR1,*,IOSTAT=IOS) DEPTH
               IF (IOS .NE. 0) THEN
                  WRITE(*,*) 'Error reading depth value, try again'
                  CYCLE
               END IF
               !If the x-axis is set to be temperature, then go ahead
               !and set the values for the chosen depth
               !Other choices for x-axis will need to be evoked explicitly
               IF (XLABEL .EQ. 'Temperature (10\u4\d K)') THEN
                  DO I=1,NMOD
                     X(I) = T(I,DEPTH)
                  END DO
                  X1(1:NMOD)=X(1:NMOD)
                  ! Now we sort X; the sorting indices are stored in X_ARGSORT;
                  ! X1 is destroyed
                  DO I=1,NMOD
                     T1=MAXVAL(X1(1:NMOD))+1.0D0
                     L=MINLOC(X1(1:NMOD),1)
                     X(I)=X1(L)
                     X_ARGSORT(I)=L
                     X1(L)=T1
                  END DO
               END IF
!               
            ELSE IF (UC(PLT_OPT) .EQ. 'XT') THEN
               XLABEL='Temperature (10\u4\d K)'
               DO I=1,NMOD
                  X(I) = T(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XOI') THEN
               XLABEL='OI Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = OI_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XO2') THEN
               XLABEL='O2 Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = O2_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XCI') THEN
               XLABEL='CI Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = CI_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XC2') THEN
               XLABEL='C2 Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = C2_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XCOI') THEN
               XLABEL='CI + OI Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = CI_POPS(I,DEPTH)+OI_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XCO2') THEN
               XLABEL='C2 + O2 Population (cm\u-3\d)'
               DO I=1,NMOD
                  X(I) = C2_POPS(I,DEPTH)+O2_POPS(I,DEPTH)
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XMOD') THEN
               XLABEL='Model'
               DO I=1,NMOD
                  X(I)=I
                  X_ARGSORT(I) = I
               END DO
!
            ELSE IF (UC(PLT_OPT) .EQ. 'XAGE') THEN
               XLABEL='Age'
               DO I=1,NMOD
                  TSTR1=TRIM(ADJUSTL(DIR_LIST(I)))//'VADAT'
                  OPEN(UNIT=22,FILE=TSTR1,ACTION='READ')
                  DO
                     READ(22,'(A)',IOSTAT=IOS) TSTR2
                     IF (IOS .NE. 0) THEN
                        WRITE(*,*) 'Error reading age from ',TSTR1
                        WRITE(*,*) 'String is: ',TSTR2
                        EXIT
                     END IF
                     IF (INDEX(TSTR2,'SN_AGE') .NE. 0) THEN
                        READ(TSTR2,*,IOSTAT=IOS) X(I), TSTR3
                        IF (IOS .NE. 0) THEN
                           WRITE(*,*) 'Error reading age from string: ',TSTR2
                        END IF
                        EXIT
                     END IF
                  END DO
               END DO
               X1(1:NMOD)=X(1:NMOD)
               ! Now we sort X; the sorting indices are stored in X_ARGSORT;
               ! X1 is destroyed
               DO I=1,NMOD
                  T1=MAXVAL(X1(1:NMOD))+1.0D0
                  L=MINLOC(X1(1:NMOD),1)
                  X(I)=X1(L)
                  X_ARGSORT(I)=L
                  X1(L)=T1
               END DO
               !
            ELSE IF (UC(PLT_OPT) .EQ. 'COCOOL') THEN
               DO I=1,NMOD
                  Y(I) = CO_COOLING(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='CO Cooling (erg s\u-1\d cm\u-3\d)'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'COPOP') THEN
               DO I=1,NMOD
                  Y(I) = COMI_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='CO Population (cm\u-3\d)'
               
            ELSE IF (UC(PLT_OPT) .EQ. 'COPOP_EVOL') THEN
               IF (XLABEL .NE. 'Age') THEN
                  WRITE(*,*) 'Must use ''Age'' x-axis with this option'
                  CYCLE
               END IF
               DO I=1,NMOD
                  Y(I) = COMI_POPS(X_ARGSORT(I),DEPTH)*(X(I)**3)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='CO Density scaled by time**3'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'CIPOP') THEN
               DO I=1,NMOD
                  Y(I) = CI_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='CI Population (cm\u-3\d)'
!
            ELSE IF (UC(PLT_OPT) .EQ. 'C2POP') THEN
               DO I=1,NMOD
                  Y(I) = C2_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='C2 Population (cm\u-3\d)'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'OIPOP') THEN
               DO I=1,NMOD
                  Y(I) = OI_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='OI Population (cm \u -3 \d )'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'O2POP') THEN
               DO I=1,NMOD
                  Y(I) = O2_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='O2 Population (cm \u -3 \d )'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'COCOOLPER') THEN
               DO I=1,NMOD
                  Y(I) = CO_COOLING(X_ARGSORT(I),DEPTH)/COMI_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='CO Cooling per Molecule (erg s\u-1\d cm\u-3\d per molecule)'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'TOTCOOL') THEN
               DO I=1,NMOD
                  Y(I) = TOT_COOLING(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Total cooling (erg s\u-1\d cm\u-3\d)'
!
            ELSE IF (UC(PLT_OPT) .EQ. 'NETCOOL') THEN
               DO I=1,NMOD
                  Y(I) = NET_COOLING(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Net Cooling (erg s\u-1\d cm\u-3\d)'
!
            ELSE IF (UC(PLT_OPT) .EQ. 'RADHEAT') THEN
               DO I=1,NMOD
                  Y(I) = RAD_HEAT(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Radioactive Decay Heating (erg s\u-1\d cm\u-3\d)'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'OTHERCOOL') THEN
               DO I=1,NMOD
                  Y(I) = TOT_COOLING(X_ARGSORT(I),DEPTH)-CO_COOLING(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Non-CO cooling (erg s\u-1\d cm\u-3\d)'
!
            ELSE IF (UC(PLT_OPT) .EQ. 'YTEMP') THEN
               DO I=1,NMOD
                  Y(I) = T(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Temperature (10\u4\d K)'
!  '   
            ELSE IF (INDEX(UC(PLT_OPT),'RXNRATE') .NE. 0) THEN
               WRITE(6,*) 'Choose a reaction: '
               DO I=1,NRXNS_DISTINCT
                  WRITE(6,'(I3,A1,2X,A)') I,':',RXNS_LIST(I)
               END DO
               CALL GEN_IN(CHOSEN_RXN,'Chosen reaction ID')
               TSTR1=RXNS_LIST(CHOSEN_RXN)
               DO I=1,NMOD
                  FILE1=TRIM(ADJUSTL(DIR_LIST(I)))//'MOLECULAR_REACTIONS_SORTED'
! '
                  RXN_RATES1(I,:)=0.0D0
                  X1(:)=0.0D0
                  OPEN(2,FILE=TRIM(ADJUSTL(FILE1)),STATUS='OLD',IOSTAT=IOS,ACTION='READ')
                  REWIND(2)
                  IF (IOS .NE. 0) THEN
                     WRITE(6,*) 'Error opening file ',FILE1
                     STOP
                  END IF
                  K=1
                  DO WHILE(K .LE. ND)
                     DO
                        READ(2,'(A)') TSTR2
                        IF (INDEX(TSTR2,'T (10^4 K)') .NE. 0) EXIT
                     END DO
                     IF (K .LE. ND-10) READ(2,'(A)') TSTR2
                     DO J=1,NRXNS
                        READ(2,'(A)') TSTR2
                        IF (INDEX(TSTR2,TRIM(ADJUSTL(TSTR1))) .NE. 0) THEN
                           READ(TSTR2,'(A28,T30,10ES15.5)') X1(1:10)
                           RXN_RATES1(I,K:K+9)=RXN_RATES1(I,K:K+9)+X1(1:10)
                        END IF
                     END DO
                     K=K+10
                  END DO
                  CLOSE(2)
               END DO
               QUERY='N'
               CALL GEN_IN(QUERY,'Plot just the rate constant?')
               IF (UC(QUERY) .EQ. 'Y') THEN
                  READ(TSTR1(1:12),'(A5,A2,A5)') RCT1,TSTR2,RCT2
                  WRITE(*,*) 'Reactants read: ',RCT1,' ',RCT2
                  DO I=1,NMOD
                     SELECT CASE (TRIM(ADJUSTL(RCT1)))
                        CASE ('CI')
                           T1 = CI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2')
                           T1 = C2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('OI')
                           T1 = OI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2')
                           T1 = O2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('COMI')
                           T1 = COMI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('COM2')
                           T1 = COM2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2MI')
                           T1 = C2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2M2')
                           T1 = C2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2MI')
                           T1 = O2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2M2')
                           T1 = O2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('CO2MI')
                           T1 = CO2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('CO2M2')
                           T1 = CO2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('ELEC')
                           T1 = ED_POPS(X_ARGSORT(I),DEPTH)
                        CASE DEFAULT
                           IF (I .EQ. 1) WRITE(*,*) 'Error reading reactants from string--will plot the total rate'
                           T1=1
                        END SELECT
!
                     SELECT CASE (TRIM(ADJUSTL(RCT2)))
                        CASE ('CI')
                           T2 = CI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2')
                           T2 = C2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('OI')
                           T2 = OI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2')
                           T2 = O2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('COMI')
                           T2 = COMI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('COM2')
                           T2 = COM2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2MI')
                           T2 = C2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('C2M2')
                           T2 = C2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2MI')
                           T2 = O2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('O2M2')
                           T2 = O2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('CO2MI')
                           T2 = CO2MI_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('CO2M2')
                           T2 = CO2M2_POPS(X_ARGSORT(I),DEPTH)
                        CASE ('ELEC')
                           T2 = ED_POPS(X_ARGSORT(I),DEPTH)
                        CASE DEFAULT
                           IF (I .EQ. 1) WRITE(*,*) 'Error reading reactants from string--will plot the total rate'
                           T2=1
                        END SELECT
!
                     Y(I) = (RXN_RATES1(X_ARGSORT(I),DEPTH)/T1)/T2
                  END DO
               ELSE 
                  CALL GEN_IN(QUERY,'Scale rates by CO abundance? (i.e., Liljegren 2020)')
                  IF (UC(QUERY) .EQ. 'Y') THEN
                     DO I=1,NMOD
                        T1 = COMI_POPS(X_ARGSORT(I),DEPTH)
                        T1 = T1 + COM2_POPS(X_ARGSORT(I),DEPTH)
                        T1 = T1 + CI_POPS(X_ARGSORT(I),DEPTH)
                        T1 = T1 + C2_POPS(X_ARGSORT(I),DEPTH)
                        T1 = T1 + OI_POPS(X_ARGSORT(I),DEPTH)
                        T1 = T1 + O2_POPS(X_ARGSORT(I),DEPTH)
                        Y(I) = RXN_RATES1(X_ARGSORT(I),DEPTH)/T1
                     END DO
                  ELSE
                     DO I=1,NMOD
                        Y(I) = RXN_RATES1(X_ARGSORT(I),DEPTH)
                     END DO
                  END IF
               END IF
               IF (INDEX(UC(PLT_OPT),'_EVOL') .NE. 0) THEN
                  IF (XLABEL .NE. 'Age') THEN
                     WRITE(*,*) 'Must use ''Age'' x-axis with this option'
                     CYCLE
                  END IF
                  DO I=1,NMOD
                     Y(I) = Y(I)*(X(I)**3)
                  END DO
               END IF
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='RXN Rate: '//TRIM(ADJUSTL(TSTR1))
               
            ELSE IF (UC(PLT_OPT) .EQ. 'ELEC') THEN
               DO I=1,NMOD
                  Y(I) = ED_POPS(X_ARGSORT(I),DEPTH)
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Electron Density (cm\u-3\d)'               
!
            ELSE IF (UC(PLT_OPT) .EQ. 'THIN') THEN
               WRITE(6,*) 'Plot optically thin CO cooling'
               THIN_COOL(:)=0.0D0
               LTE_POPS(:)=0.0D0
               PART_F=0.0D0
               DO I=1,NMOD
                  PART_F=0.0D0
                  DO L=1,NLEVELS
                     LTE_POPS(L)=COMI_LEVEL_LIST(L)%LEVEL_G*&
                          EXP(-COMI_LEVEL_LIST(L)%ENERGY_CM1*CM1TOEV/(1.0D4*T(I,DEPTH)*KBEV))
                     PART_F=PART_F+LTE_POPS(L)
                  END DO
                  LTE_POPS(:)=(COMI_POPS(I,DEPTH)/PART_F)*LTE_POPS(:)
!                  WRITE(*,*) 'Calculated LTE population total: ',SUM(LTE_POPS)
!                  WRITE(*,*) 'Population total (read from file): ',CO_POPS(I,DEPTH)
                  DO L=1,NTRANS
                     K=COMI_TRANS_LIST(L)%UPPER_LEVEL_ID
                     THIN_COOL(I)=THIN_COOL(I)+&
                     (EMIT_CONST/COMI_TRANS_LIST(L)%LAMBDA_ANG)*(COMI_TRANS_LIST(L)%EINSTEIN_A*LTE_POPS(K))
                  END DO
                  Y(I)=THIN_COOL(I)
               END DO
               Y1(1:NMOD)=Y(1:NMOD)
               DO I=1,NMOD
                  Y(I) = Y1(X_ARGSORT(I))
               END DO
               CALL DP_CURVE(NMOD,X,Y)
               YLABEL='Optically Thin CO Cooling (erg s\d-1\d cm\u-3\d)'
!
               ELSE IF (UC(PLT_OPT) .EQ. 'COION') THEN
                  DO I=1,NMOD
                     Y(I) = COM2_POPS(X_ARGSORT(I),DEPTH)/(COMI_POPS(X_ARGSORT(I),DEPTH)+COM2_POPS(X_ARGSORT(I),DEPTH))
                  END DO
                  CALL DP_CURVE(NMOD,X,Y)
                  YLABEL='CO Ion Fraction'
!                  
               ELSE IF (INDEX(UC(PLT_OPT), 'IONFRAC') .NE. 0) THEN
                  I = INDEX(UC(PLT_OPT), '_')
                  IF (I .EQ. 0) THEN
                     WRITE(*,*) 'Species not recognized'
                  ELSE
                     TSTR1 = UC(PLT_OPT(I+1:))
                  END IF
                  SELECT CASE(TRIM(ADJUSTL(TSTR1)))
                     CASE ('CI')
                        DO I=1,NMOD
                           Y(I) = CI_POPS(X_ARGSORT(I),DEPTH)/(CI_POPS(X_ARGSORT(I),DEPTH)+C2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral Carbon Fraction'
                     CASE ('C2')
                        DO I=1,NMOD
                           Y(I) = C2_POPS(X_ARGSORT(I),DEPTH)/(CI_POPS(X_ARGSORT(I),DEPTH)+C2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized Carbon Fraction'
                     CASE ('OI')
                        DO I=1,NMOD
                           Y(I) = OI_POPS(X_ARGSORT(I),DEPTH)/(OI_POPS(X_ARGSORT(I),DEPTH)+O2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral Oxygen Fraction'
                     CASE ('O2')
                        DO I=1,NMOD
                           Y(I) = O2_POPS(X_ARGSORT(I),DEPTH)/(OI_POPS(X_ARGSORT(I),DEPTH)+O2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized Oxygen Fraction'
                     CASE ('COMI')
                        DO I=1,NMOD
                           Y(I) = COMI_POPS(X_ARGSORT(I),DEPTH)/(COMI_POPS(X_ARGSORT(I),DEPTH)+COM2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral CO Fraction'
                     CASE ('COM2')
                        DO I=1,NMOD
                           Y(I) = COM2_POPS(X_ARGSORT(I),DEPTH)/(COMI_POPS(X_ARGSORT(I),DEPTH)+COM2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized CO Fraction'
                     CASE ('C2MI')
                        DO I=1,NMOD
                           Y(I) = C2MI_POPS(X_ARGSORT(I),DEPTH)/(C2MI_POPS(X_ARGSORT(I),DEPTH)+C2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral C_2 Fraction'
                     CASE ('C2M2')
                        DO I=1,NMOD
                           Y(I) = C2M2_POPS(X_ARGSORT(I),DEPTH)/(C2MI_POPS(X_ARGSORT(I),DEPTH)+C2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized C_2 Fraction'
                     CASE ('O2MI')
                        DO I=1,NMOD
                           Y(I) = O2MI_POPS(X_ARGSORT(I),DEPTH)/(O2MI_POPS(X_ARGSORT(I),DEPTH)+O2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral O_2 Fraction'
                     CASE ('O2M2')
                        DO I=1,NMOD
                           Y(I) = O2M2_POPS(X_ARGSORT(I),DEPTH)/(O2MI_POPS(X_ARGSORT(I),DEPTH)+O2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized O_2 Fraction'
                     CASE ('CO2MI')
                        DO I=1,NMOD
                           Y(I) = CO2MI_POPS(X_ARGSORT(I),DEPTH)/(CO2MI_POPS(X_ARGSORT(I),DEPTH)+CO2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Neutral CO_2 Fraction'
                     CASE ('CO2M2')
                        DO I=1,NMOD
                           Y(I) = CO2M2_POPS(X_ARGSORT(I),DEPTH)/(CO2MI_POPS(X_ARGSORT(I),DEPTH)+CO2M2_POPS(X_ARGSORT(I),DEPTH))
                        END DO
                        YLABEL='Ionized CO_2 Fraction'
                     CASE DEFAULT
                        WRITE(*,*) 'Species not recognized'
                  END SELECT
                  CALL DP_CURVE(NMOD,X,Y)
!
               ELSE IF (UC(PLT_OPT) .EQ. 'T3') THEN
                  IF (XLABEL .NE. 'Age') THEN
                     WRITE(*,*) 'Must use ''Age'' x-axis with this option'
                     CYCLE
                  END IF
                  DO I=1,NMOD
                     Y(I) = X(I)**3
                  END DO
!
               ELSE IF (UC(PLT_OPT(1:4)) .EQ. 'PEAK') THEN
                  L = INDEX(PLT_OPT,'_')
                  IF ((UC(PLT_OPT(L+1:)) .EQ. 'MASS') .OR. L .EQ. 0) THEN
                     DO I=1,ND
                        X(I) = RHO(I)
                        X1(I) = RHO(I)
                     END DO
 !                    WRITE(*,*) X(1:ND)
 !                    WRITE(*,*) RHO
                     DO I=1,ND
                        T1=MAXVAL(X1(1:ND))+1.0D0
                        J=MINLOC(X1(1:ND),1)
                        X(I)=X1(J)
                        X_ARGSORT(I)=J
                        X1(J)=T1
                     END DO
!                     WRITE(*,*) X(1:ND)
                     XLABEL='Mass Density (g cm\u-3\d)'
                  ELSE IF (UC(PLT_OPT(L+1:)) .EQ. 'DEPTH') THEN
                     DO I=1,ND
                        X(I) = I
                        X1(I) = I
                     END DO
                     DO I=1,ND
                        T1=MAXVAL(X1(1:ND))+1.0D0
                        J=MINLOC(X1(1:ND),1)
                        X(I)=X1(J)
                        X_ARGSORT(I)=J
                        X1(J)=T1
                     END DO
!                     WRITE(*,*) X(1:ND)
                     XLABEL='Depth'
                  ELSE
                     WRITE(*,*) 'Unable to parse x-axis option for plotting PEAK'
                     CYCLE
                  END IF
                  INT_ARRAY = MAXLOC(CO_COOLING,1)
!                  WRITE(*,*) INT_ARRAY
                  DO I=1,ND
                     Y1(I) = T(INT_ARRAY(I),I)
                  END DO
!                  WRITE(*,*) Y1
                  DO I=1,ND
                     Y(I) = Y1(X_ARGSORT(I))
                  END DO
                  YLABEL='Peak temp of CO Cooling Curve'
                  QUERY='Y'
                  CALL GEN_IN(QUERY,'Plot only points where a peak is resolved?')
                  IF (UC(QUERY(1:1)) .EQ. 'Y') THEN
                     L=0
                     DO I=1,ND
                        T1 = MINVAL(T(1:NMOD,I))
                        T2 = MAXVAL(T(1:NMOD,I))
                        IF (Y(I) .GT. T1 .AND. Y(I) .LT. T2) THEN
                           L=L+1
                           X1(L) = X(I)
                           Y1(L) = Y(I)
                        END IF
                     END DO
                     CALL DP_CURVE(L,X1,Y1)
                  ELSE
                     CALL DP_CURVE(ND,X,Y)
                  END IF
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'H') THEN
               WRITE(6,*) ' '
               WRITE(6,*) 'D_n: Set the current depth to n'
               WRITE(6,*) 'XT: Plot versus temperature'
               WRITE(6,*) 'XOI: Plot versus Oxygen I population-replace O with C for carbon'
               WRITE(6,*) 'XO2: Plot versus Oxygen 2 (singly ionized) population-replace O with C for carbon'
               WRITE(6,*) 'XCOI: Plot versus the sum of neutral oxygen and carbon populations'
               WRITE(6,*) 'XCO2: Plot versus sum of ionized carbon and oxygen populations'
               WRITE(6,*) 'COCOOL: Plot CO cooling across all models'
               WRITE(6,*) 'COPOP: Plot total CO population across all models'
               WRITE(6,*) 'CIPOP: Plot neutral carbon population'
               WRITE(6,*) 'OIPOP: Plot neutral Oxygen population'
               WRITE(6,*) 'COCOOLPER: Plot CO cooling per molecule (COCOOL/COPOP)'
               WRITE(6,*) 'TOTCOOL: Plot total cooling'
               WRITE(6,*) 'NETCOOL: Plot net cooling (cooling-heating)'
               WRITE(*,*) 'RADHEAT: Plot radiative decay heating rate-should be constant across models'
               WRITE(6,*) 'OTHERCOOL: Plot non-CO cooling'
               WRITE(6,*) 'RXNRATE: Plot the rate of a particular reaction'
               WRITE(6,*) '         Choose from a list of all the reactions'
               WRITE(6,*) '         Includes option to scale by reactants (i.e., the rate constant)'
               WRITE(6,*) 'THIN: Plot optically thin CO cooling rate'
               WRITE(6,*) 'COION: Plot CO ionization fraction'
               WRITE(6,*) 'PEAK_DEPTH: Plot the temperature at which the CO cooling function peaks, as a function of depth'
               WRITE(6,*) 'PEAK_MASS: Same as above, plot as a function of mass density'
               WRITE(6,*) 'P: call plotting routine'
               WRITE(6,*) 'E   : Exit'
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'P') THEN
               CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!     
            ELSE IF (UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E') THEN
               STOP
!     '
            END IF
         END DO
!
      ELSE IF (UC(QUERY) .EQ. 'S') THEN
         LU_CO=42
         INQUIRE(UNIT=LU_CO,EXIST=IS_OPEN)
         IF (IS_OPEN) THEN
            OPEN(LU_CO,STATUS='OLD',ACTION='READ')
         ELSE
            WRITE(6,*) 'fort.42 (CO cooling file) not found'
         END IF
         MIN_ITERS=1
         MAX_ITERS=1
         J=1
         DO
            READ(LU_CO,'(A)',IOSTAT=IOS) STRING
            IF (IOS .NE. 0) EXIT
            IF (LEN(TRIM(ADJUSTL(STRING))) .GE. 1 .AND. LEN(TRIM(ADJUSTL(STRING))) .LT. 6) THEN
               READ(STRING,'(I4)') I
               IF (J .LT. 2) THEN
                  MIN_ITERS=I
                  J=3
               END IF
               IF (I .GT. MAX_ITERS) MAX_ITERS=I
            END IF
         END DO
         WRITE(6,*) 'Earliest iteration: ',MIN_ITERS
         WRITE(6,*) 'Last iteration: ',MAX_ITERS
!
         OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
         IF (IOS .EQ. 0) THEN
            DO
               READ(20,'(A)',IOSTAT=IOS) STRING
               IF (IOS .NE. 0) EXIT
               IF (INDEX(STRING,'!Number of depth points') .NE. 0) THEN
                  READ(STRING,*)ND
                  WRITE(6,'(A,I4)') 'Number of depth points in the model is ',ND
! '
                  EXIT
               END IF
            END DO
         END IF
         INQUIRE(UNIT=20,OPENED=IS_OPEN)
         IF(IS_OPEN) CLOSE(20)
!
         ALLOCATE(CO_SCR_DATA(MAX_ITERS,ND,3),STAT=IOS)
         J = MAX(ND,MAX_ITERS)
         ALLOCATE(X(2*J))
         ALLOCATE(Y(2*J))
         ALLOCATE(X1(2*J))
         ALLOCATE(Y1(2*J))
         ALLOCATE(Z(2*J))
         IF (IOS .NE. 0) THEN
            WRITE(6,*) 'Error allocating space for data array CO_SCR_DATA'
! '
            STOP
         END IF
         CO_SCR_DATA(:,:,:) = 0.0D0
!
         REWIND(LU_CO)
         DO I=MIN_ITERS,MAX_ITERS
            READ(LU_CO,'(A)') STRING
            READ(LU_CO,'(I4)') J
            READ(LU_CO,'(A)') STRING
            DO J=1,ND
               READ(LU_CO,'(I4,ES12.4,ES12.4)') K,CO_SCR_DATA(I,J,2),CO_SCR_DATA(I,J,3)
            END DO
            READ(LU_CO,'(A)') STRING
         END DO
!
         DO
            WRITE(6,*) 'I: Plot CO cooling as a function of iteration (at chosen depth)'
            WRITE(6,*) 'N: Plot CO cooling as a function of depth (at chosen iteration)'
            WRITE(6,*) 'T: Plot CO cooling as a function of temperature (at chosen depth)'
            WRITE(6,*) 'E: Exit'
! '
            CO_OPT='P'
            CALL GEN_IN(CO_OPT,'CO Scratch option')
            IF (UC(CO_OPT) .EQ. 'I') THEN
               I=ND
               CALL GEN_IN(I,'Depth')
               X1(:) = 0.0D0
               Y1(:) = 0.0D0
               DO J=MIN_ITERS,MAX_ITERS
                  X1(J-MIN_ITERS+1) = J
                  Y1(J-MIN_ITERS+1) = CO_SCR_DATA(J,I,3)
               END DO
               CALL DP_CURVE(MAX_ITERS-MIN_ITERS+1,X1,Y1)
               WRITE(XLABEL,*) 'Iteration Number'
               WRITE(YLABEL,'(A,I4,A)') 'CO Cooling(depth ',I,')'
            ELSE IF(UC(CO_OPT) .EQ. 'N') THEN
               I=MIN_ITERS
               CALL GEN_IN(I,'Iteration')
               X1(:) = 0.0D0
               Y1(:) = 0.0D0
               DO J=1,ND
                  X1(J) = J
                  Y1(J) = CO_SCR_DATA(I,J,3)
               END DO
               CALL DP_CURVE(ND,X1,Y1)
               WRITE(XLABEL,*)'Depth Index'
               WRITE(YLABEL,'(A,I4,A)') 'CO Cooling(iteration ',I,')'
            ELSE IF (UC(CO_OPT) .EQ. 'T') THEN
               I=ND
               CALL GEN_IN(I,'Depth')
               X1(:) = 0.0D0
               Y1(:) = 0.0D0
               DO J=MIN_ITERS,MAX_ITERS
                  X1(J-MIN_ITERS+1) = CO_SCR_DATA(J,I,2)
                  Y1(J-MIN_ITERS+1) = CO_SCR_DATA(J,I,3)
               END DO
               K=MAX_ITERS-MIN_ITERS+1
               DO J=1,K
                  L=MINLOC(X1(J:K),1)
                  T1=X1(J)
                  X1(J)=X1(J+L-1)
                  X1(J+L-1)=T1
                  T2=Y(J)
                  Y(J)=Y(J+L-1)
                  Y(J+L-1)=T2
               END DO
               CALL DP_CURVE(K,X1,Y1)
               WRITE(XLABEL,*)'Temperature (10^4 K)'
               WRITE(YLABEL,'(A,I4,A)') 'CO Cooling (depth ',I,')'
            ELSE IF (UC(CO_OPT) .EQ. 'P') THEN
               CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
            ELSE IF (UC(CO_OPT) .EQ. 'E') THEN
               EXIT
            END IF
         END DO
!                  
      END IF
!     
      END
