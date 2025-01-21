      MODULE GRIEM_STARK_MOD
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
! Altered : 03-Dec-2023 - change continuation character to 1
!
        INTEGER, PARAMETER :: NBET=55
        INTEGER, PARAMETER :: NS=2*NBET-1
        REAL(KIND=LDP) BET(NBET)
        REAL(KIND=LDP) SS(NS)
        REAL(KIND=LDP) SX(NS)
        REAL(KIND=LDP) AS
        REAL(KIND=LDP) PS
        REAL(KIND=LDP) ODOP
!
!$OMP THREADPRIVATE(SS,SX,AS,PS,ODOP)
!
      DATA BET/ 0._LDP ,0.1_LDP ,0.2_LDP ,0.3_LDP ,0.4_LDP ,0.5_LDP , 0.6_LDP ,
     1         0.7_LDP ,0.8_LDP , 0.9_LDP ,1.0_LDP ,1.2_LDP ,1.4_LDP ,1.6_LDP ,
     1         1.8_LDP ,2.0_LDP , 2.2_LDP ,2.4_LDP ,2.6_LDP ,2.8_LDP ,3.0_LDP ,
     1         3.25_LDP ,3.5_LDP ,3.75_LDP ,4.0_LDP ,4.5_LDP ,5.0_LDP ,6.0_LDP ,
     1         8.0_LDP ,10.0_LDP ,12.5_LDP ,15._LDP ,17.5_LDP ,20.0_LDP ,
     1         25.0_LDP ,30.0_LDP ,40.0_LDP ,50.0_LDP ,60.0_LDP ,80.0_LDP ,
     1         1.0E2_LDP ,1.25E2_LDP ,1.5E2_LDP ,1.75E2_LDP ,2.0E2_LDP ,2.5E2_LDP ,3.0E2_LDP ,
     1         4.0E2_LDP ,5.0E2_LDP ,6.0E2_LDP, 8.0E2_LDP ,1.0E3_LDP ,1.5E3_LDP ,2.0E3_LDP ,2.5E3_LDP /
      END MODULE GRIEM_STARK_MOD
