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
      DATA BET/ 0.D0 ,0.1D0 ,0.2D0 ,0.3D0 ,0.4D0 ,0.5D0 , 0.6D0 ,
     1         0.7D0 ,0.8D0 , 0.9D0 ,1.0D0 ,1.2D0 ,1.4D0 ,1.6D0 ,
     1         1.8D0 ,2.0D0 , 2.2D0 ,2.4D0 ,2.6D0 ,2.8D0 ,3.0D0 ,
     1         3.25D0 ,3.5D0 ,3.75D0 ,4.0D0 ,4.5D0 ,5.0D0 ,6.0D0 ,
     1         8.0D0 ,10.0D0 ,12.5D0 ,15.D0 ,17.5D0 ,20.0D0 ,
     1         25.0D0 ,30.0D0 ,40.0D0 ,50.0D0 ,60.0D0 ,80.0D0 ,
     1         1.0D2 ,1.25D2 ,1.5D2 ,1.75D2 ,2.0D2 ,2.5D2 ,3.0D2 ,
     1         4.0D2 ,5.0D2 ,6.0D2, 8.0D2 ,1.0D3 ,1.5D3 ,2.0D3 ,2.5D3 /
      END MODULE GRIEM_STARK_MOD
