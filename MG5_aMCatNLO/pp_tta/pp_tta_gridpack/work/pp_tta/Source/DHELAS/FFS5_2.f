C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma5(2,1)
C     
      SUBROUTINE FFS5_2(F1, S3, COUP, M2, W2,F2)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(8)
      COMPLEX*16 S3(*)
      REAL*8 W2
      COMPLEX*16 P2(0:3)
      COMPLEX*16 F1(*)
      REAL*8 M2
      COMPLEX*16 DENOM
      COMPLEX*16 COUP
      F2(1) = +F1(1)+S3(1)
      F2(2) = +F1(2)+S3(2)
      F2(3) = +F1(3)+S3(3)
      F2(4) = +F1(4)+S3(4)
      P2(0) = -F2(1)
      P2(1) = -F2(2)
      P2(2) = -F2(3)
      P2(3) = -F2(4)
      DENOM = COUP/(P2(0)**2-P2(1)**2-P2(2)**2-P2(3)**2 - M2 * (M2 -CI
     $ * W2))
      F2(5)= DENOM*CI * S3(5)*(F1(7)*(P2(0)-P2(3))+(F1(8)*(+CI*(P2(2))
     $ -P2(1))-F1(5)*M2))
      F2(6)= DENOM*(-CI )* S3(5)*(F1(7)*(P2(1)+CI*(P2(2)))+(F1(8)*(
     $ -1D0)*(P2(0)+P2(3))+F1(6)*M2))
      F2(7)= DENOM*CI * S3(5)*(F1(5)*(-1D0)*(P2(0)+P2(3))+(F1(6)*(+CI
     $ *(P2(2))-P2(1))+F1(7)*M2))
      F2(8)= DENOM*(-CI )* S3(5)*(F1(5)*(P2(1)+CI*(P2(2)))+(F1(6)
     $ *(P2(0)-P2(3))-F1(8)*M2))
      END


