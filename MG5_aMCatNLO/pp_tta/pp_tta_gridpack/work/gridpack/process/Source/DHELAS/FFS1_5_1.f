C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (Identity(2,1)) + Coup(2) * (Gamma5(2,1))
C     
      SUBROUTINE FFS1_5_1(F2, S3, COUP1, COUP2, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 F2(*)
      COMPLEX*16 S3(*)
      REAL*8 M1
      REAL*8 W1
      COMPLEX*16 COUP1
      COMPLEX*16 F1(8)
      COMPLEX*16 DENOM
      COMPLEX*16 COUP2
      COMPLEX*16 P1(0:3)
      F1(1) = +F2(1)+S3(1)
      F1(2) = +F2(2)+S3(2)
      F1(3) = +F2(3)+S3(3)
      F1(4) = +F2(4)+S3(4)
      P1(0) = -F1(1)
      P1(1) = -F1(2)
      P1(2) = -F1(3)
      P1(3) = -F1(4)
      DENOM = 1D0/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI*
     $  W1))
      F1(5)= DENOM*(-CI )* S3(5)*(COUP1*(F2(7)*(P1(0)+P1(3))+(F2(8)
     $ *(P1(1)+CI*(P1(2)))-F2(5)*M1))+COUP2*(F2(7)*(P1(0)+P1(3))+(F2(8)
     $ *(P1(1)+CI*(P1(2)))+F2(5)*M1)))
      F1(6)= DENOM*CI * S3(5)*(COUP1*(F2(7)*(+CI*(P1(2))-P1(1))+(F2(8)
     $ *(P1(3)-P1(0))+F2(6)*M1))+COUP2*(F2(7)*(+CI*(P1(2))-P1(1))
     $ +(F2(8)*(P1(3)-P1(0))-F2(6)*M1)))
      F1(7)= DENOM*(-CI )* S3(5)*(COUP1*(F2(5)*(P1(0)-P1(3))+(F2(6)*(
     $ -1D0)*(P1(1)+CI*(P1(2)))-F2(7)*M1))+COUP2*(F2(5)*(P1(3)-P1(0))
     $ +(F2(6)*(P1(1)+CI*(P1(2)))-F2(7)*M1)))
      F1(8)= DENOM*(-CI )* S3(5)*(COUP1*(F2(5)*(+CI*(P1(2))-P1(1))
     $ +(F2(6)*(P1(0)+P1(3))-F2(8)*M1))+COUP2*(F2(5)*(P1(1)-CI*(P1(2)))
     $ +(F2(6)*(-1D0)*(P1(0)+P1(3))-F2(8)*M1)))
      END


