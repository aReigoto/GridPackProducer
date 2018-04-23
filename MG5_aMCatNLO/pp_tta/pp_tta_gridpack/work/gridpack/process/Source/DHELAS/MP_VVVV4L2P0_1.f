C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
C     
      SUBROUTINE MP_VVVV4L2P0_1(P2, V3, V4, COUP, M1, W1, P1, COEFF)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 TMP2
      COMPLEX*32 V3(*)
      REAL*16 M1
      INCLUDE 'coef_specs.inc'
      COMPLEX*32 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 P2(0:3)
      REAL*16 W1
      COMPLEX*32 V4(*)
      COMPLEX*32 P1(0:3)
      COMPLEX*32 COUP
      P1(0) = +P2(0)+V3(1)+V4(1)
      P1(1) = +P2(1)+V3(2)+V4(2)
      P1(2) = +P2(2)+V3(3)+V4(3)
      P1(3) = +P2(3)+V3(4)+V4(4)
      TMP2 = (V3(5)*V4(5)-V3(6)*V4(6)-V3(7)*V4(7)-V3(8)*V4(8))
      COEFF(1,0,1)= COUP*(-CI*(V3(5)*V4(5))+CI*(TMP2))
      COEFF(2,0,1)= COUP*-CI * V3(6)*V4(5)
      COEFF(3,0,1)= COUP*-CI * V3(7)*V4(5)
      COEFF(4,0,1)= COUP*-CI * V3(8)*V4(5)
      COEFF(1,0,2)= COUP*CI * V3(5)*V4(6)
      COEFF(2,0,2)= COUP*(+CI*(V3(6)*V4(6)+TMP2))
      COEFF(3,0,2)= COUP*CI * V3(7)*V4(6)
      COEFF(4,0,2)= COUP*CI * V3(8)*V4(6)
      COEFF(1,0,3)= COUP*CI * V3(5)*V4(7)
      COEFF(2,0,3)= COUP*CI * V3(6)*V4(7)
      COEFF(3,0,3)= COUP*(+CI*(V3(7)*V4(7)+TMP2))
      COEFF(4,0,3)= COUP*CI * V3(8)*V4(7)
      COEFF(1,0,4)= COUP*CI * V3(5)*V4(8)
      COEFF(2,0,4)= COUP*CI * V3(6)*V4(8)
      COEFF(3,0,4)= COUP*CI * V3(7)*V4(8)
      COEFF(4,0,4)= COUP*(+CI*(V3(8)*V4(8)+TMP2))
      END


