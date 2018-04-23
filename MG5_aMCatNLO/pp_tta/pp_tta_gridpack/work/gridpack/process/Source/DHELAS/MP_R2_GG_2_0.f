C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,1)*P(2,1)
C     
      SUBROUTINE MP_R2_GG_2_0(V1, V2, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 V2(*)
      COMPLEX*32 TMP12
      COMPLEX*32 P1(0:3)
      COMPLEX*32 VERTEX
      COMPLEX*32 COUP
      COMPLEX*32 V1(*)
      COMPLEX*32 TMP8
      P1(0) = V1(1)
      P1(1) = V1(2)
      P1(2) = V1(3)
      P1(3) = V1(4)
      TMP8 = (P1(0)*V2(5)-P1(1)*V2(6)-P1(2)*V2(7)-P1(3)*V2(8))
      TMP12 = (P1(0)*V1(5)-P1(1)*V1(6)-P1(2)*V1(7)-P1(3)*V1(8))
      VERTEX = COUP*(-CI * TMP8*TMP12)
      END


