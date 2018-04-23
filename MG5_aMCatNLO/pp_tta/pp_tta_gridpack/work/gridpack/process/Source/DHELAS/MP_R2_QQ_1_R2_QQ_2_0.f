C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
Coup(1) * (P(-1,1)*Gamma(-1,2,1)) + Coup(2) * (Identity(1,2))
C     
      SUBROUTINE MP_R2_QQ_1_R2_QQ_2_0(F1, F2, COUP1, COUP2,VERTEX)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 TMP14
      COMPLEX*32 COUP2
      COMPLEX*32 TMP17
      COMPLEX*32 F1(*)
      COMPLEX*32 COUP1
      COMPLEX*32 F2(*)
      COMPLEX*32 P1(0:3)
      COMPLEX*32 VERTEX
      P1(0) = F1(1)
      P1(1) = F1(2)
      P1(2) = F1(3)
      P1(3) = F1(4)
      TMP14 = (F2(5)*F1(5)+F2(6)*F1(6)+F2(7)*F1(7)+F2(8)*F1(8))
      TMP17 = (F1(5)*(F2(7)*(P1(0)+P1(3))+F2(8)*(P1(1)+CI*(P1(2))))
     $ +(F1(6)*(F2(7)*(P1(1)-CI*(P1(2)))+F2(8)*(P1(0)-P1(3)))+(F1(7)
     $ *(F2(5)*(P1(0)-P1(3))-F2(6)*(P1(1)+CI*(P1(2))))+F1(8)*(F2(5)*(
     $ +CI*(P1(2))-P1(1))+F2(6)*(P1(0)+P1(3))))))
      VERTEX = (-1Q0)*(+CI*(COUP1*TMP17+TMP14*COUP2))
      END


