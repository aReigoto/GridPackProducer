      SUBROUTINE LOOP_CT_CALLS_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=8)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=184, NLOOPGROUPS=89, NCTAMPS=282)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=466)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=35,NLOOPWAVEFUNCS=325)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE
     $ ,0:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.CTCALL_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     CutTools call for loop numbers 1,2,157,158,179,180,173,174,175,17
C     6,177,178
      CALL LOOP_2(6,13,DCMPLX(ZERO),DCMPLX(ZERO),2,I_SO,1)
C     CutTools call for loop numbers 3,4,5,6,159,160,181,182,183,184
      CALL LOOP_3(1,2,13,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(ZERO),3,I_SO
     $ ,2)
C     CutTools call for loop numbers 7,8
      CALL LOOP_2(6,13,DCMPLX(MDL_MB),DCMPLX(MDL_MB),2,I_SO,3)
C     CutTools call for loop numbers 9,10
      CALL LOOP_3(3,6,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,4)
C     CutTools call for loop numbers 11,12,13,14,31,32
      CALL LOOP_3(2,3,16,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,5)
C     CutTools call for loop numbers 15,16,17,18,26,29
      CALL LOOP_3(1,3,18,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,6)
C     CutTools call for loop numbers 19,20,23,24
      CALL LOOP_3(1,2,13,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,3,I_SO,7)
C     CutTools call for loop numbers 21,25
      CALL LOOP_4(1,2,15,3,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,8)
C     CutTools call for loop numbers 22,28
      CALL LOOP_4(1,2,3,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,9)
C     CutTools call for loop numbers 27,30
      CALL LOOP_4(1,3,2,15,DCMPLX(MDL_MB),DCMPLX(MDL_MB),DCMPLX(MDL_MB)
     $ ,DCMPLX(MDL_MB),4,I_SO,10)
C     CutTools call for loop numbers 33,72,124
      CALL LOOP_2(7,22,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,11)
C     CutTools call for loop numbers 34,38
      CALL LOOP_2(6,13,DCMPLX(MDL_MT),DCMPLX(MDL_MT),2,I_SO,12)
C     CutTools call for loop numbers 35
      CALL LOOP_3(5,6,7,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,13)
C     CutTools call for loop numbers 36,142,143,144
      CALL LOOP_3(5,6,7,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,14)
C     CutTools call for loop numbers 37,57,118
      CALL LOOP_2(8,23,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,15)
C     CutTools call for loop numbers 39
      CALL LOOP_3(4,6,8,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,16)
C     CutTools call for loop numbers 40,149,150,151
      CALL LOOP_3(4,6,8,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,17)
C     CutTools call for loop numbers 41,60,121
      CALL LOOP_3(3,5,23,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,18)
C     CutTools call for loop numbers 42
      CALL LOOP_4(3,5,4,6,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,19)
C     CutTools call for loop numbers 43,44
      CALL LOOP_3(3,6,15,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,20)
C     CutTools call for loop numbers 45
      CALL LOOP_4(3,4,5,6,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),3,I_SO,21)
C     CutTools call for loop numbers 46,154,155,156
      CALL LOOP_4(3,4,6,5,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,22)
C     CutTools call for loop numbers 47,75,127
      CALL LOOP_3(3,4,22,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,23)
C     CutTools call for loop numbers 48,125
      CALL LOOP_2(10,24,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,24)
C     CutTools call for loop numbers 49,58
      CALL LOOP_2(9,25,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,25)
C     CutTools call for loop numbers 50
      CALL LOOP_3(3,9,10,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,26)
C     CutTools call for loop numbers 51,130
      CALL LOOP_3(2,5,24,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,27)
C     CutTools call for loop numbers 52,56,66,69,135,136
      CALL LOOP_3(2,3,16,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,28)
C     CutTools call for loop numbers 53
      CALL LOOP_4(2,3,5,9,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,29)
C     CutTools call for loop numbers 54
      CALL LOOP_3(2,9,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,30)
C     CutTools call for loop numbers 55
      CALL LOOP_4(2,3,9,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,31)
C     CutTools call for loop numbers 59
      CALL LOOP_3(2,9,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,32)
C     CutTools call for loop numbers 61
      CALL LOOP_4(2,5,3,9,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,33)
C     CutTools call for loop numbers 62,138
      CALL LOOP_3(2,5,24,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,34)
C     CutTools call for loop numbers 63,119
      CALL LOOP_2(12,28,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,35)
C     CutTools call for loop numbers 64,73
      CALL LOOP_2(11,29,DCMPLX(ZERO),DCMPLX(MDL_MT),1,I_SO,36)
C     CutTools call for loop numbers 65
      CALL LOOP_3(3,12,11,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,37)
C     CutTools call for loop numbers 67
      CALL LOOP_3(2,7,11,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,38)
C     CutTools call for loop numbers 68
      CALL LOOP_4(2,3,4,11,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,39)
C     CutTools call for loop numbers 70,132
      CALL LOOP_3(2,4,28,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,40)
C     CutTools call for loop numbers 71
      CALL LOOP_4(2,3,11,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,41)
C     CutTools call for loop numbers 74
      CALL LOOP_3(2,7,11,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,42)
C     CutTools call for loop numbers 76
      CALL LOOP_4(2,4,3,11,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,43)
C     CutTools call for loop numbers 77,145
      CALL LOOP_3(2,4,28,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,44)
C     CutTools call for loop numbers 78,82,84,87,110,113
      CALL LOOP_3(1,3,18,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,45)
C     CutTools call for loop numbers 79
      CALL LOOP_4(1,3,5,12,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,46)
C     CutTools call for loop numbers 80,100
      CALL LOOP_3(1,5,29,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,47)
C     CutTools call for loop numbers 81
      CALL LOOP_3(1,12,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,48)
C     CutTools call for loop numbers 83
      CALL LOOP_4(1,3,12,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,49)
C     CutTools call for loop numbers 85
      CALL LOOP_3(1,7,10,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,50)
C     CutTools call for loop numbers 86
      CALL LOOP_4(1,3,4,10,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,51)
C     CutTools call for loop numbers 88
      CALL LOOP_4(1,3,10,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,52)
C     CutTools call for loop numbers 89,105
      CALL LOOP_3(1,4,25,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),2
     $ ,I_SO,53)
C     CutTools call for loop numbers 90
      CALL LOOP_5(1,2,4,5,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,54)
C     CutTools call for loop numbers 91
      CALL LOOP_5(1,2,3,4,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,55)
C     CutTools call for loop numbers 92,94,98,102
      CALL LOOP_3(1,2,13,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,3,I_SO,56)
C     CutTools call for loop numbers 93
      CALL LOOP_4(1,2,7,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,57)
C     CutTools call for loop numbers 95
      CALL LOOP_4(1,2,4,8,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,58)
C     CutTools call for loop numbers 96,109
      CALL LOOP_4(1,2,15,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,59)
C     CutTools call for loop numbers 97,112
      CALL LOOP_4(1,2,3,15,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,60)
C     CutTools call for loop numbers 99
      CALL LOOP_4(1,2,5,7,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,61)
C     CutTools call for loop numbers 101
      CALL LOOP_4(1,5,2,7,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,62)
C     CutTools call for loop numbers 103
      CALL LOOP_4(1,2,8,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,63)
C     CutTools call for loop numbers 104
      CALL LOOP_4(1,4,2,8,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),3,I_SO,64)
C     CutTools call for loop numbers 106
      CALL LOOP_5(1,2,5,4,3,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(MDL_MT),DCMPLX(MDL_MT),4,I_SO,65)
C     CutTools call for loop numbers 107
      CALL LOOP_5(1,3,2,4,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,66)
C     CutTools call for loop numbers 108
      CALL LOOP_5(1,3,4,2,5,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,67)
C     CutTools call for loop numbers 111,114
      CALL LOOP_4(1,3,2,15,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),4,I_SO,68)
C     CutTools call for loop numbers 115
      CALL LOOP_5(1,2,3,5,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,69)
C     CutTools call for loop numbers 116
      CALL LOOP_5(1,3,5,2,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT),DCMPLX(ZERO)
     $ ,DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,70)
C     CutTools call for loop numbers 117
      CALL LOOP_5(1,3,2,5,4,DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),DCMPLX(MDL_MT),4,I_SO,71)
C     CutTools call for loop numbers 120
      CALL LOOP_3(1,12,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,72)
C     CutTools call for loop numbers 122
      CALL LOOP_4(1,5,3,12,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,73)
C     CutTools call for loop numbers 123,139
      CALL LOOP_3(1,5,29,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,74)
C     CutTools call for loop numbers 126
      CALL LOOP_3(1,7,10,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,75)
C     CutTools call for loop numbers 128
      CALL LOOP_4(1,4,3,10,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,76)
C     CutTools call for loop numbers 129,146
      CALL LOOP_3(1,4,25,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(ZERO),2
     $ ,I_SO,77)
C     CutTools call for loop numbers 131
      CALL LOOP_4(1,5,2,7,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,78)
C     CutTools call for loop numbers 133
      CALL LOOP_4(1,4,2,8,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,79)
C     CutTools call for loop numbers 134
      CALL LOOP_5(1,4,3,2,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,80)
C     CutTools call for loop numbers 137
      CALL LOOP_5(1,4,2,3,5,DCMPLX(ZERO),DCMPLX(MDL_MT),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,81)
C     CutTools call for loop numbers 140
      CALL LOOP_4(1,2,7,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,82)
C     CutTools call for loop numbers 141
      CALL LOOP_4(1,2,5,7,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,83)
C     CutTools call for loop numbers 147
      CALL LOOP_4(1,2,8,4,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,84)
C     CutTools call for loop numbers 148
      CALL LOOP_4(1,2,4,8,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(ZERO),3,I_SO,85)
C     CutTools call for loop numbers 152
      CALL LOOP_5(1,2,5,3,4,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,86)
C     CutTools call for loop numbers 153
      CALL LOOP_5(1,2,4,3,5,DCMPLX(ZERO),DCMPLX(ZERO),DCMPLX(MDL_MT)
     $ ,DCMPLX(MDL_MT),DCMPLX(ZERO),4,I_SO,87)
C     CutTools call for loop numbers 161,162,163,164,165,166
      CALL LOOP_2_3(1,2,1,13,2,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,88)
C     CutTools call for loop numbers 167,168,169,170,171,172
      CALL LOOP_2_3(1,2,2,13,1,DCMPLX(ZERO),DCMPLX(ZERO),1,I_SO,89)
C     At this point, all reductions needed for (QCD=6), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 5000

      GOTO 1001
 5000 CONTINUE
      CTCALL_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

