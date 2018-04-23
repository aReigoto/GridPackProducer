      SUBROUTINE MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
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
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE
     $ ,0:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL MP_OXXXXX(P(0,4),MDL_MT,NHEL(4),+1*IC(4),W(1,4))
      CALL MP_IXXXXX(P(0,5),MDL_MT,NHEL(5),-1*IC(5),W(1,5))
      CALL MP_VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,6))
      CALL MP_FFS1_5_1(W(1,4),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,7))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),GC_5,AMP(1))
      CALL MP_FFS1_5_2(W(1,5),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,8))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),GC_5,AMP(2))
      CALL MP_FFV1_1(W(1,4),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,9))
      CALL MP_FFV1_2(W(1,5),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for born diagram with ID 3
      CALL MP_FFS1_5_0(W(1,10),W(1,9),W(1,3),GC_3010H,GC_3010A,AMP(3))
C     Amplitude(s) for born diagram with ID 4
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),GC_5,AMP(4))
      CALL MP_FFV1_2(W(1,5),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,11))
      CALL MP_FFV1_1(W(1,4),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,12))
C     Amplitude(s) for born diagram with ID 5
      CALL MP_FFS1_5_0(W(1,11),W(1,12),W(1,3),GC_3010H,GC_3010A,AMP(5))
C     Amplitude(s) for born diagram with ID 6
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),GC_5,AMP(6))
C     Amplitude(s) for born diagram with ID 7
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),GC_5,AMP(7))
C     Amplitude(s) for born diagram with ID 8
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),GC_5,AMP(8))
      CALL MP_FFV1P0_3(W(1,5),W(1,7),GC_5,ZERO,ZERO,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL MP_R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,1))
      CALL MP_R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,2))
      CALL MP_R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,3))
      CALL MP_R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,4))
      CALL MP_FFV1P0_3(W(1,8),W(1,4),GC_5,ZERO,ZERO,W(1,14))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL MP_R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,5))
      CALL MP_R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,6))
      CALL MP_R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,7))
      CALL MP_R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,8))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,9))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,10))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,11))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,12))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB,AMPL(1,13))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,14))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GT,AMPL(1,15))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,16))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GG_1EPS,AMPL(2,17))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,18))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,19))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,20))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,21))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,22))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,23))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,24))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,25))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB,AMPL(1,26))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,27))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GT,AMPL(1,28))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,29))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GG_1EPS,AMPL(2,30))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,31))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,32))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,33))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL MP_R2_GG_1_R2_GG_3_0(W(1,6),W(1,13),R2_GGQ,R2_GGB,AMPL(1,35)
     $ )
C     Counter-term amplitude(s) for loop diagram number 16
      CALL MP_R2_GG_1_R2_GG_3_0(W(1,6),W(1,14),R2_GGQ,R2_GGB,AMPL(1,36)
     $ )
      CALL MP_FFV1P0_3(W(1,5),W(1,4),GC_5,ZERO,ZERO,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL MP_VVS1_0(W(1,6),W(1,15),W(1,3),R2_GGX0_HB,AMPL(1,37))
      CALL MP_FFV1P0_3(W(1,5),W(1,9),GC_5,ZERO,ZERO,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL MP_VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGX0_HB,AMPL(1,38))
      CALL MP_FFV1P0_3(W(1,11),W(1,4),GC_5,ZERO,ZERO,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL MP_VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGX0_HB,AMPL(1,39))
      CALL MP_FFV1P0_3(W(1,5),W(1,12),GC_5,ZERO,ZERO,W(1,18))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL MP_VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGX0_HB,AMPL(1,40))
      CALL MP_FFV1P0_3(W(1,10),W(1,4),GC_5,ZERO,ZERO,W(1,19))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL MP_VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGX0_HB,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,43))
      CALL MP_VVV1P0_1(W(1,2),W(1,15),GC_4,ZERO,ZERO,W(1,20))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL MP_VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGX0_HB,AMPL(1,44))
      CALL MP_VVV1P0_1(W(1,1),W(1,15),GC_4,ZERO,ZERO,W(1,21))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL MP_VVS1_0(W(1,2),W(1,21),W(1,3),R2_GGX0_HB,AMPL(1,45))
      CALL MP_FFV1_2(W(1,5),W(1,6),GC_5,MDL_MT,MDL_WT,W(1,22))
C     Counter-term amplitude(s) for loop diagram number 41
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,22),W(1,7),R2_QQQ,R2_QQT,AMPL(1,46)
     $ )
      CALL MP_R2_QQ_2_0(W(1,22),W(1,7),UV_TMASS,AMPL(1,47))
      CALL MP_R2_QQ_2_0(W(1,22),W(1,7),UV_TMASS_1EPS,AMPL(2,48))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL MP_R2_GG_1_R2_GG_3_0(W(1,6),W(1,13),R2_GGQ,R2_GGT,AMPL(1,49)
     $ )
C     Counter-term amplitude(s) for loop diagram number 43
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,50))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,51))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,52))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,53))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQB,AMPL(1,54))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,55))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQT,AMPL(1,56))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,57))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQG_1EPS,AMPL(2,58))
      CALL MP_FFV1_0(W(1,5),W(1,7),W(1,6),R2_GQQ,AMPL(1,59))
      CALL MP_FFV1_1(W(1,4),W(1,6),GC_5,MDL_MT,MDL_WT,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,23),R2_QQQ,R2_QQT,AMPL(1,60)
     $ )
      CALL MP_R2_QQ_2_0(W(1,8),W(1,23),UV_TMASS,AMPL(1,61))
      CALL MP_R2_QQ_2_0(W(1,8),W(1,23),UV_TMASS_1EPS,AMPL(2,62))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL MP_R2_GG_1_R2_GG_3_0(W(1,6),W(1,14),R2_GGQ,R2_GGT,AMPL(1,63)
     $ )
C     Counter-term amplitude(s) for loop diagram number 47
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,64))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,65))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,66))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,67))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQB,AMPL(1,68))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,69))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQT,AMPL(1,70))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,71))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQG_1EPS,AMPL(2,72))
      CALL MP_FFV1_0(W(1,8),W(1,4),W(1,6),R2_GQQ,AMPL(1,73))
C     Counter-term amplitude(s) for loop diagram number 49
      CALL MP_FFS1_5_0(W(1,5),W(1,23),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,74))
      CALL MP_FFS1_5_0(W(1,5),W(1,23),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,75))
      CALL MP_FFS1_5_0(W(1,5),W(1,23),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,76))
C     Counter-term amplitude(s) for loop diagram number 51
      CALL MP_VVS1_0(W(1,6),W(1,15),W(1,3),R2_GGX0_HT,AMPL(1,77))
C     Counter-term amplitude(s) for loop diagram number 55
      CALL MP_FFS1_5_0(W(1,22),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,78))
      CALL MP_FFS1_5_0(W(1,22),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,79))
      CALL MP_FFS1_5_0(W(1,22),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,80))
      CALL MP_FFS1_5_1(W(1,9),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,24))
C     Counter-term amplitude(s) for loop diagram number 56
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,24),R2_QQQ,R2_QQT,AMPL(1
     $ ,81))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,24),UV_TMASS,AMPL(1,82))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,24),UV_TMASS_1EPS,AMPL(2,83))
      CALL MP_FFS1_5_2(W(1,10),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,25))
C     Counter-term amplitude(s) for loop diagram number 57
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,25),W(1,9),R2_QQQ,R2_QQT,AMPL(1,84)
     $ )
      CALL MP_R2_QQ_2_0(W(1,25),W(1,9),UV_TMASS,AMPL(1,85))
      CALL MP_R2_QQ_2_0(W(1,25),W(1,9),UV_TMASS_1EPS,AMPL(2,86))
C     Counter-term amplitude(s) for loop diagram number 58
      CALL MP_FFS1_5_0(W(1,10),W(1,9),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,87))
      CALL MP_FFS1_5_0(W(1,10),W(1,9),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,88))
      CALL MP_FFS1_5_0(W(1,10),W(1,9),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,89))
C     Counter-term amplitude(s) for loop diagram number 59
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,90))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,91))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,92))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,93))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQB,AMPL(1,94))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,95))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQT,AMPL(1,96))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,97))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQG_1EPS,AMPL(2,98))
      CALL MP_FFV1_0(W(1,5),W(1,24),W(1,2),R2_GQQ,AMPL(1,99))
C     Counter-term amplitude(s) for loop diagram number 60
      CALL MP_VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGX0_HT,AMPL(1,100))
C     Counter-term amplitude(s) for loop diagram number 62
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,101))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,102))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,103))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,104))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQB,AMPL(1,105))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,106))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQT,AMPL(1,107))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,108))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQG_1EPS,AMPL(2,109))
      CALL MP_FFV1_0(W(1,8),W(1,9),W(1,2),R2_GQQ,AMPL(1,110))
      CALL MP_FFV1_1(W(1,9),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,26))
C     Counter-term amplitude(s) for loop diagram number 65
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,26),R2_QQQ,R2_QQT,AMPL(1
     $ ,111))
      CALL MP_R2_QQ_2_0(W(1,8),W(1,26),UV_TMASS,AMPL(1,112))
      CALL MP_R2_QQ_2_0(W(1,8),W(1,26),UV_TMASS_1EPS,AMPL(2,113))
      CALL MP_FFV1_2(W(1,8),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,27))
C     Counter-term amplitude(s) for loop diagram number 66
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,27),W(1,9),R2_QQQ,R2_QQT,AMPL(1
     $ ,114))
      CALL MP_R2_QQ_2_0(W(1,27),W(1,9),UV_TMASS,AMPL(1,115))
      CALL MP_R2_QQ_2_0(W(1,27),W(1,9),UV_TMASS_1EPS,AMPL(2,116))
C     Counter-term amplitude(s) for loop diagram number 68
      CALL MP_FFS1_5_0(W(1,5),W(1,26),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,117))
      CALL MP_FFS1_5_0(W(1,5),W(1,26),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,118))
      CALL MP_FFS1_5_0(W(1,5),W(1,26),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,119))
      CALL MP_FFS1_5_2(W(1,11),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,28))
C     Counter-term amplitude(s) for loop diagram number 71
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,28),W(1,12),R2_QQQ,R2_QQT,AMPL(1
     $ ,120))
      CALL MP_R2_QQ_2_0(W(1,28),W(1,12),UV_TMASS,AMPL(1,121))
      CALL MP_R2_QQ_2_0(W(1,28),W(1,12),UV_TMASS_1EPS,AMPL(2,122))
      CALL MP_FFS1_5_1(W(1,12),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT
     $ ,W(1,29))
C     Counter-term amplitude(s) for loop diagram number 72
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,11),W(1,29),R2_QQQ,R2_QQT,AMPL(1
     $ ,123))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,29),UV_TMASS,AMPL(1,124))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,29),UV_TMASS_1EPS,AMPL(2,125))
C     Counter-term amplitude(s) for loop diagram number 73
      CALL MP_FFS1_5_0(W(1,11),W(1,12),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,126))
      CALL MP_FFS1_5_0(W(1,11),W(1,12),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,127))
      CALL MP_FFS1_5_0(W(1,11),W(1,12),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,128))
C     Counter-term amplitude(s) for loop diagram number 74
      CALL MP_VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGX0_HT,AMPL(1,129))
C     Counter-term amplitude(s) for loop diagram number 75
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,130))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,131))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,132))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,133))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQB,AMPL(1,134))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,135))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQT,AMPL(1,136))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,137))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQG_1EPS,AMPL(2,138))
      CALL MP_FFV1_0(W(1,11),W(1,7),W(1,2),R2_GQQ,AMPL(1,139))
C     Counter-term amplitude(s) for loop diagram number 78
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,140))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,141))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,142))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,143))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQB,AMPL(1,144))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,145))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQT,AMPL(1,146))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,147))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQG_1EPS,AMPL(2,148))
      CALL MP_FFV1_0(W(1,28),W(1,4),W(1,2),R2_GQQ,AMPL(1,149))
      CALL MP_FFV1_2(W(1,11),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,30))
C     Counter-term amplitude(s) for loop diagram number 80
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,30),W(1,7),R2_QQQ,R2_QQT,AMPL(1
     $ ,150))
      CALL MP_R2_QQ_2_0(W(1,30),W(1,7),UV_TMASS,AMPL(1,151))
      CALL MP_R2_QQ_2_0(W(1,30),W(1,7),UV_TMASS_1EPS,AMPL(2,152))
      CALL MP_FFV1_1(W(1,7),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,31))
C     Counter-term amplitude(s) for loop diagram number 81
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,11),W(1,31),R2_QQQ,R2_QQT,AMPL(1
     $ ,153))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,31),UV_TMASS,AMPL(1,154))
      CALL MP_R2_QQ_2_0(W(1,11),W(1,31),UV_TMASS_1EPS,AMPL(2,155))
C     Counter-term amplitude(s) for loop diagram number 83
      CALL MP_FFS1_5_0(W(1,30),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,156))
      CALL MP_FFS1_5_0(W(1,30),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,157))
      CALL MP_FFS1_5_0(W(1,30),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,158))
C     Counter-term amplitude(s) for loop diagram number 86
      CALL MP_VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGX0_HT,AMPL(1,159))
C     Counter-term amplitude(s) for loop diagram number 88
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,160))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,161))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,162))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,163))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQB,AMPL(1,164))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,165))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQT,AMPL(1,166))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,167))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQG_1EPS,AMPL(2,168))
      CALL MP_FFV1_0(W(1,5),W(1,29),W(1,1),R2_GQQ,AMPL(1,169))
C     Counter-term amplitude(s) for loop diagram number 89
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,170))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,171))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,172))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,173))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQB,AMPL(1,174))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,175))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQT,AMPL(1,176))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,177))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQG_1EPS,AMPL(2,178))
      CALL MP_FFV1_0(W(1,8),W(1,12),W(1,1),R2_GQQ,AMPL(1,179))
C     Counter-term amplitude(s) for loop diagram number 92
      CALL MP_VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGX0_HT,AMPL(1,180))
C     Counter-term amplitude(s) for loop diagram number 93
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,181))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,182))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,183))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,184))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQB,AMPL(1,185))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,186))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQT,AMPL(1,187))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,188))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQG_1EPS,AMPL(2,189))
      CALL MP_FFV1_0(W(1,10),W(1,7),W(1,1),R2_GQQ,AMPL(1,190))
C     Counter-term amplitude(s) for loop diagram number 97
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,191))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,192))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,193))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,194))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQB,AMPL(1,195))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,196))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQT,AMPL(1,197))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,198))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,199))
      CALL MP_FFV1_0(W(1,25),W(1,4),W(1,1),R2_GQQ,AMPL(1,200))
C     Counter-term amplitude(s) for loop diagram number 100
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,201))
C     Counter-term amplitude(s) for loop diagram number 102
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,202))
C     Counter-term amplitude(s) for loop diagram number 108
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,203))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,204))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,205))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,206))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQB,AMPL(1,207))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,208))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQT,AMPL(1,209))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,210))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQG_1EPS,AMPL(2,211))
      CALL MP_FFV1_0(W(1,5),W(1,31),W(1,1),R2_GQQ,AMPL(1,212))
C     Counter-term amplitude(s) for loop diagram number 113
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,213))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,214))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,215))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,216))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQB,AMPL(1,217))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,218))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQT,AMPL(1,219))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,220))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,221))
      CALL MP_FFV1_0(W(1,27),W(1,4),W(1,1),R2_GQQ,AMPL(1,222))
C     Counter-term amplitude(s) for loop diagram number 118
      CALL MP_VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGX0_HT,AMPL(1,223))
      CALL MP_FFV1_1(W(1,12),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,32))
C     Counter-term amplitude(s) for loop diagram number 126
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,32),R2_QQQ,R2_QQT,AMPL(1
     $ ,224))
      CALL MP_R2_QQ_2_0(W(1,8),W(1,32),UV_TMASS,AMPL(1,225))
      CALL MP_R2_QQ_2_0(W(1,8),W(1,32),UV_TMASS_1EPS,AMPL(2,226))
      CALL MP_FFV1_2(W(1,8),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,33))
C     Counter-term amplitude(s) for loop diagram number 127
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,33),W(1,12),R2_QQQ,R2_QQT,AMPL(1
     $ ,227))
      CALL MP_R2_QQ_2_0(W(1,33),W(1,12),UV_TMASS,AMPL(1,228))
      CALL MP_R2_QQ_2_0(W(1,33),W(1,12),UV_TMASS_1EPS,AMPL(2,229))
C     Counter-term amplitude(s) for loop diagram number 129
      CALL MP_FFS1_5_0(W(1,5),W(1,32),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,230))
      CALL MP_FFS1_5_0(W(1,5),W(1,32),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,231))
      CALL MP_FFS1_5_0(W(1,5),W(1,32),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,232))
      CALL MP_FFV1_2(W(1,10),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,34))
C     Counter-term amplitude(s) for loop diagram number 132
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,34),W(1,7),R2_QQQ,R2_QQT,AMPL(1
     $ ,233))
      CALL MP_R2_QQ_2_0(W(1,34),W(1,7),UV_TMASS,AMPL(1,234))
      CALL MP_R2_QQ_2_0(W(1,34),W(1,7),UV_TMASS_1EPS,AMPL(2,235))
      CALL MP_FFV1_1(W(1,7),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,35))
C     Counter-term amplitude(s) for loop diagram number 133
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,35),R2_QQQ,R2_QQT,AMPL(1
     $ ,236))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,35),UV_TMASS,AMPL(1,237))
      CALL MP_R2_QQ_2_0(W(1,10),W(1,35),UV_TMASS_1EPS,AMPL(2,238))
C     Counter-term amplitude(s) for loop diagram number 135
      CALL MP_FFS1_5_0(W(1,34),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A
     $ ,AMPL(1,239))
      CALL MP_FFS1_5_0(W(1,34),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,240))
      CALL MP_FFS1_5_0(W(1,34),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A
     $ ,AMPL(1,241))
C     Counter-term amplitude(s) for loop diagram number 138
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,242))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,243))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,244))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,245))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQB,AMPL(1,246))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,247))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQT,AMPL(1,248))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,249))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQG_1EPS,AMPL(2,250))
      CALL MP_FFV1_0(W(1,5),W(1,35),W(1,2),R2_GQQ,AMPL(1,251))
C     Counter-term amplitude(s) for loop diagram number 140
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,252))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,253))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,254))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,255))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQB,AMPL(1,256))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,257))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQT,AMPL(1,258))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,259))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQG_1EPS,AMPL(2,260))
      CALL MP_FFV1_0(W(1,33),W(1,4),W(1,2),R2_GQQ,AMPL(1,261))
C     Counter-term amplitude(s) for loop diagram number 143
      CALL MP_VVS1_0(W(1,2),W(1,21),W(1,3),R2_GGX0_HT,AMPL(1,262))
C     Counter-term amplitude(s) for loop diagram number 159
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,6),W(1,13),R2_GGG_1,R2_GGG_2
     $ ,AMPL(1,263))
C     Counter-term amplitude(s) for loop diagram number 160
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,6),W(1,14),R2_GGG_1,R2_GGG_2
     $ ,AMPL(1,264))
C     Counter-term amplitude(s) for loop diagram number 161
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GG,AMPL(1,265))
C     Counter-term amplitude(s) for loop diagram number 162
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GG,AMPL(1,266))
C     At this point, all CT amps needed for (QCD=6), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

