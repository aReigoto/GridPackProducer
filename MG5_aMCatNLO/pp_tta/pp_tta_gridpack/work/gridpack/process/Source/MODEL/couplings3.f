ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP3()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      R2_GGG_1 = (2.000000D+00)*MDL_COMPLEXI*MDL_G__EXP__2*MDL_NCOL
     $ /(4.800000D+01*PI**2)*(1.000000D+00/2.000000D+00+MDL_LHV)
      R2_GGG_2 = -(2.000000D+00)*MDL_COMPLEXI*MDL_G__EXP__2*MDL_NCOL
     $ /(4.800000D+01*PI**2)*MDL_LHV
      R2_QQQ = MDL_LHV*MDL_COMPLEXI*MDL_G__EXP__2*(MDL_NCOL__EXP__2
     $ -1.000000D+00)/(3.200000D+01*PI**2*MDL_NCOL)
      R2_QQT = MDL_LHV*MDL_COMPLEXI*MDL_G__EXP__2*(MDL_NCOL__EXP__2
     $ -1.000000D+00)*(2.000000D+00*MDL_MT)/(3.200000D+01*PI**2
     $ *MDL_NCOL)
      UV_3GG_1EPS = -MDL_G_UVG_1EPS_*G
      UV_3GB_1EPS = -MDL_G_UVB_1EPS_*G
      UV_GQQG_1EPS = MDL_COMPLEXI*MDL_G_UVG_1EPS_*G
      UV_GQQQ_1EPS = MDL_COMPLEXI*MDL_G_UVQ_1EPS_*G
      UV_TMASS_1EPS = MDL_TMASS_UV_1EPS_
      UVWFCT_B_0_1EPS = COND(DCMPLX(MDL_MB),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))
     $ *3.000000D+00*MDL_CF))
      UVWFCT_G_2_1EPS = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $ *4.000000D+00*MDL_TF))
      END
