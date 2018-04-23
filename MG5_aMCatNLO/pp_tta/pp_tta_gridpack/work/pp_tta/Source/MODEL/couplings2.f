ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      UV_3GB = -MDL_G_UVB_FIN_*G
      UV_3GT = -MDL_G_UVT_FIN_*G
      UV_GQQB = MDL_COMPLEXI*MDL_G_UVB_FIN_*G
      UV_GQQT = MDL_COMPLEXI*MDL_G_UVT_FIN_*G
      UV_TMASS = MDL_TMASS_UV_FIN_
      UVWFCT_T_0 = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00),DCMPLX(
     $ -((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))*MDL_CF
     $ *(4.000000D+00-3.000000D+00*REGLOG(DCMPLX(MDL_MT__EXP__2
     $ /MDL_MU_R__EXP__2)))))
      UVWFCT_G_2 = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $ *4.000000D+00*MDL_TF*REGLOG(DCMPLX(MDL_MT__EXP__2
     $ /MDL_MU_R__EXP__2))))
      UVWFCT_G_1 = COND(DCMPLX(MDL_MB),DCMPLX(0.000000D+00)
     $ ,DCMPLX(((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $ *4.000000D+00*MDL_TF*REGLOG(DCMPLX(MDL_MB__EXP__2
     $ /MDL_MU_R__EXP__2))))
      UV_X0TT_H = -((MDL_COSA*MDL_KHTT*MDL_COMPLEXI*MDL_YT)
     $ /MDL_SQRT__2)*MDL_UV_YUK_T_FIN_
      UV_X0TT_A = MDL_SINA*MDL_KATT*MDL_YT/MDL_SQRT__2
     $ *MDL_UV_YUK_T_FIN_
      R2_TTX0_H = (-((MDL_COSA*MDL_KHTT*MDL_COMPLEXI*MDL_YT)
     $ /MDL_SQRT__2))*(2.000000D+00*MDL_R2MIXEDFACTOR_FIN_)
      R2_TTX0_A = (MDL_SINA*MDL_KATT*MDL_YT/MDL_SQRT__2)*(2.000000D+00
     $ *MDL_R2MIXEDFACTOR_FIN_)
      R2_GGX0_HB = (-((MDL_COSA*MDL_KHBB*MDL_COMPLEXI*MDL_YB)
     $ /MDL_SQRT__2))*(MDL_G__EXP__2/(8.000000D+00*PI**2))*MDL_MB
      R2_GGX0_HT = (-((MDL_COSA*MDL_KHTT*MDL_COMPLEXI*MDL_YT)
     $ /MDL_SQRT__2))*(MDL_G__EXP__2/(8.000000D+00*PI**2))*MDL_MT
      UV_X0TT_H_1EPS = -((MDL_COSA*MDL_KHTT*MDL_COMPLEXI*MDL_YT)
     $ /MDL_SQRT__2)*MDL_UV_YUK_T_1EPS_
      UV_X0TT_A_1EPS = MDL_SINA*MDL_KATT*MDL_YT/MDL_SQRT__2
     $ *MDL_UV_YUK_T_1EPS_
      GC_4 = -G
      GC_5 = MDL_COMPLEXI*G
      GC_6 = MDL_COMPLEXI*MDL_G__EXP__2
      R2_3GQ = 2.000000D+00*MDL_G__EXP__3/(4.800000D+01*PI**2)
      R2_3GG = MDL_NCOL*MDL_G__EXP__3/(4.800000D+01*PI**2)*(7.000000D
     $ +00/4.000000D+00+MDL_LHV)
      R2_GQQ = -MDL_COMPLEXI*MDL_G__EXP__3/(1.600000D+01*PI**2)
     $ *((MDL_NCOL__EXP__2-1.000000D+00)/(2.000000D+00*MDL_NCOL))
     $ *(1.000000D+00+MDL_LHV)
      R2_GGQ = (2.000000D+00)*MDL_COMPLEXI*MDL_G__EXP__2/(4.800000D+01
     $ *PI**2)
      R2_GGB = (2.000000D+00)*MDL_COMPLEXI*MDL_G__EXP__2*(-6.000000D
     $ +00*MDL_MB__EXP__2)/(4.800000D+01*PI**2)
      R2_GGT = (2.000000D+00)*MDL_COMPLEXI*MDL_G__EXP__2*(-6.000000D
     $ +00*MDL_MT__EXP__2)/(4.800000D+01*PI**2)
      END
