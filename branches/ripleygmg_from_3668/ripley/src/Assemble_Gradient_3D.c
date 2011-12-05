Ripley_Assemble_Gradient_3D(Ripley_Grid *grid, escriptDataC *in, escriptDataC *out)
{
   /* GENERATOR SNIP TOP */
   if (out_data_type==RIPLEY_ELEMENTS) {
      const double tmp0_22 = -0.044658198738520451079/h1;
      const double tmp0_16 = 0.16666666666666666667/h0;
      const double tmp0_33 = -0.62200846792814621559/h1;
      const double tmp0_0 = -0.16666666666666666667/h0;
      const double tmp0_21 = -0.62200846792814621559/h1;
      const double tmp0_17 = -0.62200846792814621559/h0;
      const double tmp0_52 = 0.16666666666666666667/h2;
      const double tmp0_1 = -0.62200846792814621559/h0;
      const double tmp0_20 = 0.62200846792814621559/h1;
      const double tmp0_14 = -0.16666666666666666667/h0;
      const double tmp0_53 = -0.62200846792814621559/h2;
      const double tmp0_49 = 0.16666666666666666667/h2;
      const double tmp0_2 = 0.16666666666666666667/h0;
      const double tmp0_27 = -0.16666666666666666667/h1;
      const double tmp0_15 = -0.044658198738520451079/h0;
      const double tmp0_50 = -0.16666666666666666667/h2;
      const double tmp0_48 = 0.62200846792814621559/h2;
      const double tmp0_3 = -0.044658198738520451079/h0;
      const double tmp0_26 = 0.16666666666666666667/h1;
      const double tmp0_12 = 0.16666666666666666667/h0;
      const double tmp0_51 = 0.044658198738520451079/h2;
      const double tmp0_25 = 0.044658198738520451079/h1;
      const double tmp0_13 = 0.16666666666666666667/h0;
      const double tmp0_56 = 0.16666666666666666667/h2;
      const double tmp0_24 = 0.16666666666666666667/h1;
      const double tmp0_10 = 0.62200846792814621559/h0;
      const double tmp0_57 = 0.044658198738520451079/h2;
      const double tmp0_11 = -0.16666666666666666667/h0;
      const double tmp0_54 = -0.16666666666666666667/h2;
      const double tmp0_38 = 0.16666666666666666667/h1;
      const double tmp0_34 = 0.044658198738520451079/h1;
      const double tmp0_42 = 0.16666666666666666667/h2;
      const double tmp0_35 = -0.044658198738520451079/h1;
      const double tmp0_36 = -0.62200846792814621559/h1;
      const double tmp0_41 = -0.62200846792814621559/h2;
      const double tmp0_8 = 0.044658198738520451079/h0;
      const double tmp0_37 = -0.16666666666666666667/h1;
      const double tmp0_29 = -0.044658198738520451079/h1;
      const double tmp0_40 = -0.16666666666666666667/h2;
      const double tmp0_9 = -0.044658198738520451079/h0;
      const double tmp0_30 = 0.62200846792814621559/h1;
      const double tmp0_28 = -0.16666666666666666667/h1;
      const double tmp0_43 = 0.62200846792814621559/h2;
      const double tmp0_32 = 0.16666666666666666667/h1;
      const double tmp0_31 = 0.044658198738520451079/h1;
      const double tmp0_39 = 0.62200846792814621559/h1;
      const double tmp0_58 = -0.62200846792814621559/h2;
      const double tmp0_55 = -0.044658198738520451079/h2;
      const double tmp0_18 = 0.62200846792814621559/h0;
      const double tmp0_45 = 0.044658198738520451079/h2;
      const double tmp0_59 = 0.62200846792814621559/h2;
      const double tmp0_4 = 0.044658198738520451079/h0;
      const double tmp0_19 = 0.044658198738520451079/h0;
      const double tmp0_44 = -0.044658198738520451079/h2;
      const double tmp0_5 = 0.62200846792814621559/h0;
      const double tmp0_47 = -0.16666666666666666667/h2;
      const double tmp0_6 = -0.62200846792814621559/h0;
      const double tmp0_23 = -0.16666666666666666667/h1;
      const double tmp0_46 = -0.044658198738520451079/h2;
      const double tmp0_7 = -0.16666666666666666667/h0;
      #pragma omp parallel for private(i,k2,k1,k0)
      for (k2 =0; k2 < N2; ++k2) {
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,k2, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_1 + f_011*tmp0_3 + f_100*tmp0_5 + f_111*tmp0_4 + tmp0_0*(f_001 + f_010) + tmp0_2*(f_101 + f_110);
                  out[INDEX4(i,1,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_21 + f_010*tmp0_20 + f_101*tmp0_22 + f_111*tmp0_25 + tmp0_23*(f_001 + f_100) + tmp0_24*(f_011 + f_110);
                  out[INDEX4(i,2,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_41 + f_001*tmp0_43 + f_110*tmp0_44 + f_111*tmp0_45 + tmp0_40*(f_010 + f_100) + tmp0_42*(f_011 + f_101);
                  out[INDEX4(i,0,1,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_1 + f_011*tmp0_3 + f_100*tmp0_5 + f_111*tmp0_4 + tmp0_0*(f_001 + f_010) + tmp0_2*(f_101 + f_110);
                  out[INDEX4(i,1,1,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_27 + f_001*tmp0_29 + f_010*tmp0_26 + f_011*tmp0_31 + f_100*tmp0_33 + f_101*tmp0_28 + f_110*tmp0_30 + f_111*tmp0_32;
                  out[INDEX4(i,2,1,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_47 + f_001*tmp0_49 + f_010*tmp0_46 + f_011*tmp0_51 + f_100*tmp0_53 + f_101*tmp0_48 + f_110*tmp0_50 + f_111*tmp0_52;
                  out[INDEX4(i,0,2,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_7 + f_001*tmp0_9 + f_010*tmp0_6 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_8 + f_110*tmp0_10 + f_111*tmp0_12;
                  out[INDEX4(i,1,2,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_21 + f_010*tmp0_20 + f_101*tmp0_22 + f_111*tmp0_25 + tmp0_23*(f_001 + f_100) + tmp0_24*(f_011 + f_110);
                  out[INDEX4(i,2,2,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_47 + f_001*tmp0_49 + f_010*tmp0_53 + f_011*tmp0_48 + f_100*tmp0_46 + f_101*tmp0_51 + f_110*tmp0_50 + f_111*tmp0_52;
                  out[INDEX4(i,0,3,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_7 + f_001*tmp0_9 + f_010*tmp0_6 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_8 + f_110*tmp0_10 + f_111*tmp0_12;
                  out[INDEX4(i,1,3,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_27 + f_001*tmp0_29 + f_010*tmp0_26 + f_011*tmp0_31 + f_100*tmp0_33 + f_101*tmp0_28 + f_110*tmp0_30 + f_111*tmp0_32;
                  out[INDEX4(i,2,3,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_55 + f_001*tmp0_57 + f_110*tmp0_58 + f_111*tmp0_59 + tmp0_54*(f_010 + f_100) + tmp0_56*(f_011 + f_101);
                  out[INDEX4(i,0,4,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_7 + f_001*tmp0_6 + f_010*tmp0_9 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_10 + f_110*tmp0_8 + f_111*tmp0_12;
                  out[INDEX4(i,1,4,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_27 + f_001*tmp0_33 + f_010*tmp0_26 + f_011*tmp0_30 + f_100*tmp0_29 + f_101*tmp0_28 + f_110*tmp0_31 + f_111*tmp0_32;
                  out[INDEX4(i,2,4,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_41 + f_001*tmp0_43 + f_110*tmp0_44 + f_111*tmp0_45 + tmp0_40*(f_010 + f_100) + tmp0_42*(f_011 + f_101);
                  out[INDEX4(i,0,5,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_7 + f_001*tmp0_6 + f_010*tmp0_9 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_10 + f_110*tmp0_8 + f_111*tmp0_12;
                  out[INDEX4(i,1,5,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_35 + f_010*tmp0_34 + f_101*tmp0_36 + f_111*tmp0_39 + tmp0_37*(f_001 + f_100) + tmp0_38*(f_011 + f_110);
                  out[INDEX4(i,2,5,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_47 + f_001*tmp0_49 + f_010*tmp0_46 + f_011*tmp0_51 + f_100*tmp0_53 + f_101*tmp0_48 + f_110*tmp0_50 + f_111*tmp0_52;
                  out[INDEX4(i,0,6,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_15 + f_011*tmp0_17 + f_100*tmp0_19 + f_111*tmp0_18 + tmp0_14*(f_001 + f_010) + tmp0_16*(f_101 + f_110);
                  out[INDEX4(i,1,6,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_27 + f_001*tmp0_33 + f_010*tmp0_26 + f_011*tmp0_30 + f_100*tmp0_29 + f_101*tmp0_28 + f_110*tmp0_31 + f_111*tmp0_32;
                  out[INDEX4(i,2,6,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_47 + f_001*tmp0_49 + f_010*tmp0_53 + f_011*tmp0_48 + f_100*tmp0_46 + f_101*tmp0_51 + f_110*tmp0_50 + f_111*tmp0_52;
                  out[INDEX4(i,0,7,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_15 + f_011*tmp0_17 + f_100*tmp0_19 + f_111*tmp0_18 + tmp0_14*(f_001 + f_010) + tmp0_16*(f_101 + f_110);
                  out[INDEX4(i,1,7,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_35 + f_010*tmp0_34 + f_101*tmp0_36 + f_111*tmp0_39 + tmp0_37*(f_001 + f_100) + tmp0_38*(f_011 + f_110);
                  out[INDEX4(i,2,7,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,8)] = f_000*tmp0_55 + f_001*tmp0_57 + f_110*tmp0_58 + f_111*tmp0_59 + tmp0_54*(f_010 + f_100) + tmp0_56*(f_011 + f_101);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* close k2 loop */
   } else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {
      const double tmp0_0 = -0.25/h0;
      const double tmp0_4 = -0.25/h2;
      const double tmp0_1 = 0.25/h0;
      const double tmp0_5 = 0.25/h2;
      const double tmp0_2 = -0.25/h1;
      const double tmp0_3 = 0.25/h1;
      #pragma omp parallel for private(i,k2,k1,k0)
      for (k2 =0; k2 < N2; ++k2) {
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,k2, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_000 + f_001 + f_010 + f_011) + tmp0_1*(f_100 + f_101 + f_110 + f_111);
                  out[INDEX4(i,1,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_000 + f_001 + f_100 + f_101) + tmp0_3*(f_010 + f_011 + f_110 + f_111);
                  out[INDEX4(i,2,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_000 + f_010 + f_100 + f_110) + tmp0_5*(f_001 + f_011 + f_101 + f_111);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* close k2 loop */
   } else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_22 = -0.21132486540518711775/h1;
         const double tmp0_16 = 0.16666666666666666667/h0;
         const double tmp0_33 = -0.21132486540518711775/h2;
         const double tmp0_0 = -0.16666666666666666667/h0;
         const double tmp0_21 = -0.78867513459481288225/h1;
         const double tmp0_17 = -0.62200846792814621559/h0;
         const double tmp0_1 = -0.62200846792814621559/h0;
         const double tmp0_20 = 0.78867513459481288225/h1;
         const double tmp0_14 = -0.16666666666666666667/h0;
         const double tmp0_2 = 0.16666666666666666667/h0;
         const double tmp0_27 = 0.78867513459481288225/h1;
         const double tmp0_15 = -0.044658198738520451079/h0;
         const double tmp0_3 = -0.044658198738520451079/h0;
         const double tmp0_26 = -0.78867513459481288225/h1;
         const double tmp0_12 = 0.16666666666666666667/h0;
         const double tmp0_25 = -0.21132486540518711775/h1;
         const double tmp0_13 = 0.16666666666666666667/h0;
         const double tmp0_24 = 0.21132486540518711775/h1;
         const double tmp0_10 = 0.62200846792814621559/h0;
         const double tmp0_11 = -0.16666666666666666667/h0;
         const double tmp0_34 = 0.21132486540518711775/h2;
         const double tmp0_35 = 0.78867513459481288225/h2;
         const double tmp0_8 = 0.044658198738520451079/h0;
         const double tmp0_29 = -0.78867513459481288225/h2;
         const double tmp0_9 = -0.044658198738520451079/h0;
         const double tmp0_30 = 0.78867513459481288225/h2;
         const double tmp0_28 = -0.21132486540518711775/h2;
         const double tmp0_32 = -0.78867513459481288225/h2;
         const double tmp0_31 = 0.21132486540518711775/h2;
         const double tmp0_18 = 0.62200846792814621559/h0;
         const double tmp0_4 = 0.044658198738520451079/h0;
         const double tmp0_19 = 0.044658198738520451079/h0;
         const double tmp0_5 = 0.62200846792814621559/h0;
         const double tmp0_6 = -0.62200846792814621559/h0;
         const double tmp0_23 = 0.21132486540518711775/h1;
         const double tmp0_7 = -0.16666666666666666667/h0;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_010 = in[INDEX2(i,INDEX3(0,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(0,k1,k2, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(0,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(0,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(1,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(1,k1,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_1 + f_011*tmp0_3 + f_100*tmp0_5 + f_111*tmp0_4 + tmp0_0*(f_001 + f_010) + tmp0_2*(f_101 + f_110);
                  out[INDEX4(i,1,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_21 + f_001*tmp0_22 + f_010*tmp0_20 + f_011*tmp0_23;
                  out[INDEX4(i,2,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_29 + f_001*tmp0_30 + f_010*tmp0_28 + f_011*tmp0_31;
                  out[INDEX4(i,0,1,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_7 + f_001*tmp0_9 + f_010*tmp0_6 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_8 + f_110*tmp0_10 + f_111*tmp0_12;
                  out[INDEX4(i,1,1,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_21 + f_001*tmp0_22 + f_010*tmp0_20 + f_011*tmp0_23;
                  out[INDEX4(i,2,1,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_33 + f_001*tmp0_34 + f_010*tmp0_32 + f_011*tmp0_35;
                  out[INDEX4(i,0,2,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_7 + f_001*tmp0_6 + f_010*tmp0_9 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_10 + f_110*tmp0_8 + f_111*tmp0_12;
                  out[INDEX4(i,1,2,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_25 + f_001*tmp0_26 + f_010*tmp0_24 + f_011*tmp0_27;
                  out[INDEX4(i,2,2,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_29 + f_001*tmp0_30 + f_010*tmp0_28 + f_011*tmp0_31;
                  out[INDEX4(i,0,3,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_15 + f_011*tmp0_17 + f_100*tmp0_19 + f_111*tmp0_18 + tmp0_14*(f_001 + f_010) + tmp0_16*(f_101 + f_110);
                  out[INDEX4(i,1,3,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_25 + f_001*tmp0_26 + f_010*tmp0_24 + f_011*tmp0_27;
                  out[INDEX4(i,2,3,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_33 + f_001*tmp0_34 + f_010*tmp0_32 + f_011*tmp0_35;
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_22 = 0.21132486540518711775/h1;
         const double tmp0_16 = 0.16666666666666666667/h0;
         const double tmp0_33 = -0.78867513459481288225/h2;
         const double tmp0_0 = -0.16666666666666666667/h0;
         const double tmp0_21 = 0.78867513459481288225/h1;
         const double tmp0_17 = -0.62200846792814621559/h0;
         const double tmp0_1 = -0.62200846792814621559/h0;
         const double tmp0_20 = -0.21132486540518711775/h1;
         const double tmp0_14 = -0.16666666666666666667/h0;
         const double tmp0_2 = 0.16666666666666666667/h0;
         const double tmp0_27 = -0.21132486540518711775/h1;
         const double tmp0_15 = -0.044658198738520451079/h0;
         const double tmp0_3 = -0.044658198738520451079/h0;
         const double tmp0_26 = 0.78867513459481288225/h1;
         const double tmp0_12 = 0.16666666666666666667/h0;
         const double tmp0_25 = 0.21132486540518711775/h1;
         const double tmp0_13 = 0.16666666666666666667/h0;
         const double tmp0_24 = -0.78867513459481288225/h1;
         const double tmp0_10 = 0.62200846792814621559/h0;
         const double tmp0_11 = -0.16666666666666666667/h0;
         const double tmp0_34 = 0.78867513459481288225/h2;
         const double tmp0_35 = -0.21132486540518711775/h2;
         const double tmp0_8 = 0.044658198738520451079/h0;
         const double tmp0_29 = -0.21132486540518711775/h2;
         const double tmp0_9 = -0.044658198738520451079/h0;
         const double tmp0_30 = 0.21132486540518711775/h2;
         const double tmp0_28 = 0.78867513459481288225/h2;
         const double tmp0_32 = 0.21132486540518711775/h2;
         const double tmp0_31 = -0.78867513459481288225/h2;
         const double tmp0_18 = 0.62200846792814621559/h0;
         const double tmp0_4 = 0.044658198738520451079/h0;
         const double tmp0_19 = 0.044658198738520451079/h0;
         const double tmp0_5 = 0.62200846792814621559/h0;
         const double tmp0_6 = -0.62200846792814621559/h0;
         const double tmp0_23 = -0.78867513459481288225/h1;
         const double tmp0_7 = -0.16666666666666666667/h0;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_010 = in[INDEX2(i,INDEX3(M0-2,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(M0-2,k1,k2, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(M0-1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(M0-2,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(M0-2,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(M0-1,k1,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_1 + f_011*tmp0_3 + f_100*tmp0_5 + f_111*tmp0_4 + tmp0_0*(f_001 + f_010) + tmp0_2*(f_101 + f_110);
                  out[INDEX4(i,1,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_23 + f_101*tmp0_20 + f_110*tmp0_21 + f_111*tmp0_22;
                  out[INDEX4(i,2,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_31 + f_101*tmp0_28 + f_110*tmp0_29 + f_111*tmp0_30;
                  out[INDEX4(i,0,1,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_7 + f_001*tmp0_9 + f_010*tmp0_6 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_8 + f_110*tmp0_10 + f_111*tmp0_12;
                  out[INDEX4(i,1,1,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_23 + f_101*tmp0_20 + f_110*tmp0_21 + f_111*tmp0_22;
                  out[INDEX4(i,2,1,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_35 + f_101*tmp0_32 + f_110*tmp0_33 + f_111*tmp0_34;
                  out[INDEX4(i,0,2,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_7 + f_001*tmp0_6 + f_010*tmp0_9 + f_011*tmp0_11 + f_100*tmp0_13 + f_101*tmp0_10 + f_110*tmp0_8 + f_111*tmp0_12;
                  out[INDEX4(i,1,2,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_27 + f_101*tmp0_24 + f_110*tmp0_25 + f_111*tmp0_26;
                  out[INDEX4(i,2,2,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_31 + f_101*tmp0_28 + f_110*tmp0_29 + f_111*tmp0_30;
                  out[INDEX4(i,0,3,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_15 + f_011*tmp0_17 + f_100*tmp0_19 + f_111*tmp0_18 + tmp0_14*(f_001 + f_010) + tmp0_16*(f_101 + f_110);
                  out[INDEX4(i,1,3,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_27 + f_101*tmp0_24 + f_110*tmp0_25 + f_111*tmp0_26;
                  out[INDEX4(i,2,3,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_100*tmp0_35 + f_101*tmp0_32 + f_110*tmp0_33 + f_111*tmp0_34;
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_22 = -0.044658198738520451079/h1;
         const double tmp0_16 = -0.16666666666666666667/h1;
         const double tmp0_33 = 0.21132486540518711775/h2;
         const double tmp0_0 = -0.78867513459481288225/h0;
         const double tmp0_21 = 0.16666666666666666667/h1;
         const double tmp0_17 = -0.62200846792814621559/h1;
         const double tmp0_1 = -0.21132486540518711775/h0;
         const double tmp0_20 = 0.044658198738520451079/h1;
         const double tmp0_14 = -0.16666666666666666667/h1;
         const double tmp0_2 = 0.21132486540518711775/h0;
         const double tmp0_27 = 0.044658198738520451079/h1;
         const double tmp0_15 = -0.044658198738520451079/h1;
         const double tmp0_3 = 0.78867513459481288225/h0;
         const double tmp0_26 = 0.16666666666666666667/h1;
         const double tmp0_12 = 0.16666666666666666667/h1;
         const double tmp0_25 = 0.62200846792814621559/h1;
         const double tmp0_13 = 0.62200846792814621559/h1;
         const double tmp0_24 = -0.62200846792814621559/h1;
         const double tmp0_10 = -0.044658198738520451079/h1;
         const double tmp0_11 = 0.044658198738520451079/h1;
         const double tmp0_34 = 0.78867513459481288225/h2;
         const double tmp0_35 = -0.78867513459481288225/h2;
         const double tmp0_8 = -0.62200846792814621559/h1;
         const double tmp0_29 = 0.78867513459481288225/h2;
         const double tmp0_9 = -0.16666666666666666667/h1;
         const double tmp0_30 = 0.21132486540518711775/h2;
         const double tmp0_28 = -0.78867513459481288225/h2;
         const double tmp0_32 = -0.21132486540518711775/h2;
         const double tmp0_31 = -0.21132486540518711775/h2;
         const double tmp0_18 = 0.16666666666666666667/h1;
         const double tmp0_4 = -0.21132486540518711775/h0;
         const double tmp0_19 = 0.62200846792814621559/h1;
         const double tmp0_5 = -0.78867513459481288225/h0;
         const double tmp0_6 = 0.78867513459481288225/h0;
         const double tmp0_23 = -0.16666666666666666667/h1;
         const double tmp0_7 = 0.21132486540518711775/h0;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,0,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,0,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,0,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,0,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,1,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,1,k2+1, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,1,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_0 + f_001*tmp0_1 + f_100*tmp0_3 + f_101*tmp0_2;
                  out[INDEX4(i,1,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_8 + f_010*tmp0_13 + f_101*tmp0_10 + f_111*tmp0_11 + tmp0_12*(f_011 + f_110) + tmp0_9*(f_001 + f_100);
                  out[INDEX4(i,2,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_28 + f_001*tmp0_29 + f_100*tmp0_31 + f_101*tmp0_30;
                  out[INDEX4(i,0,1,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_0 + f_001*tmp0_1 + f_100*tmp0_3 + f_101*tmp0_2;
                  out[INDEX4(i,1,1,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_14 + f_001*tmp0_15 + f_010*tmp0_21 + f_011*tmp0_20 + f_100*tmp0_17 + f_101*tmp0_16 + f_110*tmp0_19 + f_111*tmp0_18;
                  out[INDEX4(i,2,1,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_32 + f_001*tmp0_33 + f_100*tmp0_35 + f_101*tmp0_34;
                  out[INDEX4(i,0,2,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_4 + f_001*tmp0_5 + f_100*tmp0_7 + f_101*tmp0_6;
                  out[INDEX4(i,1,2,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_14 + f_001*tmp0_17 + f_010*tmp0_21 + f_011*tmp0_19 + f_100*tmp0_15 + f_101*tmp0_16 + f_110*tmp0_20 + f_111*tmp0_18;
                  out[INDEX4(i,2,2,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_28 + f_001*tmp0_29 + f_100*tmp0_31 + f_101*tmp0_30;
                  out[INDEX4(i,0,3,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_4 + f_001*tmp0_5 + f_100*tmp0_7 + f_101*tmp0_6;
                  out[INDEX4(i,1,3,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_22 + f_010*tmp0_27 + f_101*tmp0_24 + f_111*tmp0_25 + tmp0_23*(f_001 + f_100) + tmp0_26*(f_011 + f_110);
                  out[INDEX4(i,2,3,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_32 + f_001*tmp0_33 + f_100*tmp0_35 + f_101*tmp0_34;
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_22 = 0.16666666666666666667/h1;
         const double tmp0_16 = 0.16666666666666666667/h1;
         const double tmp0_33 = 0.21132486540518711775/h2;
         const double tmp0_0 = 0.78867513459481288225/h0;
         const double tmp0_21 = -0.62200846792814621559/h1;
         const double tmp0_17 = 0.16666666666666666667/h1;
         const double tmp0_1 = -0.21132486540518711775/h0;
         const double tmp0_20 = -0.16666666666666666667/h1;
         const double tmp0_14 = 0.62200846792814621559/h1;
         const double tmp0_2 = -0.78867513459481288225/h0;
         const double tmp0_27 = -0.62200846792814621559/h1;
         const double tmp0_15 = 0.044658198738520451079/h1;
         const double tmp0_3 = 0.21132486540518711775/h0;
         const double tmp0_26 = -0.16666666666666666667/h1;
         const double tmp0_12 = -0.16666666666666666667/h1;
         const double tmp0_25 = -0.044658198738520451079/h1;
         const double tmp0_13 = -0.044658198738520451079/h1;
         const double tmp0_24 = 0.62200846792814621559/h1;
         const double tmp0_10 = 0.044658198738520451079/h1;
         const double tmp0_11 = -0.62200846792814621559/h1;
         const double tmp0_34 = -0.21132486540518711775/h2;
         const double tmp0_35 = 0.78867513459481288225/h2;
         const double tmp0_8 = 0.16666666666666666667/h1;
         const double tmp0_29 = 0.78867513459481288225/h2;
         const double tmp0_9 = 0.62200846792814621559/h1;
         const double tmp0_30 = -0.78867513459481288225/h2;
         const double tmp0_28 = -0.21132486540518711775/h2;
         const double tmp0_32 = -0.78867513459481288225/h2;
         const double tmp0_31 = 0.21132486540518711775/h2;
         const double tmp0_18 = -0.16666666666666666667/h1;
         const double tmp0_4 = 0.21132486540518711775/h0;
         const double tmp0_19 = -0.044658198738520451079/h1;
         const double tmp0_5 = -0.78867513459481288225/h0;
         const double tmp0_6 = -0.21132486540518711775/h0;
         const double tmp0_23 = 0.044658198738520451079/h1;
         const double tmp0_7 = 0.78867513459481288225/h0;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,M1-1,k2+1, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,M1-1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2+1, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,M1-2,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,M1-2,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,M1-2,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,M1-2,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_2 + f_011*tmp0_1 + f_110*tmp0_0 + f_111*tmp0_3;
                  out[INDEX4(i,1,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_11 + f_010*tmp0_9 + f_101*tmp0_13 + f_111*tmp0_10 + tmp0_12*(f_001 + f_100) + tmp0_8*(f_011 + f_110);
                  out[INDEX4(i,2,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_30 + f_011*tmp0_29 + f_110*tmp0_28 + f_111*tmp0_31;
                  out[INDEX4(i,0,1,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_2 + f_011*tmp0_1 + f_110*tmp0_0 + f_111*tmp0_3;
                  out[INDEX4(i,1,1,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_18 + f_001*tmp0_19 + f_010*tmp0_16 + f_011*tmp0_15 + f_100*tmp0_21 + f_101*tmp0_20 + f_110*tmp0_14 + f_111*tmp0_17;
                  out[INDEX4(i,2,1,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_34 + f_011*tmp0_33 + f_110*tmp0_32 + f_111*tmp0_35;
                  out[INDEX4(i,0,2,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_6 + f_011*tmp0_5 + f_110*tmp0_4 + f_111*tmp0_7;
                  out[INDEX4(i,1,2,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_18 + f_001*tmp0_21 + f_010*tmp0_16 + f_011*tmp0_14 + f_100*tmp0_19 + f_101*tmp0_20 + f_110*tmp0_15 + f_111*tmp0_17;
                  out[INDEX4(i,2,2,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_30 + f_011*tmp0_29 + f_110*tmp0_28 + f_111*tmp0_31;
                  out[INDEX4(i,0,3,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_6 + f_011*tmp0_5 + f_110*tmp0_4 + f_111*tmp0_7;
                  out[INDEX4(i,1,3,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_25 + f_010*tmp0_23 + f_101*tmp0_27 + f_111*tmp0_24 + tmp0_22*(f_011 + f_110) + tmp0_26*(f_001 + f_100);
                  out[INDEX4(i,2,3,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_010*tmp0_34 + f_011*tmp0_33 + f_110*tmp0_32 + f_111*tmp0_35;
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 3 */
      if (face_offset(4)>-1) {
         const double tmp0_22 = -0.16666666666666666667/h2;
         const double tmp0_16 = -0.62200846792814621559/h2;
         const double tmp0_33 = 0.16666666666666666667/h2;
         const double tmp0_0 = -0.78867513459481288225/h0;
         const double tmp0_21 = 0.044658198738520451079/h2;
         const double tmp0_17 = -0.044658198738520451079/h2;
         const double tmp0_1 = 0.21132486540518711775/h0;
         const double tmp0_20 = 0.62200846792814621559/h2;
         const double tmp0_14 = 0.21132486540518711775/h1;
         const double tmp0_2 = -0.21132486540518711775/h0;
         const double tmp0_27 = 0.16666666666666666667/h2;
         const double tmp0_15 = -0.78867513459481288225/h1;
         const double tmp0_3 = 0.78867513459481288225/h0;
         const double tmp0_26 = 0.62200846792814621559/h2;
         const double tmp0_12 = -0.21132486540518711775/h1;
         const double tmp0_25 = -0.62200846792814621559/h2;
         const double tmp0_13 = 0.78867513459481288225/h1;
         const double tmp0_24 = -0.044658198738520451079/h2;
         const double tmp0_10 = 0.78867513459481288225/h1;
         const double tmp0_11 = -0.21132486540518711775/h1;
         const double tmp0_34 = 0.044658198738520451079/h2;
         const double tmp0_35 = 0.62200846792814621559/h2;
         const double tmp0_8 = -0.78867513459481288225/h1;
         const double tmp0_29 = 0.044658198738520451079/h2;
         const double tmp0_9 = 0.21132486540518711775/h1;
         const double tmp0_30 = -0.044658198738520451079/h2;
         const double tmp0_28 = 0.16666666666666666667/h2;
         const double tmp0_32 = -0.16666666666666666667/h2;
         const double tmp0_31 = -0.62200846792814621559/h2;
         const double tmp0_18 = -0.16666666666666666667/h2;
         const double tmp0_4 = -0.21132486540518711775/h0;
         const double tmp0_19 = 0.16666666666666666667/h2;
         const double tmp0_5 = 0.78867513459481288225/h0;
         const double tmp0_6 = -0.78867513459481288225/h0;
         const double tmp0_23 = -0.16666666666666666667/h2;
         const double tmp0_7 = 0.21132486540518711775/h0;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,0, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,0, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,0, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,0, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,1, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_0 + f_010*tmp0_2 + f_100*tmp0_3 + f_110*tmp0_1;
                  out[INDEX4(i,1,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_8 + f_010*tmp0_10 + f_100*tmp0_11 + f_110*tmp0_9;
                  out[INDEX4(i,2,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_16 + f_001*tmp0_20 + f_110*tmp0_17 + f_111*tmp0_21 + tmp0_18*(f_010 + f_100) + tmp0_19*(f_011 + f_101);
                  out[INDEX4(i,0,1,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_0 + f_010*tmp0_2 + f_100*tmp0_3 + f_110*tmp0_1;
                  out[INDEX4(i,1,1,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_12 + f_010*tmp0_14 + f_100*tmp0_15 + f_110*tmp0_13;
                  out[INDEX4(i,2,1,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_22 + f_001*tmp0_27 + f_010*tmp0_24 + f_011*tmp0_29 + f_100*tmp0_25 + f_101*tmp0_26 + f_110*tmp0_23 + f_111*tmp0_28;
                  out[INDEX4(i,0,2,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_4 + f_010*tmp0_6 + f_100*tmp0_7 + f_110*tmp0_5;
                  out[INDEX4(i,1,2,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_8 + f_010*tmp0_10 + f_100*tmp0_11 + f_110*tmp0_9;
                  out[INDEX4(i,2,2,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_22 + f_001*tmp0_27 + f_010*tmp0_25 + f_011*tmp0_26 + f_100*tmp0_24 + f_101*tmp0_29 + f_110*tmp0_23 + f_111*tmp0_28;
                  out[INDEX4(i,0,3,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_4 + f_010*tmp0_6 + f_100*tmp0_7 + f_110*tmp0_5;
                  out[INDEX4(i,1,3,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_12 + f_010*tmp0_14 + f_100*tmp0_15 + f_110*tmp0_13;
                  out[INDEX4(i,2,3,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_30 + f_001*tmp0_34 + f_110*tmp0_31 + f_111*tmp0_35 + tmp0_32*(f_010 + f_100) + tmp0_33*(f_011 + f_101);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 4 */
      if (face_offset(5)>-1) {
         const double tmp0_22 = 0.16666666666666666667/h2;
         const double tmp0_16 = 0.62200846792814621559/h2;
         const double tmp0_33 = -0.044658198738520451079/h2;
         const double tmp0_0 = -0.78867513459481288225/h0;
         const double tmp0_21 = -0.044658198738520451079/h2;
         const double tmp0_17 = 0.16666666666666666667/h2;
         const double tmp0_1 = 0.78867513459481288225/h0;
         const double tmp0_20 = -0.16666666666666666667/h2;
         const double tmp0_14 = 0.21132486540518711775/h1;
         const double tmp0_2 = -0.21132486540518711775/h0;
         const double tmp0_27 = -0.62200846792814621559/h2;
         const double tmp0_15 = 0.78867513459481288225/h1;
         const double tmp0_3 = 0.21132486540518711775/h0;
         const double tmp0_26 = -0.16666666666666666667/h2;
         const double tmp0_12 = -0.21132486540518711775/h1;
         const double tmp0_25 = 0.16666666666666666667/h2;
         const double tmp0_13 = -0.78867513459481288225/h1;
         const double tmp0_24 = 0.044658198738520451079/h2;
         const double tmp0_10 = 0.78867513459481288225/h1;
         const double tmp0_11 = 0.21132486540518711775/h1;
         const double tmp0_34 = -0.16666666666666666667/h2;
         const double tmp0_35 = -0.62200846792814621559/h2;
         const double tmp0_8 = -0.78867513459481288225/h1;
         const double tmp0_29 = -0.044658198738520451079/h2;
         const double tmp0_9 = -0.21132486540518711775/h1;
         const double tmp0_30 = 0.044658198738520451079/h2;
         const double tmp0_28 = -0.16666666666666666667/h2;
         const double tmp0_32 = 0.62200846792814621559/h2;
         const double tmp0_31 = 0.16666666666666666667/h2;
         const double tmp0_18 = 0.044658198738520451079/h2;
         const double tmp0_4 = -0.21132486540518711775/h0;
         const double tmp0_19 = -0.62200846792814621559/h2;
         const double tmp0_5 = 0.21132486540518711775/h0;
         const double tmp0_6 = -0.78867513459481288225/h0;
         const double tmp0_23 = 0.62200846792814621559/h2;
         const double tmp0_7 = 0.78867513459481288225/h0;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,M2-1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,M2-1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,M2-1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-1, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,M2-2, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,M2-2, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,M2-2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_0 + f_011*tmp0_2 + f_101*tmp0_1 + f_111*tmp0_3;
                  out[INDEX4(i,1,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_8 + f_011*tmp0_10 + f_101*tmp0_9 + f_111*tmp0_11;
                  out[INDEX4(i,2,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_19 + f_001*tmp0_16 + f_110*tmp0_21 + f_111*tmp0_18 + tmp0_17*(f_011 + f_101) + tmp0_20*(f_010 + f_100);
                  out[INDEX4(i,0,1,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_0 + f_011*tmp0_2 + f_101*tmp0_1 + f_111*tmp0_3;
                  out[INDEX4(i,1,1,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_12 + f_011*tmp0_14 + f_101*tmp0_13 + f_111*tmp0_15;
                  out[INDEX4(i,2,1,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_26 + f_001*tmp0_22 + f_010*tmp0_29 + f_011*tmp0_24 + f_100*tmp0_27 + f_101*tmp0_23 + f_110*tmp0_28 + f_111*tmp0_25;
                  out[INDEX4(i,0,2,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_4 + f_011*tmp0_6 + f_101*tmp0_5 + f_111*tmp0_7;
                  out[INDEX4(i,1,2,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_8 + f_011*tmp0_10 + f_101*tmp0_9 + f_111*tmp0_11;
                  out[INDEX4(i,2,2,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_26 + f_001*tmp0_22 + f_010*tmp0_27 + f_011*tmp0_23 + f_100*tmp0_29 + f_101*tmp0_24 + f_110*tmp0_28 + f_111*tmp0_25;
                  out[INDEX4(i,0,3,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_4 + f_011*tmp0_6 + f_101*tmp0_5 + f_111*tmp0_7;
                  out[INDEX4(i,1,3,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_001*tmp0_12 + f_011*tmp0_14 + f_101*tmp0_13 + f_111*tmp0_15;
                  out[INDEX4(i,2,3,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,4)] = f_000*tmp0_33 + f_001*tmp0_30 + f_110*tmp0_35 + f_111*tmp0_32 + tmp0_31*(f_011 + f_101) + tmp0_34*(f_010 + f_100);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 5 */
   } else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_0 = -0.25/h0;
         const double tmp0_4 = -0.5/h2;
         const double tmp0_1 = 0.25/h0;
         const double tmp0_5 = 0.5/h2;
         const double tmp0_2 = -0.5/h1;
         const double tmp0_3 = 0.5/h1;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(0,k1,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(0,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(0,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(1,k1,k2, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(0,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_000 + f_001 + f_010 + f_011) + tmp0_1*(f_100 + f_101 + f_110 + f_111);
                  out[INDEX4(i,1,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_000 + f_001) + tmp0_3*(f_010 + f_011);
                  out[INDEX4(i,2,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_000 + f_010) + tmp0_5*(f_001 + f_011);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_0 = -0.25/h0;
         const double tmp0_4 = 0.5/h2;
         const double tmp0_1 = 0.25/h0;
         const double tmp0_5 = -0.5/h2;
         const double tmp0_2 = -0.5/h1;
         const double tmp0_3 = 0.5/h1;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(M0-2,k1,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(M0-2,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(M0-1,k1,k2+1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(M0-2,k1+1,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(M0-1,k1,k2, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(M0-2,k1+1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_000 + f_001 + f_010 + f_011) + tmp0_1*(f_100 + f_101 + f_110 + f_111);
                  out[INDEX4(i,1,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_100 + f_101) + tmp0_3*(f_110 + f_111);
                  out[INDEX4(i,2,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_101 + f_111) + tmp0_5*(f_100 + f_110);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_0 = -0.5/h0;
         const double tmp0_4 = -0.5/h2;
         const double tmp0_1 = 0.5/h0;
         const double tmp0_5 = 0.5/h2;
         const double tmp0_2 = -0.25/h1;
         const double tmp0_3 = 0.25/h1;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,0,k2, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,0,k2+1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,0,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,0,k2, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,1,k2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,1,k2+1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_000 + f_001) + tmp0_1*(f_100 + f_101);
                  out[INDEX4(i,1,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_000 + f_001 + f_100 + f_101) + tmp0_3*(f_010 + f_011 + f_110 + f_111);
                  out[INDEX4(i,2,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_000 + f_100) + tmp0_5*(f_001 + f_101);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_0 = -0.5/h0;
         const double tmp0_4 = 0.5/h2;
         const double tmp0_1 = 0.5/h0;
         const double tmp0_5 = -0.5/h2;
         const double tmp0_2 = 0.25/h1;
         const double tmp0_3 = -0.25/h1;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,M1-1,k2+1, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,M1-1,k2, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2+1, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,M1-2,k2, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,M1-2,k2+1, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,M1-2,k2+1, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,M1-2,k2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_010 + f_011) + tmp0_1*(f_110 + f_111);
                  out[INDEX4(i,1,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_010 + f_011 + f_110 + f_111) + tmp0_3*(f_000 + f_001 + f_100 + f_101);
                  out[INDEX4(i,2,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_011 + f_111) + tmp0_5*(f_010 + f_110);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 3 */
      if (face_offset(4)>-1) {
         const double tmp0_0 = -0.5/h0;
         const double tmp0_4 = -0.25/h2;
         const double tmp0_1 = 0.5/h0;
         const double tmp0_5 = 0.25/h2;
         const double tmp0_2 = -0.5/h1;
         const double tmp0_3 = 0.5/h1;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,0, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,0, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,0, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,0, M0,M1),NCOMP)];
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,1, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_000 + f_010) + tmp0_1*(f_100 + f_110);
                  out[INDEX4(i,1,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_000 + f_100) + tmp0_3*(f_010 + f_110);
                  out[INDEX4(i,2,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_000 + f_010 + f_100 + f_110) + tmp0_5*(f_001 + f_011 + f_101 + f_111);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 4 */
      if (face_offset(5)>-1) {
         const double tmp0_0 = -0.5/h0;
         const double tmp0_4 = 0.25/h2;
         const double tmp0_1 = 0.5/h0;
         const double tmp0_5 = -0.25/h2;
         const double tmp0_2 = -0.5/h1;
         const double tmp0_3 = 0.5/h1;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  const register double f_001 = in[INDEX2(i,INDEX3(k0,k1,M2-1, M0,M1),NCOMP)];
                  const register double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,M2-1, M0,M1),NCOMP)];
                  const register double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,M2-1, M0,M1),NCOMP)];
                  const register double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-1, M0,M1),NCOMP)];
                  const register double f_000 = in[INDEX2(i,INDEX3(k0,k1,M2-2, M0,M1),NCOMP)];
                  const register double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,M2-2, M0,M1),NCOMP)];
                  const register double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-2, M0,M1),NCOMP)];
                  const register double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,M2-2, M0,M1),NCOMP)];
                  out[INDEX4(i,0,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_0*(f_001 + f_011) + tmp0_1*(f_101 + f_111);
                  out[INDEX4(i,1,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_2*(f_001 + f_101) + tmp0_3*(f_011 + f_111);
                  out[INDEX4(i,2,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,3,1)] = tmp0_4*(f_001 + f_011 + f_101 + f_111) + tmp0_5*(f_000 + f_010 + f_100 + f_110);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 5 */
   } /* end of out_data_type branching */
   /* GENERATOR SNIP BOTTOM */
}
