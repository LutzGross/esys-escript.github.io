Ripley_Assemble_Interpolation_3D(Ripley_Grid *grid, Escript in, Escript out)
{
   /* GENERATOR SNIP TOP */
   if (out_data_type==RIPLEY_ELEMENTS) {
      const double tmp0_0 = 0.69856850115866751967;
      const double tmp0_4 = 0.39364916731037084426;
      const double tmp0_1 = 0.088729833462074168852;
      const double tmp0_8 = 0.028175416344814577871;
      const double tmp0_5 = 0.050000000000000000000;
      const double tmp0_2 = 0.011270166537925831148;
      const double tmp0_9 = 0.12500000000000000000;
      const double tmp0_6 = 0.0063508326896291557410;
      const double tmp0_3 = 0.0014314988413324803339;
      const double tmp0_7 = 0.22182458365518542213;
      #pragma omp parallel for private(i,k2,k1,k0)
      for (k2 =0; k2 < N2; ++k2) {
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,k1,k2, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,k2, M0,M1),NCOMP)];
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,k2+1, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_000*tmp0_0 + f_111*tmp0_3 + tmp0_1*(f_001 + f_010 + f_100) + tmp0_2*(f_011 + f_101 + f_110);
                  out[INDEX3(i,1,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_000 + f_100) + tmp0_5*(f_001 + f_010 + f_101 + f_110) + tmp0_6*(f_011 + f_111);
                  out[INDEX3(i,2,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_011*tmp0_3 + f_100*tmp0_0 + tmp0_1*(f_000 + f_101 + f_110) + tmp0_2*(f_001 + f_010 + f_111);
                  out[INDEX3(i,3,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_000 + f_010) + tmp0_5*(f_001 + f_011 + f_100 + f_110) + tmp0_6*(f_101 + f_111);
                  out[INDEX3(i,4,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_000 + f_010 + f_100 + f_110) + tmp0_8*(f_001 + f_011 + f_101 + f_111);
                  out[INDEX3(i,5,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_100 + f_110) + tmp0_5*(f_000 + f_010 + f_101 + f_111) + tmp0_6*(f_001 + f_011);
                  out[INDEX3(i,6,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_010*tmp0_0 + f_101*tmp0_3 + tmp0_1*(f_000 + f_011 + f_110) + tmp0_2*(f_001 + f_100 + f_111);
                  out[INDEX3(i,7,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_010 + f_110) + tmp0_5*(f_000 + f_011 + f_100 + f_111) + tmp0_6*(f_001 + f_101);
                  out[INDEX3(i,8,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_001*tmp0_3 + f_110*tmp0_0 + tmp0_1*(f_010 + f_100 + f_111) + tmp0_2*(f_000 + f_011 + f_101);
                  out[INDEX3(i,9,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_000 + f_001) + tmp0_5*(f_010 + f_011 + f_100 + f_101) + tmp0_6*(f_110 + f_111);
                  out[INDEX3(i,10,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_000 + f_001 + f_100 + f_101) + tmp0_8*(f_010 + f_011 + f_110 + f_111);
                  out[INDEX3(i,11,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_100 + f_101) + tmp0_5*(f_000 + f_001 + f_110 + f_111) + tmp0_6*(f_010 + f_011);
                  out[INDEX3(i,12,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_000 + f_001 + f_010 + f_011) + tmp0_8*(f_100 + f_101 + f_110 + f_111);
                  out[INDEX3(i,13,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_9*(f_000 + f_001 + f_010 + f_011 + f_100 + f_101 + f_110 + f_111);
                  out[INDEX3(i,14,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_100 + f_101 + f_110 + f_111) + tmp0_8*(f_000 + f_001 + f_010 + f_011);
                  out[INDEX3(i,15,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_010 + f_011) + tmp0_5*(f_000 + f_001 + f_110 + f_111) + tmp0_6*(f_100 + f_101);
                  out[INDEX3(i,16,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_010 + f_011 + f_110 + f_111) + tmp0_8*(f_000 + f_001 + f_100 + f_101);
                  out[INDEX3(i,17,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_110 + f_111) + tmp0_5*(f_010 + f_011 + f_100 + f_101) + tmp0_6*(f_000 + f_001);
                  out[INDEX3(i,18,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_001*tmp0_0 + f_110*tmp0_3 + tmp0_1*(f_000 + f_011 + f_101) + tmp0_2*(f_010 + f_100 + f_111);
                  out[INDEX3(i,19,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_001 + f_101) + tmp0_5*(f_000 + f_011 + f_100 + f_111) + tmp0_6*(f_010 + f_110);
                  out[INDEX3(i,20,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_010*tmp0_3 + f_101*tmp0_0 + tmp0_1*(f_001 + f_100 + f_111) + tmp0_2*(f_000 + f_011 + f_110);
                  out[INDEX3(i,21,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_001 + f_011) + tmp0_5*(f_000 + f_010 + f_101 + f_111) + tmp0_6*(f_100 + f_110);
                  out[INDEX3(i,22,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_7*(f_001 + f_011 + f_101 + f_111) + tmp0_8*(f_000 + f_010 + f_100 + f_110);
                  out[INDEX3(i,23,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_101 + f_111) + tmp0_5*(f_001 + f_011 + f_100 + f_110) + tmp0_6*(f_000 + f_010);
                  out[INDEX3(i,24,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_011*tmp0_0 + f_100*tmp0_3 + tmp0_1*(f_001 + f_010 + f_111) + tmp0_2*(f_000 + f_101 + f_110);
                  out[INDEX3(i,25,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = tmp0_4*(f_011 + f_111) + tmp0_5*(f_001 + f_010 + f_101 + f_110) + tmp0_6*(f_000 + f_100);
                  out[INDEX3(i,26,INDEX3(k0,k1,k2,N0,N1),NCOMP,27)] = f_000*tmp0_3 + f_111*tmp0_0 + tmp0_1*(f_011 + f_101 + f_110) + tmp0_2*(f_001 + f_010 + f_100);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* close k2 loop */
   } else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {
      const double tmp0_0 = 0.12500000000000000000;
      #pragma omp parallel for private(i,k2,k1,k0)
      for (k2 =0; k2 < N2; ++k2) {
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,k1,k2, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,k2+1, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,k2, M0,M1),NCOMP)];
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_000 + f_001 + f_010 + f_011 + f_100 + f_101 + f_110 + f_111);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* close k2 loop */
   } else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_0 = 0.78729833462074168852;
         const double tmp0_4 = 0.056350832689629155741;
         const double tmp0_1 = 0.10000000000000000000;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.44364916731037084426;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(0,k1,k2, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(0,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_011 = in[INDEX2(i,INDEX3(0,k1+1,k2+1, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(0,k1+1,k2, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_0 + f_011*tmp0_2 + tmp0_1*(f_001 + f_010);
                  out[INDEX3(i,1,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_010) + tmp0_4*(f_001 + f_011);
                  out[INDEX3(i,2,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_2 + f_010*tmp0_0 + tmp0_1*(f_000 + f_011);
                  out[INDEX3(i,3,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_001) + tmp0_4*(f_010 + f_011);
                  out[INDEX3(i,4,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_000 + f_001 + f_010 + f_011);
                  out[INDEX3(i,5,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_010 + f_011) + tmp0_4*(f_000 + f_001);
                  out[INDEX3(i,6,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_0 + f_010*tmp0_2 + tmp0_1*(f_000 + f_011);
                  out[INDEX3(i,7,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_001 + f_011) + tmp0_4*(f_000 + f_010);
                  out[INDEX3(i,8,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_2 + f_011*tmp0_0 + tmp0_1*(f_001 + f_010);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_0 = 0.10000000000000000000;
         const double tmp0_4 = 0.44364916731037084426;
         const double tmp0_1 = 0.78729833462074168852;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.056350832689629155741;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_101 = in[INDEX2(i,INDEX3(M0-1,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(M0-1,k1,k2, M0,M1),NCOMP)];
                  register const double f_110 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_100*tmp0_1 + f_111*tmp0_2 + tmp0_0*(f_101 + f_110);
                  out[INDEX3(i,1,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_101 + f_111) + tmp0_4*(f_100 + f_110);
                  out[INDEX3(i,2,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_101*tmp0_2 + f_110*tmp0_1 + tmp0_0*(f_100 + f_111);
                  out[INDEX3(i,3,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_110 + f_111) + tmp0_4*(f_100 + f_101);
                  out[INDEX3(i,4,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_100 + f_101 + f_110 + f_111);
                  out[INDEX3(i,5,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_100 + f_101) + tmp0_4*(f_110 + f_111);
                  out[INDEX3(i,6,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_101*tmp0_1 + f_110*tmp0_2 + tmp0_0*(f_100 + f_111);
                  out[INDEX3(i,7,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_100 + f_110) + tmp0_4*(f_101 + f_111);
                  out[INDEX3(i,8,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_100*tmp0_2 + f_111*tmp0_1 + tmp0_0*(f_101 + f_110);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_0 = 0.78729833462074168852;
         const double tmp0_4 = 0.056350832689629155741;
         const double tmp0_1 = 0.10000000000000000000;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.44364916731037084426;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,0,k2, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,0,k2+1, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,0,k2+1, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,0,k2, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_0 + f_101*tmp0_2 + tmp0_1*(f_001 + f_100);
                  out[INDEX3(i,1,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_100) + tmp0_4*(f_001 + f_101);
                  out[INDEX3(i,2,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_2 + f_100*tmp0_0 + tmp0_1*(f_000 + f_101);
                  out[INDEX3(i,3,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_001) + tmp0_4*(f_100 + f_101);
                  out[INDEX3(i,4,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_000 + f_001 + f_100 + f_101);
                  out[INDEX3(i,5,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_100 + f_101) + tmp0_4*(f_000 + f_001);
                  out[INDEX3(i,6,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_0 + f_100*tmp0_2 + tmp0_1*(f_000 + f_101);
                  out[INDEX3(i,7,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_001 + f_101) + tmp0_4*(f_000 + f_100);
                  out[INDEX3(i,8,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_2 + f_101*tmp0_0 + tmp0_1*(f_001 + f_100);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_0 = 0.10000000000000000000;
         const double tmp0_4 = 0.44364916731037084426;
         const double tmp0_1 = 0.78729833462074168852;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.056350832689629155741;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,M1-1,k2+1, M0,M1),NCOMP)];
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,M1-1,k2, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_010*tmp0_1 + f_111*tmp0_2 + tmp0_0*(f_011 + f_110);
                  out[INDEX3(i,1,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_011 + f_111) + tmp0_4*(f_010 + f_110);
                  out[INDEX3(i,2,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_011*tmp0_2 + f_110*tmp0_1 + tmp0_0*(f_010 + f_111);
                  out[INDEX3(i,3,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_110 + f_111) + tmp0_4*(f_010 + f_011);
                  out[INDEX3(i,4,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_010 + f_011 + f_110 + f_111);
                  out[INDEX3(i,5,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_010 + f_011) + tmp0_4*(f_110 + f_111);
                  out[INDEX3(i,6,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_011*tmp0_1 + f_110*tmp0_2 + tmp0_0*(f_010 + f_111);
                  out[INDEX3(i,7,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_010 + f_110) + tmp0_4*(f_011 + f_111);
                  out[INDEX3(i,8,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_010*tmp0_2 + f_111*tmp0_1 + tmp0_0*(f_011 + f_110);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 3 */
      if (face_offset(4)>-1) {
         const double tmp0_0 = 0.78729833462074168852;
         const double tmp0_4 = 0.056350832689629155741;
         const double tmp0_1 = 0.10000000000000000000;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.44364916731037084426;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,k1,0, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,0, M0,M1),NCOMP)];
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,0, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,0, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_0 + f_110*tmp0_2 + tmp0_1*(f_010 + f_100);
                  out[INDEX3(i,1,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_100) + tmp0_4*(f_010 + f_110);
                  out[INDEX3(i,2,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_010*tmp0_2 + f_100*tmp0_0 + tmp0_1*(f_000 + f_110);
                  out[INDEX3(i,3,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_000 + f_010) + tmp0_4*(f_100 + f_110);
                  out[INDEX3(i,4,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_000 + f_010 + f_100 + f_110);
                  out[INDEX3(i,5,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_100 + f_110) + tmp0_4*(f_000 + f_010);
                  out[INDEX3(i,6,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_010*tmp0_0 + f_100*tmp0_2 + tmp0_1*(f_000 + f_110);
                  out[INDEX3(i,7,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_010 + f_110) + tmp0_4*(f_000 + f_100);
                  out[INDEX3(i,8,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_000*tmp0_2 + f_110*tmp0_0 + tmp0_1*(f_010 + f_100);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 4 */
      if (face_offset(5)>-1) {
         const double tmp0_0 = 0.78729833462074168852;
         const double tmp0_4 = 0.056350832689629155741;
         const double tmp0_1 = 0.10000000000000000000;
         const double tmp0_5 = 0.25000000000000000000;
         const double tmp0_2 = 0.012701665379258311482;
         const double tmp0_3 = 0.44364916731037084426;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,k1,M2-1, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,M2-1, M0,M1),NCOMP)];
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,M2-1, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_0 + f_111*tmp0_2 + tmp0_1*(f_011 + f_101);
                  out[INDEX3(i,1,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_001 + f_101) + tmp0_4*(f_011 + f_111);
                  out[INDEX3(i,2,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_011*tmp0_2 + f_101*tmp0_0 + tmp0_1*(f_001 + f_111);
                  out[INDEX3(i,3,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_001 + f_011) + tmp0_4*(f_101 + f_111);
                  out[INDEX3(i,4,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_5*(f_001 + f_011 + f_101 + f_111);
                  out[INDEX3(i,5,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_101 + f_111) + tmp0_4*(f_001 + f_011);
                  out[INDEX3(i,6,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_011*tmp0_0 + f_101*tmp0_2 + tmp0_1*(f_001 + f_111);
                  out[INDEX3(i,7,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = tmp0_3*(f_011 + f_111) + tmp0_4*(f_001 + f_101);
                  out[INDEX3(i,8,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,9)] = f_001*tmp0_2 + f_111*tmp0_0 + tmp0_1*(f_011 + f_101);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 5 */
   } else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_011 = in[INDEX2(i,INDEX3(0,k1+1,k2+1, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(0,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(0,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_000 = in[INDEX2(i,INDEX3(0,k1,k2, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(0)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_000 + f_001 + f_010 + f_011);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k2,k1)
         for (k2 =0; k2 < N2; ++k2) {
            for (k1 =0; k1 < N1; ++k1) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_110 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(M0-1,k1,k2, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(M0-1,k1,k2+1, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(M0-1,k1+1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(1)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_100 + f_101 + f_110 + f_111);
               } /* close component loop i */
            } /* close k1 loop */
         } /* close k2 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,0,k2, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,0,k2, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,0,k2+1, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,0,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(2)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_000 + f_001 + f_100 + f_101);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k2,k0)
         for (k2 =0; k2 < N2; ++k2) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2, M0,M1),NCOMP)];
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,M1-1,k2+1, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,M1-1,k2, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,M1-1,k2+1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(3)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_010 + f_011 + f_110 + f_111);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k2 loop */
      } /* end of face 3 */
      if (face_offset(4)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_110 = in[INDEX2(i,INDEX3(k0+1,k1+1,0, M0,M1),NCOMP)];
                  register const double f_010 = in[INDEX2(i,INDEX3(k0,k1+1,0, M0,M1),NCOMP)];
                  register const double f_100 = in[INDEX2(i,INDEX3(k0+1,k1,0, M0,M1),NCOMP)];
                  register const double f_000 = in[INDEX2(i,INDEX3(k0,k1,0, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(4)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_000 + f_010 + f_100 + f_110);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 4 */
      if (face_offset(5)>-1) {
         const double tmp0_0 = 0.25000000000000000000;
         #pragma omp parallel for private(i,k1,k0)
         for (k1 =0; k1 < N1; ++k1) {
            for (k0 =0; k0 < N0; ++k0) {
               for (i =0; i < NCOMP; ++i) {
                  register const double f_011 = in[INDEX2(i,INDEX3(k0,k1+1,M2-1, M0,M1),NCOMP)];
                  register const double f_001 = in[INDEX2(i,INDEX3(k0,k1,M2-1, M0,M1),NCOMP)];
                  register const double f_101 = in[INDEX2(i,INDEX3(k0+1,k1,M2-1, M0,M1),NCOMP)];
                  register const double f_111 = in[INDEX2(i,INDEX3(k0+1,k1+1,M2-1, M0,M1),NCOMP)];
                  out[INDEX3(i,0,face_offset(5)+INDEX3(k0,k1,k2,N0,N1),NCOMP,1)] = tmp0_0*(f_001 + f_011 + f_101 + f_111);
               } /* close component loop i */
            } /* close k0 loop */
         } /* close k1 loop */
      } /* end of face 5 */
   } /* end of out_data_type branching
   /* GENERATOR SNIP BOTTOM */
}
