Ripley_Assemble_Interpolation_2D(Ripley_Grid *grid, Escript in, Escript out)
{
   /* GENERATOR SNIP TOP */
   if (out_data_type==RIPLEY_ELEMENTS) {
      const double tmp0_0 = 0.10000000000000000000;
      const double tmp0_4 = 0.056350832689629155741;
      const double tmp0_1 = 0.012701665379258311482;
      const double tmp0_5 = 0.25000000000000000000;
      const double tmp0_2 = 0.78729833462074168852;
      const double tmp0_3 = 0.44364916731037084426;
      #pragma omp parallel for private(i,k1,k0)
      for (k1 =0; k1 < N1; ++k1) {
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,k1, M0),NCOMP)];
               out[INDEX3(i,0,INDEX2(k0,k1,N0),NCOMP,9)] = f_00*tmp0_2 + f_11*tmp0_1 + tmp0_0*(f_01 + f_10);
               out[INDEX3(i,1,INDEX2(k0,k1,N0),NCOMP,9)] = tmp0_3*(f_00 + f_10) + tmp0_4*(f_01 + f_11);
               out[INDEX3(i,2,INDEX2(k0,k1,N0),NCOMP,9)] = f_01*tmp0_1 + f_10*tmp0_2 + tmp0_0*(f_00 + f_11);
               out[INDEX3(i,3,INDEX2(k0,k1,N0),NCOMP,9)] = tmp0_3*(f_00 + f_01) + tmp0_4*(f_10 + f_11);
               out[INDEX3(i,4,INDEX2(k0,k1,N0),NCOMP,9)] = tmp0_5*(f_00 + f_01 + f_10 + f_11);
               out[INDEX3(i,5,INDEX2(k0,k1,N0),NCOMP,9)] = tmp0_3*(f_10 + f_11) + tmp0_4*(f_00 + f_01);
               out[INDEX3(i,6,INDEX2(k0,k1,N0),NCOMP,9)] = f_01*tmp0_2 + f_10*tmp0_1 + tmp0_0*(f_00 + f_11);
               out[INDEX3(i,7,INDEX2(k0,k1,N0),NCOMP,9)] = tmp0_3*(f_01 + f_11) + tmp0_4*(f_00 + f_10);
               out[INDEX3(i,8,INDEX2(k0,k1,N0),NCOMP,9)] = f_00*tmp0_1 + f_11*tmp0_2 + tmp0_0*(f_01 + f_10);
            } /* close component loop i */
         } /* close k0 loop */
      } /* close k1 loop */
   } else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {
      const double tmp0_0 = 0.25000000000000000000;
      #pragma omp parallel for private(i,k1,k0)
      for (k1 =0; k1 < N1; ++k1) {
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_01 = in[INDEX2(i,INDEX2(k0,k1+1, M0),NCOMP)];
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,k1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,k1+1, M0),NCOMP)];
               out[INDEX3(i,0,INDEX2(k0,k1,N0),NCOMP,1)] = tmp0_0*(f_00 + f_01 + f_10 + f_11);
            } /* close component loop i */
         } /* close k0 loop */
      } /* close k1 loop */
   } else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_2 = 0.50000000000000000000;
         const double tmp0_1 = 0.88729833462074168852;
         const double tmp0_0 = 0.11270166537925831148;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_01 = in[INDEX2(i,INDEX2(0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(0,k1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,3)] = f_00*tmp0_1 + f_01*tmp0_0;
               out[INDEX3(i,1,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,3)] = tmp0_2*(f_00 + f_01);
               out[INDEX3(i,2,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,3)] = f_00*tmp0_0 + f_01*tmp0_1;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_2 = 0.50000000000000000000;
         const double tmp0_1 = 0.11270166537925831148;
         const double tmp0_0 = 0.88729833462074168852;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(M0-1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(M0-1,k1+1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,3)] = f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX3(i,1,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,3)] = tmp0_2*(f_10 + f_11);
               out[INDEX3(i,2,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,3)] = f_10*tmp0_1 + f_11*tmp0_0;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_2 = 0.50000000000000000000;
         const double tmp0_1 = 0.88729833462074168852;
         const double tmp0_0 = 0.11270166537925831148;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,0, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,0, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,3)] = f_00*tmp0_1 + f_10*tmp0_0;
               out[INDEX3(i,1,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,3)] = tmp0_2*(f_00 + f_10);
               out[INDEX3(i,2,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,3)] = f_00*tmp0_0 + f_10*tmp0_1;
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_2 = 0.50000000000000000000;
         const double tmp0_1 = 0.88729833462074168852;
         const double tmp0_0 = 0.11270166537925831148;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,M1-1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,M1-1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,3)] = f_01*tmp0_1 + f_11*tmp0_0;
               out[INDEX3(i,1,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,3)] = tmp0_2*(f_01 + f_11);
               out[INDEX3(i,2,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,3)] = f_01*tmp0_0 + f_11*tmp0_1;
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 3 */
   } else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_0 = 0.50000000000000000000;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_00 = in[INDEX2(i,INDEX2(0,k1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(0,k1+1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,1)] = tmp0_0*(f_00 + f_01);
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_0 = 0.50000000000000000000;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(M0-1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(M0-1,k1+1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,1)] = tmp0_0*(f_10 + f_11);
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_0 = 0.50000000000000000000;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,0, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,0, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,1)] = tmp0_0*(f_00 + f_10);
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_0 = 0.50000000000000000000;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_01 = in[INDEX2(i,INDEX2(k0,M1-1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,M1-1, M0),NCOMP)];
               out[INDEX3(i,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,1)] = tmp0_0*(f_01 + f_11);
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 3 */
   } /* end of out_data_type branching
   /* GENERATOR SNIP BOTTOM */
}
