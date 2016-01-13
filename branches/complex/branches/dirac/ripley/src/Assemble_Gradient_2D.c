Ripley_Assemble_Gradient_2D(Ripley_Grid *grid, Escript in, Escript out)
{
   /* GENERATOR SNIP TOP */
   if (out_data_type==RIPLEY_ELEMENTS) {
      const double tmp0_16 = -0.88729833462074168852/h1;
      const double tmp0_18 = 0.11270166537925831148/h1;
      const double tmp0_13 = -0.88729833462074168852/h1;
      const double tmp0_0 = 0.88729833462074168852/h0;
      const double tmp0_17 = 0.88729833462074168852/h1;
      const double tmp0_4 = 0.5/h0;
      const double tmp0_19 = -0.11270166537925831148/h1;
      const double tmp0_10 = -0.11270166537925831148/h1;
      const double tmp0_1 = 0.11270166537925831148/h0;
      const double tmp0_8 = -0.88729833462074168852/h0;
      const double tmp0_14 = -0.5/h1;
      const double tmp0_5 = -0.5/h0;
      const double tmp0_11 = 0.11270166537925831148/h1;
      const double tmp0_2 = -0.11270166537925831148/h0;
      const double tmp0_9 = -0.11270166537925831148/h0;
      const double tmp0_15 = 0.5/h1;
      const double tmp0_6 = 0.11270166537925831148/h0;
      const double tmp0_3 = -0.88729833462074168852/h0;
      const double tmp0_12 = 0.88729833462074168852/h1;
      const double tmp0_7 = 0.88729833462074168852/h0;
      #pragma omp parallel for private(i,k1,k0)
      for (k1 =0; k1 < N1; ++k1) {
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_3 + f_01*tmp0_2 + f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX4(i,1,0,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_13 + f_01*tmp0_12 + f_10*tmp0_10 + f_11*tmp0_11;
               out[INDEX4(i,0,1,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_3 + f_01*tmp0_2 + f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX4(i,1,1,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_14*(f_00 + f_10) + tmp0_15*(f_01 + f_11);
               out[INDEX4(i,0,2,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_3 + f_01*tmp0_2 + f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX4(i,1,2,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_19 + f_01*tmp0_18 + f_10*tmp0_16 + f_11*tmp0_17;
               out[INDEX4(i,0,3,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_4*(f_10 + f_11) + tmp0_5*(f_00 + f_01);
               out[INDEX4(i,1,3,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_13 + f_01*tmp0_12 + f_10*tmp0_10 + f_11*tmp0_11;
               out[INDEX4(i,0,4,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_4*(f_10 + f_11) + tmp0_5*(f_00 + f_01);
               out[INDEX4(i,1,4,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_14*(f_00 + f_10) + tmp0_15*(f_01 + f_11);
               out[INDEX4(i,0,5,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_4*(f_10 + f_11) + tmp0_5*(f_00 + f_01);
               out[INDEX4(i,1,5,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_19 + f_01*tmp0_18 + f_10*tmp0_16 + f_11*tmp0_17;
               out[INDEX4(i,0,6,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_9 + f_01*tmp0_8 + f_10*tmp0_6 + f_11*tmp0_7;
               out[INDEX4(i,1,6,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_13 + f_01*tmp0_12 + f_10*tmp0_10 + f_11*tmp0_11;
               out[INDEX4(i,0,7,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_9 + f_01*tmp0_8 + f_10*tmp0_6 + f_11*tmp0_7;
               out[INDEX4(i,1,7,INDEX2(k0,k1,N0),NCOMP,2,9)] = tmp0_14*(f_00 + f_10) + tmp0_15*(f_01 + f_11);
               out[INDEX4(i,0,8,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_9 + f_01*tmp0_8 + f_10*tmp0_6 + f_11*tmp0_7;
               out[INDEX4(i,1,8,INDEX2(k0,k1,N0),NCOMP,2,9)] = f_00*tmp0_19 + f_01*tmp0_18 + f_10*tmp0_16 + f_11*tmp0_17;
            } /* close component loop i */
         } /* close k0 loop */
      } /* close k1 loop */
   } else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {
      const double tmp0_3 = 0.5/h1;
      const double tmp0_2 = -0.5/h1;
      const double tmp0_1 = -0.5/h0;
      const double tmp0_0 = 0.5/h0;
      #pragma omp parallel for private(i,k1,k0)
      for (k1 =0; k1 < N1; ++k1) {
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_0*(f_10 + f_11) + tmp0_1*(f_00 + f_01);
               out[INDEX4(i,1,0,INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_2*(f_00 + f_10) + tmp0_3*(f_01 + f_11);
            } /* close component loop i */
         } /* close k0 loop */
      } /* close k1 loop */
   } else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_0 = 0.88729833462074168852/h0;
         const double tmp0_4 = 0.5/h0;
         const double tmp0_10 = 1.0/h1;
         const double tmp0_1 = 0.11270166537925831148/h0;
         const double tmp0_8 = -0.88729833462074168852/h0;
         const double tmp0_5 = -0.5/h0;
         const double tmp0_11 = -1/h1;
         const double tmp0_2 = -0.11270166537925831148/h0;
         const double tmp0_9 = -0.11270166537925831148/h0;
         const double tmp0_6 = 0.11270166537925831148/h0;
         const double tmp0_3 = -0.88729833462074168852/h0;
         const double tmp0_7 = 0.88729833462074168852/h0;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(0,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_3 + f_01*tmp0_2 + f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX4(i,1,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_11 + f_01*tmp0_10;
               out[INDEX4(i,0,1,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = tmp0_4*(f_10 + f_11) + tmp0_5*(f_00 + f_01);
               out[INDEX4(i,1,1,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_11 + f_01*tmp0_10;
               out[INDEX4(i,0,2,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_9 + f_01*tmp0_8 + f_10*tmp0_6 + f_11*tmp0_7;
               out[INDEX4(i,1,2,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_11 + f_01*tmp0_10;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_0 = 0.88729833462074168852/h0;
         const double tmp0_4 = 0.5/h0;
         const double tmp0_10 = -1/h1;
         const double tmp0_1 = 0.11270166537925831148/h0;
         const double tmp0_8 = -0.88729833462074168852/h0;
         const double tmp0_5 = -0.5/h0;
         const double tmp0_11 = 1.0/h1;
         const double tmp0_2 = -0.11270166537925831148/h0;
         const double tmp0_9 = -0.11270166537925831148/h0;
         const double tmp0_6 = 0.11270166537925831148/h0;
         const double tmp0_3 = -0.88729833462074168852/h0;
         const double tmp0_7 = 0.88729833462074168852/h0;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(M0-1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(M0-1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(M0-2,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(M0-2,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_3 + f_01*tmp0_2 + f_10*tmp0_0 + f_11*tmp0_1;
               out[INDEX4(i,1,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_10*tmp0_10 + f_11*tmp0_11;
               out[INDEX4(i,0,1,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = tmp0_4*(f_10 + f_11) + tmp0_5*(f_00 + f_01);
               out[INDEX4(i,1,1,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_10*tmp0_10 + f_11*tmp0_11;
               out[INDEX4(i,0,2,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_9 + f_01*tmp0_8 + f_10*tmp0_6 + f_11*tmp0_7;
               out[INDEX4(i,1,2,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_10*tmp0_10 + f_11*tmp0_11;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_0 = 1.0/h0;
         const double tmp0_4 = 0.11270166537925831148/h1;
         const double tmp0_10 = 0.88729833462074168852/h1;
         const double tmp0_1 = -1/h0;
         const double tmp0_8 = -0.88729833462074168852/h1;
         const double tmp0_5 = 0.88729833462074168852/h1;
         const double tmp0_11 = 0.11270166537925831148/h1;
         const double tmp0_2 = -0.11270166537925831148/h1;
         const double tmp0_9 = -0.11270166537925831148/h1;
         const double tmp0_6 = -0.5/h1;
         const double tmp0_3 = -0.88729833462074168852/h1;
         const double tmp0_7 = 0.5/h1;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,0, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,0, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_1 + f_10*tmp0_0;
               out[INDEX4(i,1,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_3 + f_01*tmp0_5 + f_10*tmp0_2 + f_11*tmp0_4;
               out[INDEX4(i,0,1,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_1 + f_10*tmp0_0;
               out[INDEX4(i,1,1,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = tmp0_6*(f_00 + f_10) + tmp0_7*(f_01 + f_11);
               out[INDEX4(i,0,2,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_1 + f_10*tmp0_0;
               out[INDEX4(i,1,2,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_9 + f_01*tmp0_11 + f_10*tmp0_8 + f_11*tmp0_10;
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_0 = 1.0/h0;
         const double tmp0_4 = -0.11270166537925831148/h1;
         const double tmp0_10 = -0.88729833462074168852/h1;
         const double tmp0_1 = -1/h0;
         const double tmp0_8 = 0.88729833462074168852/h1;
         const double tmp0_5 = -0.88729833462074168852/h1;
         const double tmp0_11 = -0.11270166537925831148/h1;
         const double tmp0_2 = 0.11270166537925831148/h1;
         const double tmp0_9 = 0.11270166537925831148/h1;
         const double tmp0_6 = 0.5/h1;
         const double tmp0_3 = 0.88729833462074168852/h1;
         const double tmp0_7 = -0.5/h1;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,M1-1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,M1-1, M0),NCOMP)];
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,M1-2, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,M1-2, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_01*tmp0_1 + f_11*tmp0_0;
               out[INDEX4(i,1,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_5 + f_01*tmp0_3 + f_10*tmp0_4 + f_11*tmp0_2;
               out[INDEX4(i,0,1,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_01*tmp0_1 + f_11*tmp0_0;
               out[INDEX4(i,1,1,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = tmp0_6*(f_01 + f_11) + tmp0_7*(f_00 + f_10);
               out[INDEX4(i,0,2,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_01*tmp0_1 + f_11*tmp0_0;
               out[INDEX4(i,1,2,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,3)] = f_00*tmp0_11 + f_01*tmp0_9 + f_10*tmp0_10 + f_11*tmp0_8;
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 3 */
   } else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {
      if (face_offset(0)>-1) {
         const double tmp0_3 = -1/h1;
         const double tmp0_2 = 1.0/h1;
         const double tmp0_1 = -0.5/h0;
         const double tmp0_0 = 0.5/h0;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(0,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(0,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_0*(f_10 + f_11) + tmp0_1*(f_00 + f_01);
               out[INDEX4(i,1,0,face_offset(0)+INDEX2(k0,k1,N0),NCOMP,2,1)] = f_00*tmp0_3 + f_01*tmp0_2;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 0 */
      if (face_offset(1)>-1) {
         const double tmp0_3 = 1.0/h1;
         const double tmp0_2 = -1/h1;
         const double tmp0_1 = -0.5/h0;
         const double tmp0_0 = 0.5/h0;
         #pragma omp parallel for private(i,k1)
         for (k1 =0; k1 < N1; ++k1) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(M0-1,k1, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(M0-1,k1+1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(M0-2,k1+1, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(M0-2,k1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_0*(f_10 + f_11) + tmp0_1*(f_00 + f_01);
               out[INDEX4(i,1,0,face_offset(1)+INDEX2(k0,k1,N0),NCOMP,2,1)] = f_10*tmp0_2 + f_11*tmp0_3;
            } /* close component loop i */
         } /* close k1 loop */
      } /* end of face 1 */
      if (face_offset(2)>-1) {
         const double tmp0_3 = 0.5/h1;
         const double tmp0_2 = -0.5/h1;
         const double tmp0_1 = -1/h0;
         const double tmp0_0 = 1.0/h0;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,0, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,0, M0),NCOMP)];
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,1, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,1)] = f_00*tmp0_1 + f_10*tmp0_0;
               out[INDEX4(i,1,0,face_offset(2)+INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_2*(f_00 + f_10) + tmp0_3*(f_01 + f_11);
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 2 */
      if (face_offset(3)>-1) {
         const double tmp0_3 = -0.5/h1;
         const double tmp0_2 = 0.5/h1;
         const double tmp0_1 = -1/h0;
         const double tmp0_0 = 1.0/h0;
         #pragma omp parallel for private(i,k0)
         for (k0 =0; k0 < N0; ++k0) {
            for (i =0; i < NCOMP; ++i) {
               register const double f_11 = in[INDEX2(i,INDEX2(k0+1,M1-1, M0),NCOMP)];
               register const double f_01 = in[INDEX2(i,INDEX2(k0,M1-1, M0),NCOMP)];
               register const double f_10 = in[INDEX2(i,INDEX2(k0+1,M1-2, M0),NCOMP)];
               register const double f_00 = in[INDEX2(i,INDEX2(k0,M1-2, M0),NCOMP)];
               out[INDEX4(i,0,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,1)] = f_01*tmp0_1 + f_11*tmp0_0;
               out[INDEX4(i,1,0,face_offset(3)+INDEX2(k0,k1,N0),NCOMP,2,1)] = tmp0_2*(f_01 + f_11) + tmp0_3*(f_00 + f_10);
            } /* close component loop i */
         } /* close k0 loop */
      } /* end of face 3 */
   } /* end of out_data_type branching
   /* GENERATOR SNIP BOTTOM */
}
