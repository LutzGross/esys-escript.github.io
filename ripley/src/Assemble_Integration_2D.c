Ripley_Assemble_Integration_2D(Ripley_Grid *grid, Escript in, double *out)
{
   /* GENERATOR SNIP TOP */
   if (out_data_type==RIPLEY_ELEMENTS) {
      const double w_0 = h0*h1/4;
      for (k1_0 = 0; k1_0 <2; k1_0++) { /* coloring */
          #pragma omp parallel for private(k0,k1)
          for (k1 = k1_0; k1< N1; k1=k1+2) {
                for (k0 = 0; k0< N0; ++k0)  {
                      const index_t e = k0 + M0 * k1;
                      for (i =0; i < NCOMP; ++i) {
			      const register double f_0 = in[INDEX3(i,0,e, NCOMP,4)];
                              const register double f_1 = in[INDEX3(i,1,e, NCOMP,4)];
                              const register double f_2 = in[INDEX3(i,2,e, NCOMP,4)];
                              const register double f_3 = in[INDEX3(i,3,e, NCOMP,4)];
			      out[i]+=(f_0+f_1+f_2+f_3)*w_0;
                      }  /* close component loop i */
		} /* close k0 loop */
	  } /* close k1 loop */
      } /* close coloring k1 loop */
   } else if (out_data_type==RIPLEY_REDUCED_ELEMENTS) {
      const double w_0 = h0*h1;
      for (k1_0 = 0; k1_0 <2; k1_0++) { /* coloring */
          #pragma omp parallel for private(k0,k1)
          for (k1 = k1_0; k1< N1; k1=k1+2) {
                for (k0 = 0; k0< N0; ++k0)  {
                      const index_t e = k0 + M0 * k1;
                      for (i =0; i < NCOMP; ++i) {
			      const register double f_0 = in[INDEX3(i,0,e, NCOMP,1)];
			      out[i]+=f_0*w_0;
                      }  /* close component loop i */
		} /* close k0 loop */
	  } /* close k1 loop */
      } /* close coloring k1 loop */
   } else if (out_data_type==RIPLEY_BOUNDARY_ELEMENTS) {
      #pragma omp parallel private(k0, k0_0, k1, k1_0)
      { 
	 const double w_0 = h0/2.;
	 const double w_1 = h1/2.;
         if (face_offset(0)>-1) {
	     for (k0_0=0; k0_0<2; k0_0++) {
                #pragma omp for
	        for (k0 = k0_0; k0< N0; k0=k0+2) {
                    const index_t e = k0 + face_offset(0);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,2)];
                        const register double f_1 = in[INDEX3(i,1,e, NCOMP,2)];
       	                out[i]+=(f_0+f_1)*w_0;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 } 
         if (face_offset(1)>-1) {
	    for (k1_0=0; k1_0<2; k1_0++) {
                #pragma omp for
	        for (k1 = k1_0; k1< N1; k1=k1+2) {
                    const index_t e = k1 + face_offset(1);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,2)];
                        const register double f_1 = in[INDEX3(i,1,e, NCOMP,2)];
       	                out[i]+=(f_0+f_1)*w_1;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 } 
         if (face_offset(2)>-1) {
	    for (k0_0=0; k0_0<2; k0_0++) {
                #pragma omp for
	        for (k0 = k0_0; k0< N0; k0=k0+2) {
                    const index_t e = k0 + face_offset(2);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,2)];
                        const register double f_1 = in[INDEX3(i,1,e, NCOMP,2)];
       	                out[i]+=(f_0+f_1)*w_0;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 }
         if (face_offset(3)>-1) {
	    for (k1_0=0; k1_0<2; k1_0++) {
                #pragma omp for
	        for (k1 = k1_0; k1< N1; k1=k1+2) {
                    const index_t e = k1 + face_offset(3);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,2)];
                        const register double f_1 = in[INDEX3(i,1,e, NCOMP,2)];
       	                out[i]+=(f_0+f_1)*w_0;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 }
      } /* end parallel region */
   } else if (out_data_type==RIPLEY_REDUCED_BOUNDARY_ELEMENTS) {
      #pragma omp parallel private(k0, k0_0, k1, k1_0)
      { 
         if (face_offset(0)>-1) {
	     for (k0_0=0; k0_0<2; k0_0++) {
                #pragma omp for
	        for (k0 = k0_0; k0< N0; k0=k0+2) {
                    const index_t e = k0 + face_offset(0);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,1)];
       	                out[i]+=f_0*h_0;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 } 
         if (face_offset(1)>-1) {
	    for (k1_0=0; k1_0<2; k1_0++) {
                #pragma omp for
	        for (k1 = k1_0; k1< N1; k1=k1+2) {
                    const index_t e = k1 + face_offset(1);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,1)];
       	                out[i]+=f_0*h1;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 } 
         if (face_offset(2)>-1) {
	    for (k0_0=0; k0_0<2; k0_0++) {
                #pragma omp for
	        for (k0 = k0_0; k0< N0; k0=k0+2) {
                    const index_t e = k0 + face_offset(2);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,1)];
       	                out[i]+=f_0*h0;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 }
         if (face_offset(3)>-1) {
	    for (k1_0=0; k1_0<2; k1_0++) {
                #pragma omp for
	        for (k1 = k1_0; k1< N1; k1=k1+2) {
                    const index_t e = k1 + face_offset(3);
                    for (i =0; i < NCOMP; ++i) {
	                const register double f_0 = in[INDEX3(i,0,e, NCOMP,1)];
       	                out[i]+=f_0*w_1;
                    }  /* close component loop i */
                } /* close k0 loop */
	     }
	 }
      } /* end parallel region */
   } /* end of out_data_type branching */
   /* GENERATOR SNIP BOTTOM */
}
