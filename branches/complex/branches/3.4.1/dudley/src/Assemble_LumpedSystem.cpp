
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

/************************************************************************************/

/*    assembles the mass matrix in lumped form                */

/*    The coefficient D has to be defined on the integration points or not present. */

/*    lumpedMat has to be initialized before the routine is called. */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ShapeTable.h"

/************************************************************************************/

void Dudley_Assemble_LumpedSystem(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escriptDataC * lumpedMat,
				  escriptDataC * D, const bool useHRZ)
{

    bool reducedIntegrationOrder = FALSE, expandedD;
    char error_msg[LenErrorMsg_MAX];
    Dudley_Assemble_Parameters p;
    dim_t dimensions[ESCRIPT_MAX_DATA_RANK], k, e, len_EM_lumpedMat, q, s;
    type_t funcspace;
    index_t color, *row_index = NULL;
    __const double *D_p = NULL;
    const double *S = NULL;
    double *EM_lumpedMat = NULL, *lumpedMat_p = NULL;
    register double rtmp;
    register double m_t = 0., diagS = 0.;


    Dudley_resetError();

    if (nodes == NULL || elements == NULL)
	return;
    if (isEmpty(lumpedMat) || isEmpty(D))
	return;
    if (isEmpty(lumpedMat) && !isEmpty(D))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_LumpedSystem: coefficients are non-zero but no lumped matrix is given.");
	return;
    }
    funcspace = getFunctionSpaceType(D);
    /* check if all function spaces are the same */
    if (funcspace == DUDLEY_ELEMENTS)
    {
	reducedIntegrationOrder = FALSE;
    }
    else if (funcspace == DUDLEY_FACE_ELEMENTS)
    {
	reducedIntegrationOrder = FALSE;
    }
    else if (funcspace == DUDLEY_REDUCED_ELEMENTS)
    {
	reducedIntegrationOrder = TRUE;
    }
    else if (funcspace == DUDLEY_REDUCED_FACE_ELEMENTS)
    {
	reducedIntegrationOrder = TRUE;
    }
    else
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_LumpedSystem: assemblage failed because of illegal function space.");
    }
    if (!Dudley_noError())
	return;

    /* set all parameters in p */
    Dudley_Assemble_getAssembleParameters(nodes, elements, NULL, lumpedMat, reducedIntegrationOrder, &p);
    if (!Dudley_noError())
	return;

    /* check if all function spaces are the same */
    if (!numSamplesEqual(D, p.numQuad, elements->numElements))
    {
	sprintf(error_msg, "Dudley_Assemble_LumpedSystem: sample points of coefficient D don't match (%d,%d)", p.numQuad,
		elements->numElements);
	Dudley_setError(TYPE_ERROR, error_msg);
    }

    /*  check the dimensions: */
    if (p.numEqu == 1)
    {
	if (!isEmpty(D))
	{
	    if (!isDataPointShapeEqual(D, 0, dimensions))
	    {
		Dudley_setError(TYPE_ERROR, "Dudley_Assemble_LumpedSystem: coefficient D, rank 0 expected.");
	    }

	}
    }
    else
    {
	if (!isEmpty(D))
	{
	    dimensions[0] = p.numEqu;
	    if (!isDataPointShapeEqual(D, 1, dimensions))
	    {
		sprintf(error_msg, "Dudley_Assemble_LumpedSystem: coefficient D, expected shape (%d,)", dimensions[0]);
		Dudley_setError(TYPE_ERROR, error_msg);
	    }
	}
    }
    if (Dudley_noError())
    {
	requireWrite(lumpedMat);
	lumpedMat_p = getSampleDataRW(lumpedMat, 0);
	
	if (funcspace==DUDLEY_POINTS) {
              #pragma omp parallel private(color, D_p)
	      {
		    for (color=elements->minColor;color<=elements->maxColor;color++) {
		      /*  open loop over all elements: */
		      #pragma omp for private(e) schedule(static)
		      for(e=0;e<elements->numElements;e++){
			  if (elements->Color[e]==color) {
			    D_p=getSampleDataRO(D, e);
			    if (NULL!=D_p)  Dudley_Util_AddScatter(1,
                                                        &(p.row_DOF[elements->Nodes[INDEX2(0,e,p.NN)]]),
                                                        p.numEqu,
                                                        D_p,
                                                        lumpedMat_p, 
                                                        p.row_DOF_UpperBound);
			  } /* end color check */
		      } /* end element loop */
		  } /* end color loop */
	    } /* end parallel region */
        } else {  
	      
	      len_EM_lumpedMat = p.numShapes * p.numEqu;

	      expandedD = isExpanded(D);
	      if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &S))
	      {
		  Dudley_setError(TYPE_ERROR, "Dudley_Assemble_LumpedSystem: Unable to locate shape function.");
	      }
	      #pragma omp parallel private(color, EM_lumpedMat, row_index, D_p, s, q, k, rtmp, diagS, m_t)
	      {
		  EM_lumpedMat = new double[len_EM_lumpedMat];
		  row_index = new index_t[p.numShapes];
		  if (!Dudley_checkPtr(EM_lumpedMat) && !Dudley_checkPtr(row_index))
		  {
		      if (p.numEqu == 1)
		      {		/* single equation */
			  if (expandedD)
			  {		/* with expanded D */
			      for (color = elements->minColor; color <= elements->maxColor; color++)
			      {
				  /*  open loop over all elements: */
      #pragma omp for private(e) schedule(static)
				  for (e = 0; e < elements->numElements; e++)
				  {
				      if (elements->Color[e] == color)
				      {
					  double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
					  D_p = getSampleDataRO(D, e);
					  if (useHRZ)	{
					    m_t = 0;	/* mass of the element: m_t */
					    for (q = 0; q < p.numQuad; q++)
						m_t += vol * D_p[INDEX2(q, 0, p.numQuad)];
					    diagS = 0;	/* diagonal sum: S */
					    for (s = 0; s < p.numShapes; s++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						  rtmp +=
							vol * D_p[INDEX2(q, 0, p.numQuad)] * S[INDEX2(s, q, p.numShapes)] *
							S[INDEX2(s, q, p.numShapes)];
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
						diagS += rtmp;
					    }
					    /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					    rtmp = m_t / diagS;
					    for (s = 0; s < p.numShapes; s++)
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] *= rtmp;

					  } else {/* row-sum lumping */
					    for (s = 0; s < p.numShapes; s++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						  rtmp += vol * S[INDEX2(s, q, p.numShapes)] * D_p[INDEX2(q, 0, p.numQuad)];
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
					    }
					  }
					  for (q = 0; q < p.numShapes; q++)
					  {
					      row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
					  }
					  Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_lumpedMat, lumpedMat_p,
								p.row_DOF_UpperBound);
				      }	/* end color check */
				  }	/* end element loop */
			      }	/* end color loop */
			  }
			  else
			  {		/* with constant D */

			      for (color = elements->minColor; color <= elements->maxColor; color++)
			      {
				  /*  open loop over all elements: */
      #pragma omp for private(e) schedule(static)
				  for (e = 0; e < elements->numElements; e++)
				  {
				      if (elements->Color[e] == color)
				      {
					  double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
					  D_p = getSampleDataRO(D, e);
					  if (useHRZ)	{	/* HRZ lumping */
					    m_t = 0;	/* mass of the element: m_t */
					    for (q = 0; q < p.numQuad; q++)
						m_t += vol;
					    diagS = 0;	/* diagonal sum: S */
					    for (s = 0; s < p.numShapes; s++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						{
						  rtmp += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
						}
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp;
						diagS += rtmp;
					    }
					    /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					    rtmp = m_t / diagS * D_p[0];
					    for (s = 0; s < p.numShapes; s++)
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] *= rtmp;
					  } else {			/* row-sum lumping */
					    for (s = 0; s < p.numShapes; s++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						  rtmp += vol * S[INDEX2(s, q, p.numShapes)];
						EM_lumpedMat[INDEX2(0, s, p.numEqu)] = rtmp * D_p[0];
					    }
					  }
					  for (q = 0; q < p.numShapes; q++)
					      row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
					  Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_lumpedMat, lumpedMat_p,
								p.row_DOF_UpperBound);
				      }	/* end color check */
				  }	/* end element loop */
			      }	/* end color loop */

			  }
		      }
		      else
		      {		/* system of  equation */
			  if (expandedD)
			  {		/* with expanded D */
			      for (color = elements->minColor; color <= elements->maxColor; color++)
			      {
				  /*  open loop over all elements: */
      #pragma omp for private(e) schedule(static)
				  for (e = 0; e < elements->numElements; e++)
				  {
				      if (elements->Color[e] == color)
				      {
					  double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
					  D_p = getSampleDataRO(D, e);

					  if (useHRZ)	{	/* HRZ lumping */
					    for (k = 0; k < p.numEqu; k++)
					    {
						m_t = 0;	/* mass of the element: m_t */
						for (q = 0; q < p.numQuad; q++)
						  m_t += vol * D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)];

						diagS = 0;	/* diagonal sum: S */
						for (s = 0; s < p.numShapes; s++)
						{
						  rtmp = 0;
						  for (q = 0; q < p.numQuad; q++)
							rtmp +=
							    vol * D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)] *
							    S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
						  diagS += rtmp;
						}
						/* rescale diagonals by m_t/diagS to ensure consistent mass over element */
						rtmp = m_t / diagS;
						for (s = 0; s < p.numShapes; s++)
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] *= rtmp;
					    }
					  } else {				/* row-sum lumping */
					    for (s = 0; s < p.numShapes; s++)
					    {
						for (k = 0; k < p.numEqu; k++)
						{
						  rtmp = 0.;
						  for (q = 0; q < p.numQuad; q++)
							rtmp +=
							    vol * S[INDEX2(s, q, p.numShapes)] *
							    D_p[INDEX3(k, q, 0, p.numEqu, p.numQuad)];
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
						}
					    }
					  }
					  for (q = 0; q < p.numShapes; q++)
					      row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
					  Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_lumpedMat, lumpedMat_p,
								p.row_DOF_UpperBound);
				      }	/* end color check */
				  }	/* end element loop */
			      }	/* end color loop */
			  }
			  else
			  {		/* with constant D */
			      for (color = elements->minColor; color <= elements->maxColor; color++)
			      {
				  /*  open loop over all elements: */
      #pragma omp for private(e) schedule(static)
				  for (e = 0; e < elements->numElements; e++)
				  {
				      if (elements->Color[e] == color)
				      {
					  double vol = p.row_jac->absD[e] * p.row_jac->quadweight;
					  D_p = getSampleDataRO(D, e);

					  if (useHRZ)		{ /* HRZ lumping */
					    m_t = 0;	/* mass of the element: m_t */
					    for (q = 0; q < p.numQuad; q++)
						m_t += vol;
					    diagS = 0;	/* diagonal sum: S */
					    for (s = 0; s < p.numShapes; s++)
					    {
						rtmp = 0;
						for (q = 0; q < p.numQuad; q++)
						  rtmp += vol * S[INDEX2(s, q, p.numShapes)] * S[INDEX2(s, q, p.numShapes)];
						for (k = 0; k < p.numEqu; k++)
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp;
						diagS += rtmp;
					    }
					    /* rescale diagonals by m_t/diagS to ensure consistent mass over element */
					    rtmp = m_t / diagS;
					    for (s = 0; s < p.numShapes; s++)
					    {
						for (k = 0; k < p.numEqu; k++)
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] *= rtmp * D_p[k];
					    }
					  } else {				/* row-sum lumping */
					    for (s = 0; s < p.numShapes; s++)
					    {
						for (k = 0; k < p.numEqu; k++)
						{
						  rtmp = 0.;
						  for (q = 0; q < p.numQuad; q++)
						      rtmp += vol * S[INDEX2(s, q, p.numShapes)];
						  EM_lumpedMat[INDEX2(k, s, p.numEqu)] = rtmp * D_p[k];
						}
					    }
					  }
					  for (q = 0; q < p.numShapes; q++)
					      row_index[q] = p.row_DOF[elements->Nodes[INDEX2(q, e, p.NN)]];
					  Dudley_Util_AddScatter(p.numShapes, row_index, p.numEqu, EM_lumpedMat, lumpedMat_p,
								p.row_DOF_UpperBound);
				      }	/* end color check */
				  }	/* end element loop */
			      }	/* end color loop */
			  }
		      }
		  }			/* end of pointer check */
		  delete[] EM_lumpedMat;
		  delete[] row_index;
	      }			/* end parallel region */
	}
    }
}
