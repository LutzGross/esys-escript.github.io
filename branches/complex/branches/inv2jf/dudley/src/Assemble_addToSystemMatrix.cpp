
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

/* Dudley: SystemMatrix and SystemVector */

/*  adds the matrix array[Equa,Sol,NN,NN] onto the matrix in. */
/* the rows/columns are given by */
/*  i_Equa+Equa*Nodes_Equa[Nodes[j_Equa]] (i_Equa=0:Equa; j_Equa=0:NN_Equa). */
/*  the routine has to be called from a parallel region                        */

/*  This routine assumes that in->Equa=in->Sol=1, i.e. */
/*  array is fully packed. */
/* TODO: the case in->Equa!=1  */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"

/************************************************************************************/

void Dudley_Assemble_addToSystemMatrix(paso::SystemMatrix_ptr in, const dim_t NN_Equa, const index_t * Nodes_Equa, const dim_t num_Equa,
				       const dim_t NN_Sol, const index_t * Nodes_Sol, const dim_t num_Sol, const double *array)
{
    index_t index_offset = (in->type & MATRIX_FORMAT_OFFSET1 ? 1 : 0);
    dim_t k_Equa, j_Equa, j_Sol, k_Sol, i_Equa, i_Sol, l_col, l_row, ic, ir, k, i_row, i_col;
    index_t *mainBlock_ptr, *mainBlock_index, *col_coupleBlock_ptr, *col_coupleBlock_index, *row_coupleBlock_ptr,
	*row_coupleBlock_index;
    double *mainBlock_val, *row_coupleBlock_val, *col_coupleBlock_val;
    dim_t row_block_size = in->row_block_size;
    dim_t col_block_size = in->col_block_size;
    dim_t block_size = in->block_size;
    dim_t num_subblocks_Equa = num_Equa / row_block_size;
    dim_t num_subblocks_Sol = num_Sol / col_block_size;
    dim_t numMyCols = in->pattern->mainPattern->numInput;
    dim_t numMyRows = in->pattern->mainPattern->numOutput;

    if (in->type & MATRIX_FORMAT_CSC)
    {
	/* MATRIX_FORMAT_CSC does not support MPI !!!!! */
	mainBlock_ptr = in->mainBlock->pattern->ptr;
	mainBlock_index = in->mainBlock->pattern->index;
	mainBlock_val = in->mainBlock->val;
	col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
	col_coupleBlock_index = in->col_coupleBlock->pattern->index;
	col_coupleBlock_val = in->col_coupleBlock->val;
	row_coupleBlock_ptr = in->row_coupleBlock->pattern->ptr;
	row_coupleBlock_index = in->row_coupleBlock->pattern->index;
	row_coupleBlock_val = in->row_coupleBlock->val;

	for (k_Sol = 0; k_Sol < NN_Sol; ++k_Sol)
	{			/* Down columns of array */
	    j_Sol = Nodes_Sol[k_Sol];
	    for (l_col = 0; l_col < num_subblocks_Sol; ++l_col)
	    {
		i_col = j_Sol * num_subblocks_Sol + l_col;
		if (i_col < numMyCols)
		{
		    for (k_Equa = 0; k_Equa < NN_Equa; ++k_Equa)
		    {		/* Across cols of array */
			j_Equa = Nodes_Equa[k_Equa];
			for (l_row = 0; l_row < num_subblocks_Equa; ++l_row)
			{
			    i_row = j_Equa * num_subblocks_Equa + index_offset + l_row;
			    if (i_row < numMyRows + index_offset)
			    {
				for (k = mainBlock_ptr[i_col] - index_offset;
				     k < mainBlock_ptr[i_col + 1] - index_offset; ++k)
				{
				    if (mainBlock_index[k] == i_row)
				    {
					/* Entry array(k_Equa, j_Sol) is a block (col_block_size x col_block_size) */
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						mainBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
					    }
					}
					break;
				    }
				}
			    }
			    else
			    {
				for (k = col_coupleBlock_ptr[i_col] - index_offset;
				     k < col_coupleBlock_ptr[i_col + 1] - index_offset; ++k)
				{
				    if (row_coupleBlock_index[k] == i_row - numMyRows)
				    {
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						row_coupleBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];;
					    }
					}
					break;
				    }
				}
			    }
			}
		    }
		}
		else
		{
		    for (k_Equa = 0; k_Equa < NN_Equa; ++k_Equa)
		    {		/* Across rows of array */
			j_Equa = Nodes_Equa[k_Equa];
			for (l_row = 0; l_row < num_subblocks_Equa; ++l_row)
			{
			    i_row = j_Equa * num_subblocks_Equa + index_offset + l_row;
			    if (i_row < numMyRows + index_offset)
			    {
				for (k = col_coupleBlock_ptr[i_col - numMyCols] - index_offset;
				     k < col_coupleBlock_ptr[i_col - numMyCols + 1] - index_offset; ++k)
				{
				    if (col_coupleBlock_index[k] == i_row)
				    {
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						col_coupleBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
					    }
					}
					break;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    else if (in->type & MATRIX_FORMAT_TRILINOS_CRS)
    {
	/* this needs to be modified */
#ifdef TRILINOS
	for (k_Equa = 0; k_Equa < NN_Equa; ++k_Equa)
	{			/* Down columns of array */
	    j_Equa = Nodes_Equa[k_Equa];
	    if (j_Equa < in->mainBlock->pattern->output_node_distribution->numLocal)
	    {
		for (k_Sol = 0; k_Sol < NN_Sol; ++k_Sol)
		{		/* Across rows of array */
		    j_Sol = Nodes_Sol[k_Sol];
		    for (l_row = 0; l_row < num_subblocks_Equa; ++l_row)
		    {
			irow = j_Equa * row_block_size + l_row;
			for (l_col = 0; l_col < col_block_size; ++l_col)
			{
			    icol = j_Sol * col_block_size + index_offset + l_col;
			    /* irow is local and icol is global */
			    Trilinos_SumIntoMyValues(in->trilinos_data, irow, icol,
						     array[INDEX4
							   (l_row, l_col, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)]);
			}
		    }
		}
	    }
	}
#endif
    }
    else
    {
	mainBlock_ptr = in->mainBlock->pattern->ptr;
	mainBlock_index = in->mainBlock->pattern->index;
	mainBlock_val = in->mainBlock->val;
	col_coupleBlock_ptr = in->col_coupleBlock->pattern->ptr;
	col_coupleBlock_index = in->col_coupleBlock->pattern->index;
	col_coupleBlock_val = in->col_coupleBlock->val;
	row_coupleBlock_ptr = in->row_coupleBlock->pattern->ptr;
	row_coupleBlock_index = in->row_coupleBlock->pattern->index;
	row_coupleBlock_val = in->row_coupleBlock->val;

	for (k_Equa = 0; k_Equa < NN_Equa; ++k_Equa)
	{			/* Down columns of array */
	    j_Equa = Nodes_Equa[k_Equa];
	    for (l_row = 0; l_row < num_subblocks_Equa; ++l_row)
	    {
		i_row = j_Equa * num_subblocks_Equa + l_row;
		/* only look at the matrix rows stored on this processor */
		if (i_row < numMyRows)
		{
		    for (k_Sol = 0; k_Sol < NN_Sol; ++k_Sol)
		    {		/* Across rows of array */
			j_Sol = Nodes_Sol[k_Sol];
			for (l_col = 0; l_col < num_subblocks_Sol; ++l_col)
			{
			    /* only look at the matrix rows stored on this processor */
			    i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
			    if (i_col < numMyCols + index_offset)
			    {
				for (k = mainBlock_ptr[i_row] - index_offset;
				     k < mainBlock_ptr[i_row + 1] - index_offset; ++k)
				{
				    if (mainBlock_index[k] == i_col)
				    {
					/* Entry array(k_Sol, j_Equa) is a block (row_block_size x col_block_size) */
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						mainBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
					    }
					}
					break;
				    }
				}
			    }
			    else
			    {
				for (k = col_coupleBlock_ptr[i_row] - index_offset;
				     k < col_coupleBlock_ptr[i_row + 1] - index_offset; ++k)
				{
				    if (col_coupleBlock_index[k] == i_col - numMyCols)
				    {
					/* Entry array(k_Sol, j_Equa) is a block (row_block_size x col_block_size) */
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						col_coupleBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
					    }
					}
					break;
				    }
				}
			    }
			}
		    }
		}
		else
		{
		    for (k_Sol = 0; k_Sol < NN_Sol; ++k_Sol)
		    {		/* Across rows of array */
			j_Sol = Nodes_Sol[k_Sol];
			for (l_col = 0; l_col < num_subblocks_Sol; ++l_col)
			{
			    i_col = j_Sol * num_subblocks_Sol + index_offset + l_col;
			    if (i_col < numMyCols + index_offset)
			    {
				for (k = row_coupleBlock_ptr[i_row - numMyRows] - index_offset;
				     k < row_coupleBlock_ptr[i_row - numMyRows + 1] - index_offset; ++k)
				{
				    if (row_coupleBlock_index[k] == i_col)
				    {
					/* Entry array(k_Sol, j_Equa) is a block (row_block_size x col_block_size) */
					for (ic = 0; ic < col_block_size; ++ic)
					{
					    i_Sol = ic + col_block_size * l_col;
					    for (ir = 0; ir < row_block_size; ++ir)
					    {
						i_Equa = ir + row_block_size * l_row;
						row_coupleBlock_val[k * block_size + ir + row_block_size * ic] +=
						    array[INDEX4
							  (i_Equa, i_Sol, k_Equa, k_Sol, num_Equa, num_Sol, NN_Equa)];
					    }
					}
					break;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
