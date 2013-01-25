

#include "PatternBuilder.h"

extern "C" {
#include <paso/Pattern.h>
}

namespace paso
{

PatternBuilder* makePB(Esys_MPIInfo* inf, size_t numrows, unsigned char maxneighbours)
{
    return new PatternBuilder(numrows, maxneighbours, inf);
}

void deallocPB(PatternBuilder* pb)
{
    delete pb;
}

Paso_Pattern* createPasoPattern(int type, dim_t numOutput, dim_t numInput, index_t numrows, index_t total_elts)
{
    Paso_Pattern* out=(Paso_Pattern*)malloc(sizeof(Paso_Pattern));
    out->type=type;
    out->reference_counter=1;
    out->numOutput=numOutput;
    out->numInput=numInput;
    out->ptr=new index_t[numrows];
    out->index=new index_t[total_elts];
    out->main_iptr = NULL;
    out->coloring = NULL;
    out->numColors=-1;
    out->len=total_elts;
    return out;
}

PatternBuilder::PatternBuilder(size_t numrows, unsigned char maxneighbours, Esys_MPIInfo* info)
{
    this->numrows=numrows;
    this->maxneighbours=maxneighbours;
    values=new double[numrows*maxneighbours];
    pos=new int[numrows*maxneighbours];
    psums=new unsigned char[numrows*2];		// one total for local links and one for remote (coupled) links
    int i;
    mpi_info=Esys_MPIInfo_getReference(info);
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numrows;++i)
    {
        for (unsigned char j=0;j<maxneighbours;++j)
	{
	    values[i*maxneighbours+j]=0;
	    pos[i*maxneighbours+j]=-1;
	}
	psums[2*i]=0;
	psums[2*i+1]=0;
    }
    input_d=output_d=0;
}

PatternBuilder::~PatternBuilder()
{
    delete[] values;  
    delete[] pos;
    delete[] psums;
    Esys_MPIInfo_free(mpi_info);
    Paso_Distribution_free(input_d);
    Paso_Distribution_free(output_d);
}


void PatternBuilder::setDistribution(index_t* ranges, index_t m, index_t b, bool chooseoutput)
{
    Paso_Distribution* di=Paso_Distribution_alloc(mpi_info, ranges, m, b);
    if (chooseoutput)
    {
        output_d=di;
    }
    else
    {
        input_d=di;
    }
}

//low and high are inclusive
Paso_SystemMatrixPattern* PatternBuilder::generatePattern(size_t local_low, size_t local_high, Paso_Connector* c)
{
    const int mtype=MATRIX_FORMAT_DEFAULT;	// hard coding this for now
  
  
    // compute totals
    int i;
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numrows;++i)
    {
        for (unsigned char j=0;j<maxneighbours;++j)
	{
	    if (pos[i*maxneighbours+j]!=-1)
	    {
	        if ((pos[i*maxneighbours+j]>local_high) || (pos[i*maxneighbours+j]<local_low))
		{
		    psums[2*i+1]++;	// non-local
		}
		else
		{
		    psums[2*i]++;	// local
		}
	    }
	}
    }
    // now we do partial sums
    // not parallelising this yet
    for (i=1;i<numrows;++i)
    {
        psums[i]+=psums[i-1];
    }
    // now we have enough info to populate a pattern
    Paso_Pattern* mainpat=createPasoPattern(mtype, numrows, numrows, numrows, psums[numrows]);
    
    // These patterns are dummies right now.
    Paso_Pattern* colpat=createPasoPattern(mtype, numrows, numrows, numrows, psums[numrows]);
    Paso_Pattern* rowpat=createPasoPattern(mtype, numrows, numrows, numrows, psums[numrows]);

    mainpat->ptr[0]=0;
    rowpat->ptr[0]=0;
    colpat->ptr[0]=0;

    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<numrows;++i)
    {
        mainpat->ptr[i+1]=psums[2*i];
	rowpat->ptr[i+1]=psums[2*i+1];
	
	// dodgy hack - colpat should not be a clone of rowpat
	colpat->ptr[i+1]=psums[2*i+1];
	
	int l=0, r=0;	//local and remote counters
	for (int j=0;j<maxneighbours;++j)
	{
	    if (pos[i*maxneighbours+j]!=-1)
	    {
	        if ((pos[i*maxneighbours+j]>local_high) || (pos[i*maxneighbours+j]<local_low))
		{
		    rowpat->index[psums[2*i+1]+r]=pos[i*maxneighbours+j];	// non-local
		    
		    // dodgy hack - colpat should not be a clone of rowpat
		    colpat->index[psums[2*i+1]+r]=pos[i*maxneighbours+j];	// non-local
		    
		    r++;
		}
		else
		{
		    mainpat->index[psums[2*i]+l]=pos[i*maxneighbours+j];	// local
		    l++;
		}
	    }	  
	}
    }
    
    // at the moment we are only considering one rank
    // so what are m and b
    //Paso_Distribution_alloc( Esys_MPIInfo *mpi_info, index_t* first_component, index_t m, index_t b);


    Paso_Connector* row_connector=c;
    Paso_Connector* col_connector=Paso_Connector_copy(c);
    
    // the type field must match in all the pieces but we are creating all of them now so that should be ok
    
    return Paso_SystemMatrixPattern_alloc( mtype, output_d, input_d, mainpat, colpat, rowpat, col_connector, row_connector);    
}
}
