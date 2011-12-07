
extern "C"
{
#include <paso/SystemMatrixPattern.h>
}

namespace paso
{

// The idea is to create one of these, fill in the fields and then have it generate a pattern  
// NOTE!: unlike finley this will not check to see that you have these in the correct order!
class PatternBuilder
{
  public:
    PatternBuilder(size_t numrows, unsigned char maxneighbours, Esys_MPIInfo* inf);
    ~PatternBuilder();

    double* values;
    int* pos;

    void setDistribution(index_t* ranges, index_t m, index_t b, bool chooseoutput=false);
    
    Paso_SystemMatrixPattern* generatePattern(size_t local_low, size_t local_high, Paso_Connector* c);
    
    
  private:
    size_t numrows;
    unsigned char maxneighbours;
    unsigned char* psums;
    Esys_MPIInfo* mpi_info;
    Paso_Distribution* input_d;
    Paso_Distribution* output_d;
};



PatternBuilder* makePB(Esys_MPIInfo* inf, size_t numrows, unsigned char maxneighbours);
void deallocPB(PatternBuilder* pb);
Paso_Pattern* createPasoPattern(int type, dim_t numOutput, dim_t numInput, index_t numrows, index_t total_elts);

}