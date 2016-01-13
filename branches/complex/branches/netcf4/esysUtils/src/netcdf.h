

// Pull in the netcdf header
#ifdef NETCDF_CPPV4
#include <ncFile.h>
#include <ncInt.h>
#include <ncDim.h>
#include <ncVar.h>
// This still needs to be tested.
// does putAtt return null object on failure?
#define ENCDF_ATTR(X,Y,Z) ((X).putAtt(Y,netCDF::NcInt(),Z)==NcGroupAtt())
#define ENCDF_DIM(X,Y,Z) ((X).addDim(Y,Z))
// Yes dropping the D param is deliberate
#define ENCDF_VAR(A,B,C,D,E) ((A).addVar(B,C,E))
#define ENCDF_VARS(A,B,C,D) ((A).addVar(B,C,D))
// Yes dropping the C param is deliberate
#define ENCDF_PUT(A,B,C,D) try {(A).putVar(B);} catch(netCDF::exceptions::NcException e){throw DataException(D);}
#define ENCDF_REPMODE NcFile::replace
#define ENCDF_ROMODE NcFile::readonly
#define ENCDF_ERRSETUP

#define ENCDF_DIMARR(X,Y) std::vector<NcDim> X(Y);
#define ENCDF_BADDIM(X) ((X)==NcDim())
#define ENCDF_BADVAR(X) ((X)==NcVar())
typedef netCDF::NcVar ENCDF_VART;


#define ENCDF_GATT(X,Y) ((X).getAtt(Y))

typedef netCDF::NcGroupAtt ENCDF_ATTT;
#define ENCDF_BADATT(X) ((X)==NcGroupAtt())


// This macro is not safe if there is more than one value to read
#define ENCDF_GETINT(A,X,Y) ((X).getValues(&A))

#define ENCDF_FREE_ATT(X)

#define ENCDF_DIMCOUNT(X) (X).getDimCount()
#define ENCDF_GETDIM(X,Y) (X).getDim(Y)

typedef netCDF::NcDim ENCDF_DIMTYPE;

#define ENCDF_GETVAR(X,Y) (X).getVar(Y)

// This call is not safe if the buffer is not big enough
// needs to be replaced with a better call before release
#define ENCDF_VARGET(X,Y,Z,Q) try {(X).getVar(Y);} catch (netCDF::exceptions::NcException ne) {throw DataException(Q);}

#define ENCDF_SIZE(X) (X).getSize()

#else
#include <netcdfcpp.h>
namespace netCDF
{
}

#define ENCDF_ATTR(X,Y,Z) ((X).add_attr(Y,Z))
#define ENCDF_DIM(X,Y,Z) ((X).add_dim(Y,Z))
#define ENCDF_VAR(A,B,C,D,E) ((A).add_var(B,C,D,E)
#define ENCDF_VARS(A,B,C,D) ((A).add_var(B,C,D)
#define ENCDF_PUT(A,B,C,D) if (!(A)->put(B,C)) {throw DataException(D);}
#define ENCDF_REPMODE NcFile::Replace
#define ENCDF_ROMODE NcFile::ReadOnly
#define ENCDF_ERRSETUP NcError err(NcError::verbose_nonfatal);

#define ENCDF_DIMARR(X,Y) const NcDim* X[Y];
#define ENCDF_BADDIM(X) (X)
#define ENCDF_BADVAR(X) (X)
typedef netCDF::NcVar* ENCDF_VART;

#define ENCDF_GATT(X,Y) ((X).get_att(Y))

typedef netCDF::NcAtt* ENCDF_ATTT;

#define ENCDF_BADATT(X) (X)
#define ENCDF_GETINT(A,X,Y) (A=(X)->as_int(Y))


#define ENCDF_FREE_ATT(X) delete (X)

#define ENCDF_DIMCOUNT(X) (X).num_dims()
#define ENCDF_GETDIM(X,Y) (X).get_dim(Y)
#define ENCDF_BADDIM(X) (X)

typedef netCDF::NcDim* ENCDF_DIMTYPE;

#define ENCDF_GETVAR(X,Y) (X).get_var(Y)

#define ENCDF_VARGET(X,Y,Z,Q) if (! ((X).get((Y),(Z)))){throw DataException(Q);}

#define ENCDF_SIZE(X) (X)->size()
#endif

netCDF::NcFile* createNcFile(const char* name, bool readonly);
