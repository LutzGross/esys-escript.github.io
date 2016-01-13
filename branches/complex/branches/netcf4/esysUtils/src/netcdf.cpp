
#include "netcdf.h"

netCDF::NcFile* createNcFile(const char* name, bool readonly)
{
#ifdef NETCDF_CPPV4
    try
    {
#endif
        netCDF::NcFile* res=new NcFile(name, readonly? ENCDF_ROMODE:ENCDF_REPMODE);
#ifndef NETCDF_CPPV4
        if (!res->is_valid())
        {
            return 0;
        }
#endif
        return res;
#ifdef NETCDF_CPPV4
    } catch(netCDF::exceptions::NcException e)
    {
        return 0;
    }
#endif
}
