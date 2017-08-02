
/*****************************************************************************
*
* Copyright (c) 2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <cstring>
#include <fstream>
#include <cstring>
#ifdef NETCDF4
  #include "NCHelper.h"
#endif


namespace escript
{
char NcFType(const std::string& name)
{
    // since we don't have a parameter path for specifying file type
    // we'll look at the magic numbers
    std::ifstream f(name.c_str());
    if (!f)
    {
        return '?';
    }
    char buff[10];
    f.read(buff, 9);
    if (!f)
    {
        return '?';
    }
    buff[9]=0;
    if (strncmp(buff, "CDF\x01",4)==0)
    {
        return 'c';
    }
    else if (strncmp(buff, "CDF\x02",4)==0)
    {
        return 'C';
    }
    else if (strncmp(buff, "\x89HDF\r\n\x1a\n", 8)==0)
    {        
#ifdef NETCDF4                   
        return '4';
#endif    
        return 'u';
    }
    else
    {
        return '?';
    }
}


#ifdef NETCDF4

using namespace netCDF;

bool openNcFile(netCDF::NcFile& ncf, const std::string& name)
{
    try
    {
        char type=NcFType(name);
        NcFile::FileFormat fm=NcFile::FileFormat::classic;        
        switch (type)
        {
            case 'c':    fm=NcFile::FileFormat::classic; break;
            case 'C':    fm=NcFile::FileFormat::classic64; break;
            case '4':    fm=NcFile::FileFormat::nc4; break;    // using this rather than nc4classic since we don't intend to write into the file
            default:
                return false;
        }
        ncf.open(name.c_str(), NcFile::FileMode::read, fm);
    }
    catch (exceptions::NcException e)
    {
        return false;
    }
    return true;
}



#endif
}
