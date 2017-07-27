
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

#include <fstream>
#ifdef NETCDF4
  #include "NCHelper.h"
#endif

char NcFType(const std::string& name)
{
    // since we don't have a parameter path for specifying file type
    // we'll look at the magic numbers
    std::ifstream f(name.c_str());
    if (!f)
    {
        return '?';
    }
    char buff[5];
    f.read(buff, 4);
    if (!f)
    {
        return false;
    }
    buff[4]=0;
    if (strcmp(buff, "CDF\x01")==0)
    {
        return 'c';
    }
    else if (strcmp(buff, "CDF\x02")==0)
    {
        return 'C';
    }
#ifdef NETCDF4          // if you don't support v4, we won't report it   
    else if (strncmp(buff, "HD5", 3)==0)
    {
        return '4';
    }
#endif    
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
