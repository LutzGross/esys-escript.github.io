
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

#ifdef NETCDF4
#include "NCHelper.h"
#include <fstream>

using namespace netCDF;

bool openNcFile(netCDF::NcFile& ncf, const std::string& name)
{
    try
    {
        // since we don't have a parameter path for specifying file type
        // we'll look at the magic numbers
        std::ifstream f(name.c_str());
        if (!f)
        {
            return false;
        }
        char buff[5];
        f.read(buff, 4);
        if (!f)
        {
            return false;
        }
        buff[4]=0;
        NcFile::FileFormat fm=NcFile::FileFormat::classic;
        if (strcmp(buff, "CDF\x01")==0)
        {
            fm=NcFile::FileFormat::classic;
        }
        else if (strcmp(buff, "CDF\x02")==0)
        {
            fm=NcFile::FileFormat::classic64;
        }
        else if (strncmp(buff, "HD5", 3)==0)
        {
            fm=NcFile::FileFormat::nc4;     // using this rather than nc4classic since we don't intend to write into the file
        }
        else
        {
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
