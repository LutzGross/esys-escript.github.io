
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

#ifndef __WEIPA_FILEWRITER_H__
#define __WEIPA_FILEWRITER_H__

#include <weipa/weipa.h>

#include <fstream>

namespace weipa {

class FileWriter
{
public:
    FileWriter();

#if HAVE_MPI
    FileWriter(MPI_Comm comm);
#endif

    bool openFile(std::string filename);
    bool writeOrdered(std::ostringstream& oss);
    bool writeShared(std::ostringstream& oss);
    void close();

private:
    int mpiRank, mpiSize;
#if HAVE_MPI
    MPI_Comm mpiComm;
    MPI_File fileHandle;
#else
    void* mpiComm;
#endif
    std::ofstream ofs;
};

} // namespace weipa

#endif // __WEIPA_FILEWRITER_H__

