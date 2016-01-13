
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <weipa/FileWriter.h>
#include <sstream>
#include <iostream>


using namespace std;

namespace weipa {

FileWriter::FileWriter() :
    mpiRank(0),
    mpiSize(1)
{
}

#if HAVE_MPI
FileWriter::FileWriter(MPI_Comm comm) :
    mpiComm(comm)
{
    MPI_Comm_rank(mpiComm, &mpiRank);
    MPI_Comm_size(mpiComm, &mpiSize);
}
#endif

bool FileWriter::openFile(string filename)
{
    bool success=false;

    if (mpiSize>1) {
#if HAVE_MPI
        // remove file first if it exists
        int error = 0;
        int mpiErr;
        if (mpiRank == 0) {
            ifstream f(filename.c_str());
            if (f.is_open()) {
                f.close();
                if (remove(filename.c_str())) {
                    error=1;
                }
            }
        }
        MPI_Allreduce(&error, &mpiErr, 1, MPI_INT, MPI_MAX, mpiComm);
        if (mpiErr != 0) {
            cerr << "Error removing " << filename << "!" << endl;
            return false;
        }

        MPI_Info mpiInfo = MPI_INFO_NULL;
        int amode = MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_UNIQUE_OPEN;
        mpiErr = MPI_File_open(mpiComm, const_cast<char*>(filename.c_str()),
                amode, mpiInfo, &fileHandle);
        if (mpiErr == MPI_SUCCESS) {
            mpiErr = MPI_File_set_view(fileHandle, 0, MPI_CHAR, MPI_CHAR,
                    const_cast<char*>("native"), mpiInfo);
        }
        if (mpiErr != MPI_SUCCESS) {
            cerr << "Error opening " << filename << " for parallel writing!" << endl;
        } else {
            success=true;
        }
#endif
    } else {
        ofs.open(filename.c_str());
        success = !ofs.fail();
    }
    return success;
}

bool FileWriter::writeShared(ostringstream& oss)
{
    bool success=false;
    if (mpiSize>1) {
#if HAVE_MPI
        MPI_Status mpiStatus;
        string contents = oss.str();
        int mpiErr = MPI_File_write_shared(
            fileHandle, const_cast<char*>(contents.c_str()),
            contents.length(), MPI_CHAR, &mpiStatus);
        oss.str("");
        success=(mpiErr==0);
#endif
    } else {
        ofs << oss.str();
        oss.str("");
        success=!ofs.fail();
    }
    return success;
}

bool FileWriter::writeOrdered(ostringstream& oss)
{
    bool success=false;
    if (mpiSize>1) {
#if HAVE_MPI
        MPI_Status mpiStatus;
        string contents = oss.str();
        int mpiErr = MPI_File_write_ordered(
            fileHandle, const_cast<char*>(contents.c_str()),
            contents.length(), MPI_CHAR, &mpiStatus);
        oss.str("");
        success=(mpiErr==0);
#endif
    } else {
        ofs << oss.str();
        oss.str("");
        success=!ofs.fail();
    }
    return success;
}

void FileWriter::close()
{
    if (mpiSize>1) {
#if HAVE_MPI
        MPI_File_close(&fileHandle);
#endif
    } else {
        ofs.close();
    }
}

} // namespace weipa

