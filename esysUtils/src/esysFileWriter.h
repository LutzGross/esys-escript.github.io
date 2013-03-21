
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

#ifndef __ESYS_FILEWRITER_H__
#define __ESYS_FILEWRITER_H__

#include <fstream>
#include <iostream>
#include <sstream>

#ifdef ESYS_MPI
#include <mpi.h>
#endif

namespace esysUtils {

class FileWriter
{
public:
    FileWriter() : mpiRank(0), mpiSize(1) {}

#ifdef ESYS_MPI
    FileWriter(MPI_Comm comm) : mpiComm(comm)
    {
        MPI_Comm_rank(mpiComm, &mpiRank);
        MPI_Comm_size(mpiComm, &mpiSize);
    }
#endif

    bool openFile(std::string filename, size_t initialSize=0)
    {
        bool success=false;

        if (mpiSize>1) {
#ifdef ESYS_MPI
            // remove file first if it exists
            int error = 0;
            int mpiErr;
            if (mpiRank == 0) {
                std::ifstream f(filename.c_str());
                if (f.is_open()) {
                    f.close();
                    if (std::remove(filename.c_str())) {
                        error=1;
                    }
                }
            }
            MPI_Allreduce(&error, &mpiErr, 1, MPI_INT, MPI_MAX, mpiComm);
            if (mpiErr != 0) {
                std::cerr << "Error removing " << filename << "!" << std::endl;
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
            if (mpiErr == MPI_SUCCESS) {
                mpiErr = MPI_File_set_size(fileHandle, initialSize);
            }
            if (mpiErr != MPI_SUCCESS) {
                std::cerr << "Error opening " << filename << " for parallel writing!" << std::endl;
            } else {
                success=true;
            }
#endif
        } else {
            std::ios_base::openmode mode = std::ios_base::binary;
            ofs.open(filename.c_str(), mode);
            success = !ofs.fail();
            if (success && initialSize>0) {
                ofs.seekp(initialSize-1, ofs.beg).put(0).seekp(0, ofs.beg);
                success = !ofs.fail();
            }
        }
        return success;
    }

    bool writeOrdered(std::ostringstream& oss)
    {
        bool success=false;
        if (mpiSize>1) {
#ifdef ESYS_MPI
            MPI_Status mpiStatus;
            std::string contents = oss.str();
            int mpiErr = MPI_File_write_ordered(
                fileHandle, const_cast<char*>(contents.c_str()),
                contents.length(), MPI_CHAR, &mpiStatus);
            oss.str(std::string());
            success=(mpiErr==0);
#endif
        } else {
            ofs << oss.str();
            oss.str(std::string());
            success=!ofs.fail();
        }
        return success;
    }

    bool writeShared(std::ostringstream& oss)
    {
        bool success=false;
        if (mpiSize>1) {
#ifdef ESYS_MPI
            MPI_Status mpiStatus;
            std::string contents = oss.str();
            int mpiErr = MPI_File_write_shared(
                fileHandle, const_cast<char*>(contents.c_str()),
                contents.length(), MPI_CHAR, &mpiStatus);
            oss.str(std::string());
            success=(mpiErr==0);
#endif
        } else {
            ofs << oss.str();
            oss.str(std::string());
            success=!ofs.fail();
        }
        return success;
    }

    bool writeAt(std::ostringstream& oss, long offset)
    {
        bool success=false;
        if (mpiSize>1) {
#ifdef ESYS_MPI
            MPI_Status mpiStatus;
            std::string contents = oss.str();
            int mpiErr = MPI_File_write_at(
                fileHandle, offset, const_cast<char*>(contents.c_str()),
                contents.length(), MPI_CHAR, &mpiStatus);
            oss.str(std::string());
            success=(mpiErr==0);
#endif
        } else {
            ofs.seekp(offset);
            ofs << oss.str();
            oss.str(std::string());
            success=!ofs.fail();
        }
        return success;
    }

    void close()
    {
        if (mpiSize>1) {
#ifdef ESYS_MPI
            MPI_File_close(&fileHandle);
#endif
        } else {
            ofs.close();
        }
    }

private:
    int mpiRank, mpiSize;
#ifdef ESYS_MPI
    MPI_Comm mpiComm;
    MPI_File fileHandle;
#else
    void* mpiComm;
#endif
    std::ofstream ofs;
};


} // namespace esysUtils

#endif //  __ESYS_FILEWRITER_H__

