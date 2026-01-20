
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_FILEWRITER_H__
#define __ESCRIPT_FILEWRITER_H__

#include <escript/EsysMPI.h>

#include <fstream>
#include <iostream>
#include <sstream>

namespace escript {

class FileWriter
{
public:
    FileWriter(MPI_Comm comm = MPI_COMM_NULL) :
        mpiComm(comm), mpiRank(0), mpiSize(1), m_open(false)
    {
#ifdef ESYS_MPI
        if (comm != MPI_COMM_NULL) {
            MPI_Comm_rank(mpiComm, &mpiRank);
            MPI_Comm_size(mpiComm, &mpiSize);
        }
#endif
    }

    ~FileWriter()
    {
        if (m_open)
            close();
    }

    bool openFile(std::string filename, size_t initialSize=0,
                  bool binary=false, bool append=false)
    {
        // close any open file first
        if (m_open)
            close();

        bool success=false;

        if (mpiSize > 1) {
#ifdef ESYS_MPI
            int mpiErr;
            if (!append) {
                // remove file first if it exists
                int error = 0;
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
                    std::cerr << "Error removing " << filename << "!"
                              << std::endl;
                    return false;
                }
            }

            MPI_Info mpiInfo = MPI_INFO_NULL;
            int amode = MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_UNIQUE_OPEN;
            if (append)
                amode |= MPI_MODE_APPEND;

            mpiErr = MPI_File_open(mpiComm, const_cast<char*>(filename.c_str()),
                    amode, mpiInfo, &fileHandle);
            if (mpiErr == MPI_SUCCESS) {
                mpiErr = MPI_File_set_view(fileHandle, 0, MPI_CHAR, MPI_CHAR,
                        const_cast<char*>("native"), mpiInfo);
            }
            if (mpiErr == MPI_SUCCESS) {
                if (append) {
                    mpiErr = MPI_File_seek_shared(fileHandle, 0, MPI_SEEK_END);
                } else {
                    mpiErr = MPI_File_set_size(fileHandle, initialSize);
                }
            }
            if (mpiErr != MPI_SUCCESS) {
                char errorstr[MPI_MAX_ERROR_STRING];
                int len;
                MPI_Error_string(mpiErr, errorstr, &len);
                std::cerr << "Error opening " << filename
                          << " for parallel writing: " << errorstr << std::endl;
            } else {
                success=true;
            }
#endif
        } else {
            std::ios_base::openmode mode =
                        (binary ? std::ios_base::binary : std::ios_base::out);
            if (append)
                mode |= std::ios_base::app;

            ofs.open(filename.c_str(), mode);
            success = !ofs.fail();
            if (success && initialSize>0 && !append) {
                ofs.seekp(initialSize-1, ofs.beg).put(0).seekp(0, ofs.beg);
                success = !ofs.fail();
            }
        }
        m_open = success;
        return success;
    }

    bool writeOrdered(std::ostringstream& oss)
    {
        if (!m_open)
            return false;

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
        if (!m_open)
            return false;

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
        if (!m_open)
            return false;

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
        if (!m_open)
            return;

        if (mpiSize>1) {
#ifdef ESYS_MPI
            MPI_File_close(&fileHandle);
#endif
        } else {
            ofs.close();
        }
        m_open = false;
    }

private:
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-private-field"
    MPI_Comm mpiComm;
    int mpiRank;
#pragma clang diagnostic pop
    int mpiSize;
    bool m_open;
#ifdef ESYS_MPI
    MPI_File fileHandle;
#endif
    std::ofstream ofs;
};


} // namespace escript

#endif // __ESCRIPT_FILEWRITER_H__

