
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#include <Python.h>
#include <iostream>
#include <stdexcept>

#include <escript/EsysMPI.h>

#ifdef ESYS_MPI

int main(int argc, char **argv)
{
    int status = 0;
    int provided;
    try
    {
        /*
         * Initialise MPI
         */
        /* status = MPI_Init(&argc, &argv); */
        status = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided );
        if (status != MPI_SUCCESS) {
            std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
            return status;
        }
	int rank=0;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)!=MPI_SUCCESS)
	{
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}

        if( rank )
        {
            char fname[256];
            sprintf( fname, "stdout_%04d.out", rank );
            if (freopen( fname, "w+", stdout ) == NULL) {
                std::cerr << "Warning: Unable to redirect stdout." << std::endl;
            }
            sprintf( fname, "stdout_%04d.err", rank );
            if (freopen(fname, "w+", stderr) == NULL) {
                std::cerr << "Warning: Unable to redirect stderr." << std::endl;
            }
        }

        /*
         * Start the python parser
         */
#ifdef ESPYTHON3
        wchar_t** wargv=new wchar_t*[argc+1];
        for (int i=0;i<argc;++i)
        {
            int len=strlen(argv[i]);
            wargv[i]=new wchar_t[len+1];
            for (int j=0;j<len;++j)
            {
                wargv[i][j]=wchar_t(argv[i][j]);
            }
            wargv[i][len]=0;
        }
        wargv[argc]=0;
        status = Py_Main(argc, wargv);
#else
        status = Py_Main(argc, argv);
#endif

        /*
         * Close down MPI.
         * status==1 : uncaught python exception
         * status==2 : invalid python cmd line
         * status>2 : supposed to be param of sys.exit()
         *            sys.exit doesn't return as it should.
         *
         * I have made an exception for 2 because calling MPI_Abort
         * can display pretty ugly messages for not typing params
         * properly.
         * Yes that means you probably shouldn't use 2 as an exit code
         * but we don't really recommend sys.exit anyway.
         */
        if ((status!=0) && (status!=2))
        {
            MPI_Abort(MPI_COMM_WORLD,status);
        }  
        else
        {
            /*
             * Finalise MPI for a clean exit.
             */
            MPI_Finalize();
        }
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "EXCEPTION: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        throw;
    }
    catch (char *e)
    {
        std::cerr << "EXCEPTION: " << e << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        throw;
    }
    catch (...)
    {
        std::cerr << "EXCEPTION: " << "UNKNOWN." << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        throw;
    }

    return status;
}

#else

int main(int argc, char **argv)
{
    std::cout << "Escript must be compiled with ESYS_MPI defined to make the MPI version available" << std::endl;
    return 0;
}

#endif // ESYS_MPI

