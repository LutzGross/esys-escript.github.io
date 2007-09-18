
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifdef PASO_MPI
#include <mpi.h>
#endif
#include "MeshAdapterFactory.h"
#include "FinleyError.h"

#include <boost/python/extract.hpp>

#include <sstream>

using namespace std;
using namespace escript;

namespace finley {

  AbstractContinuousDomain* loadMesh(const std::string& fileName) 
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

    fMesh=Finley_Mesh_load(fName);
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    return temp;
  }

  AbstractContinuousDomain* readMesh(const std::string& fileName,
  				     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

    fMesh=Finley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    return temp;
  }

  AbstractContinuousDomain* readGmsh(const std::string& fileName,
                                     int numDim,
                                     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     int optimize)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

    fMesh=Finley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimize ? TRUE : FALSE));
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    
    /* win32 refactor */
    TMPMEMFREE(fName);
    
    return temp;
  }

  AbstractContinuousDomain* brick(int n0,int n1,int n2,int order,
		    double l0,double l1,double l2,
		    int periodic0,int periodic1,
		    int periodic2,
		    int integrationOrder,
                    int reducedIntegrationOrder,
		    int useElementsOnFace,
		    int useFullElementOrder,
                    int optimize)
  {
    int numElements[]={n0,n1,n2};
    double length[]={l0,l1,l2};
    int periodic[]={periodic0, periodic1, periodic2};

    //
    // linearInterpolation
    Finley_Mesh* fMesh=NULL;

    if (order==1) {
      fMesh=Finley_RectangularMesh_Hex8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE)) ;
    } 
		else if (order==2) {
      fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					 useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE)) ;
    } else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }
  AbstractContinuousDomain*  rectangle(int n0,int n1,int order,
			double l0, double l1,
			int periodic0,int periodic1,
			int integrationOrder,
                        int reducedIntegrationOrder,
			int useElementsOnFace,
		        int useFullElementOrder,
                        int optimize)
  {
    int numElements[]={n0,n1};
    double length[]={l0,l1};
    int periodic[]={periodic0, periodic1};

    Finley_Mesh* fMesh=0;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Rec4(numElements, length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE));
    }
    else if (order==2) {
      fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder,(optimize ? TRUE : FALSE));
    }
    else {
      stringstream temp;
      temp << "Illegal interpolation order: " << order;
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }

  AbstractContinuousDomain* meshMerge(const boost::python::list& meshList)
  {
    Finley_Mesh* fMesh=0;
    //
    // extract the meshes from meshList
    int numMsh=boost::python::extract<int>(meshList.attr("__len__")());
    Finley_Mesh **mshes = (numMsh) ? TMPMEMALLOC(numMsh,Finley_Mesh*) : (Finley_Mesh**)NULL;
    for (int i=0;i<numMsh;++i) {
         AbstractContinuousDomain& meshListMember=boost::python::extract<AbstractContinuousDomain&>(meshList[i]);
         const MeshAdapter* finley_meshListMember=static_cast<const MeshAdapter*>(&meshListMember);
         mshes[i]=finley_meshListMember->getFinley_Mesh();
    }
    //
    // merge the meshes:
    fMesh=Finley_Mesh_merge(numMsh,mshes);
	  TMPMEMFREE(mshes);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);

    return temp;
  }
  AbstractContinuousDomain*  glueFaces(const boost::python::list& meshList,
                 	               double safety_factor, 
			               double tolerance,
                                       int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_glueFaces(fMesh,safety_factor,tolerance,(optimize ? TRUE : FALSE));

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }
  AbstractContinuousDomain*  joinFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance,
                        int optimize)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_joinFaces(fMesh,safety_factor,tolerance, (optimize ? TRUE : FALSE));
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }

  // end of namespace

}
