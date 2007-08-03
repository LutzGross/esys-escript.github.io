/* $Id$ */
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

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

  AbstractContinuousDomain* readMesh(const std::string& fileName,
  				     int integrationOrder,
                                     int reducedIntegrationOrder,
                                     bool optimizeLabeling)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

#ifndef PASO_MPI
    fMesh=Finley_Mesh_read(fName,integrationOrder, reducedIntegrationOrder, (optimizeLabeling ? TRUE : FALSE));
#else
    {
      stringstream temp;
      temp << "Unable to read meshes from file under MPI yet...";
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
#endif
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
                                     bool optimizeLabeling)
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    Finley_Mesh* fMesh=0;
    // Win32 refactor
    char *fName = ((fileName.size()+1)>0) ? TMPMEMALLOC((fileName.size()+1),char) : (char*)NULL;
    strcpy(fName,fileName.c_str());

#ifndef PASO_MPI
    fMesh=Finley_Mesh_readGmsh(fName, numDim, integrationOrder, reducedIntegrationOrder, (optimizeLabeling ? TRUE : FALSE));
#else
    {
      stringstream temp;
      temp << "Unable to read gmsh meshes from file under MPI yet...";
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
#endif
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
		    int useFullElementOrder) 
  {
//     cout << "n0=" << n0 << " n1=" << n1 << " n2=" << n2
// 	 << " order=" << order 
// 	 << " l0=" << l0 << " l1=" << l1 << " l2=" << l2
// 	 << " periodic0=" << periodic0 
// 	 << " periodic1=" << periodic1 
// 	 << " periodic2=" << periodic2
// 	 << " integerationOrder=" << integrationOrder
// 	 << " useElementsOnFace=" << useElementsOnFace << endl;
      
    int numElements[]={n0,n1,n2};
    double length[]={l0,l1,l2};
    int periodic[]={periodic0, periodic1, periodic2};

    //
    // linearInterpolation
    Finley_Mesh* fMesh=NULL;

    if (order==1) {
      fMesh=Finley_RectangularMesh_Hex8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder) ;
    } 
		else if (order==2) {
      fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					 useElementsOnFace,useFullElementOrder) ;
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
		        int useFullElementOrder) 
  {
    int numElements[]={n0,n1};
    double length[]={l0,l1};
    int periodic[]={periodic0, periodic1};

    Finley_Mesh* fMesh=0;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Rec4(numElements, length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder);
    }
    else if (order==2) {
      fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,reducedIntegrationOrder,
					useElementsOnFace,useFullElementOrder);
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
#ifndef PASO_MPI
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
#else
    {
      stringstream temp;
      temp << "meshMerge() not available in MPI yet...";
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
#endif
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);

    return temp;
  }
  AbstractContinuousDomain*  glueFaces(const boost::python::list& meshList,
                 		double safety_factor, 
			double tolerance,
                        bool optimizeLabeling)
  {
    Finley_Mesh* fMesh=0;
#ifndef PASO_MPI
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_glueFaces(fMesh,safety_factor,tolerance,(optimizeLabeling ? TRUE : FALSE));

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
#else
    {
      stringstream temp;
      temp << "glueFaces() not available in MPI yet...";
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }

    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return (AbstractContinuousDomain*)0;
#endif

  }
  AbstractContinuousDomain*  joinFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance,
                        bool optimizeLabeling)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
#ifndef PASO_MPI
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_joinFaces(fMesh,safety_factor,tolerance, (optimizeLabeling ? TRUE : FALSE));
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
#else
    {
      stringstream temp;
      temp << "joinFaces() not available in MPI yet...";
      setFinleyError(VALUE_ERROR,temp.str().c_str());
    }
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return (AbstractContinuousDomain*)0;

#endif
  }

  // end of namespace

}
