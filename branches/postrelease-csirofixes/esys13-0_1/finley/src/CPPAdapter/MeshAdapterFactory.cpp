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

#include "MeshAdapterFactory.h"
#include "FinleyError.h"

#include <boost/python/extract.hpp>

#include <sstream>

using namespace std;
using namespace escript;

namespace finley {

  AbstractContinuousDomain* readMesh(const std::string& fileName,
				     int integrationOrder) 
  {
    //
    // create a copy of the filename to overcome the non-constness of call
    // to Finley_Mesh_read
    char fName[fileName.size()+1];
    strcpy(fName,fileName.c_str());
    Finley_Mesh* fMesh=Finley_Mesh_read(fName,integrationOrder);
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }

  AbstractContinuousDomain* brick(int n0,int n1,int n2,int order,
		    double l0,double l1,double l2,
		    int periodic0,int periodic1,
		    int periodic2,
		    int integrationOrder,
		    int useElementsOnFace) 
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
    Finley_Mesh* fMesh;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Hex8(numElements,length,periodic,integrationOrder,
					useElementsOnFace) ;
    } else if (order==2) {
      fMesh=Finley_RectangularMesh_Hex20(numElements,length,periodic,integrationOrder,
					 useElementsOnFace) ;
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
			int useElementsOnFace) 
  {
    int numElements[]={n0,n1};
    double length[]={l0,l1};
    int periodic[]={periodic0, periodic1};

    Finley_Mesh* fMesh;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Rec4(numElements, length,periodic,integrationOrder,
					useElementsOnFace);
    } else if (order==2) {
      fMesh=Finley_RectangularMesh_Rec8(numElements,length,periodic,integrationOrder,
					useElementsOnFace);
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
  AbstractContinuousDomain*  interval(int n0,int order,double l0,int periodic0,
		       int integrationOrder,
		       int useElementsOnFace) 
  {
    int numElements[]={n0};
    double length[]={l0};
    int periodic[]={periodic0};
    Finley_Mesh* fMesh;
    if (order==1) {
      fMesh=Finley_RectangularMesh_Line2(numElements, length,periodic,integrationOrder,
					 useElementsOnFace);
    } else if (order==2) {
      fMesh=Finley_RectangularMesh_Line3(numElements,length,periodic,integrationOrder,
					 useElementsOnFace);
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
  AbstractContinuousDomain* meshMerge(const boost::python::list& meshList)
  {
    Finley_Mesh* fMesh=0;
    //
    // extract the meshes from meshList
    int numMsh=boost::python::extract<int>(meshList.attr("__len__")());
    Finley_Mesh* mshes[numMsh];
    for (int i=0;i<numMsh;++i) {
         AbstractContinuousDomain& meshListMember=boost::python::extract<AbstractContinuousDomain&>(meshList[i]);
         const MeshAdapter* finley_meshListMember=static_cast<const MeshAdapter*>(&meshListMember);
         mshes[i]=finley_meshListMember->getFinley_Mesh();
    }
    //
    // merge the meshes:
    fMesh=Finley_Mesh_merge(numMsh,mshes);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    AbstractContinuousDomain* temp=new MeshAdapter(fMesh);
    return temp;
  }
  AbstractContinuousDomain*  glueFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // glue the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_glueFaces(fMesh,safety_factor,tolerance);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }
  AbstractContinuousDomain*  joinFaces(const boost::python::list& meshList,
			double safety_factor, 
			double tolerance)
  {
    Finley_Mesh* fMesh=0;
    //
    // merge the meshes:
    AbstractContinuousDomain* merged_meshes=meshMerge(meshList);
    //
    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes);
    fMesh=merged_finley_meshes->getFinley_Mesh();
    Finley_Mesh_joinFaces(fMesh,safety_factor,tolerance);
    //
    // Convert any finley errors into a C++ exception
    checkFinleyError();
    return merged_meshes;
  }

}  // end of namespace
