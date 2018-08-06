
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include <finley/DomainFactory.h>

#include <escript/index.h>
#include <escript/SubWorld.h>

#ifdef ESYS_HAVE_NETCDF
 #ifdef NETCDF4
  #include <ncDim.h>
  #include <ncVar.h>
  #include <ncFile.h>
  
 #include <escript/NCHelper.h>  
  
 #else
  #include <netcdfcpp.h>
 #endif
#endif


#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>

using namespace std;
using namespace escript;
namespace bp = boost::python;

#ifdef NETCDF4
using namespace netCDF;
#endif

namespace finley {

#ifdef ESYS_HAVE_NETCDF
#ifdef NETCDF4

// A convenience method to retrieve an integer attribute from a NetCDF file
template<typename T>
T ncReadAtt(NcFile& dataFile, const string& fName, const string& attrName)
{
    NcGroupAtt attr = dataFile.getAtt(attrName.c_str());
    if (attr.isNull()) {
        stringstream msg;
        msg << "loadMesh: Error retrieving integer attribute '" << attrName
            << "' from NetCDF file '" << fName << "'";
        throw IOError(msg.str());
    }
    T value;
    attr.getValues(&value);
    return value;
}
#else    
    
// A convenience method to retrieve an integer attribute from a NetCDF file
template<typename T>
T ncReadAtt(NcFile* dataFile, const string& fName, const string& attrName)
{
    NcAtt* attr = dataFile->get_att(attrName.c_str());
    if (!attr) {
        stringstream msg;
        msg << "loadMesh: Error retrieving integer attribute '" << attrName
            << "' from NetCDF file '" << fName << "'";
        throw IOError(msg.str());
    }
    T value = (sizeof(T) > 4 ? attr->as_long(0) : attr->as_int(0));
    delete attr;
    return value;
}

#endif
#endif

inline void cleanupAndThrow(FinleyDomain* dom, string msg)
{
    delete dom;
    string msgPrefix("loadMesh: NetCDF operation failed - ");
    throw IOError(msgPrefix+msg);
}

#ifdef NETCDF4

Domain_ptr FinleyDomain::load(const string& fileName)
{
#ifdef ESYS_HAVE_NETCDF
    JMPI mpiInfo = makeInfo(MPI_COMM_WORLD);
    const string fName(mpiInfo->appendRankToFileName(fileName));

    // Open NetCDF file for reading
    NcGroupAtt attr;
    NcVar nc_var_temp;
    NcFile dataFile;
    // Create the NetCDF file.
    if (!openNcFile(dataFile, fileName))
    {
        throw IOError("load: opening of netCDF file for input failed.");
    }    
    // Read NetCDF integer attributes

    // index_size was only introduced with 64-bit index support so fall back
    // to 32 bits if not found.
    int index_size;
    try {
        index_size = ncReadAtt<int>(dataFile, fName, "index_size");
    } catch (IOError& e) {
        index_size = 4;
    }
    // technically we could cast if reading 32-bit data on 64-bit escript
    // but cost-benefit analysis clearly favours this implementation for now
    if (sizeof(index_t) != index_size) {
        throw IOError("loadMesh: size of index types at runtime differ from dump file");
    }

    int mpi_size = ncReadAtt<int>(dataFile, fName, "mpi_size");
    int mpi_rank = ncReadAtt<int>(dataFile, fName, "mpi_rank");
    int numDim = ncReadAtt<int>(dataFile, fName, "numDim");
    int order = ncReadAtt<int>(dataFile, fName, "order");
    int reduced_order = ncReadAtt<int>(dataFile, fName, "reduced_order");
    dim_t numNodes = ncReadAtt<dim_t>(dataFile, fName, "numNodes");
    dim_t num_Elements = ncReadAtt<dim_t>(dataFile, fName, "num_Elements");
    dim_t num_FaceElements = ncReadAtt<dim_t>(dataFile, fName, "num_FaceElements");
    dim_t num_ContactElements = ncReadAtt<dim_t>(dataFile, fName, "num_ContactElements");
    dim_t num_Points = ncReadAtt<dim_t>(dataFile, fName, "num_Points");
    int num_Elements_numNodes = ncReadAtt<int>(dataFile, fName, "num_Elements_numNodes");
    int Elements_TypeId = ncReadAtt<int>(dataFile, fName, "Elements_TypeId");
    int num_FaceElements_numNodes = ncReadAtt<int>(dataFile, fName, "num_FaceElements_numNodes");
    int FaceElements_TypeId = ncReadAtt<int>(dataFile, fName, "FaceElements_TypeId");
    int num_ContactElements_numNodes = ncReadAtt<int>(dataFile, fName, "num_ContactElements_numNodes");
    int ContactElements_TypeId = ncReadAtt<int>(dataFile, fName, "ContactElements_TypeId");
    int Points_TypeId = ncReadAtt<int>(dataFile, fName, "Points_TypeId");
    int num_Tags = ncReadAtt<int>(dataFile, fName, "num_Tags");

    // Verify size and rank
    if (mpiInfo->size != mpi_size) {
        stringstream msg;
        msg << "loadMesh: The NetCDF file '" << fName
            << "' can only be read on " << mpi_size
            << " CPUs. Currently running: " << mpiInfo->size;
        throw FinleyException(msg.str());
    }
    if (mpiInfo->rank != mpi_rank) {
        stringstream msg;
        msg << "loadMesh: The NetCDF file '" << fName
            << "' should be read on CPU #" << mpi_rank
            << " and NOT on #" << mpiInfo->rank;
        throw FinleyException(msg.str());
    }

    // Read mesh name
    if ((attr=dataFile.getAtt("Name")).isNull() ) {
        stringstream msg;
        msg << "loadMesh: Error retrieving mesh name from NetCDF file '"
            << fName << "'";
        throw IOError(msg.str());
    }
    string name;
    attr.getValues(name);

    // allocate mesh
    FinleyDomain* dom = new FinleyDomain(name, numDim, mpiInfo);

    // read nodes
    NodeFile* nodes = dom->getNodes();
    nodes->allocTable(numNodes);

    try
    {
        // Nodes_Id
        nc_var_temp = dataFile.getVar("Nodes_Id");
        nc_var_temp.getVar(&nodes->Id[0]); // numNodes values read
        // Nodes_Tag
        nc_var_temp = dataFile.getVar("Nodes_Tag");
        nc_var_temp.getVar(&nodes->Tag[0]);   // numNodes
        // Nodes_gDOF
        nc_var_temp = dataFile.getVar("Nodes_gDOF");
        nc_var_temp.getVar(&nodes->globalDegreesOfFreedom[0]);   //numNodes
        // Nodes_gNI
        nc_var_temp = dataFile.getVar("Nodes_gNI");
        nc_var_temp.getVar(&nodes->globalNodesIndex[0]); // numNodes
        // Nodes_grDfI
        nc_var_temp = dataFile.getVar("Nodes_grDfI");
        nc_var_temp.getVar(&nodes->globalReducedDOFIndex[0]); // numNodes
        // Nodes_grNI
        nc_var_temp = dataFile.getVar("Nodes_grNI");
        nc_var_temp.getVar(&nodes->globalReducedNodesIndex[0]);  // numNodes
        // Nodes_Coordinates
        nc_var_temp = dataFile.getVar("Nodes_Coordinates");
        nc_var_temp.getVar(&nodes->Coordinates[0]); // (numNodes, numDim) 
    }
    catch (exceptions::NcException& e)
    {
        cleanupAndThrow(dom, "Read vars from file");
    }
    nodes->updateTagList();

    // read elements
    const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(
                (ElementTypeId)Elements_TypeId, order, reduced_order));
    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    dom->setElements(elements);
    elements->allocTable(num_Elements);
    elements->minColor = 0;
    elements->maxColor = num_Elements-1;
    if (num_Elements > 0) {
        try
        {
            // Elements_Id
            nc_var_temp = dataFile.getVar("Elements_Id");
            nc_var_temp.getVar(&elements->Id[0]);    // num_Elements
            // Elements_Tag
            nc_var_temp = dataFile.getVar("Elements_Tag");
            nc_var_temp.getVar(&elements->Tag[0]); // num_Elements
            // Elements_Owner
            nc_var_temp = dataFile.getVar("Elements_Owner");
            nc_var_temp.getVar(&elements->Owner[0]);  // num_Elements
            // Elements_Color
            nc_var_temp = dataFile.getVar("Elements_Color");
            nc_var_temp.getVar(&elements->Color[0]); // num_Elements
        }
        catch (exceptions::NcException& e)
        {
            cleanupAndThrow(dom, "Readig element vars");            
        }
       // Now we need to adjust maxColor
       index_t mc = elements->Color[0];
       for (index_t i = 1; i < num_Elements; ++i) {
           if (mc < elements->Color[i]) {
               mc = elements->Color[i];
           }
       }
       elements->maxColor = mc;
       // Elements_Nodes
       int* Elements_Nodes = new int[num_Elements*num_Elements_numNodes];       
       try
       {
            nc_var_temp = dataFile.getVar("Elements_Nodes");
            nc_var_temp.getVar(&Elements_Nodes[0]);    // (num_Elements, num_Elements_numNodes) )             
       }
       catch (exceptions::NcException& e)
       {
           delete[] Elements_Nodes;
           cleanupAndThrow(dom, "get_var(Elements_Nodes)");
       }
       // Copy temp array into elements->Nodes
       for (index_t i = 0; i < num_Elements; i++) {
           for (int j = 0; j < num_Elements_numNodes; j++) {
               elements->Nodes[INDEX2(j,i,num_Elements_numNodes)]
                    = Elements_Nodes[INDEX2(j,i,num_Elements_numNodes)];
           }
       }
       delete[] Elements_Nodes;
    } // num_Elements > 0
    elements->updateTagList();

    // get the face elements
    const_ReferenceElementSet_ptr refFaceElements(
            new ReferenceElementSet((ElementTypeId)FaceElements_TypeId,
                order, reduced_order));
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    dom->setFaceElements(faces);
    faces->allocTable(num_FaceElements);
    faces->minColor = 0;
    faces->maxColor = num_FaceElements-1;
    if (num_FaceElements > 0) {
        try
        {
            // FaceElements_Id
            nc_var_temp = dataFile.getVar("FaceElements_Id");
            nc_var_temp.getVar(&faces->Id[0]); // num_FaceElements
            // FaceElements_Tag
            nc_var_temp = dataFile.getVar("FaceElements_Tag");
            nc_var_temp.getVar(&faces->Tag[0]); // num_FaceElements
            // FaceElements_Owner
            nc_var_temp = dataFile.getVar("FaceElements_Owner");
            nc_var_temp.getVar(&faces->Owner[0]);   // num_FaceElements
            // FaceElements_Color
            nc_var_temp = dataFile.getVar("FaceElements_Color");
            nc_var_temp.getVar(&faces->Color[0]); // num_FaceElements
        }
        catch (exceptions::NcException& e)
        {
            cleanupAndThrow(dom, "read face variables");
        }
        // Now we need to adjust maxColor
        index_t mc = faces->Color[0];
        for (index_t i = 1; i < num_FaceElements; ++i) {
            if (mc < faces->Color[i]) {
                mc = faces->Color[i];
            }
        }
        faces->maxColor = mc;
        // FaceElements_Nodes
        int* FaceElements_Nodes = new int[num_FaceElements*num_FaceElements_numNodes];
        try
        {
            nc_var_temp = dataFile.getVar("FaceElements_Nodes");
            nc_var_temp.getVar(&(FaceElements_Nodes[0])); // num_FaceElements, num_FaceElements_numNodes) ) 
        }
        catch (exceptions::NcException& e)
        {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(dom, "read face elements");
        }
        if ((nc_var_temp = dataFile.getVar("FaceElements_Nodes")), nc_var_temp.isNull()) {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(dom, "get_var(FaceElements_Nodes)");
        }
        nc_var_temp.getVar(&(FaceElements_Nodes[0]));  // num_FaceElements, num_FaceElements_numNodes
        // Copy temp array into faces->Nodes
        for (index_t i = 0; i < num_FaceElements; i++) {
            for (int j = 0; j < num_FaceElements_numNodes; j++) {
                faces->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
            }
        }
        delete[] FaceElements_Nodes;
    } // num_FaceElements > 0
    faces->updateTagList();

    // get the Contact elements
    const_ReferenceElementSet_ptr refContactElements(
         new ReferenceElementSet((ElementTypeId)ContactElements_TypeId,
             order, reduced_order));
    ElementFile* contacts = new ElementFile(refContactElements, mpiInfo);
    dom->setContactElements(contacts);
    contacts->allocTable(num_ContactElements);
    contacts->minColor = 0;
    contacts->maxColor = num_ContactElements-1;
    if (num_ContactElements > 0) {
        // ContactElements_Id
        if (( nc_var_temp = dataFile.getVar("ContactElements_Id")), nc_var_temp.isNull() )
            cleanupAndThrow(dom, "get_var(ContactElements_Id)");
        nc_var_temp.getVar(&contacts->Id[0]);       // num_ContactElements
        // ContactElements_Tag
        if (( nc_var_temp = dataFile.getVar("ContactElements_Tag")), nc_var_temp.isNull() )
            cleanupAndThrow(dom, "get_var(ContactElements_Tag)");
        nc_var_temp.getVar(&contacts->Tag[0]);  //, num_ContactElements
        // ContactElements_Owner
        if (( nc_var_temp = dataFile.getVar("ContactElements_Owner")), nc_var_temp.isNull() )
            cleanupAndThrow(dom, "get_var(ContactElements_Owner)");
        nc_var_temp.getVar(&contacts->Owner[0]);    //, num_ContactElements)
        // ContactElements_Color
        if (( nc_var_temp = dataFile.getVar("ContactElements_Color")), nc_var_temp.isNull() )
            cleanupAndThrow(dom, "get_var(ContactElements_Color)");
        nc_var_temp.getVar(&contacts->Color[0]);  //, num_ContactElements
        // Now we need to adjust maxColor
        index_t mc = contacts->Color[0];
        for (index_t i = 1; i < num_ContactElements; ++i) {
            if (mc < contacts->Color[i]) {
                mc = contacts->Color[i];
            }
        }
        contacts->maxColor = mc;
        // ContactElements_Nodes
        int* ContactElements_Nodes = new int[num_ContactElements*num_ContactElements_numNodes];
        if ((nc_var_temp = dataFile.getVar("ContactElements_Nodes")), nc_var_temp.isNull()) {
            delete[] ContactElements_Nodes;
            cleanupAndThrow(dom, "get_var(ContactElements_Nodes)");
        }
        nc_var_temp.getVar(&ContactElements_Nodes[0]);    // num_ContactElements, num_ContactElements_numNodes
        // Copy temp array into contacts->Nodes
        for (index_t i = 0; i < num_ContactElements; i++) {
            for (int j = 0; j < num_ContactElements_numNodes; j++) {
                contacts->Nodes[INDEX2(j,i,num_ContactElements_numNodes)] = ContactElements_Nodes[INDEX2(j,i,num_ContactElements_numNodes)];
            }
        }
        delete[] ContactElements_Nodes;
    } // num_ContactElements > 0
    contacts->updateTagList();

    // get the Points (nodal elements)
    const_ReferenceElementSet_ptr refPoints(new ReferenceElementSet(
                (ElementTypeId)Points_TypeId, order, reduced_order));
    ElementFile* points = new ElementFile(refPoints, mpiInfo);
    dom->setPoints(points);
    points->allocTable(num_Points);
    points->minColor = 0;
    points->maxColor = num_Points-1;
    if (num_Points > 0) {
        // Points_Id
        if (( nc_var_temp = dataFile.getVar("Points_Id")), nc_var_temp.isNull())
            cleanupAndThrow(dom, "get_var(Points_Id)");
        nc_var_temp.getVar(&points->Id[0]); // num_Points
        // Points_Tag
        if (( nc_var_temp = dataFile.getVar("Points_Tag")), nc_var_temp.isNull())
            cleanupAndThrow(dom, "get_var(Points_Tag)");
        nc_var_temp.getVar(&points->Tag[0]);    // num_Points
        // Points_Owner
        if (( nc_var_temp = dataFile.getVar("Points_Owner")), nc_var_temp.isNull())
            cleanupAndThrow(dom, "get_var(Points_Owner)");
        nc_var_temp.getVar(&points->Owner[0]);  // num_Points
        // Points_Color
        if (( nc_var_temp = dataFile.getVar("Points_Color")), nc_var_temp.isNull())
            cleanupAndThrow(dom, "get_var(Points_Color)");
        nc_var_temp.getVar(&points->Color[0]);  // num_Points
        // Now we need to adjust maxColor
        index_t mc = points->Color[0];
        for (index_t i = 1; i < num_Points; ++i) {
            if (mc < points->Color[i]) {
                mc = points->Color[i];
            }
        }
        points->maxColor = mc;
        // Points_Nodes
        int* Points_Nodes = new int[num_Points];
        if ((nc_var_temp = dataFile.getVar("Points_Nodes")), nc_var_temp.isNull()) {
            delete[] Points_Nodes;
            cleanupAndThrow(dom, "get_var(Points_Nodes)");
        }
        nc_var_temp.getVar(&Points_Nodes[0]);  // num_Points
        // Copy temp array into points->Nodes
        for (index_t i = 0; i < num_Points; i++) {
            points->Id[points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
        }
        delete[] Points_Nodes;
    } // num_Points > 0
    points->updateTagList();

    // get the tags
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        int *Tags_keys = new int[num_Tags];
        char name_temp[4096];
        int i;

        // Tags_keys
        if (( nc_var_temp = dataFile.getVar("Tags_keys")), nc_var_temp.isNull() ) {
            delete[] Tags_keys;
            cleanupAndThrow(dom, "get_var(Tags_keys)");
        }
        nc_var_temp.getVar(&Tags_keys[0]); // num_Tags
        for (i=0; i<num_Tags; i++) {
          // Retrieve tag name
          sprintf(name_temp, "Tags_name_%d", i);
          if ((attr=dataFile.getAtt(name_temp)), attr.isNull() ) {
              delete[] Tags_keys;
              stringstream msg;
              msg << "get_att(" << name_temp << ")";
              cleanupAndThrow(dom, msg.str());
          }
          std::string name;
          attr.getValues(name);
//           boost::scoped_array<char> name(attr->as_string(0));
//           delete attr;
//           dom->setTagMap(name.get(), Tags_keys[i]);
          dom->setTagMap(name.c_str(), Tags_keys[i]);
        }
        delete[] Tags_keys;
    }

    // Nodes_DofDistribution
    IndexVector first_DofComponent(mpi_size+1);
    if ((nc_var_temp = dataFile.getVar("Nodes_DofDistribution")), nc_var_temp.isNull() ) {
        cleanupAndThrow(dom, "get_var(Nodes_DofDistribution)");
    }
    nc_var_temp.getVar(&first_DofComponent[0]); // mpi_size+1

    // Nodes_NodeDistribution
    IndexVector first_NodeComponent(mpi_size+1);
    if ((nc_var_temp = dataFile.getVar("Nodes_NodeDistribution")), nc_var_temp.isNull() ) {
        cleanupAndThrow(dom, "get_var(Nodes_NodeDistribution)");
    }
    nc_var_temp.getVar(&first_NodeComponent[0]);    // mpi_size+1
    dom->createMappings(first_DofComponent, first_NodeComponent);

    return dom->getPtr();
#else
    throw FinleyException("loadMesh: not compiled with NetCDF. Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}


#else

Domain_ptr FinleyDomain::load(const string& fileName)
{
#ifdef ESYS_HAVE_NETCDF
    JMPI mpiInfo = makeInfo(MPI_COMM_WORLD);
    const string fName(mpiInfo->appendRankToFileName(fileName));

    // Open NetCDF file for reading
    NcAtt *attr;
    NcVar *nc_var_temp;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    // Create the NetCDF file.
    NcFile dataFile(fName.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid()) {
        stringstream msg;
        msg << "loadMesh: Opening NetCDF file '" << fName << "' for reading failed.";
        throw IOError(msg.str());
    }

    // Read NetCDF integer attributes

    // index_size was only introduced with 64-bit index support so fall back
    // to 32 bits if not found.
    int index_size;
    try {
        index_size = ncReadAtt<int>(&dataFile, fName, "index_size");
    } catch (IOError& e) {
        index_size = 4;
    }
    // technically we could cast if reading 32-bit data on 64-bit escript
    // but cost-benefit analysis clearly favours this implementation for now
    if (sizeof(index_t) != index_size) {
        throw IOError("loadMesh: size of index types at runtime differ from dump file");
    }

    int mpi_size = ncReadAtt<int>(&dataFile, fName, "mpi_size");
    int mpi_rank = ncReadAtt<int>(&dataFile, fName, "mpi_rank");
    int numDim = ncReadAtt<int>(&dataFile, fName, "numDim");
    int order = ncReadAtt<int>(&dataFile, fName, "order");
    int reduced_order = ncReadAtt<int>(&dataFile, fName, "reduced_order");
    dim_t numNodes = ncReadAtt<dim_t>(&dataFile, fName, "numNodes");
    dim_t num_Elements = ncReadAtt<dim_t>(&dataFile, fName, "num_Elements");
    dim_t num_FaceElements = ncReadAtt<dim_t>(&dataFile, fName, "num_FaceElements");
    dim_t num_ContactElements = ncReadAtt<dim_t>(&dataFile, fName, "num_ContactElements");
    dim_t num_Points = ncReadAtt<dim_t>(&dataFile, fName, "num_Points");
    int num_Elements_numNodes = ncReadAtt<int>(&dataFile, fName, "num_Elements_numNodes");
    int Elements_TypeId = ncReadAtt<int>(&dataFile, fName, "Elements_TypeId");
    int num_FaceElements_numNodes = ncReadAtt<int>(&dataFile, fName, "num_FaceElements_numNodes");
    int FaceElements_TypeId = ncReadAtt<int>(&dataFile, fName, "FaceElements_TypeId");
    int num_ContactElements_numNodes = ncReadAtt<int>(&dataFile, fName, "num_ContactElements_numNodes");
    int ContactElements_TypeId = ncReadAtt<int>(&dataFile, fName, "ContactElements_TypeId");
    int Points_TypeId = ncReadAtt<int>(&dataFile, fName, "Points_TypeId");
    int num_Tags = ncReadAtt<int>(&dataFile, fName, "num_Tags");

    // Verify size and rank
    if (mpiInfo->size != mpi_size) {
        stringstream msg;
        msg << "loadMesh: The NetCDF file '" << fName
            << "' can only be read on " << mpi_size
            << " CPUs. Currently running: " << mpiInfo->size;
        throw FinleyException(msg.str());
    }
    if (mpiInfo->rank != mpi_rank) {
        stringstream msg;
        msg << "loadMesh: The NetCDF file '" << fName
            << "' should be read on CPU #" << mpi_rank
            << " and NOT on #" << mpiInfo->rank;
        throw FinleyException(msg.str());
    }

    // Read mesh name
    if (! (attr=dataFile.get_att("Name")) ) {
        stringstream msg;
        msg << "loadMesh: Error retrieving mesh name from NetCDF file '"
            << fName << "'";
        throw IOError(msg.str());
    }
    boost::scoped_array<char> name(attr->as_string(0));
    delete attr;

    // allocate mesh
    FinleyDomain* dom = new FinleyDomain(name.get(), numDim, mpiInfo);

    // read nodes
    NodeFile* nodes = dom->getNodes();
    nodes->allocTable(numNodes);
    // Nodes_Id
    if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
        cleanupAndThrow(dom, "get_var(Nodes_Id)");
    if (! nc_var_temp->get(&nodes->Id[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_Id)");
    // Nodes_Tag
    if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
        cleanupAndThrow(dom, "get_var(Nodes_Tag)");
    if (! nc_var_temp->get(&nodes->Tag[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_Tag)");
    // Nodes_gDOF
    if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
        cleanupAndThrow(dom, "get_var(Nodes_gDOF)");
    if (! nc_var_temp->get(&nodes->globalDegreesOfFreedom[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_gDOF)");
    // Nodes_gNI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
        cleanupAndThrow(dom, "get_var(Nodes_gNI)");
    if (! nc_var_temp->get(&nodes->globalNodesIndex[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_gNI)");
    // Nodes_grDfI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
        cleanupAndThrow(dom, "get_var(Nodes_grDfI)");
    if (! nc_var_temp->get(&nodes->globalReducedDOFIndex[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_grDfI)");
    // Nodes_grNI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
        cleanupAndThrow(dom, "get_var(Nodes_grNI)");
    if (! nc_var_temp->get(&nodes->globalReducedNodesIndex[0], numNodes) )
        cleanupAndThrow(dom, "get(Nodes_grNI)");
    // Nodes_Coordinates
    if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates")))
        cleanupAndThrow(dom, "get_var(Nodes_Coordinates)");
    if (! nc_var_temp->get(&nodes->Coordinates[0], numNodes, numDim) )
        cleanupAndThrow(dom, "get(Nodes_Coordinates)");

    nodes->updateTagList();

    // read elements
    const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(
                (ElementTypeId)Elements_TypeId, order, reduced_order));
    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    dom->setElements(elements);
    elements->allocTable(num_Elements);
    elements->minColor = 0;
    elements->maxColor = num_Elements-1;
    if (num_Elements > 0) {
       // Elements_Id
       if (! ( nc_var_temp = dataFile.get_var("Elements_Id")) )
           cleanupAndThrow(dom, "get_var(Elements_Id)");
       if (! nc_var_temp->get(&elements->Id[0], num_Elements) )
           cleanupAndThrow(dom, "get(Elements_Id)");
       // Elements_Tag
       if (! ( nc_var_temp = dataFile.get_var("Elements_Tag")) )
           cleanupAndThrow(dom, "get_var(Elements_Tag)");
       if (! nc_var_temp->get(&elements->Tag[0], num_Elements) )
           cleanupAndThrow(dom, "get(Elements_Tag)");
       // Elements_Owner
       if (! ( nc_var_temp = dataFile.get_var("Elements_Owner")) )
           cleanupAndThrow(dom, "get_var(Elements_Owner)");
       if (! nc_var_temp->get(&elements->Owner[0], num_Elements) )
           cleanupAndThrow(dom, "get(Elements_Owner)");
       // Elements_Color
       if (! ( nc_var_temp = dataFile.get_var("Elements_Color")) )
           cleanupAndThrow(dom, "get_var(Elements_Color)");
       if (! nc_var_temp->get(&elements->Color[0], num_Elements) )
           cleanupAndThrow(dom, "get(Elements_Color)");
       // Now we need to adjust maxColor
       index_t mc = elements->Color[0];
       for (index_t i = 1; i < num_Elements; ++i) {
           if (mc < elements->Color[i]) {
               mc = elements->Color[i];
           }
       }
       elements->maxColor = mc;
       // Elements_Nodes
       int* Elements_Nodes = new int[num_Elements*num_Elements_numNodes];
       if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
           delete[] Elements_Nodes;
           cleanupAndThrow(dom, "get_var(Elements_Nodes)");
       }
       if (! nc_var_temp->get(&Elements_Nodes[0], num_Elements, num_Elements_numNodes) ) {
           delete[] Elements_Nodes;
           cleanupAndThrow(dom, "get(Elements_Nodes)");
       }

       // Copy temp array into elements->Nodes
       for (index_t i = 0; i < num_Elements; i++) {
           for (int j = 0; j < num_Elements_numNodes; j++) {
               elements->Nodes[INDEX2(j,i,num_Elements_numNodes)]
                    = Elements_Nodes[INDEX2(j,i,num_Elements_numNodes)];
           }
       }
       delete[] Elements_Nodes;
    } // num_Elements > 0
    elements->updateTagList();

    // get the face elements
    const_ReferenceElementSet_ptr refFaceElements(
            new ReferenceElementSet((ElementTypeId)FaceElements_TypeId,
                order, reduced_order));
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    dom->setFaceElements(faces);
    faces->allocTable(num_FaceElements);
    faces->minColor = 0;
    faces->maxColor = num_FaceElements-1;
    if (num_FaceElements > 0) {
        // FaceElements_Id
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Id")) )
            cleanupAndThrow(dom, "get_var(FaceElements_Id)");
        if (! nc_var_temp->get(&faces->Id[0], num_FaceElements) )
            cleanupAndThrow(dom, "get(FaceElements_Id)");
        // FaceElements_Tag
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Tag")) )
            cleanupAndThrow(dom, "get_var(FaceElements_Tag)");
        if (! nc_var_temp->get(&faces->Tag[0], num_FaceElements) )
            cleanupAndThrow(dom, "get(FaceElements_Tag)");
        // FaceElements_Owner
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Owner")) )
            cleanupAndThrow(dom, "get_var(FaceElements_Owner)");
        if (! nc_var_temp->get(&faces->Owner[0], num_FaceElements) )
            cleanupAndThrow(dom, "get(FaceElements_Owner)");
        // FaceElements_Color
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Color")) )
            cleanupAndThrow(dom, "get_var(FaceElements_Color)");
        if (! nc_var_temp->get(&faces->Color[0], num_FaceElements) )
            cleanupAndThrow(dom, "get(FaceElements_Color)");
        // Now we need to adjust maxColor
        index_t mc = faces->Color[0];
        for (index_t i = 1; i < num_FaceElements; ++i) {
            if (mc < faces->Color[i]) {
                mc = faces->Color[i];
            }
        }
        faces->maxColor = mc;
        // FaceElements_Nodes
        int* FaceElements_Nodes = new int[num_FaceElements*num_FaceElements_numNodes];
        if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(dom, "get_var(FaceElements_Nodes)");
        }
        if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(dom, "get(FaceElements_Nodes)");
        }
        // Copy temp array into faces->Nodes
        for (index_t i = 0; i < num_FaceElements; i++) {
            for (int j = 0; j < num_FaceElements_numNodes; j++) {
                faces->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
            }
        }
        delete[] FaceElements_Nodes;
    } // num_FaceElements > 0
    faces->updateTagList();

    // get the Contact elements
    const_ReferenceElementSet_ptr refContactElements(
         new ReferenceElementSet((ElementTypeId)ContactElements_TypeId,
             order, reduced_order));
    ElementFile* contacts = new ElementFile(refContactElements, mpiInfo);
    dom->setContactElements(contacts);
    contacts->allocTable(num_ContactElements);
    contacts->minColor = 0;
    contacts->maxColor = num_ContactElements-1;
    if (num_ContactElements > 0) {
        // ContactElements_Id
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Id")) )
            cleanupAndThrow(dom, "get_var(ContactElements_Id)");
        if (! nc_var_temp->get(&contacts->Id[0], num_ContactElements) )
            cleanupAndThrow(dom, "get(ContactElements_Id)");
        // ContactElements_Tag
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Tag")) )
            cleanupAndThrow(dom, "get_var(ContactElements_Tag)");
        if (! nc_var_temp->get(&contacts->Tag[0], num_ContactElements) )
            cleanupAndThrow(dom, "get(ContactElements_Tag)");
        // ContactElements_Owner
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Owner")) )
            cleanupAndThrow(dom, "get_var(ContactElements_Owner)");
        if (! nc_var_temp->get(&contacts->Owner[0], num_ContactElements) )
            cleanupAndThrow(dom, "get(ContactElements_Owner)");
        // ContactElements_Color
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Color")) )
            cleanupAndThrow(dom, "get_var(ContactElements_Color)");
        if (! nc_var_temp->get(&contacts->Color[0], num_ContactElements) )
            cleanupAndThrow(dom, "get(ContactElements_Color)");
        // Now we need to adjust maxColor
        index_t mc = contacts->Color[0];
        for (index_t i = 1; i < num_ContactElements; ++i) {
            if (mc < contacts->Color[i]) {
                mc = contacts->Color[i];
            }
        }
        contacts->maxColor = mc;
        // ContactElements_Nodes
        int* ContactElements_Nodes = new int[num_ContactElements*num_ContactElements_numNodes];
        if (!(nc_var_temp = dataFile.get_var("ContactElements_Nodes"))) {
            delete[] ContactElements_Nodes;
            cleanupAndThrow(dom, "get_var(ContactElements_Nodes)");
        }
        if (! nc_var_temp->get(&ContactElements_Nodes[0], num_ContactElements, num_ContactElements_numNodes) ) {
            delete[] ContactElements_Nodes;
            cleanupAndThrow(dom, "get(ContactElements_Nodes)");
        }
        // Copy temp array into contacts->Nodes
        for (index_t i = 0; i < num_ContactElements; i++) {
            for (int j = 0; j < num_ContactElements_numNodes; j++) {
                contacts->Nodes[INDEX2(j,i,num_ContactElements_numNodes)] = ContactElements_Nodes[INDEX2(j,i,num_ContactElements_numNodes)];
            }
        }
        delete[] ContactElements_Nodes;
    } // num_ContactElements > 0
    contacts->updateTagList();

    // get the Points (nodal elements)
    const_ReferenceElementSet_ptr refPoints(new ReferenceElementSet(
                (ElementTypeId)Points_TypeId, order, reduced_order));
    ElementFile* points = new ElementFile(refPoints, mpiInfo);
    dom->setPoints(points);
    points->allocTable(num_Points);
    points->minColor = 0;
    points->maxColor = num_Points-1;
    if (num_Points > 0) {
        // Points_Id
        if (! ( nc_var_temp = dataFile.get_var("Points_Id")))
            cleanupAndThrow(dom, "get_var(Points_Id)");
        if (! nc_var_temp->get(&points->Id[0], num_Points))
            cleanupAndThrow(dom, "get(Points_Id)");
        // Points_Tag
        if (! ( nc_var_temp = dataFile.get_var("Points_Tag")))
            cleanupAndThrow(dom, "get_var(Points_Tag)");
        if (! nc_var_temp->get(&points->Tag[0], num_Points))
            cleanupAndThrow(dom, "get(Points_Tag)");
        // Points_Owner
        if (! ( nc_var_temp = dataFile.get_var("Points_Owner")))
            cleanupAndThrow(dom, "get_var(Points_Owner)");
        if (!nc_var_temp->get(&points->Owner[0], num_Points))
            cleanupAndThrow(dom, "get(Points_Owner)");
        // Points_Color
        if (! ( nc_var_temp = dataFile.get_var("Points_Color")))
            cleanupAndThrow(dom, "get_var(Points_Color)");
        if (!nc_var_temp->get(&points->Color[0], num_Points))
            cleanupAndThrow(dom, "get(Points_Color)");
        // Now we need to adjust maxColor
        index_t mc = points->Color[0];
        for (index_t i = 1; i < num_Points; ++i) {
            if (mc < points->Color[i]) {
                mc = points->Color[i];
            }
        }
        points->maxColor = mc;
        // Points_Nodes
        int* Points_Nodes = new int[num_Points];
        if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
            delete[] Points_Nodes;
            cleanupAndThrow(dom, "get_var(Points_Nodes)");
        }
        if (! nc_var_temp->get(&Points_Nodes[0], num_Points) ) {
            delete[] Points_Nodes;
            cleanupAndThrow(dom, "get(Points_Nodes)");
        }
        // Copy temp array into points->Nodes
        for (index_t i = 0; i < num_Points; i++) {
            points->Id[points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
        }
        delete[] Points_Nodes;
    } // num_Points > 0
    points->updateTagList();

    // get the tags
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        int *Tags_keys = new int[num_Tags];
        char name_temp[4096];
        int i;

        // Tags_keys
        if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) ) {
            delete[] Tags_keys;
            cleanupAndThrow(dom, "get_var(Tags_keys)");
        }
        if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
            delete[] Tags_keys;
            cleanupAndThrow(dom, "get(Tags_keys)");
        }
        for (i=0; i<num_Tags; i++) {
          // Retrieve tag name
          sprintf(name_temp, "Tags_name_%d", i);
          if (! (attr=dataFile.get_att(name_temp)) ) {
              delete[] Tags_keys;
              stringstream msg;
              msg << "get_att(" << name_temp << ")";
              cleanupAndThrow(dom, msg.str());
          }
          boost::scoped_array<char> name(attr->as_string(0));
          delete attr;
          dom->setTagMap(name.get(), Tags_keys[i]);
        }
        delete[] Tags_keys;
    }

    // Nodes_DofDistribution
    IndexVector first_DofComponent(mpi_size+1);
    if (! (nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) ) {
        cleanupAndThrow(dom, "get_var(Nodes_DofDistribution)");
    }
    if (!nc_var_temp->get(&first_DofComponent[0], mpi_size+1)) {
        cleanupAndThrow(dom, "get(Nodes_DofDistribution)");
    }

    // Nodes_NodeDistribution
    IndexVector first_NodeComponent(mpi_size+1);
    if (! (nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) ) {
        cleanupAndThrow(dom, "get_var(Nodes_NodeDistribution)");
    }
    if (!nc_var_temp->get(&first_NodeComponent[0], mpi_size+1)) {
        cleanupAndThrow(dom, "get(Nodes_NodeDistribution)");
    }
    dom->createMappings(first_DofComponent, first_NodeComponent);

    return dom->getPtr();
#else
    throw FinleyException("loadMesh: not compiled with NetCDF. Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}

#endif

Domain_ptr readMesh_driver(const bp::list& args)
{
    int l = len(args);
    if (l < 7) {
        throw ValueError("Insufficient arguments to readMesh_driver");
    }
    string fileName = bp::extract<string>(args[0])();
    int integrationOrder = bp::extract<int>(args[1])();
    int reducedIntegrationOrder = bp::extract<int>(args[2])();
    bool optimize = bp::extract<bool>(args[3])();
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    bp::list pypoints = bp::extract<bp::list>(args[4]);
    bp::list pytags = bp::extract<bp::list>(args[5]);
    int numpts = bp::extract<int>(pypoints.attr("__len__")());
    int numtags = bp::extract<int>(pytags.attr("__len__")());

    bp::object pworld = args[6];
    JMPI info;
    if (!pworld.is_none()) {
        bp::extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = makeInfo(MPI_COMM_WORLD);
    }
    Domain_ptr dom(FinleyDomain::read(info, fileName, integrationOrder,
                                      reducedIntegrationOrder, optimize));

    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = pypoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
              points.push_back(bp::extract<double>(temp[k]));
        }
    }
    // bricks use up to 200 but the existing tag check will find that
    int curmax = 40;
    const TagMap& tagmap = fd->getTagMap();
    // first we work out what tags are already in use
    for (TagMap::const_iterator it = tagmap.begin(); it != tagmap.end(); ++it) {
        if (it->second > curmax) {
            curmax = it->second+1;
        }
    }

    tags.resize(numtags, -1);
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(pytags[i]);
        bp::extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i] = ex_int();
            if (tags[i] >= curmax) {
                curmax = tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::const_iterator it = tagmap.find(s);
            if (it != tagmap.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                fd->setTagMap(s, curmax);
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Unable to extract tag value.");
        }
    }
    // now we need to add the dirac points
    fd->addDiracPoints(points, tags);
    return dom;
}

Domain_ptr readGmsh_driver(const bp::list& args)
{
    int l = len(args);
    if (l < 7) {
        throw ValueError("Insufficient arguments to readMesh_driver");
    }
    string fileName = bp::extract<string>(args[0])();
    int numDim = bp::extract<int>(args[1])();
    int integrationOrder = bp::extract<int>(args[2])();
    int reducedIntegrationOrder = bp::extract<int>(args[3])();
    bool optimize = bp::extract<bool>(args[4])();
    bool useMacroElements = bp::extract<bool>(args[5])();
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    bp::list pypoints = bp::extract<bp::list>(args[6]);
    bp::list pytags = bp::extract<bp::list>(args[7]);
    int numpts = bp::extract<int>(pypoints.attr("__len__")());
    int numtags = bp::extract<int>(pytags.attr("__len__")());
    bp::object pworld = args[8];
    JMPI info;
    if (!pworld.is_none()) {
        bp::extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = makeInfo(MPI_COMM_WORLD);
    }
    Domain_ptr dom(FinleyDomain::readGmsh(info, fileName, numDim,
                                     integrationOrder, reducedIntegrationOrder,
                                     optimize, useMacroElements));
    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = pypoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(bp::extract<double>(temp[k]));
        }
    }
    int curmax = 40; // bricks use up to 30
    const TagMap& tagmap = fd->getTagMap();
    // first we work out what tags are already in use
    for (TagMap::const_iterator it = tagmap.begin(); it != tagmap.end(); ++it) {
        if (it->second > curmax) {
            curmax = it->second+1;
        }
    }

    tags.resize(numtags, -1);
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(pytags[i]);
        bp::extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i] = ex_int();
            if (tags[i] >= curmax) {
                curmax = tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::const_iterator it = tagmap.find(s);
            if (it != tagmap.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                fd->setTagMap(s, curmax);
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Unable to extract tag value");
        }
    }
    // now we need to add the dirac points
    fd->addDiracPoints(points, tags);
    return dom;
}

Domain_ptr brick(JMPI info, dim_t n0, dim_t n1, dim_t n2, int order,
                 double l0, double l1, double l2,
                 bool periodic0, bool periodic1, bool periodic2,
                 int integrationOrder, int reducedIntegrationOrder,
                 bool useElementsOnFace, bool useFullElementOrder,
                 bool optimize, const std::vector<double>& points,
                 const std::vector<int>& tags,
                 const std::map<std::string, int>& tagNamesToNums)
{
    Domain_ptr dom;
    if (order == 1) {
        dom = FinleyDomain::createHex8(n0, n1, n2, l0, l1, l2, periodic0,
                   periodic1, periodic2, integrationOrder,
                   reducedIntegrationOrder, useElementsOnFace, optimize, info);
    } else if (order == 2) {
        dom = FinleyDomain::createHex20(n0, n1, n2, l0, l1, l2, periodic0,
                                   periodic1, periodic2, integrationOrder,
                                   reducedIntegrationOrder, useElementsOnFace,
                                   useFullElementOrder, false, optimize, info);
    } else if (order == -1) {
        dom = FinleyDomain::createHex20(n0, n1, n2, l0, l1, l2, periodic0,
                                   periodic1, periodic2, integrationOrder,
                                   reducedIntegrationOrder, useElementsOnFace,
                                   useFullElementOrder, true, optimize, info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw ValueError(message.str());
    }

    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());
    fd->addDiracPoints(points, tags);
    for (TagMap::const_iterator it = tagNamesToNums.begin(); it != tagNamesToNums.end(); ++it) {
        fd->setTagMap(it->first, it->second);
    }
    fd->getPoints()->updateTagList();
    return dom;
}

Domain_ptr brick_driver(const bp::list& args)
{
    // we need to convert lists to stl vectors
    bp::list pypoints = bp::extract<bp::list>(args[15]);
    bp::list pytags = bp::extract<bp::list>(args[16]);
    int numpts = bp::extract<int>(pypoints.attr("__len__")());
    int numtags = bp::extract<int>(pytags.attr("__len__")());
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);
    for (int i = 0; i < numpts; ++i) {
        bp::object temp = pypoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(bp::extract<double>(temp[k]));
        }
    }
    map<string, int> namestonums;
    int curmax = 40; // bricks use up to 30
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(pytags[i]);
        bp::extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i] = ex_int();
            if (tags[i] >= curmax) {
                curmax = tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::iterator it = namestonums.find(s);
            if (it != namestonums.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                namestonums[s] = curmax;
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Unable to extract tag value.");
        }
    }
    bp::object pworld = args[17];
    JMPI info;
    if (!pworld.is_none()) {
        bp::extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = makeInfo(MPI_COMM_WORLD);
    }
    return brick(info, static_cast<dim_t>(bp::extract<float>(args[0])),
                 static_cast<dim_t>(bp::extract<float>(args[1])),
                 static_cast<dim_t>(bp::extract<float>(args[2])),
                 bp::extract<int>(args[3]), bp::extract<double>(args[4]),
                 bp::extract<double>(args[5]), bp::extract<double>(args[6]),
                 bp::extract<int>(args[7]), bp::extract<int>(args[8]),
                 bp::extract<int>(args[9]), bp::extract<int>(args[10]),
                 bp::extract<int>(args[11]), bp::extract<int>(args[12]),
                 bp::extract<int>(args[13]), bp::extract<int>(args[14]),
                 points, tags, namestonums);
}

Domain_ptr rectangle(JMPI info, dim_t n0, dim_t n1, int order,
                     double l0, double l1, bool periodic0, bool periodic1,
                     int integrationOrder, int reducedIntegrationOrder,
                     bool useElementsOnFace, bool useFullElementOrder,
                     bool optimize, const vector<double>& points,
                     const vector<int>& tags,
                     const std::map<std::string, int>& tagNamesToNums)
{
    Domain_ptr dom;
    if (order == 1) {
        dom = FinleyDomain::createRec4(n0, n1, l0, l1, periodic0, periodic1,
                                     integrationOrder, reducedIntegrationOrder,
                                     useElementsOnFace, optimize, info);
    } else if (order == 2) {
        dom = FinleyDomain::createRec8(n0, n1, l0, l1, periodic0, periodic1,
                 integrationOrder, reducedIntegrationOrder,
                 useElementsOnFace,useFullElementOrder, false, optimize, info);
    } else if (order == -1) {
        dom = FinleyDomain::createRec8(n0, n1, l0, l1, periodic0, periodic1,
                 integrationOrder, reducedIntegrationOrder,
                 useElementsOnFace, useFullElementOrder, true, optimize, info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw ValueError(message.str());
    }

    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());
    fd->addDiracPoints(points, tags);
    for (TagMap::const_iterator it = tagNamesToNums.begin(); it != tagNamesToNums.end(); ++it)
    {
        fd->setTagMap(it->first, it->second);
    }
    fd->getPoints()->updateTagList();
    return dom;
}

Domain_ptr rectangle_driver(const bp::list& args)
{
    // we need to convert lists to stl vectors
    bp::list pypoints = bp::extract<bp::list>(args[12]);
    bp::list pytags = bp::extract<bp::list>(args[13]);
    int numpts = bp::extract<int>(pypoints.attr("__len__")());
    int numtags = bp::extract<int>(pytags.attr("__len__")());
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);
    for (int i = 0; i < numpts; ++i) {
        bp::object temp = pypoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(bp::extract<double>(temp[k]));
        }
    }
    TagMap tagstonames;
    int curmax = 40;
    // but which order to assign tags to names?????
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(pytags[i]);
        bp::extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i] = ex_int();
            if (tags[i] >= curmax) {
                curmax = tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::iterator it = tagstonames.find(s);
            if (it != tagstonames.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                tagstonames[s] = curmax;
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Unable to extract tag value.");
        }
    }
    bp::object pworld = args[14];
    JMPI info;
    if (!pworld.is_none()) {
        bp::extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = makeInfo(MPI_COMM_WORLD);
    }

    return rectangle(info, static_cast<dim_t>(bp::extract<float>(args[0])),
                     static_cast<dim_t>(bp::extract<float>(args[1])),
                     bp::extract<int>(args[2]), bp::extract<double>(args[3]),
                     bp::extract<double>(args[4]), bp::extract<int>(args[5]),
                     bp::extract<int>(args[6]), bp::extract<int>(args[7]),
                     bp::extract<int>(args[8]), bp::extract<int>(args[9]),
                     bp::extract<int>(args[10]), bp::extract<int>(args[11]),
                     points, tags, tagstonames);
}

Domain_ptr meshMerge(const bp::list& meshList)
{
    // extract the meshes from meshList
    int num = bp::extract<int>(meshList.attr("__len__")());
    vector<const FinleyDomain*> meshes(num);
    for (int i = 0; i < num; ++i) {
        AbstractContinuousDomain& meshListMember = bp::extract<AbstractContinuousDomain&>(meshList[i]);
        meshes[i] = dynamic_cast<const FinleyDomain*>(&meshListMember);
    }

    // merge the meshes
    FinleyDomain* dom = FinleyDomain::merge(meshes);

    return dom->getPtr();
}

Domain_ptr glueFaces(const bp::list& meshList, double safetyFactor,
                     double tolerance, bool optimize)
{
    // merge the meshes
    Domain_ptr merged_meshes = meshMerge(meshList);

    // glue the faces
    FinleyDomain* merged = dynamic_cast<FinleyDomain*>(merged_meshes.get());
    merged->glueFaces(safetyFactor, tolerance, optimize);
    return merged_meshes;
}

Domain_ptr joinFaces(const bp::list& meshList, double safetyFactor,
                     double tolerance, bool optimize)
{
    // merge the meshes
    Domain_ptr merged_meshes = meshMerge(meshList);

    // join the faces
    FinleyDomain* merged = dynamic_cast<FinleyDomain*>(merged_meshes.get());
    merged->joinFaces(safetyFactor, tolerance, optimize);
    return merged_meshes;
}

} // namespace finley
