
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "MeshAdapterFactory.h"
#include <finley/FinleyException.h>

#include <escript/index.h>

#ifdef ESYS_HAVE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>

using namespace std;
using namespace escript;

namespace finley {

#ifdef ESYS_HAVE_NETCDF
// A convenience method to retrieve an integer attribute from a NetCDF file
template<typename T>
T ncReadAtt(NcFile* dataFile, const string& fName, const string& attrName)
{
    NcAtt* attr = dataFile->get_att(attrName.c_str());
    if (!attr) {
        stringstream msg;
        msg << "loadMesh: Error retrieving integer attribute '" << attrName
            << "' from NetCDF file '" << fName << "'";
        throw escript::IOError(msg.str());
    }
    T value = (sizeof(T) > 4 ? attr->as_long(0) : attr->as_int(0));
    delete attr;
    return value;
}
#endif

inline void cleanupAndThrow(Mesh* mesh, string msg)
{
    delete mesh;
    string msgPrefix("loadMesh: NetCDF operation failed - ");
    throw escript::IOError(msgPrefix+msg);
}

Domain_ptr loadMesh(const std::string& fileName)
{
#ifdef ESYS_HAVE_NETCDF
    escript::JMPI mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
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
        throw escript::IOError(msg.str());
    }

    // Read NetCDF integer attributes

    // index_size was only introduced with 64-bit index support so fall back
    // to 32 bits if not found.
    int index_size;
    try {
        index_size = ncReadAtt<int>(&dataFile, fName, "index_size");
    } catch (escript::IOError& e) {
        index_size = 4;
    }
    // technically we could cast if reading 32-bit data on 64-bit escript
    // but cost-benefit analysis clearly favours this implementation for now
    if (sizeof(index_t) != index_size) {
        throw escript::IOError("loadMesh: size of index types at runtime differ from dump file");
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
        throw escript::IOError(msg.str());
    }
    boost::scoped_array<char> name(attr->as_string(0));
    delete attr;

    // allocate mesh
    Mesh* mesh = new Mesh(name.get(), numDim, mpiInfo);

    // read nodes
    mesh->Nodes->allocTable(numNodes);
    // Nodes_Id
    if (! ( nc_var_temp = dataFile.get_var("Nodes_Id")) )
        cleanupAndThrow(mesh, "get_var(Nodes_Id)");
    if (! nc_var_temp->get(&mesh->Nodes->Id[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_Id)");
    // Nodes_Tag
    if (! ( nc_var_temp = dataFile.get_var("Nodes_Tag")) )
        cleanupAndThrow(mesh, "get_var(Nodes_Tag)");
    if (! nc_var_temp->get(&mesh->Nodes->Tag[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_Tag)");
    // Nodes_gDOF
    if (! ( nc_var_temp = dataFile.get_var("Nodes_gDOF")) )
        cleanupAndThrow(mesh, "get_var(Nodes_gDOF)");
    if (! nc_var_temp->get(&mesh->Nodes->globalDegreesOfFreedom[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_gDOF)");
    // Nodes_gNI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_gNI")) )
        cleanupAndThrow(mesh, "get_var(Nodes_gNI)");
    if (! nc_var_temp->get(&mesh->Nodes->globalNodesIndex[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_gNI)");
    // Nodes_grDfI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_grDfI")) )
        cleanupAndThrow(mesh, "get_var(Nodes_grDfI)");
    if (! nc_var_temp->get(&mesh->Nodes->globalReducedDOFIndex[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_grDfI)");
    // Nodes_grNI
    if (! ( nc_var_temp = dataFile.get_var("Nodes_grNI")) )
        cleanupAndThrow(mesh, "get_var(Nodes_grNI)");
    if (! nc_var_temp->get(&mesh->Nodes->globalReducedNodesIndex[0], numNodes) )
        cleanupAndThrow(mesh, "get(Nodes_grNI)");
    // Nodes_Coordinates
    if (!(nc_var_temp = dataFile.get_var("Nodes_Coordinates")))
        cleanupAndThrow(mesh, "get_var(Nodes_Coordinates)");
    if (! nc_var_temp->get(&mesh->Nodes->Coordinates[0], numNodes, numDim) )
        cleanupAndThrow(mesh, "get(Nodes_Coordinates)");

    mesh->Nodes->updateTagList();

    // read elements
    const_ReferenceElementSet_ptr refElements(new ReferenceElementSet(
                (ElementTypeId)Elements_TypeId, order, reduced_order));
    mesh->Elements=new ElementFile(refElements, mpiInfo);
    mesh->Elements->allocTable(num_Elements);
    mesh->Elements->minColor = 0;
    mesh->Elements->maxColor = num_Elements-1;
    if (num_Elements > 0) {
       // Elements_Id
       if (! ( nc_var_temp = dataFile.get_var("Elements_Id")) )
           cleanupAndThrow(mesh, "get_var(Elements_Id)");
       if (! nc_var_temp->get(&mesh->Elements->Id[0], num_Elements) )
           cleanupAndThrow(mesh, "get(Elements_Id)");
       // Elements_Tag
       if (! ( nc_var_temp = dataFile.get_var("Elements_Tag")) )
           cleanupAndThrow(mesh, "get_var(Elements_Tag)");
       if (! nc_var_temp->get(&mesh->Elements->Tag[0], num_Elements) )
           cleanupAndThrow(mesh, "get(Elements_Tag)");
       // Elements_Owner
       if (! ( nc_var_temp = dataFile.get_var("Elements_Owner")) )
           cleanupAndThrow(mesh, "get_var(Elements_Owner)");
       if (! nc_var_temp->get(&mesh->Elements->Owner[0], num_Elements) )
           cleanupAndThrow(mesh, "get(Elements_Owner)");
       // Elements_Color
       if (! ( nc_var_temp = dataFile.get_var("Elements_Color")) )
           cleanupAndThrow(mesh, "get_var(Elements_Color)");
       if (! nc_var_temp->get(&mesh->Elements->Color[0], num_Elements) )
           cleanupAndThrow(mesh, "get(Elements_Color)");
       // Now we need to adjust maxColor
       index_t mc = mesh->Elements->Color[0];
       for (index_t i = 1; i < num_Elements; ++i) {
           if (mc < mesh->Elements->Color[i]) {
               mc = mesh->Elements->Color[i];
           }
       }
       mesh->Elements->maxColor = mc;
       // Elements_Nodes
       int* Elements_Nodes = new int[num_Elements*num_Elements_numNodes];
       if (!(nc_var_temp = dataFile.get_var("Elements_Nodes"))) {
           delete[] Elements_Nodes;
           cleanupAndThrow(mesh, "get_var(Elements_Nodes)");
       }
       if (! nc_var_temp->get(&Elements_Nodes[0], num_Elements, num_Elements_numNodes) ) {
           delete[] Elements_Nodes;
           cleanupAndThrow(mesh, "get(Elements_Nodes)");
       }

       // Copy temp array into mesh->Elements->Nodes
       for (index_t i = 0; i < num_Elements; i++) {
           for (int j = 0; j < num_Elements_numNodes; j++) {
               mesh->Elements->Nodes[INDEX2(j,i,num_Elements_numNodes)]
                    = Elements_Nodes[INDEX2(j,i,num_Elements_numNodes)];
           }
       }
       delete[] Elements_Nodes;
    } // num_Elements > 0
    mesh->Elements->updateTagList();

    // get the face elements
    const_ReferenceElementSet_ptr refFaceElements(
            new ReferenceElementSet((ElementTypeId)FaceElements_TypeId,
                order, reduced_order));
    mesh->FaceElements = new ElementFile(refFaceElements, mpiInfo);
    mesh->FaceElements->allocTable(num_FaceElements);
    mesh->FaceElements->minColor = 0;
    mesh->FaceElements->maxColor = num_FaceElements-1;
    if (num_FaceElements > 0) {
        // FaceElements_Id
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Id")) )
            cleanupAndThrow(mesh, "get_var(FaceElements_Id)");
        if (! nc_var_temp->get(&mesh->FaceElements->Id[0], num_FaceElements) )
            cleanupAndThrow(mesh, "get(FaceElements_Id)");
        // FaceElements_Tag
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Tag")) )
            cleanupAndThrow(mesh, "get_var(FaceElements_Tag)");
        if (! nc_var_temp->get(&mesh->FaceElements->Tag[0], num_FaceElements) )
            cleanupAndThrow(mesh, "get(FaceElements_Tag)");
        // FaceElements_Owner
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Owner")) )
            cleanupAndThrow(mesh, "get_var(FaceElements_Owner)");
        if (! nc_var_temp->get(&mesh->FaceElements->Owner[0], num_FaceElements) )
            cleanupAndThrow(mesh, "get(FaceElements_Owner)");
        // FaceElements_Color
        if (! ( nc_var_temp = dataFile.get_var("FaceElements_Color")) )
            cleanupAndThrow(mesh, "get_var(FaceElements_Color)");
        if (! nc_var_temp->get(&mesh->FaceElements->Color[0], num_FaceElements) )
            cleanupAndThrow(mesh, "get(FaceElements_Color)");
        // Now we need to adjust maxColor
        index_t mc = mesh->FaceElements->Color[0];
        for (index_t i = 1; i < num_FaceElements; ++i) {
            if (mc < mesh->FaceElements->Color[i]) {
                mc = mesh->FaceElements->Color[i];
            }
        }
        mesh->FaceElements->maxColor = mc;
        // FaceElements_Nodes
        int* FaceElements_Nodes = new int[num_FaceElements*num_FaceElements_numNodes];
        if (!(nc_var_temp = dataFile.get_var("FaceElements_Nodes"))) {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(mesh, "get_var(FaceElements_Nodes)");
        }
        if (! nc_var_temp->get(&(FaceElements_Nodes[0]), num_FaceElements, num_FaceElements_numNodes) ) {
            delete[] FaceElements_Nodes;
            cleanupAndThrow(mesh, "get(FaceElements_Nodes)");
        }
        // Copy temp array into mesh->FaceElements->Nodes
        for (index_t i = 0; i < num_FaceElements; i++) {
            for (int j = 0; j < num_FaceElements_numNodes; j++) {
                mesh->FaceElements->Nodes[INDEX2(j,i,num_FaceElements_numNodes)] = FaceElements_Nodes[INDEX2(j,i,num_FaceElements_numNodes)];
            }
        }
        delete[] FaceElements_Nodes;
    } // num_FaceElements > 0
    mesh->FaceElements->updateTagList();

    // get the Contact elements
    const_ReferenceElementSet_ptr refContactElements(
         new ReferenceElementSet((ElementTypeId)ContactElements_TypeId,
             order, reduced_order));
    mesh->ContactElements = new ElementFile(refContactElements, mpiInfo);
    mesh->ContactElements->allocTable(num_ContactElements);
    mesh->ContactElements->minColor = 0;
    mesh->ContactElements->maxColor = num_ContactElements-1;
    if (num_ContactElements > 0) {
        // ContactElements_Id
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Id")) )
            cleanupAndThrow(mesh, "get_var(ContactElements_Id)");
        if (! nc_var_temp->get(&mesh->ContactElements->Id[0], num_ContactElements) )
            cleanupAndThrow(mesh, "get(ContactElements_Id)");
        // ContactElements_Tag
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Tag")) )
            cleanupAndThrow(mesh, "get_var(ContactElements_Tag)");
        if (! nc_var_temp->get(&mesh->ContactElements->Tag[0], num_ContactElements) )
            cleanupAndThrow(mesh, "get(ContactElements_Tag)");
        // ContactElements_Owner
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Owner")) )
            cleanupAndThrow(mesh, "get_var(ContactElements_Owner)");
        if (! nc_var_temp->get(&mesh->ContactElements->Owner[0], num_ContactElements) )
            cleanupAndThrow(mesh, "get(ContactElements_Owner)");
        // ContactElements_Color
        if (! ( nc_var_temp = dataFile.get_var("ContactElements_Color")) )
            cleanupAndThrow(mesh, "get_var(ContactElements_Color)");
        if (! nc_var_temp->get(&mesh->ContactElements->Color[0], num_ContactElements) )
            cleanupAndThrow(mesh, "get(ContactElements_Color)");
        // Now we need to adjust maxColor
        index_t mc = mesh->ContactElements->Color[0];
        for (index_t i = 1; i < num_ContactElements; ++i) {
            if (mc < mesh->ContactElements->Color[i]) {
                mc = mesh->ContactElements->Color[i];
            }
        }
        mesh->ContactElements->maxColor = mc;
        // ContactElements_Nodes
        int* ContactElements_Nodes = new int[num_ContactElements*num_ContactElements_numNodes];
        if (!(nc_var_temp = dataFile.get_var("ContactElements_Nodes"))) {
            delete[] ContactElements_Nodes;
            cleanupAndThrow(mesh, "get_var(ContactElements_Nodes)");
        }
        if (! nc_var_temp->get(&ContactElements_Nodes[0], num_ContactElements, num_ContactElements_numNodes) ) {
            delete[] ContactElements_Nodes;
            cleanupAndThrow(mesh, "get(ContactElements_Nodes)");
        }
        // Copy temp array into mesh->ContactElements->Nodes
        for (index_t i = 0; i < num_ContactElements; i++) {
            for (int j = 0; j < num_ContactElements_numNodes; j++) {
                mesh->ContactElements->Nodes[INDEX2(j,i,num_ContactElements_numNodes)] = ContactElements_Nodes[INDEX2(j,i,num_ContactElements_numNodes)];
            }
        }
        delete[] ContactElements_Nodes;
    } // num_ContactElements > 0
    mesh->ContactElements->updateTagList();

    // get the Points (nodal elements)
    const_ReferenceElementSet_ptr refPoints(new ReferenceElementSet(
                (ElementTypeId)Points_TypeId, order, reduced_order));
    mesh->Points = new ElementFile(refPoints, mpiInfo);
    mesh->Points->allocTable(num_Points);
    mesh->Points->minColor = 0;
    mesh->Points->maxColor = num_Points-1;
    if (num_Points > 0) {
        // Points_Id
        if (! ( nc_var_temp = dataFile.get_var("Points_Id")))
            cleanupAndThrow(mesh, "get_var(Points_Id)");
        if (! nc_var_temp->get(&mesh->Points->Id[0], num_Points))
            cleanupAndThrow(mesh, "get(Points_Id)");
        // Points_Tag
        if (! ( nc_var_temp = dataFile.get_var("Points_Tag")))
            cleanupAndThrow(mesh, "get_var(Points_Tag)");
        if (! nc_var_temp->get(&mesh->Points->Tag[0], num_Points))
            cleanupAndThrow(mesh, "get(Points_Tag)");
        // Points_Owner
        if (! ( nc_var_temp = dataFile.get_var("Points_Owner")))
            cleanupAndThrow(mesh, "get_var(Points_Owner)");
        if (!nc_var_temp->get(&mesh->Points->Owner[0], num_Points))
            cleanupAndThrow(mesh, "get(Points_Owner)");
        // Points_Color
        if (! ( nc_var_temp = dataFile.get_var("Points_Color")))
            cleanupAndThrow(mesh, "get_var(Points_Color)");
        if (!nc_var_temp->get(&mesh->Points->Color[0], num_Points))
            cleanupAndThrow(mesh, "get(Points_Color)");
        // Now we need to adjust maxColor
        index_t mc = mesh->Points->Color[0];
        for (index_t i = 1; i < num_Points; ++i) {
            if (mc < mesh->Points->Color[i]) {
                mc = mesh->Points->Color[i];
            }
        }
        mesh->Points->maxColor = mc;
        // Points_Nodes
        int* Points_Nodes = new int[num_Points];
        if (!(nc_var_temp = dataFile.get_var("Points_Nodes"))) {
            delete[] Points_Nodes;
            cleanupAndThrow(mesh, "get_var(Points_Nodes)");
        }
        if (! nc_var_temp->get(&Points_Nodes[0], num_Points) ) {
            delete[] Points_Nodes;
            cleanupAndThrow(mesh, "get(Points_Nodes)");
        }
        // Copy temp array into mesh->Points->Nodes
        for (index_t i = 0; i < num_Points; i++) {
            mesh->Points->Id[mesh->Points->Nodes[INDEX2(0,i,1)]] = Points_Nodes[i];
        }
        delete[] Points_Nodes;
    } // num_Points > 0
    mesh->Points->updateTagList();

    // get the tags
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        int *Tags_keys = new int[num_Tags];
        char name_temp[4096];
        int i;

        // Tags_keys
        if (! ( nc_var_temp = dataFile.get_var("Tags_keys")) ) {
            delete[] Tags_keys;
            cleanupAndThrow(mesh, "get_var(Tags_keys)");
        }
        if (! nc_var_temp->get(&Tags_keys[0], num_Tags) ) {
            delete[] Tags_keys;
            cleanupAndThrow(mesh, "get(Tags_keys)");
        }
        for (i=0; i<num_Tags; i++) {
          // Retrieve tag name
          sprintf(name_temp, "Tags_name_%d", i);
          if (! (attr=dataFile.get_att(name_temp)) ) {
              delete[] Tags_keys;
              stringstream msg;
              msg << "get_att(" << name_temp << ")";
              cleanupAndThrow(mesh, msg.str());
          }
          boost::scoped_array<char> name(attr->as_string(0));
          delete attr;
          mesh->addTagMap(name.get(), Tags_keys[i]);
        }
        delete[] Tags_keys;
    }

    // Nodes_DofDistribution
    std::vector<index_t> first_DofComponent(mpi_size+1);
    if (! (nc_var_temp = dataFile.get_var("Nodes_DofDistribution")) ) {
        cleanupAndThrow(mesh, "get_var(Nodes_DofDistribution)");
    }
    if (!nc_var_temp->get(&first_DofComponent[0], mpi_size+1)) {
        cleanupAndThrow(mesh, "get(Nodes_DofDistribution)");
    }

    // Nodes_NodeDistribution
    std::vector<index_t> first_NodeComponent(mpi_size+1);
    if (! (nc_var_temp = dataFile.get_var("Nodes_NodeDistribution")) ) {
        cleanupAndThrow(mesh, "get_var(Nodes_NodeDistribution)");
    }
    if (!nc_var_temp->get(&first_NodeComponent[0], mpi_size+1)) {
        cleanupAndThrow(mesh, "get(Nodes_NodeDistribution)");
    }
    mesh->createMappings(first_DofComponent, first_NodeComponent);

    AbstractContinuousDomain* dom(new MeshAdapter(mesh));
    return dom->getPtr();
#else
    throw FinleyException("loadMesh: not compiled with NetCDF. Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}

Domain_ptr readMesh(escript::JMPI info, const std::string& fileName,
                    int integrationOrder, int reducedIntegrationOrder,
                    bool optimize, const std::vector<double>& points,
                    const std::vector<int>& tags)
{
    if (fileName.size() == 0 )
        throw escript::ValueError("Null file name!");

    Mesh* fMesh = Mesh::read(info, fileName, integrationOrder, reducedIntegrationOrder, optimize);
    MeshAdapter* ma = new MeshAdapter(fMesh);
    ma->addDiracPoints(points, tags);
    return Domain_ptr(ma);
}

Domain_ptr readMesh_driver(const boost::python::list& args)
{
    using boost::python::extract;
    int l = len(args);
    if (l < 7) {
        throw escript::ValueError("Insufficient arguments to readMesh_driver");
    }
    string fileName = extract<string>(args[0])();
    int integrationOrder = extract<int>(args[1])();
    int reducedIntegrationOrder = extract<int>(args[2])();
    bool optimize = extract<bool>(args[3])();
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    boost::python::list pypoints = extract<boost::python::list>(args[4]);
    boost::python::list pytags = extract<boost::python::list>(args[5]);
    int numpts = extract<int>(pypoints.attr("__len__")());
    int numtags = extract<int>(pytags.attr("__len__")());

    boost::python::object pworld = args[6];
    escript::JMPI info;
    if (!pworld.is_none()) {
        extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw escript::ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = escript::makeInfo(MPI_COMM_WORLD);
    }
    Domain_ptr result = readMesh(info, fileName, integrationOrder,
                                 reducedIntegrationOrder, optimize, points, tags);

    for (int i = 0; i < numpts; ++i) {
        boost::python::object temp = pypoints[i];
        int l = extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
              points.push_back(extract<double>(temp[k]));
        }
    }
    // bricks use up to 200 but the existing tag check will find that
    int curmax = 40;
    TagMap& tagmap = dynamic_cast<MeshAdapter*>(result.get())->getMesh()->tagMap;
    // first we work out what tags are already in use
    for (TagMap::iterator it = tagmap.begin(); it != tagmap.end(); ++it) {
        if (it->second > curmax) {
            curmax = it->second+1;
        }
    }

    tags.resize(numtags, -1);
    for (int i = 0; i < numtags; ++i) {
        extract<int> ex_int(pytags[i]);
        extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i] = ex_int();
            if (tags[i] >= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::iterator it = tagmap.find(s);
            if (it != tagmap.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                result->setTagMap(s,curmax);
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Unable to extract tag value.");
        }
    }
    // now we need to add the dirac points
    dynamic_cast<MeshAdapter*>(result.get())->addDiracPoints(points, tags);
    return result;
}

Domain_ptr readGmsh(escript::JMPI info, const std::string& fileName,
                    int numDim, int integrationOrder,
                    int reducedIntegrationOrder, bool optimize,
                    bool useMacroElements, const std::vector<double>& points,
                    const std::vector<int>& tags)
{
    if (fileName.size() == 0 )
        throw escript::ValueError("Null file name!");

    Mesh* fMesh = Mesh::readGmsh(info, fileName, numDim, integrationOrder, reducedIntegrationOrder, optimize, useMacroElements);
    MeshAdapter* ma = new MeshAdapter(fMesh);
    ma->addDiracPoints(points, tags);
    return Domain_ptr(ma);
}

Domain_ptr readGmsh_driver(const boost::python::list& args)
{
    using boost::python::extract;
    int l=len(args);
    if (l<7) {
        throw escript::ValueError("Insufficient arguments to readMesh_driver");
    }
    string fileName=extract<string>(args[0])();
    int numDim=extract<int>(args[1])();
    int integrationOrder=extract<int>(args[2])();
    int reducedIntegrationOrder=extract<int>(args[3])();
    bool optimize=extract<bool>(args[4])();
    bool useMacroElements=extract<bool>(args[5])();
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    boost::python::list pypoints=extract<boost::python::list>(args[6]);
    boost::python::list pytags=extract<boost::python::list>(args[7]);
    int numpts=extract<int>(pypoints.attr("__len__")());
    int numtags=extract<int>(pytags.attr("__len__")());
    boost::python::object pworld=args[8];
    escript::JMPI info;
    if (!pworld.is_none()) {
        extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw escript::ValueError("Invalid escriptWorld parameter.");
        }
        info=ex()->getMPI();
    } else {
        info=escript::makeInfo(MPI_COMM_WORLD);
    }
    Domain_ptr result = readGmsh(info, fileName, numDim, integrationOrder,
                                 reducedIntegrationOrder, optimize,
                                 useMacroElements, points, tags);

    for (int i=0;i<numpts;++i) {
        boost::python::object temp=pypoints[i];
        int l=extract<int>(temp.attr("__len__")());
        for (int k=0;k<l;++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    int curmax=40; // bricks use up to 30
    TagMap& tagmap=dynamic_cast<MeshAdapter*>(result.get())->getMesh()->tagMap;
    // first we work out what tags are already in use
    for (TagMap::iterator it=tagmap.begin(); it!=tagmap.end(); ++it) {
        if (it->second>curmax) {
            curmax=it->second+1;
        }
    }

    tags.resize(numtags, -1);
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i]>= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s=ex_str();
            map<string, int>::iterator it=tagmap.find(s);
            if (it!=tagmap.end()) {
                // we have the tag already so look it up
                tags[i]=it->second;
            } else {
                result->setTagMap(s,curmax);
                tags[i]=curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Error - Unable to extract tag value");
        }
    }
    // now we need to add the dirac points
    dynamic_cast<MeshAdapter*>(result.get())->addDiracPoints(points, tags);
    return result;
}

Domain_ptr brick(escript::JMPI info, dim_t n0, dim_t n1, dim_t n2, int order,
                 double l0, double l1, double l2,
                 bool periodic0, bool periodic1, bool periodic2,
                 int integrationOrder, int reducedIntegrationOrder,
                 bool useElementsOnFace, bool useFullElementOrder,
                 bool optimize, const std::vector<double>& points,
                 const std::vector<int>& tags,
                 const std::map<std::string, int>& tagNamesToNums)
{
    const dim_t numElements[] = {n0, n1, n2};
    const double length[] = {l0, l1, l2};
    const bool periodic[] = {periodic0, periodic1, periodic2};

    Mesh* fMesh = NULL;
    if (order==1) {
        fMesh=RectangularMesh_Hex8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, optimize, info);
    } else if (order==2) {
        fMesh=RectangularMesh_Hex20(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, false, optimize, info);
    } else if (order==-1) {
        fMesh=RectangularMesh_Hex20(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, true, optimize, info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw escript::ValueError(message.str());
    }

    MeshAdapter* dom = new MeshAdapter(fMesh);
    dom->addDiracPoints(points, tags);
    Mesh* out = dom->getMesh();
    for (TagMap::const_iterator it = tagNamesToNums.begin(); it != tagNamesToNums.end(); ++it) {
        out->addTagMap(it->first, it->second);
    }
    out->Points->updateTagList();
    return Domain_ptr(dom);
}

Domain_ptr brick_driver(const boost::python::list& args)
{
    using boost::python::extract;

    // we need to convert lists to stl vectors
    boost::python::list pypoints = extract<boost::python::list>(args[15]);
    boost::python::list pytags = extract<boost::python::list>(args[16]);
    int numpts = extract<int>(pypoints.attr("__len__")());
    int numtags = extract<int>(pytags.attr("__len__")());
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);
    for (int i = 0; i < numpts; ++i) {
        boost::python::object temp = pypoints[i];
        int l = extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    map<string, int> namestonums;
    int curmax=40; // bricks use up to 30
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i]>= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            string s = ex_str();
            TagMap::iterator it=namestonums.find(s);
            if (it != namestonums.end()) {
                // we have the tag already so look it up
                tags[i] = it->second;
            } else {
                namestonums[s] = curmax;
                tags[i] = curmax;
                curmax++;
            }
        } else {
            throw FinleyException("Error - Unable to extract tag value.");
        }
    }
    boost::python::object pworld=args[17];
    escript::JMPI info;
    if (!pworld.is_none()) {
        extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw escript::ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = escript::makeInfo(MPI_COMM_WORLD);
    }
    return brick(info, static_cast<dim_t>(extract<float>(args[0])),
                   static_cast<dim_t>(extract<float>(args[1])),
                   static_cast<dim_t>(extract<float>(args[2])),
                   extract<int>(args[3]), extract<double>(args[4]),
                   extract<double>(args[5]), extract<double>(args[6]),
                   extract<int>(args[7]), extract<int>(args[8]),
                   extract<int>(args[9]), extract<int>(args[10]),
                   extract<int>(args[11]), extract<int>(args[12]),
                   extract<int>(args[13]), extract<int>(args[14]),
                   points, tags, namestonums);
}

Domain_ptr rectangle(escript::JMPI info, dim_t n0, dim_t n1, int order,
                     double l0, double l1, bool periodic0, bool periodic1,
                     int integrationOrder, int reducedIntegrationOrder,
                     bool useElementsOnFace, bool useFullElementOrder,
                     bool optimize, const vector<double>& points,
                     const vector<int>& tags,
                     const std::map<std::string, int>& tagNamesToNums)
{
    const dim_t numElements[] = {n0, n1};
    const double length[] = {l0, l1};
    const bool periodic[] = {periodic0, periodic1};

    Mesh* fMesh = NULL;
    if (order==1) {
        fMesh=RectangularMesh_Rec4(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, optimize, info);
    } else if (order==2) {
        fMesh=RectangularMesh_Rec8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace,useFullElementOrder, false, optimize, info);
    } else if (order==-1) {
        fMesh=RectangularMesh_Rec8(numElements, length, periodic,
                integrationOrder, reducedIntegrationOrder,
                useElementsOnFace, useFullElementOrder, true, optimize, info);
    } else {
        stringstream message;
        message << "Illegal interpolation order " << order;
        throw escript::ValueError(message.str());
    }

    MeshAdapter* dom = new MeshAdapter(fMesh);
    dom->addDiracPoints(points, tags);
    Mesh* out = dom->getMesh();
    for (TagMap::const_iterator it = tagNamesToNums.begin(); it != tagNamesToNums.end(); ++it)
    {
        out->addTagMap(it->first, it->second);
    }
    out->Points->updateTagList();
    return Domain_ptr(dom);
}

Domain_ptr meshMerge(const boost::python::list& meshList)
{
    // extract the meshes from meshList
    int num = boost::python::extract<int>(meshList.attr("__len__")());
    vector<Mesh*> meshes(num);
    for (int i = 0; i < num; ++i) {
        AbstractContinuousDomain& meshListMember = boost::python::extract<AbstractContinuousDomain&>(meshList[i]);
        const MeshAdapter* finley_meshListMember = static_cast<const MeshAdapter*>(&meshListMember);
        meshes[i] = finley_meshListMember->getMesh();
    }

    // merge the meshes
    Mesh* fMesh = Mesh_merge(meshes);

    return Domain_ptr(new MeshAdapter(fMesh));
}

Domain_ptr rectangle_driver(const boost::python::list& args)
{
    using boost::python::extract;

    // we need to convert lists to stl vectors
    boost::python::list pypoints=extract<boost::python::list>(args[12]);
    boost::python::list pytags=extract<boost::python::list>(args[13]);
    int numpts = extract<int>(pypoints.attr("__len__")());
    int numtags = extract<int>(pytags.attr("__len__")());
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);
    for (int i = 0; i < numpts; ++i) {
        boost::python::object temp = pypoints[i];
        int l = extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    TagMap tagstonames;
    int curmax = 40;
    // but which order to assign tags to names?????
    for (int i = 0; i < numtags; ++i) {
        extract<int> ex_int(pytags[i]);
        extract<string> ex_str(pytags[i]);
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
    boost::python::object pworld = args[14];
    escript::JMPI info;
    if (!pworld.is_none()) {
        extract<SubWorld_ptr> ex(pworld);
        if (!ex.check()) {
            throw escript::ValueError("Invalid escriptWorld parameter.");
        }
        info = ex()->getMPI();
    } else {
        info = escript::makeInfo(MPI_COMM_WORLD);
    }

    return rectangle(info, static_cast<dim_t>(extract<float>(args[0])),
                       static_cast<dim_t>(extract<float>(args[1])),
                       extract<int>(args[2]), extract<double>(args[3]),
                       extract<double>(args[4]), extract<int>(args[5]),
                       extract<int>(args[6]), extract<int>(args[7]),
                       extract<int>(args[8]), extract<int>(args[9]),
                       extract<int>(args[10]), extract<int>(args[11]),
                       points, tags, tagstonames);
}

Domain_ptr glueFaces(const boost::python::list& meshList, double safety_factor,
                     double tolerance, bool optimize)
{
    // merge the meshes:
    Domain_ptr merged_meshes = meshMerge(meshList);

    // glue the faces:
    const MeshAdapter* merged_finley_meshes = dynamic_cast<const MeshAdapter*>(merged_meshes.get());
    Mesh* fMesh = merged_finley_meshes->getMesh();
    fMesh->glueFaces(safety_factor, tolerance, optimize);

    return merged_meshes;
}

Domain_ptr joinFaces(const boost::python::list& meshList, double safety_factor,
                     double tolerance, bool optimize)
{
    // merge the meshes:
    Domain_ptr merged_meshes = meshMerge(meshList);

    // join the faces:
    const MeshAdapter* merged_finley_meshes=static_cast<const MeshAdapter*>(merged_meshes.get());
    Mesh* fMesh=merged_finley_meshes->getMesh();
    fMesh->joinFaces(safety_factor, tolerance, optimize);

    return merged_meshes;
}

} // end of namespace

