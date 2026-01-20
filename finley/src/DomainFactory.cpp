
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

#include <finley/DomainFactory.h>

#include <escript/index.h>

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <sstream>

using namespace std;
using namespace escript;


namespace finley {

Domain_ptr FinleyDomain::load(const string& fileName, const bp::object& comm)
{
    int error = 0;
    std::string msg;
#ifdef ESYS_HAVE_HDF5
    #ifdef ESYS_INDEXTYPE_LONG
        H5::DataType h5_type_index = H5::PredType::NATIVE_LONG;
    #else
        H5::DataType h5_type_index = H5::PredType::NATIVE_INT;
    #endif
    JMPI mpiInfo = makeInfoFromPyComm(comm);
    const string fName(mpiInfo->appendRankToFileName(fileName));
    try
    {
            string name;
            int numDim = -1;
            dim_t numNodes = -1;
            // ... Open file ....
            H5::H5File h5_file(fName, H5F_ACC_RDONLY);
            // .... read meta data ...
            H5::Group h5_grp_meta=h5_file.openGroup("Meta");
            // .... get MPI information ...........................
            long long h5_values_mpi[2];
            H5::Attribute h5_attr_mpi(h5_grp_meta.openAttribute("mpi"));
            H5::DataType h5_type_mpi(h5_attr_mpi.getDataType());
            if ( h5_type_mpi != H5::PredType::NATIVE_LLONG ) {
                 throw FinleyException("Error - finley.load: illegal data type for MPI informations in HDF5 file.");
            }
            if ( h5_attr_mpi.getStorageSize() != 2 * h5_type_mpi.getSize() ) {
                 throw FinleyException("Error - finley.load: MPI information in HDF5 file needs to be a two value.");
            }
            h5_attr_mpi.read(h5_type_mpi, &h5_values_mpi[0]);
            // Verify size and rank
            if (mpiInfo->size != h5_values_mpi[0]) {
                stringstream msg;
                msg << "finley.load: The HDF5 file '" << fName << "' can only be read on " << h5_values_mpi[0] << " MPI ranks. Currently running: " << mpiInfo->size;
                throw FinleyException(msg.str());
            }
            if (mpiInfo->rank != h5_values_mpi[1]) {
                stringstream msg;
                msg << "finley.load: The HDF5 file '" << fName << "' should be read on MPI rank #" << h5_values_mpi[1] << " and NOT on #" << mpiInfo->rank;
                throw FinleyException(msg.str());
            }
            // ... retrieve name ....
            H5::Attribute h5_attr_name(h5_grp_meta.openAttribute("name"));
            H5::DataType h5_type_name(h5_attr_name.getDataType());
            if (  h5_type_name != H5::StrType(0, H5T_VARIABLE)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for name in HDF5 file.");
            }
            //DataSpace dataspace(H5S_SCALAR);
            h5_attr_name.read(h5_type_name, name);

            H5::Group h5_grp_nodes = h5_file.openGroup("Nodes");
            // ... retrieve dimension ....
            uint h5_numDim;
            H5::Attribute h5_attr_dim(h5_grp_nodes.openAttribute("numDim"));
            H5::DataType h5_type_dim(h5_attr_dim.getDataType());
            if (  h5_type_dim != H5::PredType::NATIVE_UINT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for numDim in HDF5 file.");
            }
            if ( h5_attr_dim.getStorageSize() !=  1 * h5_type_dim.getSize() ) {
                 throw FinleyException("Error - finley.load: numDim  in HDF5 file needs to be a single value.");
            }
            h5_attr_dim.read(h5_type_dim, &h5_numDim);
            numDim = h5_numDim;
            // allocate domain:
            FinleyDomain* dom = new FinleyDomain(name, numDim, mpiInfo);
            NodeFile* nodes = dom->getNodes();
            // .......... get nodes .................................................
            long h5_numNodes;
            H5::Attribute h5_attr_nn(h5_grp_nodes.openAttribute("numNodes"));
            H5::DataType h5_type_nn(h5_attr_nn.getDataType());
            if (  h5_type_nn != H5::PredType::NATIVE_LONG  ) {
                 throw FinleyException("Error - finley.load: illegal data type for numNodes in HDF5 file.");
            }
            if ( h5_attr_nn.getStorageSize() !=  1 * h5_type_nn.getSize() ) {
                 throw FinleyException("Error - finley.load: numNodes in HDF5 file needs to be a single value.");
            }
            h5_attr_nn.read(h5_type_dim, &h5_numNodes);
            numNodes = h5_numNodes;
            nodes->allocTable(numNodes);

            H5::DataSet h5_ds_nodes =h5_grp_nodes.openDataSet("Coordinates");
            H5::DataType h5_type_nodes(h5_ds_nodes.getDataType());
            if (  h5_type_nodes != H5::PredType::NATIVE_DOUBLE  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node coordinates in HDF5 file.");
            }
            if ( h5_ds_nodes.getStorageSize() !=  numNodes * numDim * h5_type_nodes.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node coordinates in HDF5 file is incorrect.");
            }
            h5_ds_nodes.read(&nodes->Coordinates[0], h5_type_nodes);

            H5::DataSet h5_ds_ids =h5_grp_nodes.openDataSet("Ids");
            H5::DataType h5_type_ids(h5_ds_ids.getDataType());
            if (  h5_type_ids != H5::PredType::NATIVE_INT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node ids in HDF5 file.");
            }
            if ( h5_ds_ids.getStorageSize() !=  numNodes * h5_type_ids.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node ids in HDF5 file is incorrect.");
            }
            h5_ds_ids.read(&nodes->Id[0], h5_type_ids);

            H5::DataSet h5_ds_tags =h5_grp_nodes.openDataSet("Tags");
            H5::DataType h5_type_tags(h5_ds_tags.getDataType());
            if (  h5_type_tags != H5::PredType::NATIVE_INT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node tags in HDF5 file.");
            }
            if ( h5_ds_tags.getStorageSize() !=  numNodes * h5_type_tags.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node tags in HDF5 file is incorrect.");
            }
            h5_ds_tags.read(&nodes->Tag[0], h5_type_tags);

            H5::DataSet h5_ds_gDOF =h5_grp_nodes.openDataSet("globalDegreesOfFreedom");
            H5::DataType h5_type_gDOF(h5_ds_gDOF.getDataType());
            if (  h5_type_gDOF != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node global DOF index in HDF5 file.");
            }
            if ( h5_ds_gDOF.getStorageSize() !=  numNodes * h5_type_gDOF.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node global DOF index in HDF5 file is incorrect.");
            }
            h5_ds_gDOF.read(&nodes->globalDegreesOfFreedom[0], h5_type_gDOF);

            H5::DataSet h5_ds_gNI =h5_grp_nodes.openDataSet("globalNodesIndex");
            H5::DataType h5_type_gNI(h5_ds_gNI.getDataType());
            if (  h5_type_gNI != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node index in HDF5 file.");
            }
            if ( h5_ds_gNI.getStorageSize() !=  numNodes * h5_type_gNI.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node index in HDF5 file is incorrect.");
            }
            h5_ds_gNI.read(&nodes->globalNodesIndex[0], h5_type_gNI);


            H5::DataSet h5_ds_gRDOFI =h5_grp_nodes.openDataSet("globalReducedDOFIndex");
            H5::DataType h5_type_gRDOFI(h5_ds_gRDOFI.getDataType());
            if (  h5_type_gRDOFI != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for global reduced DOF index in HDF5 file.");
            }
            if ( h5_ds_gRDOFI.getStorageSize() !=  numNodes * h5_type_gRDOFI.getSize() ) {
                 throw FinleyException("Error - finley.load: number of global reduces DOF index in HDF5 file is incorrect.");
            }
            h5_ds_gRDOFI.read(&nodes->globalReducedDOFIndex[0], h5_type_gRDOFI);

            H5::DataSet h5_ds_gRNI =h5_grp_nodes.openDataSet("globalReducedNodesIndex");
            H5::DataType h5_type_gRNI(h5_ds_gRNI.getDataType());
            if (  h5_type_gRNI != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for global reduced node index in HDF5 file.");
            }
            if ( h5_ds_gRNI.getStorageSize() !=  numNodes * h5_type_gRNI.getSize() ) {
                 throw FinleyException("Error - finley.load: number of global reduced node index in HDF5 file is incorrect.");
            }
            h5_ds_gRNI.read(&nodes->globalReducedNodesIndex[0], h5_type_gRNI);
            nodes->updateTagList();
            // ... end nodes ..................
            // .... read elements:
            int integrationOrder=-1, reducedIntegrationOrder=-1;
            H5::Attribute h5_attr_iO(h5_grp_meta.openAttribute("integrationOrder"));
            H5::DataType h5_type_iO(h5_attr_iO.getDataType());
            if (  h5_type_iO != H5::PredType::NATIVE_INT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for integrationOrder in HDF5 file.");
            }
            if ( h5_attr_iO.getStorageSize() !=  1 * h5_type_iO.getSize() ) {
                 throw FinleyException("Error - finley.load: integrationOrder  in HDF5 file needs to be a single value.");
            }
            h5_attr_iO.read(h5_type_iO, &integrationOrder);

            H5::Attribute h5_attr_riO(h5_grp_meta.openAttribute("reducedIntegrationOrder"));
            H5::DataType h5_type_riO(h5_attr_riO.getDataType());
            if (  h5_type_riO != H5::PredType::NATIVE_INT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for reducedIntegrationOrder in HDF5 file.");
            }
            if ( h5_attr_riO.getStorageSize() !=  1 * h5_type_riO.getSize() ) {
                 throw FinleyException("Error - finley.load: reducedIntegrationOrder  in HDF5 file needs to be a single value.");
            }
            h5_attr_riO.read(h5_type_riO, &reducedIntegrationOrder);

            // ... load elements ....
            dom->setElements(loadElements_hdf5(h5_file.openGroup("Elements"), integrationOrder, reducedIntegrationOrder, mpiInfo));
            dom->setFaceElements(loadElements_hdf5(h5_file.openGroup("FaceElements"),  integrationOrder, reducedIntegrationOrder, mpiInfo));
            dom->setContactElements(loadElements_hdf5(h5_file.openGroup("ContactElements"),  integrationOrder, reducedIntegrationOrder, mpiInfo));
            dom->setPoints(loadElements_hdf5(h5_file.openGroup("Points"),  integrationOrder, reducedIntegrationOrder, mpiInfo));
            // .... end read elements
            // ...  read tag map ..... :
            H5::Group h5_grp_tags=h5_file.openGroup("Tags");
            uint num_Tags=0;
            H5::Attribute h5_attr_ntags(h5_grp_tags.openAttribute("numTags"));
            H5::DataType h5_type_ntags(h5_attr_ntags.getDataType());
            if (  h5_type_ntags != H5::PredType::NATIVE_UINT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for number of tags in HDF5 file.");
            }
            if ( h5_attr_ntags.getStorageSize() !=  1 * h5_type_ntags.getSize() ) {
                 throw FinleyException("Error - finley.load: number of tags in HDF5 file needs to be a single value.");
            }
            h5_attr_ntags.read(h5_type_ntags, &num_Tags);

            int *Tag_keys = new int[num_Tags];
            H5::DataSet h5_ds_tagkys =h5_grp_tags.openDataSet("TagIds");
            H5::DataType h5_type_tagkys(h5_ds_tagkys.getDataType());
            if (  h5_type_tagkys != H5::PredType::NATIVE_INT  ) {
                 throw FinleyException("Error - finley.load: illegal data type for tag keys in HDF5 file.");
            }
            if ( h5_ds_tagkys.getStorageSize() !=  num_Tags * h5_type_tagkys.getSize() ) {
                 throw FinleyException("Error - finley.load: number of  tag keys  in HDF5 file is incorrect.");
            }
            h5_ds_tagkys.read(&Tag_keys[0], h5_type_tagkys);

            vector<const char*> Tag_names;
            Tag_names.resize(num_Tags);

            H5::DataSet h5_ds_tagnames =h5_grp_tags.openDataSet("TagNames");
           H5::DataType h5_type_tagnames(h5_ds_tagnames.getDataType());
           if (  h5_type_tagnames != H5::StrType(0, H5T_VARIABLE)  ) {
                         throw FinleyException("Error - finley.load: illegal data type for tage names in HDF5 file.");
           }
           h5_ds_tagnames.read(&Tag_names[0], h5_type_tagnames);

            for (uint i = 0; i < num_Tags; i++) {
                dom->setTagMap(Tag_names[i], Tag_keys[i]);
            }
            // ... node & DOF distributions:
            IndexVector first_DofComponent(mpiInfo->size + 1);
            IndexVector first_NodeComponent(mpiInfo->size + 1);
            H5::DataSet h5_ds_ndist =h5_grp_nodes.openDataSet("NodeDistribution");
            H5::DataType h5_type_ndist(h5_ds_ndist.getDataType());
            if (  h5_type_ndist != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for node distribution in HDF5 file.");
            }
            if ( h5_ds_ndist.getStorageSize() !=  ( mpiInfo->size + 1 ) * h5_type_ndist.getSize() ) {
                 throw FinleyException("Error - finley.load: number of node distribution entries in HDF5 file is incorrect.");
            }
            h5_ds_ndist.read(&first_DofComponent[0], h5_type_ndist);

            H5::DataSet h5_ds_dofdist =h5_grp_nodes.openDataSet("DoFDistribution");
            H5::DataType h5_type_dofdist(h5_ds_dofdist.getDataType());
            if (  h5_type_dofdist != H5::DataType(h5_type_index)  ) {
                 throw FinleyException("Error - finley.load: illegal data type for DOF distribution in HDF5 file.");
            }
            if ( h5_ds_dofdist.getStorageSize() !=  ( mpiInfo->size + 1 ) * h5_type_dofdist.getSize() ) {
                 throw FinleyException("Error - finley.load: number of DOF distribution entries in HDF5 file is incorrect.");
            }
            h5_ds_dofdist.read(&first_NodeComponent[0], h5_type_dofdist);
            dom->createMappings(first_DofComponent, first_NodeComponent);
            // ...  all done
            return dom->getPtr();
    }

    // catch failure caused by the H5File operations
    catch (H5::Exception& e)
    {
        error=1;
        e.printErrorStack();
        msg=e.getCDetailMsg();
    }
    catch (FinleyException& e) {
        error=1;
        msg=e.what();
    }
    int gerror = error;
    checkResult(error, gerror, mpiInfo);
    if (gerror > 0) {
        char* gmsg;
        shipString(msg.c_str(), &gmsg, mpiInfo->comm);
        //delete dom;
        throw FinleyException(gmsg);
    }
#else
    throw FinleyException("loadMesh: not compiled with HDF5. Please contact your installation manager.");
#endif // ESYS_HAVE_HDF5
    return NULL;
}

Domain_ptr readMesh_driver(const string& fileName,
                           int integrationOrder,
                           int reducedIntegrationOrder,
                           bool optimize,
                           const bp::list& diracPoints,
                           const bp::list& diracTags,
                           const bp::object& comm)
{
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    int numpts = bp::len(diracPoints);
    int numtags = bp::len(diracTags);

    // Handle optional MPI communicator (defaults to MPI_COMM_WORLD if None)
    JMPI info = makeInfoFromPyComm(comm);

    Domain_ptr dom(FinleyDomain::read(info, fileName, integrationOrder,
                                      reducedIntegrationOrder, optimize));

    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = diracPoints[i];
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
        bp::extract<int> ex_int(diracTags[i]);
        bp::extract<string> ex_str(diracTags[i]);
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

Domain_ptr readGmsh_driver(const string& fileName,
                            int numDim,
                            int integrationOrder,
                            int reducedIntegrationOrder,
                            bool optimize,
                            bool useMacroElements,
                            const bp::list& diracPoints,
                            const bp::list& diracTags,
                            const bp::object& comm)
{
    vector<double> points;
    vector<int> tags;

    // we need to convert lists to stl vectors
    int numpts = bp::len(diracPoints);
    int numtags = bp::len(diracTags);

    // Handle optional MPI communicator (defaults to MPI_COMM_WORLD if None)
    JMPI info = makeInfoFromPyComm(comm);

    Domain_ptr dom(FinleyDomain::readGmsh(info, fileName, numDim,
                                     integrationOrder, reducedIntegrationOrder,
                                     optimize, useMacroElements));
    FinleyDomain* fd = dynamic_cast<FinleyDomain*>(dom.get());

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = diracPoints[i];
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
        bp::extract<int> ex_int(diracTags[i]);
        bp::extract<string> ex_str(diracTags[i]);
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

Domain_ptr brick_driver(const bp::tuple& n_tuple,
                        int order,
                        const bp::tuple& l_tuple,
                        const bp::tuple& periodic_tuple,
                        int integrationOrder, int reducedIntegrationOrder,
                        bool useElementsOnFace, bool useFullElementOrder,
                        bool optimize,
                        const bp::list& diracPoints,
                        const bp::list& diracTags,
                        const bp::object& comm)
{
    // Extract dimensions from tuple
    dim_t n0 = bp::extract<dim_t>(n_tuple[0]);
    dim_t n1 = bp::extract<dim_t>(n_tuple[1]);
    dim_t n2 = bp::extract<dim_t>(n_tuple[2]);

    // Extract lengths from tuple
    double l0 = bp::extract<double>(l_tuple[0]);
    double l1 = bp::extract<double>(l_tuple[1]);
    double l2 = bp::extract<double>(l_tuple[2]);

    // Extract periodic flags from tuple
    bool periodic0 = bp::extract<bool>(periodic_tuple[0]);
    bool periodic1 = bp::extract<bool>(periodic_tuple[1]);
    bool periodic2 = bp::extract<bool>(periodic_tuple[2]);

    // Convert lists to stl vectors
    int numpts = bp::len(diracPoints);
    int numtags = bp::len(diracTags);
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = diracPoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(bp::extract<double>(temp[k]));
        }
    }

    map<string, int> namestonums;
    int curmax = 40; // bricks use up to 30
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(diracTags[i]);
        bp::extract<string> ex_str(diracTags[i]);
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

    // Handle optional MPI communicator (defaults to MPI_COMM_WORLD if None)
    JMPI info = makeInfoFromPyComm(comm);

    return brick(info, n0, n1, n2, order, l0, l1, l2,
                 periodic0, periodic1, periodic2,
                 integrationOrder, reducedIntegrationOrder,
                 useElementsOnFace, useFullElementOrder, optimize,
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

Domain_ptr rectangle_driver(const bp::tuple& n_tuple,
                            int order,
                            const bp::tuple& l_tuple,
                            const bp::tuple& periodic_tuple,
                            int integrationOrder, int reducedIntegrationOrder,
                            bool useElementsOnFace, bool useFullElementOrder,
                            bool optimize,
                            const bp::list& diracPoints,
                            const bp::list& diracTags,
                            const bp::object& comm)
{
    // Extract dimensions from tuple
    dim_t n0 = bp::extract<dim_t>(n_tuple[0]);
    dim_t n1 = bp::extract<dim_t>(n_tuple[1]);

    // Extract lengths from tuple
    double l0 = bp::extract<double>(l_tuple[0]);
    double l1 = bp::extract<double>(l_tuple[1]);

    // Extract periodic flags from tuple
    bool periodic0 = bp::extract<bool>(periodic_tuple[0]);
    bool periodic1 = bp::extract<bool>(periodic_tuple[1]);

    // Convert lists to stl vectors
    int numpts = bp::len(diracPoints);
    int numtags = bp::len(diracTags);
    vector<double> points;
    vector<int> tags;
    tags.resize(numtags, -1);

    for (int i = 0; i < numpts; ++i) {
        bp::object temp = diracPoints[i];
        int l = bp::extract<int>(temp.attr("__len__")());
        for (int k = 0; k < l; ++k) {
            points.push_back(bp::extract<double>(temp[k]));
        }
    }

    TagMap tagstonames;
    int curmax = 40;
    for (int i = 0; i < numtags; ++i) {
        bp::extract<int> ex_int(diracTags[i]);
        bp::extract<string> ex_str(diracTags[i]);
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

    // Handle optional MPI communicator (defaults to MPI_COMM_WORLD if None)
    JMPI info = makeInfoFromPyComm(comm);

    return rectangle(info, n0, n1, order, l0, l1,
                     periodic0, periodic1,
                     integrationOrder, reducedIntegrationOrder,
                     useElementsOnFace, useFullElementOrder, optimize,
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
