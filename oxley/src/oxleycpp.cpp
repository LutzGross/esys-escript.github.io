/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <oxley/Brick.h>
#include <oxley/FinleyConverter.h>
#include <oxley/Rectangle.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyDomain.h>
#include <oxley/RefinementQueue.h>

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#endif

// using namespace boost::python::numpy;

namespace oxley {

// converts the diagnostic's counts to a python list
boost::python::list _lnodesDegree2Report(const oxley::OxleyDomain& dom)
{
    const std::vector<long> r = oxley::lnodesDegree2Report(dom);
    boost::python::list out;
    for (size_t i = 0; i < r.size(); ++i)
        out.append(r[i]);
    return out;
}

// one python tuple per degree-2 lnodes slot: (element, slot, gid, x, y, z,
// faceCode, isCorner)
boost::python::list _lnodesDegree2Slots(const oxley::OxleyDomain& dom)
{
    const std::vector<oxley::SlotRecord> recs = oxley::lnodesDegree2Slots(dom);
    boost::python::list out;
    for (size_t i = 0; i < recs.size(); ++i) {
        const oxley::SlotRecord& r = recs[i];
        out.append(boost::python::make_tuple(r.element, r.slot, r.gid,
                                            r.x, r.y, r.z, r.faceCode, r.corner));
    }
    return out;
}

// Convert a Python refine_level (an int, or a flat sequence of ints of length
// n0*n1[*n2] in row-major block order) into the per-block vector the domain
// constructors expect. A scalar becomes a size-1 vector (uniform, broadcast in
// the constructor).
static std::vector<int> extractRefineLevel(const object& refine_level, dim_t numBlocks)
{
    std::vector<int> levels;
    extract<int> as_int(refine_level);
    if (as_int.check()) {
        levels.push_back(as_int());
        return levels;
    }
    extract<boost::python::list> as_list(refine_level);
    if (!as_list.check())
        throw OxleyException("refine_level must be an int or a list of ints");
    boost::python::list pylist = as_list();
    int n = extract<int>(pylist.attr("__len__")());
    if (n != (int) numBlocks)
        throw OxleyException("refine_level list length must equal the number of blocks (n0*n1[*n2])");
    levels.reserve(n);
    for (int i = 0; i < n; ++i)
        levels.push_back(extract<int>(pylist[i]));
    return levels;
}

escript::Domain_ptr _rectangle(double _n0, double _n1,
                        const object& l0, const object& l1,
                        const object& objpoints, const object& objtags,
                        const object& py_comm, const object& refine_level)
{
    // The assembler always uses a fixed 2-point Gauss rule, which is exact to
    // cubic, i.e. integration order 3. m_order records this actual order.
    const int order = 3;

    // Number of nodes in each direction
    dim_t n0=static_cast<dim_t>(_n0), n1=static_cast<dim_t>(_n1);
    double x0=0., x1=1., y0=0., y1=1.;

    // Length of the domain in each direction
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw OxleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw OxleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw OxleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw OxleyException("Argument l1 must be a float or 2-tuple");

    // process tags and points
    boost::python::list pypoints=extract<boost::python::list>(objpoints);
    boost::python::list pytags=extract<boost::python::list>(objtags);
    int numpts=extract<int>(pypoints.attr("__len__")());
    int numtags=extract<int>(pytags.attr("__len__")());
    std::vector<double> points;
    std::vector<int> tags;
    tags.resize(numtags, -1);
    for (int i=0;i<numpts;++i) {
        tuple temp = extract<tuple>(pypoints[i]);
        int l=extract<int>(temp.attr("__len__")());
        if (l != 2)
            throw OxleyException("Number of coordinates for each dirac point must match dimensions.");
        for (int k=0;k<l;++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    std::map<std::string, int> tagstonames;
    int curmax=40;
    // but which order to assign tags to names?????
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<std::string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i] >= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            std::string s=ex_str();
            std::map<std::string, int>::iterator it=tagstonames.find(s);
            if (it!=tagstonames.end()) {
                // we have the tag already so look it up
                tags[i]=it->second;
            } else {
                tagstonames[s]=curmax;
                tags[i]=curmax;
                curmax++;
            }
        } else {
            throw OxleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw OxleyException("Number of tags does not match number of points.");

    // Handle optional MPI communicator
    escript::JMPI jmpi = escript::makeInfoFromPyComm(py_comm);

    std::vector<int> levels = extractRefineLevel(refine_level, n0*n1);

    return escript::Domain_ptr(new Rectangle(jmpi, order, n0,n1, x0,y0, x1,y1,
                                points, tags, tagstonames, levels));
}


escript::Domain_ptr _brick(double _n0, double _n1, double _n2,
                        const object& l0, const object& l1, const object& l2,
                        const object& objpoints, const object& objtags,
                        const object& py_comm, const object& refine_level)
{
    // The assembler always uses a fixed 2-point Gauss rule, which is exact to
    // cubic, i.e. integration order 3. m_order records this actual order.
    const int order = 3;

    // Number of nodes in each direction
    dim_t n0=static_cast<dim_t>(_n0), n1=static_cast<dim_t>(_n1), n2=static_cast<dim_t>(_n2);;
    double x0=0., x1=1., y0=0., y1=1., z0=0., z1=1.;

    // Length of the domain in each direction
    if (extract<tuple>(l0).check()) {
        tuple x=extract<tuple>(l0);
        if (len(x)==2) {
            x0=extract<double>(x[0]);
            x1=extract<double>(x[1]);
        } else
            throw OxleyException("Argument l0 must be a float or 2-tuple");
    } else if (extract<double>(l0).check()) {
        x1=extract<double>(l0);
    } else
        throw OxleyException("Argument l0 must be a float or 2-tuple");

    if (extract<tuple>(l1).check()) {
        tuple y=extract<tuple>(l1);
        if (len(y)==2) {
            y0=extract<double>(y[0]);
            y1=extract<double>(y[1]);
        } else
            throw OxleyException("Argument l1 must be a float or 2-tuple");
    } else if (extract<double>(l1).check()) {
        y1=extract<double>(l1);
    } else
        throw OxleyException("Argument l1 must be a float or 2-tuple");

    if (extract<tuple>(l2).check()) {
        tuple z=extract<tuple>(l2);
        if (len(z)==2) {
            z0=extract<double>(z[0]);
            z1=extract<double>(z[1]);
        } else
            throw OxleyException("Argument l2 must be a float or 2-tuple");
    } else if (extract<double>(l2).check()) {
        z1=extract<double>(l2);
    } else
        throw OxleyException("Argument l2 must be a float or 2-tuple");

    // process tags and points
    boost::python::list pypoints=extract<boost::python::list>(objpoints);
    boost::python::list pytags=extract<boost::python::list>(objtags);
    int numpts=extract<int>(pypoints.attr("__len__")());
    int numtags=extract<int>(pytags.attr("__len__")());
    std::vector<double> points;
    std::vector<int> tags;
    tags.resize(numtags, -1);
    for (int i=0;i<numpts;++i) {
        tuple temp = extract<tuple>(pypoints[i]);
        int l=extract<int>(temp.attr("__len__")());
        if (l != 3)
            throw OxleyException("Number of coordinates for each dirac point must match dimensions.");
        for (int k=0;k<l;++k) {
            points.push_back(extract<double>(temp[k]));
        }
    }
    std::map<std::string, int> tagstonames;
    int curmax=40;
    // but which order to assign tags to names?????
    for (int i=0;i<numtags;++i) {
        extract<int> ex_int(pytags[i]);
        extract<std::string> ex_str(pytags[i]);
        if (ex_int.check()) {
            tags[i]=ex_int();
            if (tags[i]>= curmax) {
                curmax=tags[i]+1;
            }
        } else if (ex_str.check()) {
            std::string s=ex_str();
            std::map<std::string, int>::iterator it=tagstonames.find(s);
            if (it!=tagstonames.end()) {
                // we have the tag already so look it up
                tags[i]=it->second;
            } else {
                tagstonames[s]=curmax;
                tags[i]=curmax;
                curmax++;
            }
        } else {
            throw OxleyException("Error - Unable to extract tag value.");
        }
    }
    if (numtags != numpts)
        throw OxleyException("Number of tags does not match number of points.");

    // Handle optional MPI communicator
    escript::JMPI jmpi = escript::makeInfoFromPyComm(py_comm);

    std::vector<int> levels = extractRefineLevel(refine_level, n0*n1*n2);

    return escript::Domain_ptr(new Brick(jmpi, order, n0,n1,n2, x0,y0,z0, x1,y1,z1,
                                points, tags, tagstonames, levels));
}

// //tmp
// oxley::RefinementQueue_Ptr _refinementQueue()
// {
//     return oxley::RefinementQueue_Ptr(new RefinementQueue());
//     // return oxley::RefinementQueue2D_Ptr(new RefinementQueue2D());
// }


oxley::RefinementQueue2D_Ptr _refinementQueue2D()
{
    return oxley::RefinementQueue2D_Ptr(new RefinementQueue2D());
}

oxley::RefinementQueue3D_Ptr _refinementQueue3D()
{
    return oxley::RefinementQueue3D_Ptr(new RefinementQueue3D());
}

BOOST_PYTHON_MODULE(oxleycpp)
{

    def("Rectangle", oxley::_rectangle, (
    arg("n0"),arg("n1"),
    arg("l0")=1.0,arg("l1")=1.0,
    arg("diracPoints")=list(), arg("diracTags")=list(),
    arg("comm")=object(), arg("refine_level")=0),
    "Creates a rectangular p4est mesh of n0 x n1 blocks over the rectangle [0,l0] x [0,l1],\n"
    "each block subdivided refine_level times.\n\n"
    ":param n0: number of blocks in direction 0\n:type n0: ``int``\n"
    ":param n1: number of blocks in direction 1\n:type n1: ``int``\n"
    ":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
    ":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
    ":param refine_level: refinement level applied to every block, or a flat list\n"
    "    of one level per block in row-major block order (differing levels create\n"
    "    hanging nodes at the block seams)\n:type refine_level: ``int`` or ``list`` of ``int``\n"
    ":param comm: MPI communicator (optional, from mpi4py)\n:type comm: ``mpi4py.MPI.Comm``");

    def("Brick", oxley::_brick, (
    arg("n0"),arg("n1"),arg("n2"),
    arg("l0")=1.0,arg("l1")=1.0,arg("l2")=1.0,
    arg("diracPoints")=list(), arg("diracTags")=list(),
    arg("comm")=object(), arg("refine_level")=0),
    "Creates a brick p8est mesh of n0 x n1 x n2 blocks over [0,l0] x [0,l1] x [0,l2],\n"
    "each block subdivided refine_level times.\n\n"
    ":param n0: number of blocks in direction 0\n:type n0: ``int``\n"
    ":param n1: number of blocks in direction 1\n:type n1: ``int``\n"
    ":param n2: number of blocks in direction 2\n:type n2: ``int``\n"
    ":param l0: length of side 0 or coordinate range of side 0\n:type l0: ``float`` or ``tuple``\n"
    ":param l1: length of side 1 or coordinate range of side 1\n:type l1: ``float`` or ``tuple``\n"
    ":param l2: length of side 2 or coordinate range of side 1\n:type l2: ``float`` or ``tuple``\n"
    ":param refine_level: refinement level applied to every block, or a flat list\n"
    "    of one level per block in row-major block order (differing levels create\n"
    "    hanging nodes at the block seams)\n:type refine_level: ``int`` or ``list`` of ``int``\n"
    ":param comm: MPI communicator (optional, from mpi4py)\n:type comm: ``mpi4py.MPI.Comm``");

    // def("RefinementQueue", oxley::_refinementQueue, 
    //     "Creates a refinement zone parent class.\n\n"
    //     );

    def("RefinementQueue2D", oxley::_refinementQueue2D, 
        "Creates a refinement zone of dimension 2.\n\n"
        );

    def("RefinementQueue3D", oxley::_refinementQueue3D, 
        "Creates a refinement zone of dimension 3.\n\n"
        );

// #ifdef ESYS_HAVE_BOOST_NUMPY
//     // def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
//     def("addSurface", oxley::_addCurve, (arg("domain"),arg("x"),arg("y")));
//     def("addSurface", oxley::_addSurface, (arg("domain"),arg("x"),arg("y"),arg("z")));
// #endif

    class_<oxley::OxleyDomain, bases<escript::AbstractContinuousDomain>, boost::noncopyable >
        ("OxleyDomain", "", no_init)
        .def("toFinley", &oxley::toFinley,
            (arg("self"), arg("order")=-1, arg("reducedOrder")=-1,
             arg("optimize")=false, arg("simplices")=true),
            "returns a finley domain describing the same mesh\n\n"
            "Each octant is split into simplices, Tri3 in 2D and Tet4 in 3D. "
            "The forest supplies the geometry; finley owns the degrees of "
            "freedom and the parallel overlap. Only conforming forests can be "
            "converted so far.\n\n"
            ":param order: integration order, -1 for the default\n"
            ":type order: ``int``\n"
            ":param reducedOrder: reduced integration order, -1 for the default\n"
            ":type reducedOrder: ``int``\n"
            ":param optimize: whether to let finley repartition with ParMETIS\n"
            ":type optimize: ``bool``\n"
            ":param simplices: when False emit one Rec4/Hex8 per octant "
            "instead of splitting; for debugging only\n"
            ":type simplices: ``bool``\n"
            ":rtype: `Domain`")
        .def("lnodesDegree2Report", &oxley::_lnodesDegree2Report,
            (arg("self")),
            "diagnostic: does a degree-2 lnodes number the hanging positions?\n\n"
            "returns [octants, cornerSlots, distinctIds, idsAtSeveralPositions, "
            "hangingOctants, localNodes, distinctPositions, positionsWithSeveralIds]. "
            "idsAtSeveralPositions and positionsWithSeveralIds must both be 0.\n\n"
            ":rtype: ``list``")
        .def("lnodesDegree2Slots", &oxley::_lnodesDegree2Slots,
            (arg("self")),
            "diagnostic: every degree-2 lnodes slot as a tuple\n"
            "(element, slot, globalId, x, y, z, faceCode, isCorner), where x,y,z is "
            "the position OF THE SLOT and globalId is what the slot holds\n\n"
            ":rtype: ``list``")
        .def("isConforming", &oxley::OxleyDomain::isConforming,
            "returns True when no element in the forest has a hanging node\n\n"
            ":rtype: ``bool``")
        .def("addToRHS",&oxley::OxleyDomain::addToRHSFromPython,
            args("rhs", "data"),
            "adds a PDE onto the stiffness matrix mat and a rhs, "
            "results depends on domain\n\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\n:type data: `list`\n")
        .def("addToSystem",&oxley::OxleyDomain::addToSystemFromPython,
            args("mat", "rhs", "data"),
            "adds a PDE to the system, results depend on domain\n\n"
            ":param mat:\n:type mat: `OperatorAdapter`\n"
            ":param rhs:\n:type rhs: `Data`\n"
            ":param data:\n:type data: `list`\n")
        .def("createAssembler", &oxley::OxleyDomain::createAssemblerFromPython,
            args("typename", "options"),
            "request from the domain an assembler of the specified type, if "
            "supported, using the supplied options (if provided)"
            ":param typename:\n:type typename: `string`\n"
            ":param options:\n:type options: `list`\n")
        .def("dump", &oxley::OxleyDomain::dump, args("filename"),
            "Dumps the mesh to a silo file with the name `filename`."
            ":param filename:\n:type typename: `string`\n")
        .def("getDataShape", &oxley::OxleyDomain::getDataShape, args("functionSpaceCode"),
            ":return: a pair (dps, ns) where dps is the number of data points per sample, and ns is the number of samples\n"
            ":rtype: ``tuple``")        
        .def("getDescription", &oxley::OxleyDomain::getDescription,
                "Prints out a description of the mesh.")
        .def("getDim", &oxley::OxleyDomain::getDim, ":rtype: ``int``")
        .def("getNormal",&oxley::OxleyDomain::getNormal,
            ":return: boundary normals at the quadrature point on the face elements\n"
            ":rtype: `Data`")
        .def("getNumVertices", &oxley::OxleyDomain::getNumVertices,
            "Returns the number of corners in the mesh.\n"
            ":rtype: ``int``")
        .def("getSystemMatrixTypeId",&oxley::OxleyDomain::getSystemMatrixTypeId,
            args("options"),
            ":return: the identifier of the matrix type to be used for the global stiffness matrix when particular solver options are used.\n"
            ":rtype: ``int``\n"
            ":param options:\n:type options: `SolverBuddy`\n")
        .def("setTagMap",&oxley::OxleyDomain::setTagMap,args("name","tag"),
            "Give a tag number a name.\n\n:param name: Name for the tag\n:type name: ``string``\n"
            ":param tag: numeric id\n:type tag: ``int``\n:note: Tag names must be unique within a domain")
        .def("getTransportTypeId",&oxley::OxleyDomain::getTransportTypeId,
            args("solver", "preconditioner", "package", "symmetry"),
            ":return: the identifier of the transport problem type to be used when a particular solver, preconditioner, package and symmetric matrix is used.\n"
            ":rtype: ``int``\n"
            ":param solver:\n:type solver: ``int``\n"
            ":param preconditioner:\n:type preconditioner: ``int``\n"
            ":param package:\n:type package: ``int``\n"
            ":param symmetry:\n:type symmetry: ``int``")
        .def("getX",&oxley::OxleyDomain::getX, ":return: locations in the FEM nodes\n\n"
            ":rtype: `Data`")
        #ifdef ESYS_HAVE_BOOST_NUMPY
        .def("getMeshInfo",&oxley::OxleyDomain::getMeshInfo,
            (arg("materializeHanging")=false),
            "Returns an lnodes-based view of the mesh (node coordinates, global ids,\n"
            "element-to-node connectivity and element tags) as a dict of numpy arrays.\n\n"
            ":param materializeHanging: give every hanging position a node of its own,\n"
            "    appended after the lnodes nodes and described by the returned\n"
            "    numRealNodes/constrainedNodes/constraintMasters/constraintWeights.\n"
            "    Without it an element whose corner hangs lists a master there, a node\n"
            "    that lies outside the element.\n"
            ":type materializeHanging: ``bool``\n"
            ":rtype: ``dict``")
        #endif
        #ifdef ESYS_HAVE_TRILINOS
        #endif
        .def("saveFsType",&oxley::OxleyDomain::saveFsType, arg("rhs"), "saves the fs type")
        .def("getOrigFsType",&oxley::OxleyDomain::getOrigFsType, "returns the fs type")
        .def("loadMesh", &oxley::OxleyDomain::loadMesh, (arg("filename")),
                "Loads a mesh (in p4est format)\n"
                ":param filename: The name of the file to load\n")
        .def("newOperator",&oxley::OxleyDomain::newSystemMatrix,
            args("row_blocksize", "row_functionspace", "column_blocksize", "column_functionspace", "type"),
            "creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros\n\n"
            ":param row_blocksize:\n:type row_blocksize: ``int``\n"
            ":param row_functionspace:\n:type row_functionspace: `FunctionSpace`\n"
            ":param column_blocksize:\n:type column_blocksize: ``int``\n"
            ":param column_functionspace:\n:type column_functionspace: `FunctionSpace`\n"
            ":param type:\n:type type: ``int``")
        .def("saveMesh", &oxley::OxleyDomain::saveMesh, (arg("filename")),
                "Saves the mesh to file using p4est format\n"
                ":param filename: The name of the output file\n")
        .def("showTagNames",&oxley::OxleyDomain::showTagNames,
                ":return: A space separated list of tag names\n:rtype: ``string``")
        .def("writeToVTK", &oxley::OxleyDomain::writeToVTK, (arg("filename"), arg("writeMesh")=false),
                "Writes the mesh to a VTK file.\n"
                ":param filename: The name of the output file\n"
                ":param writeMesh: Boolean: Only writes the mesh to file")

        ;

    // These two class exports are necessary to ensure that the extra methods
    // added by oxley make it to python.
    //
    // The refinement methods that used to live here are gone: refinement now
    // goes through RefinementQueue2D/3D, which applies to a domain and hands
    // back a NEW one. They mutated the domain in place, which left every Data
    // already defined over it silently stale.
    class_<oxley::Brick, bases<oxley::OxleyDomain> >("OxleyBrick", "", no_init)
        ;

    class_<oxley::Rectangle, bases<oxley::OxleyDomain> >("OxleyRectangle", "", no_init)
        .def("interpolate", &oxley::Rectangle::interpolateAcross, (arg("target"),arg("source")),
                "Interpolates source to target\n"
                ":param source:\n:type Data: The source Data object. \n"
                ":param target:\n:type Data: The target Data object. \n")
        ;

    class_<oxley::RefinementQueue>("RefinementQueue", "")

    ;

    class_<oxley::RefinementQueue2D, bases<oxley::RefinementQueue>>("RefinementQueue2D")
        .def("refineUniform", &oxley::RefinementQueue2D::refineUniform, (arg("level")=-1),
                "Queues a refinement of EVERY element, the old refineMesh(\"uniform\").\n"
                ":param level: levels of refinement, default the queue's own\n"
                ":type level: ``int``")
        .def("apply", &oxley::RefinementQueue2D::apply, (arg("domain")),
                "Applies the queued refinements to a domain and returns the RESULT as\n"
                "a new domain. The domain passed in is not modified, so the caller\n"
                "keeps a usable handle on the coarser mesh and on any Data over it.\n\n"
                ":param domain: the domain to refine\n"
                ":type domain: `Domain`\n"
                ":return: the refined domain\n:rtype: `Domain`")
        .def("setRefinementLevel", &oxley::RefinementQueue2D::setRefinementLevel, (arg("level")),
                "Sets the level of refinement\n"
                ":param level:\n:type int: the level of the refinement.\n")
        .def("print", &oxley::RefinementQueue2D::print,
                "Prints the current queue to console\n")
        .def("remove", &oxley::RefinementQueue2D::deleteFromQueue, (arg("n")),
                "Removes the n^th item from the queue\n"
                ":param n:\n:type int: the refinement to remove.\n")
        .def("refinePoint", &oxley::RefinementQueue2D::refinePoint, (arg("x0"),arg("y0"),arg("level")=-1),
                "Refines the mesh around the point (x0,y0) to the level of refinement"
                "set by setRefinementLevel \n"
                ":param x0:\n:type float: x coordinate of the point to be refined.\n"
                ":param y0:\n:type float: y coordinate of the point to be refined.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineRegion", &oxley::RefinementQueue2D::refineRegion, (arg("x0"),arg("y0"),arg("x1"),arg("y1"),arg("level")=-1),
                "Refines the mesh around a rectangular region bound by the points (x0,y0)"
                "and (x1,y1) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the upper left coordinate.\n"
                ":param y0:\n:type float: y coordinate of the upper left coordinate.\n"
                ":param x1:\n:type float: x coordinate of the lower right coordinate.\n"
                ":param y1:\n:type float: y coordinate of the lower right coordinate.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineCircle", &oxley::RefinementQueue2D::refineCircle, (arg("x0"),arg("y0"),arg("r"),arg("level")=-1),
                "Refines the mesh around a circular region with radius r and center"
                "and (x0,y0) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the center of the circle.\n"
                ":param y0:\n:type float: y coordinate of the center of the circle.\n"
                ":param r :\n:type float: the radius of the circle.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineBorder", static_cast<void (oxley::RefinementQueue2D::*)(std::string, float, int)>(&oxley::RefinementQueue2D::refineBorder), (arg("border"),arg("dx"),arg("level")=-1),
                "Refines the border of the mesh to depth dx to the level of refinement"
                "set by setRefinementLevel\n"
                ":param Border:\n:type string: The border to refine (top,bottom,right,left).\n"
                ":param dx:\n:type float: the depth of the refinement.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineMask", &oxley::RefinementQueue2D::refineMask, (args("mask")),
                "Refines the mesh in regions defined by a mask\n"
                ":param mask:\n:type Data: a mask.\n")
        ;

    class_<oxley::RefinementQueue3D, bases<oxley::RefinementQueue>>("RefinementQueue3D")
        .def("refineUniform", &oxley::RefinementQueue3D::refineUniform, (arg("level")=-1),
                "Queues a refinement of EVERY element, the old refineMesh(\"uniform\").\n"
                ":param level: levels of refinement, default the queue's own\n"
                ":type level: ``int``")
        .def("apply", &oxley::RefinementQueue3D::apply, (arg("domain")),
                "Applies the queued refinements to a domain and returns the RESULT as\n"
                "a new domain. The domain passed in is not modified, so the caller\n"
                "keeps a usable handle on the coarser mesh and on any Data over it.\n\n"
                ":param domain: the domain to refine\n"
                ":type domain: `Domain`\n"
                ":return: the refined domain\n:rtype: `Domain`")
        .def("setRefinementLevel", &oxley::RefinementQueue3D::setRefinementLevel, (args("level")),
                "Sets the level of refinement\n"
                ":param level:\n:type int: the level of the refinement.\n")
        .def("print", &oxley::RefinementQueue3D::print, (arg("level")),
                "Prints the current queue to console\n")
        .def("remove", &oxley::RefinementQueue3D::deleteFromQueue, (arg("n")),
                "Removes the n^th item from the queue\n"
                ":param n:\n:type int: the refinement to remove.\n")
        .def("refinePoint", &oxley::RefinementQueue3D::refinePoint, (arg("x0"),arg("y0"),arg("z0"),arg("level")=-1),
                "Refines the mesh around the point (x0,y0,z0) to the level of refinement"
                "set by setRefinementLevel \n"
                ":param x0:\n:type float: x coordinate of the point to be refined.\n"
                ":param y0:\n:type float: y coordinate of the point to be refined.\n"
                ":param z0:\n:type float: z coordinate of the point to be refined.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineRegion", &oxley::RefinementQueue3D::refineRegion, (arg("x0"),arg("y0"),arg("z0"),arg("x1"),arg("y1"),arg("z1"),arg("level")=-1),
                "Refines the mesh around a rectangular region bound by the points (x0,y0,z0)"
                "and (x1,y1,z1) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the upper left coordinate.\n"
                ":param y0:\n:type float: y coordinate of the upper left coordinate.\n"
                ":param z0:\n:type float: z coordinate of the upper left coordinate.\n"
                ":param x1:\n:type float: x coordinate of the lower right coordinate.\n"
                ":param y1:\n:type float: y coordinate of the lower right coordinate.\n"
                ":param z1:\n:type float: z coordinate of the lower right coordinate.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineSphere", &oxley::RefinementQueue3D::refineSphere, (arg("x0"),arg("y0"),arg("z0"),arg("r"),arg("level")=-1),
                "Refines the mesh around a spherical region with radius r and center"
                "and (x0,y0,z0) to the level of refinement set by setRefinementLevel\n"
                ":param x0:\n:type float: x coordinate of the center of the circle.\n"
                ":param y0:\n:type float: y coordinate of the center of the circle.\n"
                ":param z0:\n:type float: z coordinate of the center of the circle.\n"
                ":param r :\n:type float: the radius of the circle.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineBorder", static_cast<void (oxley::RefinementQueue3D::*)(std::string, float, int)>(&oxley::RefinementQueue3D::refineBorder), (arg("border"),arg("dx"),arg("level")=-1),
                "Refines the border of the mesh to depth dx to the level of refinement"
                "set by setRefinementLevel\n"
                ":param Border:\n:type string: The border to refine (top,bottom,right,left).\n"
                ":param dx:\n:type float: the depth of the refinement.\n"
                ":param level:\n:type float: the level of refinement.\n")
        .def("refineMask", &oxley::RefinementQueue3D::refineMask, (arg("mask"),arg("level")=-1),
                "Refines the mesh in regions defined by a mask\n"
                ":param mask:\n:type Data: a mask.\n")
        ;

    class_<oxley::AbstractAssembler, oxley::Assembler_ptr, boost::noncopyable >  ("AbstractAssembler", "", no_init);

    // register_ptr_to_python<boost::shared_ptr<RefinementQueue>>();
    // register_ptr_to_python<boost::shared_ptr<RefinementQueue2D>>();
    // register_ptr_to_python<boost::shared_ptr<RefinementQueue3D>>();

}

} //namespace oxley
