// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   Intrepid2_OrientationToolsDefMatrixData.hpp
    \brief  Definition file for matrix data in the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  template<typename DT>
  template<typename BasisHostType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createCoeffMatrixInternal(const BasisHostType* basis, const bool inverse) {
    const std::string name(basis->getName());
    CoeffMatrixDataViewType matData;

    const auto cellTopo = basis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    ordinal_type matDim = 0, matDim1 = 0, matDim2 = 0, numOrts = 0, numSubCells;
    for(ordinal_type i=0; i<numEdges; ++i) {
      matDim1 = std::max(matDim1, basis->getDofCount(1,i));
      numOrts = std::max(numOrts,2);
    }
    for(ordinal_type i=0; i<numFaces; ++i) {
      matDim2 = std::max(matDim2, basis->getDofCount(2,i));
      numOrts = std::max(numOrts,2*ordinal_type(cellTopo.getSideCount(2,i)));
    }
    matDim = std::max(matDim1,matDim2);
    numSubCells = (matDim1>0)*numEdges + (matDim2>0)*numFaces;


    matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::"+name,
                                      numSubCells,
                                      numOrts,
                                      matDim,
                                      matDim);

    if(basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD) {
      init_HGRAD(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL) {
      init_HCURL(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV) {
      init_HDIV(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HVOL) {
      init_HVOL(matData, basis, inverse);
    }
    return matData;
  }

  //
  // HGRAD elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HGRAD(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
                                 BasisHostType const *cellBasis,
                                 const bool inverse) {

    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    basisPtr;
    BasisHostType const *subcellBasis;
    
    { //edges
      subcellBasis = cellBasis; // if (dim==1)
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 2) continue;
        if(cellTopo.getDimension()!=1) {
          basisPtr = cellBasis->getSubCellRefBasis(1,edgeId);
          subcellBasis = basisPtr.get();
        }

        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
              edgeId, edgeOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HGRAD
            (mat,
              *subcellBasis, *cellBasis,
              edgeId, edgeOrt, inverse);
        }
      }
    }
    { //faces
      subcellBasis = cellBasis; // if(dim==2)
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadrilaterals (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        if(cellTopo.getDimension()!=2) {
          basisPtr = cellBasis->getSubCellRefBasis(2,faceId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
              numEdges+faceId, faceOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HGRAD
            (mat,
             *subcellBasis, *cellBasis,
             faceId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HCURL elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HCURL(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
      BasisHostType const *cellBasis, const bool inverse) {
    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    basisPtr;
    BasisHostType const* subcellBasis;

    { // edges
      subcellBasis = cellBasis; // if (dim==1)
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 1) continue;
        if(cellTopo.getDimension()!=1) {
          basisPtr = cellBasis->getSubCellRefBasis(1,edgeId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
                                     edgeId, edgeOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
              *subcellBasis, *cellBasis,
              edgeId, edgeOrt, inverse);
        }
      }
    }
    { //faces
      subcellBasis = cellBasis; // if (dim==2)
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadratures (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        if(cellTopo.getDimension()!=2) {
          basisPtr = cellBasis->getSubCellRefBasis(2,faceId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
                                     numEdges+faceId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL
            (mat,
             *subcellBasis, *cellBasis,
             faceId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HDIV elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HDIV(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
      BasisHostType const *cellBasis, const bool inverse) {
    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numSides = cellTopo.getSideCount();
    const ordinal_type sideDim = cellTopo.getDimension()-1;
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    subcellBasisPtr;

    {
      for (ordinal_type sideId=0;sideId<numSides;++sideId) {
        if(cellBasis->getDofCount(sideDim, sideId) < 1) continue;
        const ordinal_type numOrt = (sideDim == 1) ? 2 : 2*cellTopo.getSideCount(sideDim,sideId);
        subcellBasisPtr = cellBasis->getSubCellRefBasis(sideDim,sideId);
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData, 
                                     sideId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
              *subcellBasisPtr, *cellBasis,
              sideId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HVOL elements (used for 2D and 1D side cells)
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HVOL(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
                                 BasisHostType const *cellBasis, const bool inverse) {

    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = (cellTopo.getDimension()==1);
    const ordinal_type numFaces = (cellTopo.getDimension()==2);

    {
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 1) continue;
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
              edgeId, edgeOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HVOL
            (mat, *cellBasis, edgeOrt, inverse);
        }
      }
    }
    {
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadratures (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
              numEdges+faceId, faceOrt,
              Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HVOL
              (mat, *cellBasis, faceOrt, inverse);
        }
      }
    }
  }

  template<typename DT>
  template<typename BasisType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createCoeffMatrix(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      ortCoeffData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::CoeffMatrixDataViewType>();
    });

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto found = ortCoeffData.find(key);
    
    CoeffMatrixDataViewType matData;
    if (found == ortCoeffData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        matData = createCoeffMatrixInternal(basis_host.getRawPtr());
        ortCoeffData.insert(std::make_pair(key, matData));
      }
    } else {
      matData = found->second;
    }
    
    return matData;
  }

  template<typename DT>
  template<typename BasisType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createInvCoeffMatrix(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      ortInvCoeffData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::CoeffMatrixDataViewType>();
    });

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto found = ortInvCoeffData.find(key);
    
    CoeffMatrixDataViewType matData;
    if (found == ortInvCoeffData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        matData = createCoeffMatrixInternal(basis_host.getRawPtr(),true);
        ortInvCoeffData.insert(std::make_pair(key, matData));
      }
    } else {
      matData = found->second;
    }
    
    return matData;
  }
  
  template<typename DT>
  void OrientationTools<DT>::clearCoeffMatrix() {
    ortCoeffData.clear();
    ortInvCoeffData.clear();
  }
}

#endif
