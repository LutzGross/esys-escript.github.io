#ifndef NgpFieldTestUtils_hpp
#define NgpFieldTestUtils_hpp

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace ngp_field_test_utils {

template<typename T>
void check_field_data_on_device(stk::mesh::NgpMesh& ngpMesh,
                                stk::mesh::NgpField<T>& ngpField,
                                const stk::mesh::Selector& selector,
                                T expectedValue)
{
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      const int numComponents = ngpField.get_num_components_per_entity(entity);
      for(int i=0; i<numComponents; ++i) {
        STK_NGP_ThrowRequire(ngpField(entity, i) == expectedValue);
      }
    }
  );
}

template<typename T>
void check_field_data_on_host(const stk::mesh::HostMesh& stkMesh,
                              const stk::mesh::FieldBase& stkField,
                              const stk::mesh::Selector& selector,
                              T expectedValue)
{
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::FastMeshIndex& fastMeshIndex) {
      stk::mesh::Entity entity = stkMesh.get_entity(stkField.entity_rank(), fastMeshIndex);
      const int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      const T* fieldData = reinterpret_cast<const T*>(stk::mesh::field_data(stkField, entity));
      for(int i=0; i<numComponents; ++i) {
        STK_ThrowRequire(fieldData[i] == expectedValue);
      }
    }
  );
}

template<typename T>
void check_field_data_on_host(const stk::mesh::BulkData& stkMesh,
                              const stk::mesh::FieldBase& stkField,
                              const stk::mesh::Selector& selector,
                              T expectedValue)
{
  stk::mesh::for_each_entity_run(stkMesh, stkField.entity_rank(), selector,
    [&](const stk::mesh::BulkData& bulk, const stk::mesh::Entity entity) {
      const int numComponents = stk::mesh::field_scalars_per_entity(stkField, entity);
      const T* fieldData = reinterpret_cast<const T*>(stk::mesh::field_data(stkField, entity));
      for(int i=0; i<numComponents; ++i) {
        STK_ThrowRequire(fieldData[i] == expectedValue);
      }
    }
  );
}

} // ngp_field_test_utils

#endif

