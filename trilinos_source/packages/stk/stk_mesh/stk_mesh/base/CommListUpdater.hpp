// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_COMMLISTUPDATER_HPP
#define STK_COMMLISTUPDATER_HPP

#include <stk_mesh/base/EntityCommListInfo.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

namespace stk {
namespace mesh {

class CommListUpdater  : public CommMapChangeListener {
public:
    CommListUpdater(EntityCommListInfoVector& comm_list,
                    std::vector<EntityComm*>& entity_comms)
    : m_comm_list(comm_list), m_entity_comms(entity_comms)
    {}
    virtual ~CommListUpdater(){}

    virtual void removedKey(const EntityKey& key) {
        EntityCommListInfoVector::iterator iter =
                std::lower_bound(m_comm_list.begin(), m_comm_list.end(), key);
        if (iter != m_comm_list.end() && iter->key == key) {
            iter->entity_comm = nullptr;
            m_entity_comms[iter->entity.local_offset()] = nullptr;
        }
    }

private:
  EntityCommListInfoVector& m_comm_list;
  std::vector<EntityComm*>& m_entity_comms;
};

}
}

#endif
