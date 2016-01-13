
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

/**************************************************************/

/*   Finley: Mesh: dump and read mesh to/from NetCDF file     */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/


void Finley_Mesh_dump(Finley_Mesh *in,char* fname) {
/*
dump:
  size, rank
  in->Name;                           
  in->order;                       
  in->reduced_order;                
  in->Nodes;              
          Id;    
          Tag;  
          globalDegreesOfFreedom  
          globalReducedDOFIndex
          globalReducedNodesIndex
          globalNodesIndex    
          Coordinates             
          degreesOfFreedomDistribution->distribution
  in->Elements;        
          ReferenceElement->ElementTypeId
          Owner;
          Id;          
          Tag;        
          Nodes;     
          Color;    
  in->FaceElements;    
          ReferenceElement->ElementTypeId
          Owner;
          Id;          
          Tag;        
          Nodes;     
          Color;    
  in->ContactElements; 
          ReferenceElement->ElementTypeId
          Owner;
          Id;          
          Tag;        
          Nodes;     
          Color;    
  in->Points;          
          ReferenceElement->ElementTypeId
          Owner;
          Id;          
          Tag;        
          Nodes;     
          Color;    
  in->TagMap;              
          name
          tag_key
          next;
             ...
*/
    Finley_setError(IO_ERROR,"Mesh_dump: not implemented yet.");
    return;
}
Finley_Mesh* Finley_Mesh_load(char* fname) {
  Finley_Mesh* out=NULL;
  index_t *distribution=NULL;
  Finley_setError(IO_ERROR,"Mesh_load: not implemented yet.");
  /* call Finley_Mesh_createMappings(out,distribution); */
  return out;
}
