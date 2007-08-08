/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: dump and read mesh to/from NetCDF file     */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:*/

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
          globalDegreesOfFreedom;  
          Coordinates;             
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

  Finley_setError(IO_ERROR,"Mesh_load: not implemented yet.");
  return NULL;
}
