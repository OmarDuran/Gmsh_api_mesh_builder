//
//  TGeometryBuilder.h
//  GmshMeshBuilder
//
//  Created by Omar Durán on 3/7/19.
//

#ifndef TGeometryBuilder_h
#define TGeometryBuilder_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <thread>
#include "gmsh.h"

class TGeometryBuilder {
    
    /// Default constructor
    TGeometryBuilder();
    
    /// Copy constructor
    TGeometryBuilder(const TGeometryBuilder & other);
    
    /// Assignement constructor
    const TGeometryBuilder & operator=(const TGeometryBuilder & other);
    
    /// Destructor constructor
    ~TGeometryBuilder();
    
};


#endif /* TGeometryBuilder_h */
