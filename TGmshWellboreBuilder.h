//
//  TGmshWellboreBuilder.hpp
//  GmshMeshBuilder
//
//  Created by Omar Durán on 2/13/19.
//

#ifndef TGmshWellboreBuilder_h
#define TGmshWellboreBuilder_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <thread>
#include "gmsh.h"

class TGmshWellboreBuilder{
  
private:
    
    /// Stands for the wellbore trajectory points in (x,y,z) format
    std::vector<std::vector<double>> m_wellbore_trajectory;
    
    /// Stands for wellbore radius
    double m_wellbore_radius;
    
    gmsh::vectorpair m_base_dim_tags;
    
public:
    
    /// Default constructor
    TGmshWellboreBuilder();
    
    /// Constructor based on wellbore trajectory points
    TGmshWellboreBuilder(double & wellbore_radius, std::vector<std::vector<double>> & wellbore_trajectory);
    
    /// Copy constructor
    TGmshWellboreBuilder(const TGmshWellboreBuilder & other);
    
    /// Assignement constructor
    const TGmshWellboreBuilder & operator=(const TGmshWellboreBuilder & other);
    
    /// Destructor constructor
    ~TGmshWellboreBuilder();
    
    /// Set the wellbore radius
    void SetWellRadius(double wellbore_radius);
    
    /// Get the wellbore radius
    double WellRadius();
    
    /// Create the wellbore volume by sections and retun th entities dim_tag array
    gmsh::vectorpair DrawWellboreBySections();
    
    /// Create the wellbore volume by shell and retun th entities dim_tag array
    gmsh::vectorpair DrawWellboreByShell();
    
public: //  convert later to private
    
    /// Create Ring with normal pointing in z direction
    void CreateRingBase();
    
    /// Create a copy of the ring base
    gmsh::vectorpair CopyRingBase();
    
    /// Apply a rotation over the referred ring
    void RotateRing(gmsh::vectorpair ring_dim_tags, std::vector<double> & axis);
    
    /// Apply a translation over the referred ring
    void TranslateRing(gmsh::vectorpair ring_dim_tags, std::vector<double> & point);
    
    /// Remove Ring base
    void RemoveRingBase();
    
public:
    
    static double dot(std::vector<double> & u, std::vector<double> & v);
    
    static std::vector<double> cross(std::vector<double> & u, std::vector<double> & v);
    
    static double norm(std::vector<double> & v);
    
    static void normalize(std::vector<double> & v);
    
    static double angle(std::vector<double> & u, std::vector<double> & v);
    
};

#endif /* TGmshWellboreBuilder_h */
