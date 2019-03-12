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
#include <fstream>
#include <sstream> 
#include <map>
#include <cmath>
#include <thread>
#include "gmsh.h"


class TGeometryBuilder {
    
public:
    
    /// Vector of points of the internal wire
    std::vector<std::vector<double>> m_internal_wire_pts;
    
    /// Vector of points of the external wire
    std::vector<std::vector<double>> m_external_wire_pts;
    
    /// Vector of internal points
    std::vector<std::vector<double>> m_internal_pts;
    
    /// Vector of DFN points
    std::vector<std::vector<double>> m_fracture_pts;
    
    /// Vector of internal wire point tags
    std::vector<int> m_internal_wire_point_tags;
    
    /// Vector of internal wire curve tags
    std::vector<int> m_internal_wire_curve_tags;
    
    /// Internal wire loop tag
    int m_internal_wire_curve_loop;
    
    /// Vector of external wire point tags
    std::vector<int> m_external_wire_point_tags;
    
    /// Vector of external wire curve tags
    std::vector<int> m_external_wire_curve_tags;
    
    /// External wire loop tag
    int m_external_wire_curve_loop;
    
    /// Vector of internal point tags
    std::vector<int> m_internal_point_tags;
    
    /// Vector of base fracture points
    std::vector<int> m_base_fracture_point_tags;
    
    /// Vector of fracture points
    std::vector<int> m_fracture_point_tags;
    
    /// Vector of base fracture curves
    std::vector<int> m_base_fracture_curve_tags;
    
    /// Map of base fracture curve tag to subfracture curve tags
    std::map<int, std::vector<int> > m_fracture_curve_tags;
    
    /// Entity tag for the wellbore region and boundaries
    std::vector<int>  m_wellbore_region_tags;
    
    /// Physical tag for the wellbore region and boundaries
    std::vector<int>  m_wellbore_region_physical_tags;
    
    /// Load the user-defined points (x,y,z)
    std::vector<std::vector<double>> LoadPoints(std::string & file_name, int n_data);

    /// Construct the a wire based on provided points. Uses gmsh as engine.
    void BuildWire(std::vector<std::vector<double>> & points, std::vector<int> & ignored_points, std::vector<int> & point_tags, std::vector<int> & curve_tags, int & curve_loop);
    
    /// Construct the internal wire. Uses gmsh as engine.
    void BuildInternalWire();
    
    /// Construct the exertnal wire. Uses gmsh as engine.
    void BuildExternalWire();
    
    /// Construct internal points. Uses gmsh as engine.
    void BuildInternalPoints();
    
    /// Construct the DFN (including intersections fracture-frature or fracture-boundary). Uses gmsh as engine.
    void BuildDFN();
    
    /// Draw a circle points with radius r and center x_center
    void DrawCirclePoints(double r, std::vector<double> x_center, std::vector<std::vector<double>> & points);

    /// Draw a rectangle points with with x_mix coordinate and x_max coordinate
    void DrawRectanglePoints(std::vector<double> x_mix, std::vector<double> x_max, std::vector<std::vector<double>> & points);
    
    /// Draw the DFN add points and 1D fractures
    void DrawDFN();
    
    /// Computes the intersection between two fractures and return false when the intersection is empty.
    bool IntersectFractures(int object_tag, int tool_tag, std::vector<gmsh::vectorpair> & map_dim_tags);
    
    bool ComputeFracturesIntersections(std::map<int,std::vector<int>> & objects, std::map<int,std::vector<int>>& tools, std::map<int,std::vector<int>>  & fractures);
    
    std::vector<int> ComputeAssociatedMicroFractures(std::pair<int,std::vector<int>> & fracture_data,std::map<int,std::vector<int>>  & fractures);
    
public:
    
    /// Default constructor
    TGeometryBuilder();
    
    /// Copy constructor
    TGeometryBuilder(const TGeometryBuilder & other);
    
    /// Assignement constructor
    const TGeometryBuilder & operator=(const TGeometryBuilder & other);
    
    /// Destructor constructor
    ~TGeometryBuilder();
    
    /// Load the user-defined internal wire points (x,y,z). The wire is closed.
    void LoadInternalWire(std::string & file_name, int n_data);
    
    /// Load the user-defined external wire points (x,y,z). The wire is closed.
    void LoadExternalWire(std::string & file_name, int n_data);
    
    /// Load the user-defined points (x,y,z). They are used for fixing rigid modes.
    void LoadInternalPoints(std::string & file_name, int n_data);
    
    /// Load a DFN points (x,y,z). It requires an even number of lines.
    void LoadDiscreteFractureNetwork(std::string & file_name, int n_data);
    
    /// Draws a internal circle with radius r and center x_center
    void DrawInternalCricle(double r, std::vector<double> x_center);
    
    /// Draws a external circle with radius r and center x_center
    void DrawExternalCricle(double r, std::vector<double> x_center);
    
    /// Draws a external rectangle with x_mix coordinate and x_max coordinate
    void DrawExternalRectangle(std::vector<double> x_mix, std::vector<double> x_max);
    
    /// Draw a wellbore region domain based on the internal and external wires
    void DrawWellboreRegion();
    
    /// Compute reservoir physical tags (domain and boundaries)
    void ComputeReservoirPhysicalTags();
    
};


#endif /* TGeometryBuilder_h */
