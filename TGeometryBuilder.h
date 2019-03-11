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
#include <cmath>
#include <thread>
#include "gmsh.h"

class TGeometryBuilder {
    
private:
    
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
    
    
};


#endif /* TGeometryBuilder_h */
