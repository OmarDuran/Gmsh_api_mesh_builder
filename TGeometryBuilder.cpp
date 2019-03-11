//
//  TGeometryBuilder.cpp
//  GmshMeshBuilder
//
//  Created by Omar Durán on 3/7/19.
//

#include "TGeometryBuilder.h"


TGeometryBuilder::TGeometryBuilder(){
    
}

TGeometryBuilder::TGeometryBuilder(const TGeometryBuilder & other){
    
}

const TGeometryBuilder & TGeometryBuilder::operator=(const TGeometryBuilder & other){
    
}

TGeometryBuilder::~TGeometryBuilder(){
    
}

std::vector<std::vector<double>> TGeometryBuilder::LoadPoints(std::string & file_name, int n_data){
    
    std::ifstream in(file_name.c_str());
    std::vector<std::vector<double>> points(n_data);
    double x,y,z;
    int c = 0;
    while(in)
    {
        in >> x;
        in >> y;
        in >> z;
        std::vector<double> p = {x,y,z};
        points[c] = p;
        c++;
        if (c == n_data) {
            break;
        }
    }
    
    return points;
}

void TGeometryBuilder::LoadInternalWire(std::string & file_name, int n_data){
    m_internal_wire_pts = LoadPoints(file_name, n_data);
}

void TGeometryBuilder::LoadExternalWire(std::string & file_name, int n_data){
    m_external_wire_pts = LoadPoints(file_name, n_data);
}

void TGeometryBuilder::LoadInternalPoints(std::string & file_name, int n_data){
    m_internal_pts = LoadPoints(file_name, n_data);
}

void TGeometryBuilder::LoadDiscreteFractureNetwork(std::string & file_name, int n_data){
    m_fracture_pts = LoadPoints(file_name, n_data);
}

void TGeometryBuilder::DrawInternalCricle(double r, std::vector<double> x_center){
    
}

void TGeometryBuilder::DrawExternalCricle(double r, std::vector<double> x_center){
    
}

void TGeometryBuilder::DrawExternalRectangle(std::vector<double> x_mix, std::vector<double> x_max){
    
}

void TGeometryBuilder::BuildWire(std::vector<std::vector<double>> & points, std::vector<int> & ignored_points, std::vector<int> & point_tags, std::vector<int> & curve_tags, int & curve_loop){
    
    int n_point = points.size();
    for (int ip = 0; ip < n_point; ip++) {
        
        bool ignore_point_Q = false;
        for (auto ignored_point: ignored_points) {
            ignore_point_Q = ip == ignored_point;
            break;
        }
        
        if (ignore_point_Q) {
            continue;
        }
        
        std::vector<double> point = points[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        point_tags.push_back(point_tag);
    }
    
    int n_point_tag = point_tags.size();
    for (int i_tag = 0; i_tag < n_point_tag - 1; i_tag++) {
        int ini_tag = point_tags[i_tag];
        int end_tag = point_tags[i_tag+1];
        int curve_tag = gmsh::model::occ::addLine(ini_tag, end_tag);
        curve_tags.push_back(curve_tag);
    }

    curve_loop = gmsh::model::occ::addWire(curve_tags);
}

void TGeometryBuilder::BuildInternalWire(){
    std::vector<int> ignored_points;
    BuildWire(m_internal_wire_pts, ignored_points, m_internal_wire_point_tags, m_internal_wire_curve_tags, m_internal_wire_curve_loop);
}

void TGeometryBuilder::BuildExternalWire(){
    std::vector<int> ignored_points;
    BuildWire(m_external_wire_pts, ignored_points, m_external_wire_point_tags, m_external_wire_curve_tags, m_external_wire_curve_loop);
}

void TGeometryBuilder::BuildInternalPoints(){
    
    int n_point = m_internal_pts.size();
    for (int ip = 0; ip < n_point; ip++) {
        
        std::vector<double> point = m_internal_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_internal_point_tags.push_back(point_tag);
    }
    
}

void TGeometryBuilder::BuildDFN(){
    
}
