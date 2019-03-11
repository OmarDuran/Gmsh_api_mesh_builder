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
    
    int n_point = m_fracture_pts.size();
    for (int ip = 0; ip < n_point; ip++) {
        
        std::vector<double> point = m_internal_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_base_fracture_point_tags.push_back(point_tag);
    }
    
}

void TGeometryBuilder::DrawCirclePoints(double r, std::vector<double> x_center, std::vector<std::vector<double>> & points){
    
    std::vector<double> point(x_center);
    // p_center
    points.push_back(point);
    
    // p1
    point[0] = r;
    point[1] = 0;
    points.push_back(point);
    
    // p2
    point[0] = 0;
    point[1] = r;
    points.push_back(point);
    
    // p3
    point[0] = -r;
    point[1] = 0;
    points.push_back(point);
    
    // p4
    point[0] = 0;
    point[1] = -r;
    points.push_back(point);
    
}

void TGeometryBuilder::DrawRectanglePoints(std::vector<double> x_mix, std::vector<double> x_max, std::vector<std::vector<double>> & points){
    
    double z = 0.0;
    if (x_mix[2]==x_max[2]) {
        z = x_mix[2];
    }
    
    std::vector<double> point(x_mix);
    // p1
    point[0] = x_mix[0];
    point[1] = x_mix[1];
    point[2] = z;
    points.push_back(point);
    
    // p2
    point[0] = x_max[0];
    point[1] = x_mix[1];
    point[2] = z;
    points.push_back(point);
    
    // p2
    point[0] = x_max[0];
    point[1] = x_max[1];
    point[2] = z;
    points.push_back(point);
    
    // p4
    point[0] = x_mix[0];
    point[1] = x_max[1];
    point[2] = z;
    points.push_back(point);
    
}

void TGeometryBuilder::DrawInternalCricle(double r, std::vector<double> x_center){
    
    DrawCirclePoints(r,x_center,m_internal_wire_pts);
    int n_point = m_internal_wire_pts.size();
    for (int ip = 0; ip < n_point; ip++) {
    
        std::vector<double> point = m_internal_wire_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_internal_point_tags.push_back(point_tag);
    }
    
    int pc = m_internal_point_tags[0];
    int p1 = m_internal_point_tags[1];
    int p2 = m_internal_point_tags[2];
    int p3 = m_internal_point_tags[3];
    int p4 = m_internal_point_tags[4];
    
    int c1 = gmsh::model::occ::addEllipseArc(p1, pc, p2);
    int c2 = gmsh::model::occ::addEllipseArc(p2, pc, p3);
    int c3 = gmsh::model::occ::addEllipseArc(p3, pc, p4);
    int c4 = gmsh::model::occ::addEllipseArc(p4, pc, p1);
    
    m_internal_wire_curve_tags.push_back(c1);
    m_internal_wire_curve_tags.push_back(c2);
    m_internal_wire_curve_tags.push_back(c3);
    m_internal_wire_curve_tags.push_back(c4);
    m_internal_wire_curve_loop = gmsh::model::occ::addWire(m_internal_wire_curve_tags);
    
}

void TGeometryBuilder::DrawExternalCricle(double r, std::vector<double> x_center){
    
    DrawCirclePoints(r,x_center,m_external_wire_pts);
    int n_point = m_external_wire_pts.size();
    for (int ip = 0; ip < n_point; ip++) {
        
        std::vector<double> point = m_external_wire_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_external_wire_point_tags.push_back(point_tag);
    }
    
    
    int pc = m_external_wire_point_tags[0];
    int p1 = m_external_wire_point_tags[1];
    int p2 = m_external_wire_point_tags[2];
    int p3 = m_external_wire_point_tags[3];
    int p4 = m_external_wire_point_tags[4];
    
    int c1 = gmsh::model::occ::addEllipseArc(p1, pc, p2);
    int c2 = gmsh::model::occ::addEllipseArc(p2, pc, p3);
    int c3 = gmsh::model::occ::addEllipseArc(p3, pc, p4);
    int c4 = gmsh::model::occ::addEllipseArc(p4, pc, p1);
    
    m_external_wire_curve_tags.push_back(c1);
    m_external_wire_curve_tags.push_back(c2);
    m_external_wire_curve_tags.push_back(c3);
    m_external_wire_curve_tags.push_back(c4);
    m_external_wire_curve_loop = gmsh::model::occ::addWire(m_external_wire_curve_tags);
    
}

void TGeometryBuilder::DrawExternalRectangle(std::vector<double> x_mix, std::vector<double> x_max){
    
    DrawRectanglePoints(x_mix, x_max, m_external_wire_pts);
    int n_point = m_external_wire_pts.size();
    for (int ip = 0; ip < n_point; ip++) {
        
        std::vector<double> point = m_external_wire_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_external_wire_point_tags.push_back(point_tag);
    }
    
    int p1 = m_external_wire_point_tags[1];
    int p2 = m_external_wire_point_tags[2];
    int p3 = m_external_wire_point_tags[3];
    int p4 = m_external_wire_point_tags[4];
    
    int c1 = gmsh::model::occ::addLine(p1, p2);
    int c2 = gmsh::model::occ::addLine(p2, p3);
    int c3 = gmsh::model::occ::addLine(p3, p4);
    int c4 = gmsh::model::occ::addLine(p4, p1);
    
    m_external_wire_curve_tags.push_back(c1);
    m_external_wire_curve_tags.push_back(c2);
    m_external_wire_curve_tags.push_back(c3);
    m_external_wire_curve_tags.push_back(c4);
    m_external_wire_curve_loop = gmsh::model::occ::addWire(m_external_wire_curve_tags);
    
}
