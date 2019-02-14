//
//  TGmshWellboreBuilder.cpp
//  GmshMeshBuilder
//
//  Created by Omar Durán on 2/13/19.
//

#include "TGmshWellboreBuilder.h"

/// Default constructor
TGmshWellboreBuilder::TGmshWellboreBuilder(){
    
}

TGmshWellboreBuilder::TGmshWellboreBuilder(double & wellbore_radius, std::vector<std::vector<double>> & wellbore_trajectory){
    m_wellbore_radius       = wellbore_radius;
    m_wellbore_trajectory   = wellbore_trajectory;
}

TGmshWellboreBuilder::TGmshWellboreBuilder(const TGmshWellboreBuilder & other){
    m_wellbore_radius       = other.m_wellbore_radius;
    m_wellbore_trajectory   = other.m_wellbore_trajectory;
}

const TGmshWellboreBuilder & TGmshWellboreBuilder::operator=(const TGmshWellboreBuilder & other){
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_wellbore_radius       = other.m_wellbore_radius;
    m_wellbore_trajectory   = other.m_wellbore_trajectory;
    return *this;
}

TGmshWellboreBuilder::~TGmshWellboreBuilder(){
    
}

void TGmshWellboreBuilder::SetWellRadius(double wellbore_radius){
    m_wellbore_radius = wellbore_radius;
}

double TGmshWellboreBuilder::WellRadius(){
    return m_wellbore_radius;
}

gmsh::vectorpair TGmshWellboreBuilder::DrawWellbore(){
    
    gmsh::vectorpair wb_dim_tags;
    std::vector<gmsh::vectorpair> rings;
    CreateRingBase();
    
    int n_points = m_wellbore_trajectory.size();
    for (int i = 0; i < n_points - 1; i++) {
        std::vector<double> p_i = m_wellbore_trajectory[i];
        std::vector<double> p_e = m_wellbore_trajectory[i+1];
        std::vector<double> axial_dir = {p_e[0]-p_i[0],p_e[1]-p_i[1],p_e[2]-p_i[2]};
        normalize(axial_dir);
        gmsh::vectorpair ring_dim_tags = CopyRingBase();
        RotateRing(ring_dim_tags, axial_dir);
        TranslateRing(ring_dim_tags, p_i);
        rings.push_back(ring_dim_tags);
    }
    
    
    { /// Ending point case
        std::vector<double> p_i = m_wellbore_trajectory[n_points-2];
        std::vector<double> p_e = m_wellbore_trajectory[n_points-1];
        std::vector<double> axial_dir = {p_e[0]-p_i[0],p_e[1]-p_i[1],p_e[2]-p_i[2]};
        normalize(axial_dir);
        gmsh::vectorpair ring_dim_tags = CopyRingBase();
        RotateRing(ring_dim_tags, axial_dir);
        TranslateRing(ring_dim_tags, p_e);
        rings.push_back(ring_dim_tags);
    }
    
    
    for (int i = 0; i < n_points - 1; i++)
    {
//        int i = 0;
        std::vector<int> surface_tags;
        {
            gmsh::vectorpair ring_dim_tags_i = rings[i];
            gmsh::vectorpair ring_dim_tags_e = rings[i+1];
            
            std::vector<int> axial_lines_ring_dim_tags;
            for (int k = 1; k < 5; k++)
            { /// Adding lines
                int p_i = ring_dim_tags_i[k].second;
                int p_e = ring_dim_tags_e[k].second;
                int line_tag = gmsh::model::occ::addLine(p_i, p_e);
                axial_lines_ring_dim_tags.push_back(line_tag);
                
            }
            
            std::vector<int> plane_ring_dim_tags_i;
            std::vector<int> plane_ring_dim_tags_e;
            for (int k = 5; k < 9; k++)
            { /// Adding Plane surfaces
                plane_ring_dim_tags_i.push_back(ring_dim_tags_i[k].second);
                plane_ring_dim_tags_e.push_back(ring_dim_tags_e[k].second);
            }
            int wired_tag = gmsh::model::occ::addWire(plane_ring_dim_tags_i);
            int plane_surface_tag = gmsh::model::occ::addSurfaceFilling(wired_tag);
            surface_tags.push_back(plane_surface_tag);
            std::vector<int> axial={0,1,2,3,0};
            std::vector<int> planes={0,1,2,3};
            for (int k = 0; k < 4; k++)
            {
                std::vector<int> curve_ring_dim_tags;
                int a_i = axial[k];
                int a_e = axial[k+1];
                int p = planes[k];
                curve_ring_dim_tags.push_back(axial_lines_ring_dim_tags[a_i]);
                curve_ring_dim_tags.push_back(plane_ring_dim_tags_e[p]);
                curve_ring_dim_tags.push_back(axial_lines_ring_dim_tags[a_e]);
                curve_ring_dim_tags.push_back(plane_ring_dim_tags_i[p]);
                int curve_wired_tag = gmsh::model::occ::addWire(curve_ring_dim_tags);
                int curve_surface_tag = gmsh::model::occ::addSurfaceFilling(curve_wired_tag);
                surface_tags.push_back(curve_surface_tag);
            }
 
        }
        
    }
    
    
    /// Deleting ring base
//    RemoveRingBase();
    return wb_dim_tags;
}

void TGmshWellboreBuilder::CreateRingBase(){
    
    m_base_dim_tags.resize(9);
    for (int i = 0; i < 5; i++) {
        m_base_dim_tags[i].first = 0;
    }
    int pc_id = gmsh::model::occ::addPoint(0,0,0);
    m_base_dim_tags[0].second = pc_id;
    int p1_id = gmsh::model::occ::addPoint(m_wellbore_radius, 0, 0);
    m_base_dim_tags[1].second = p1_id;
    int p2_id = gmsh::model::occ::addPoint(0, m_wellbore_radius, 0);
    m_base_dim_tags[2].second = p2_id;
    int p3_id = gmsh::model::occ::addPoint(-m_wellbore_radius, 0, 0);
    m_base_dim_tags[3].second = p3_id;
    int p4_id = gmsh::model::occ::addPoint(0, -m_wellbore_radius, 0);
    m_base_dim_tags[4].second = p4_id;
    
    for (int i = 5; i < 9; i++) {
        m_base_dim_tags[i].first = 1;
    }
    m_base_dim_tags[5].second = gmsh::model::occ::addCircleArc(p1_id, pc_id, p2_id);
    m_base_dim_tags[6].second = gmsh::model::occ::addCircleArc(p2_id, pc_id, p3_id);
    m_base_dim_tags[7].second = gmsh::model::occ::addCircleArc(p3_id, pc_id, p4_id);
    m_base_dim_tags[8].second = gmsh::model::occ::addCircleArc(p4_id, pc_id, p1_id);
}


gmsh::vectorpair TGmshWellboreBuilder::CopyRingBase(){
    gmsh::vectorpair copy_dim_tags;
    gmsh::model::occ::copy(m_base_dim_tags, copy_dim_tags);
    return copy_dim_tags;
}

void TGmshWellboreBuilder::RotateRing(gmsh::vectorpair ring_dim_tags, std::vector<double> & axis){
    std::vector<double> z_unit = {0,0,1};
    double omega = angle(z_unit, axis);
    std::vector<double> dir = cross(z_unit, axis);
    gmsh::model::occ::rotate(ring_dim_tags, 0, 0, 0, dir[0], dir[1], dir[2], omega);
}

void TGmshWellboreBuilder::TranslateRing(gmsh::vectorpair ring_dim_tags, std::vector<double> & point){
    gmsh::model::occ::translate(ring_dim_tags, point[0], point[1], point[2]);
}

void TGmshWellboreBuilder::RemoveRingBase(){
    gmsh::model::occ::remove(m_base_dim_tags);
}

double TGmshWellboreBuilder::dot(std::vector<double> & u, std::vector<double> & v){
    double u_dot_v = 0;
    for (int i = 0; i < 3; i++) {
        u_dot_v += u[i]*v[i];
    }
    return u_dot_v;
}

std::vector<double> TGmshWellboreBuilder::cross(std::vector<double> & u, std::vector<double> & v){
    std::vector<double> cross = {-u[2]*v[1]+u[1]*v[2],u[2]*v[0]-u[0]*v[2],-u[1]*v[0]+u[0]*v[1]};
    return cross;
}

double TGmshWellboreBuilder::norm(std::vector<double> & v){
    if (v.size()!=3) {
        std::cout << "TGmshWellboreBuilder:: required a three entries." << std::endl;
    }
    double norm_v = dot(v, v);
    norm_v = sqrt(norm_v);
    return norm_v;
}

void TGmshWellboreBuilder::normalize(std::vector<double> & v){
    double norm_v = norm(v);
    for (int i = 0; i < 3; i++) {
        v[i] /= norm_v;
    }
}

double TGmshWellboreBuilder::angle(std::vector<double> & u, std::vector<double> & v){
    double norm_u = norm(u);
    double norm_v = norm(v);
    double u_dot_v = dot(u, v);
    double angle = acos(u_dot_v/(norm_u*norm_v));
    return angle;
}


