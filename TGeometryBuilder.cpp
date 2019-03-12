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

void TGeometryBuilder::DrawBaseFractures(){
    
    int n_point = m_fracture_pts.size();
    
    if (n_point <= 1) { /// No DFN is generated
        return;
    }
    
    bool even_number_Q = n_point % 2 == 0;
    if (!even_number_Q){
        std::cout << "TGeometryBuilder:: Odd fracture points, the last point is ignored." <<std::endl;
        m_fracture_pts.resize(n_point-1);
        n_point = m_fracture_pts.size();
    }
    
    for (int ip = 0; ip < n_point; ip++) {
        
        std::vector<double> point = m_fracture_pts[ip];
        double x = point[0];
        double y = point[1];
        double z = point[2];
        int point_tag = gmsh::model::occ::addPoint(x, y, z);
        m_base_fracture_point_tags.push_back(point_tag);
    }
    
    int n_point_tag = m_base_fracture_point_tags.size();
    for (int i_tag = 0; i_tag < n_point_tag - 1; i_tag += 2) {
        int ini_tag = m_base_fracture_point_tags[i_tag];
        int end_tag = m_base_fracture_point_tags[i_tag+1];
        int curve_tag = gmsh::model::occ::addLine(ini_tag, end_tag);
        m_base_fracture_curve_tags.push_back(curve_tag);
    }
    
}

void TGeometryBuilder::DrawDFN(){
    
    DrawBaseFractures();
    gmsh::model::occ::synchronize();
    int n_point = m_fracture_pts.size();
    if (n_point <= 1) { /// No DFN is generated
        return;
    }
    
    std::map<int,std::vector<int>> fractures;
    std::map<int,std::vector<int>> objects;
    std::map<int,std::vector<int>> tools;
    
    std::pair<int, std::vector<int> > chunk_tag_sub_tag;
    chunk_tag_sub_tag.second.resize(1);
    for (auto f: m_base_fracture_curve_tags) {
        chunk_tag_sub_tag.first = f;
        chunk_tag_sub_tag.second[0] = f;
        fractures.insert(chunk_tag_sub_tag);
        objects.insert(chunk_tag_sub_tag);
        tools.insert(chunk_tag_sub_tag);
        m_fracture_curve_tags.insert(chunk_tag_sub_tag);
    }
    
    bool there_are_intersections_Q = true;
    while (there_are_intersections_Q) {
        there_are_intersections_Q = ComputeFracturesIntersections(objects, tools, fractures);
    }
    
    chunk_tag_sub_tag.second.resize(0);
    /// Computing physical tag for fractures by group
    std::map<int,std::vector<int>> fractures_to_micro;
    for (auto f : m_fracture_curve_tags) {
        /// compute micro fractures
        chunk_tag_sub_tag = f;
        chunk_tag_sub_tag.second = ComputeAssociatedMicroFractures(chunk_tag_sub_tag, fractures);
        chunk_tag_sub_tag.first = f.first;
        fractures_to_micro.insert(chunk_tag_sub_tag);
        
    }
    m_fracture_curve_tags = fractures_to_micro; /// updating the micro fractures
    
}

/// Compute DFN physical tags (fractures and end points)
void TGeometryBuilder::ComputeDFNPhysicalTags(){
    
    int c_p_tag;
    int domain_tags = m_wellbore_region_tags.size();
    if (domain_tags == 0) {
        c_p_tag = 1;
    }else{
        c_p_tag = m_wellbore_region_tags[2]+1;
    }
    
    /// Functional physical tag for fractures
    int c = 1;
    int dim = 1;
    for (auto f : m_fracture_curve_tags) {
        std::vector<int> tags;
        for (auto micro_f : f.second) {
            tags.push_back(micro_f);
        }
        gmsh::model::addPhysicalGroup(dim, tags);
        std::stringstream f_name;
        f_name << "fracture_" << c;
        std::string name = f_name.str();
        gmsh::model::setPhysicalName(dim, c_p_tag, name);
        c++;
        c_p_tag++;
    }
}

void TGeometryBuilder::EmbedDFNInsideReservoir(){
    
    int domain_tags = m_wellbore_region_tags.size();
    if (domain_tags == 0) {
        return;
    }
    int domain_tag = m_wellbore_region_tags[0];
    std::vector<int> fractures_tags;
    for (auto f : m_fracture_curve_tags) {
        for (auto micro_f : f.second) {
            fractures_tags.push_back(micro_f);
        }
    }
    gmsh::model::mesh::embed(1,fractures_tags,2,domain_tag);
}

std::vector<int> TGeometryBuilder::ComputeAssociatedMicroFractures(std::pair<int,std::vector<int>> & fracture_data,std::map<int,std::vector<int>>  & fractures){

    std::vector<int> micro_fractures;
    std::pair<int,std::vector<int>> fracture = fracture_data;
    std::map<int,std::vector<int>>::iterator it;
    
    for (auto tag: fracture.second) {
        
        it = fractures.find(tag);
        bool is_member_Q = it != fractures.end();
        int n_micro_fractures = it->second.size();
        bool has_micro_fractures_Q = n_micro_fractures > 1;
        if (is_member_Q && has_micro_fractures_Q) {
            std::pair<int,std::vector<int>> fracture;
            fracture.first = it->first;
            fracture.second = it->second;
            std::vector<int> sub_micro_fractures = ComputeAssociatedMicroFractures(fracture, fractures);
            for (auto k : sub_micro_fractures) {
                micro_fractures.push_back(k);
            }
            
        }
        else{
            int micro_fracture_tag = it->second[0];
            micro_fractures.push_back(micro_fracture_tag);
        }
    }

    return micro_fractures;
}

bool TGeometryBuilder::ComputeFracturesIntersections(std::map<int,std::vector<int>> & objects, std::map<int,std::vector<int>>& tools, std::map<int,std::vector<int>>  & fractures){
    
    bool there_is_intersection_Q = false;
    std::pair<int, std::vector<int> > chunk_tag_sub_tag;
    for (auto object : objects) {
        for (auto tool: tools) {
            std::vector<gmsh::vectorpair> map_dim_tags;
            there_is_intersection_Q = IntersectFractures(object.first, tool.first, map_dim_tags);
            if (there_is_intersection_Q && (object.first != tool.first)) { /// required deletion and expansion for both vectors
                
                /// dropout object and tool on objects vector
                objects.erase(object.first);
                objects.erase(tool.first);
                /// dropout object and tool on tools vector
                tools.erase(tool.first);
                tools.erase(object.first);
                
                /// adjusting objects structure
                {

                    
                    /// Appending new micro fractures associated to the object
                    chunk_tag_sub_tag.first = object.first;
                    chunk_tag_sub_tag.second.resize(0);
                    for (auto micro_f: map_dim_tags[0]) { /// Associated to the object
                        chunk_tag_sub_tag.second.push_back(micro_f.second);
                    }
                    fractures.erase(object.first);
                    fractures.insert(chunk_tag_sub_tag);
                    
                    /// Appending new fractures
                    for (auto micro_f: map_dim_tags[0]) { /// Asspciated to the object
                        chunk_tag_sub_tag.second.resize(0);
                        chunk_tag_sub_tag.first = micro_f.second;
                        chunk_tag_sub_tag.second.push_back(micro_f.second);
                        fractures.insert(chunk_tag_sub_tag);
                        objects.insert(chunk_tag_sub_tag);
                    }
                    
                }
                
                /// adjusting tools structure
                {
                    
                    /// Appending new micro fractures associated to the object
                    chunk_tag_sub_tag.first = tool.first;
                    chunk_tag_sub_tag.second.resize(0);
                    for (auto micro_f: map_dim_tags[1]) { /// Associated to the tool
                        chunk_tag_sub_tag.second.push_back(micro_f.second);
                    }
                    fractures.erase(tool.first);
                    fractures.insert(chunk_tag_sub_tag);
                    
                    /// Appending new fractures
                    for (auto micro_f: map_dim_tags[1]) { /// Asspciated to the object
                        chunk_tag_sub_tag.second.resize(0);
                        chunk_tag_sub_tag.first = micro_f.second;
                        chunk_tag_sub_tag.second.push_back(micro_f.second);
                        fractures.insert(chunk_tag_sub_tag);
                        tools.insert(chunk_tag_sub_tag);
                    }
                    
                }
                return there_is_intersection_Q;
                break;
            }
            else{
                int nothing_to_do = 0;
            }
        }
    }
    
    return there_is_intersection_Q;
    
}

bool TGeometryBuilder::IntersectFractures(int object_tag, int tool_tag, std::vector<gmsh::vectorpair> & map_dim_tags){
    
    if (object_tag == tool_tag ) { // there is nothing to do for this case.
        return true;
    }
    
    gmsh::vectorpair dim_tag_object;
    dim_tag_object.push_back(std::make_pair(1, object_tag));
    gmsh::vectorpair dim_tag_tool;
    dim_tag_tool.push_back(std::make_pair(1, tool_tag));
    gmsh::vectorpair outDimTags;
    gmsh::model::occ::fragment(dim_tag_object, dim_tag_tool, outDimTags, map_dim_tags);
    gmsh::model::occ::synchronize(); /// expensive
    int n_data = outDimTags.size();
    if (n_data == 2) {
        return false;
    }else{
        return true;
    }
    
}

void TGeometryBuilder::DrawCirclePoints(double r, std::vector<double> x_center, std::vector<std::vector<double>> & points){
    
    std::vector<double> point(x_center);
    // p_center
    points.push_back(point);
    
    // p1
    point[0] = r + x_center[0];
    point[1] = 0 + x_center[1];
    points.push_back(point);
    
    // p2
    point[0] = 0 + x_center[0];
    point[1] = r + x_center[1];
    points.push_back(point);
    
    // p3
    point[0] = -r + x_center[0];
    point[1] = +0 + x_center[1];
    points.push_back(point);
    
    // p4
    point[0] = +0 + x_center[0];
    point[1] = -r + x_center[1];
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
    
    // p3
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
    
    int p1 = m_external_wire_point_tags[0];
    int p2 = m_external_wire_point_tags[1];
    int p3 = m_external_wire_point_tags[2];
    int p4 = m_external_wire_point_tags[3];
    
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

void TGeometryBuilder::DrawWellboreRegion(){
 
    /// construction by already defined wires    
    std::vector<int> wire_tags;
    wire_tags.push_back(m_internal_wire_curve_loop);
    wire_tags.push_back(m_external_wire_curve_loop);
    int wellbore_region_tag = gmsh::model::occ::addPlaneSurface(wire_tags);
    m_wellbore_region_tags.push_back(wellbore_region_tag);
}

void TGeometryBuilder::ComputeReservoirPhysicalTags(){
    
    int wellbore_tag = m_wellbore_region_tags[0];
    int c_p_tag = 1;
    int dim = 2;
    std::vector<int> tags;
    tags.push_back(wellbore_tag);
    gmsh::model::addPhysicalGroup(dim, tags);
    std::stringstream f_name;
    f_name << "reservoir";
    std::string name = f_name.str();
    gmsh::model::setPhysicalName(dim, c_p_tag, name);
    m_wellbore_region_physical_tags.push_back(c_p_tag);
    c_p_tag++;
    
    /// Physical tag for computational boundaries
    dim--;
    { ///  internal
        gmsh::model::addPhysicalGroup(dim, m_internal_wire_curve_tags);
        std::stringstream s_bc_name;
        s_bc_name << "bc_internal";
        std::string bc_name = s_bc_name.str();
        gmsh::model::setPhysicalName(dim, c_p_tag, bc_name);
        m_wellbore_region_physical_tags.push_back(c_p_tag);
        c_p_tag++;
    }
    
    { ///  external
        gmsh::model::addPhysicalGroup(dim, m_external_wire_curve_tags);
        std::stringstream s_bc_name;
        s_bc_name << "bc_external";
        std::string bc_name = s_bc_name.str();
        gmsh::model::setPhysicalName(dim, c_p_tag, bc_name);
        m_wellbore_region_physical_tags.push_back(c_p_tag);
        c_p_tag++;
    }
}


/// Divide wellbore boundary into n_points + 1 elements
void TGeometryBuilder::RefineWellboreElements(int n_points){
    for (auto bc: m_internal_wire_curve_tags) {
        gmsh::model::mesh::setTransfiniteCurve(bc, n_points);
    }
}

void TGeometryBuilder::RefineDFN(double omega, double size_ratio){
    
    const int dim = 1;
    /// Computing minimun fracture size
    double min_size = 183729;
    double max_size = 0.0;
    for (auto f : m_fracture_curve_tags) {
        for (auto micro_f : f.second) {
            double xmin,ymin,zmin;
            double xmax,ymax,zmax;
            gmsh::model::getBoundingBox(dim, micro_f, xmin, ymin, zmin, xmax, ymax, zmax);
            double squared_norm = (xmin-xmax)*(xmin-xmax)+(ymin-ymax)*(ymin-ymax)+(zmin-zmax)*(zmin-zmax);
            double norm = sqrt(squared_norm);
            if (norm <= min_size) {
                min_size = norm;
            }
            if (norm > max_size) {
                max_size = norm;
            }
        }
    }
    
    double avg_scale = (1.0 - omega) * max_size - omega * min_size;
    
    std::map<int,int> fracture_tag_n_points;
    for (auto f : m_fracture_curve_tags) {
        for (auto micro_f : f.second) {
            double xmin,ymin,zmin;
            double xmax,ymax,zmax;
            gmsh::model::getBoundingBox(dim, micro_f, xmin, ymin, zmin, xmax, ymax, zmax);
            double squared_norm = (xmin-xmax)*(xmin-xmax)+(ymin-ymax)*(ymin-ymax)+(zmin-zmax)*(zmin-zmax);
            double norm = sqrt(squared_norm);
            double ratio = norm / avg_scale;
            int n_points = floor(ratio);
            if (n_points == 0) {
                n_points = 2;
            }
            double max_size_ratio = norm / max_size;
            if (size_ratio > max_size_ratio) {
                fracture_tag_n_points.insert(std::make_pair(micro_f, n_points));
            }
        }
    }
    
    for (auto f_data: fracture_tag_n_points) {
        int fracture_tag = f_data.first;
        int n_points = f_data.second;
        gmsh::model::mesh::setTransfiniteCurve(fracture_tag, n_points);
    }
    
}

