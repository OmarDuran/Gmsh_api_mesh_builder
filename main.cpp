
#include <iostream>
#include <fstream>

#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshBluider.h"


#include <cmath>
#include <thread>
#include "gmsh.h"

/// Wellbore builder
#include "TGmshWellboreBuilder.h"
#include "TPZGmshReader.h"

/// Welbore 2D with fractures intersection automatically computed
#include "TGeometryBuilder.h"


/// Geometry construction from: http://en.wikipedia.org/wiki/Constructive_solid_geometry
void Constructive_solid_geometry();

/// Geometry that represents a 3D wellbore inside a irregular reservoir
void Wellbore_trajectory_3D();

/// Geometry that represents a 2D wellbore inside a irregular reservoir with line fractures
void Wellbore_2D_with_factures();

/// Read a wellbore trajectory defined with xyz data.
std::vector<std::vector<double>> ReadWellTrajectory(std::string & file_name, int n_data);

int main()
{
    Wellbore_2D_with_factures();
//    Wellbore_trajectory_3D();
//    Constructive_solid_geometry();
    return 0;
}

void Wellbore_2D_with_factures(){
    
    gmsh::initialize(); /// Mandatory
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("boolean");
    gmsh::option::setNumber("Mesh.Algorithm", 6);
    
    int f1pt_i = gmsh::model::occ::addPoint(0.25,0,0);
    int f1pt_e = gmsh::model::occ::addPoint(1.5,1,0);

    int f2pt_i = gmsh::model::occ::addPoint(1.0,1,0);
    int f2pt_e = gmsh::model::occ::addPoint(1.25,-0.5,0);

    int f1 = gmsh::model::occ::addLine(f1pt_i, f1pt_e);
    int f2 = gmsh::model::occ::addLine(f2pt_i, f2pt_e);
    
    double x,y,z;
    x = y = z = 0.0;
    double r = 2.0;
    double r_w = 0.127;
    int circle = gmsh::model::occ::addCircle(x, y, z, r);
    int circle_rw = gmsh::model::occ::addCircle(x, y, z, r_w);
//    gmsh::model::occ::synchronize();
    std::vector<int> curve_tags;
    curve_tags.push_back(circle);
    curve_tags.push_back(circle_rw);
    int res_loop = gmsh::model::geo::addCurveLoop(curve_tags);
    std::vector<int> curve_loop;
    curve_loop.push_back(res_loop);
//    int wire = gmsh::model::occ::addWire(curve_tags);
    int area = gmsh::model::geo::addSurfaceFilling(curve_loop);
//    int area = gmsh::model::occ::addSurfaceFilling(res_loop);
//     int area = gmsh::model::occ::addPlaneSurface(curve_tags);
//    gmsh::model::occ::addSurfaceFilling
    gmsh::vectorpair dim_tag_object;
    dim_tag_object.push_back(std::make_pair(1, f1));
    gmsh::vectorpair dim_tag_tool;
    dim_tag_tool.push_back(std::make_pair(1, f2));
    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;
//    gmsh::model::occ::intersect(dim_tag_object, dim_tag_tool, outDimTags, outDimTagsMap);
    gmsh::model::occ::fragment(dim_tag_object, dim_tag_tool, outDimTags, outDimTagsMap);
    
    gmsh::model::occ::synchronize();
    std::vector<int> fractures;
    for (auto i: outDimTags) {
        fractures.push_back(i.second);
    }
    gmsh::model::mesh::embed(1,fractures,2,area);
    
    gmsh::model::mesh::field::setAsBoundaryLayer(circle);
//    std::string f1_name="fracture_1";
//    gmsh::model::setPhysicalName(1, f1, f1_name);
//    std::string f2_name="fracture_2";
//    gmsh::model::setPhysicalName(1, f2, f2_name);
    
    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(2);
//    gmsh::model::mesh::setRecombine(2, area);
//    gmsh::model::mesh::generate(2);
    //    gmsh::model::mesh::setOrder(2);
    //gmsh::model::mesh::partition(4)
    gmsh::write("wellbore_2D.msh");
    
    gmsh::finalize();
}


void InsertTheElements(TPZGeoMesh * gmesh);

void Wellbore_trajectory_3D(){
    
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    
    gmsh::model::add("boolean");
    
    
    gmsh::option::setNumber("Mesh.Algorithm", 6);

    
    {/// Functional but expensive
        
        /// wellbore radius
        REAL r_wb = 15.0;
        REAL characteristic_length = 1.0*r_wb;
        
        gmsh::option::setNumber("Mesh.CharacteristicLengthMin", characteristic_length);
        gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 50.0);
        
        int n_data = 50;
        std::string file_name = "producer_trajectory.txt";
        std::vector<std::vector<double>> wb_trajectory = ReadWellTrajectory(file_name,n_data);
        TGmshWellboreBuilder wb_builder(r_wb,wb_trajectory);
        if(1){ // volume option
            wb_builder.SetCharacteristicLength(characteristic_length);
            gmsh::vectorpair wb_dim_tags = wb_builder.DrawWellboreByShell();
            
            std::vector<double> x_min = {-250, -250., -50.};
            std::vector<double> dx = {990, 638.413, 110.653};
            int box_tag = gmsh::model::occ::addBox(x_min[0],x_min[1],x_min[2], dx[0],dx[1],dx[2]);
            
            std::vector<std::pair<int, int> > ov;
            std::vector<std::pair<int, int> > wellbore_volume;
            std::vector<std::vector<std::pair<int, int> > > ovv;
            int n_volumes = wb_dim_tags.size();
            std::vector<std::pair<int, int> > sector_volume;
            sector_volume.push_back(wb_dim_tags[0]);
            for (int i =1 ; i < n_volumes; i++) {
                std::vector<std::pair<int, int> > next_sector_volume;
                next_sector_volume.push_back(wb_dim_tags[i]);
                gmsh::model::occ::fuse(sector_volume, next_sector_volume, wellbore_volume, ovv);
                sector_volume = wellbore_volume;
            }
            if (n_volumes == 1) {
                wellbore_volume = wb_dim_tags;
            }

            int tag = 0;
            bool removeObject = true;
            bool removeTool = false;
            /// computing wellbore - reservoir intersection volume and keep the tool (the reservoir volume)
            std::vector<std::pair<int, int> > wellbore_reservoir_intersection;
            gmsh::model::occ::intersect(wellbore_volume,{{3, box_tag}}, wellbore_reservoir_intersection, ovv, tag, removeObject, removeTool);
            removeObject = true;
            removeTool = true;
            /// computing holled reservoir by cutting wellbore - reservoir intersection volume and delete the tool.
            gmsh::model::occ::cut({{3, box_tag}}, wellbore_reservoir_intersection, ov, ovv, tag, removeObject, removeTool);
        }

    }
    
    gmsh::model::occ::synchronize();
    
    /// physiscal taggging
    // Get all elementary entities in the model
    std::vector<std::pair<int, int> > entities_3d;
    gmsh::model::getEntities(entities_3d,3);
    int n_tags_3d = entities_3d.size();
    std::vector<int> reservoir_3d(n_tags_3d);
    int c=0;
    for (auto i: entities_3d) {
        reservoir_3d[c] = i.second;
        c++;
    }
    
    /// Physical for the reservoir
    int physical_tag_3d;
    physical_tag_3d = gmsh::model::addPhysicalGroup(3, reservoir_3d);
    std::string reservoir("reservoir");
    gmsh::model::setPhysicalName(3, reservoir_3d[0], reservoir);
    
    /// Computing the reservoir boundaries
    gmsh::vectorpair bc_dim_tags;
    bool combined = false;
    bool oriented = false;
    bool recursive = false;
    gmsh::model::getBoundary(entities_3d, bc_dim_tags, combined, oriented, recursive);
    
    c = 0; //  The first 6 bcs are the external ones
    std::vector<int> bc_reservoir_external_tags(6);
    for (int i = 0;  i < 6; i++) {
        bc_reservoir_external_tags[c] = bc_dim_tags[i].second;
        c++;
    }
    std::vector<int> bc_W(1),bc_S(1),bc_T(1),bc_N(1),bc_B(1),bc_E(1);
    bc_W[0] = bc_reservoir_external_tags[0]; // W -> 0
    gmsh::model::addPhysicalGroup(2, bc_W);
    gmsh::model::setPhysicalName(2, 2, "BCWest");
    
    bc_S[0] = bc_reservoir_external_tags[1]; // S -> 1
    gmsh::model::addPhysicalGroup(2, bc_S);
    gmsh::model::setPhysicalName(2, 3, "BCSouth");
    
    bc_T[0] = bc_reservoir_external_tags[2]; // T -> 2
    gmsh::model::addPhysicalGroup(2, bc_T);
    gmsh::model::setPhysicalName(2, 4, "BCTop");
    
    bc_N[0] = bc_reservoir_external_tags[3]; // N -> 3
    gmsh::model::addPhysicalGroup(2, bc_N);
    gmsh::model::setPhysicalName(2, 5, "BCNorth");
    
    bc_B[0] = bc_reservoir_external_tags[4]; // B -> 4
    gmsh::model::addPhysicalGroup(2, bc_B);
    gmsh::model::setPhysicalName(2, 6, "BCBottom");
    
    bc_E[0] = bc_reservoir_external_tags[5]; // E -> 5
    gmsh::model::addPhysicalGroup(2, bc_E);
    gmsh::model::setPhysicalName(2, 7, "BCEast");
    
    int n_bc = bc_dim_tags.size();
    
    //  Cased Hole section
    int n_cased_bc = 28;
    c = 0;
    std::vector<int> bc_cased_wellbore_tags(n_cased_bc);
    for (int i = 6;  i < n_cased_bc + 6; i++) {
        bc_cased_wellbore_tags[c] = bc_dim_tags[i].second;
        c++;
    }
    gmsh::model::addPhysicalGroup(2, bc_cased_wellbore_tags);
    gmsh::model::setPhysicalName(2, 8, "BCCasedHole");
    
    //  Open Hole section
    int n_open_bc = n_bc - n_cased_bc - 6 - 1;
    c = 0;
    std::vector<int> bc_open_wellbore_tags(n_open_bc);
    for (int i = n_cased_bc + 6;  i < n_bc; i++) {
        bc_open_wellbore_tags[c] = bc_dim_tags[i].second;
        c++;
    }
    gmsh::model::addPhysicalGroup(2, bc_open_wellbore_tags);
    gmsh::model::setPhysicalName(2, 9, "BCOpenHole");
    
    //  The wellbore lid
    int n_lid_bc = n_bc - n_open_bc - n_cased_bc - 6;
    c = 0;
    std::vector<int> bc_lid_wellbore_tags(n_lid_bc);
    for (int i = n_open_bc + n_cased_bc + 6;  i < n_bc; i++) {
        bc_lid_wellbore_tags[c] = bc_dim_tags[i].second;
        c++;
    }
    gmsh::model::addPhysicalGroup(2, bc_lid_wellbore_tags);
    gmsh::model::setPhysicalName(2, 10, "BCWellboreLid");
    
    /// Meshing directives
    gmsh::model::mesh::generate(3);
//    gmsh::model::mesh::refine();
    //    gmsh::model::mesh::setOrder(2);
    //gmsh::model::mesh::partition(4)
    gmsh::write("wellbore_geo.msh");
    
    std::string geometry_file = "wellbore_geo.msh";
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.0");
    TPZGeoMesh * gmesh = Geometry.GeometricGmshMesh(geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
    std::string vtk_file = "wellbore_geo";
    TPZGeoMeshBluider::PrintGeometry(gmesh,vtk_file);
    
//    /// Building the geomesh object
//    TPZGeoMesh * gmesh = new TPZGeoMesh;
//    int mesh_dimension  = gmsh::model::getDimension();
//    gmesh->SetDimension(mesh_dimension);
//
//    //// Gathering required information for constuction of a TPZGeoMesh object
//    {
//        std::vector<int> node_identifiers;
//        std::vector<double> coord;
//        std::vector<double> parametricCoord;
//        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
//        TPZGeoMeshBluider::InsertNodes(gmesh, node_identifiers, coord);
//        InsertTheElements(gmesh);
//    }
//
//
//    gmesh->BuildConnectivity();
//    std::string name = "wellbore";
//    TPZGeoMeshBluider::PrintGeometry(gmesh,name);
    
    gmsh::finalize();
    
}

std::vector<std::vector<double>> ReadWellTrajectory(std::string & file_name, int n_data){
    
    std::ifstream in(file_name.c_str());
    std::vector<std::vector<double>> trajectory(n_data);
    double x,y,z;
    int c = 0;
    while(in)
    {
        in >> x;
        in >> y;
        in >> z;
        std::vector<double> p = {x,y,z};
        trajectory[c] = p;
        c++;
        if (c == n_data) {
            break;
        }
    }
    return trajectory;
}

void InsertTheElements(TPZGeoMesh * gmesh){
    
    
    std::vector<std::pair<int, int> > dim_to_physical_groups;
    gmsh::model::getPhysicalGroups(dim_to_physical_groups);
    
    std::vector<std::pair<int, int> > entities_0d;
    std::vector<std::pair<int, int> > entities_1d;
    std::vector<std::pair<int, int> > entities_2d;
    std::vector<std::pair<int, int> > entities_3d;
    gmsh::model::getEntities(entities_0d,0);
    gmsh::model::getEntities(entities_1d,1);
    gmsh::model::getEntities(entities_2d,2);
    gmsh::model::getEntities(entities_3d,3);
    
    /// inserting the elements
    {
        
        for (auto group: dim_to_physical_groups) {
            
            int dim = group.first;
            int physical_identifier = group.second;
            
            std::vector<std::pair<int, int> > dim_to_entities;
            switch (dim) {
                case 0:
                    dim_to_entities = entities_0d;
                    break;
                case 1:
                    dim_to_entities = entities_1d;
                    break;
                case 2:
                    dim_to_entities = entities_2d;
                    break;
                case 3:
                    dim_to_entities = entities_3d;
                    break;
                default:
                    break;
            }
            
            for (auto entity: dim_to_entities) {
                
                int tag = entity.second;
                std::vector<int> group_element_types;
                std::vector<std::vector<int> > group_element_identifiers;
                std::vector<std::vector<int> > group_node_identifiers;
                gmsh::model::mesh::getElements(group_element_types, group_element_identifiers, group_node_identifiers, dim, tag);
                int n_types = group_element_types.size();
                for (int itype = 0; itype < n_types; itype++){
                    int el_type = group_element_types[itype];
                    int n_nodes = TPZGeoMeshBluider::GetNumberofNodes(el_type);
                    std::vector<int> node_identifiers(n_nodes);
                    int n_elements = group_element_identifiers[itype].size();
                    for (int iel = 0; iel < n_elements; iel++) {
                        int el_identifier = group_element_identifiers[itype][iel]-1;
                        for (int inode = 0; inode < n_nodes; inode++) {
                            node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
                        }
                        TPZGeoMeshBluider::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
                    }
                    
                }
            }
            
        }
        
    }
    
}


void Constructive_solid_geometry(){
    
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    
    gmsh::model::add("boolean");
    
    
    gmsh::option::setNumber("Mesh.Algorithm", 6);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.25);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.25);
    
    double R = 1.4, Rs = R*.7, Rt = R*1.25;
    
    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    gmsh::model::occ::addBox(-R,-R,-R, 2*R,2*R,2*R, 1);
    gmsh::model::occ::addSphere(0,0,0,Rt, 2);
    gmsh::model::occ::intersect({{3, 1}}, {{3, 2}}, ov, ovv, 3);
    gmsh::model::occ::addCylinder(-2*R,0,0, 4*R,0,0, Rs, 4);
    gmsh::model::occ::addCylinder(0,-2*R,0, 0,4*R,0, Rs, 5);
    gmsh::model::occ::addCylinder(0,0,-2*R, 0,0,4*R, Rs, 6);
    gmsh::model::occ::fuse({{3, 4}, {3, 5}}, {{3, 6}}, ov, ovv, 7);
    gmsh::model::occ::cut({{3, 3}}, {{3, 7}}, ov, ovv, 8);
    
    gmsh::model::occ::synchronize();
    
    /// physiscal taggging
    {
        // Get all elementary entities in the model
        std::vector<std::pair<int, int> > entities_0d;
        std::vector<std::pair<int, int> > entities_1d;
        std::vector<std::pair<int, int> > entities_2d;
        std::vector<std::pair<int, int> > entities_3d;
        gmsh::model::getEntities(entities_0d,0);
        gmsh::model::getEntities(entities_1d,1);
        gmsh::model::getEntities(entities_2d,2);
        gmsh::model::getEntities(entities_3d,3);
        
        
        int n_tags_2d = entities_2d.size();
        int n_tags_3d = entities_3d.size();
        std::vector<int> tags_2d(n_tags_2d);
        std::vector<int> tags_3d(n_tags_3d);
        int c=0;
        for (auto i: entities_2d) {
            tags_2d[c] = i.second;
            c++;
        }
        
        c=0;
        for (auto i: entities_3d) {
            tags_3d[c] = i.second;
            c++;
        }
        
        /// tag just boundary and volumetric mesh
        int physical_tag_2d,physical_tag_3d;
        physical_tag_3d = gmsh::model::addPhysicalGroup(3, tags_3d);
        physical_tag_2d = gmsh::model::addPhysicalGroup(2, tags_2d);
        
    }
    
    /// Meshing directives
    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::refine();
    gmsh::model::mesh::refine();
    //    gmsh::model::mesh::setOrder(2);
    //gmsh::model::mesh::partition(4)
    gmsh::write("geometry.msh");
    
    
    
    /// Building the geomesh object
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    int mesh_dimension  = gmsh::model::getDimension();
    gmesh->SetDimension(mesh_dimension);
    
    //// Gathering required information for constuction of a TPZGeoMesh object
    {
        std::vector<int> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
        TPZGeoMeshBluider::InsertNodes(gmesh, node_identifiers, coord);
        
        InsertTheElements(gmesh);
    }
    
    
    gmesh->BuildConnectivity();
    std::string name = "geometry";
    TPZGeoMeshBluider::PrintGeometry(gmesh,name);
    
    gmsh::finalize();
}




