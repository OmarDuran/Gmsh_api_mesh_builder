
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


/// Geometry construction from: http://en.wikipedia.org/wiki/Constructive_solid_geometry
void Constructive_solid_geometry();

/// Geometry that represents a 3D wellbore inside a irregular reservoir
void Wellbore_trajectory_3D();

int main()
{
    
    Wellbore_trajectory_3D();
//    Constructive_solid_geometry();
    return 0;
}

void InsertTheElements(TPZGeoMesh * gmesh);

void Wellbore_trajectory_3D(){
    
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    
    gmsh::model::add("boolean");
    
    
    gmsh::option::setNumber("Mesh.Algorithm", 6);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 1.0);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.0);
    
    {/// Functional but expensive
        
        /// wellbore radius
        REAL r_wb = 0.5;
        /// Defining the wire
        std::vector<double> p0 = {0., 0., 1.};
        std::vector<double> p1 = {0.3141592653589793, 0.23344536385590542, 0.9723699203976766};
        std::vector<double> p2 = {0.6283185307179586, 0.4539904997395468, 0.8910065241883679};
        std::vector<double> p3 = {0.9424777960769379, 0.6494480483301835, 0.760405965600031};
        std::vector<double> p4 = {1.2566370614359172, 0.8090169943749475, 0.587785252292473};
        std::vector<double> p5 = {1.5707963267948966, 0.9238795325112867, 0.38268343236508984};
        std::vector<double> p6 = {1.8849555921538759, 0.9876883405951377, 0.15643446504023115};
        std::vector<double> p7 = {2.199114857512855, 0.996917333733128, -0.07845909572784467};
        std::vector<double> p8 = {2.5132741228718345, 0.9510565162951535, -0.30901699437494756};
        std::vector<double> p9 = {2.827433388230814, 0.8526401643540923, -0.5224985647159488};
        std::vector<double> p10 = {3.141592653589793, 0.7071067811865476, -0.7071067811865475};
        std::vector<double> p11 = {3.455751918948773, 0.5224985647159489, -0.8526401643540922};
        std::vector<double> p12 = {3.7699111843077517, 0.30901699437494795, -0.9510565162951534};
        std::vector<double> p13 = {4.084070449666731, 0.07845909572784462, -0.996917333733128};
        std::vector<double> p14 = {4.39822971502571, -0.1564344650402303, -0.9876883405951378};
        std::vector<double> p15 = {4.71238898038469, -0.38268343236508967, -0.9238795325112868};
        std::vector<double> p16 = {5.026548245743669, -0.5877852522924734, -0.8090169943749472};
        std::vector<double> p17 = {5.340707511102648, -0.7604059656000306, -0.6494480483301841};
        std::vector<double> p18 = {5.654866776461628, -0.8910065241883678, -0.4539904997395469};
        std::vector<double> p19 = {5.969026041820607, -0.9723699203976764, -0.233445363855906};
        std::vector<double> p20 = {6.283185307179586, -1., -1.8369701987210297e-16};
        
        std::vector<std::vector<double>> wb_trajectory;
        wb_trajectory.push_back(p0);
        wb_trajectory.push_back(p1);
        wb_trajectory.push_back(p2);
        wb_trajectory.push_back(p3);
        wb_trajectory.push_back(p4);
        wb_trajectory.push_back(p5);
        wb_trajectory.push_back(p6);
        wb_trajectory.push_back(p7);
        wb_trajectory.push_back(p8);
        wb_trajectory.push_back(p9);
        wb_trajectory.push_back(p10);
        wb_trajectory.push_back(p11);
        wb_trajectory.push_back(p12);
        wb_trajectory.push_back(p13);
        wb_trajectory.push_back(p14);
        wb_trajectory.push_back(p15);
        wb_trajectory.push_back(p16);
        wb_trajectory.push_back(p17);
        wb_trajectory.push_back(p18);
        wb_trajectory.push_back(p19);
        wb_trajectory.push_back(p20);
    
        TGmshWellboreBuilder wb_builder(r_wb,wb_trajectory);
        if(1){ // volume option
//            gmsh::vectorpair wb_dim_tags = wb_builder.DrawWellboreBySections();
            gmsh::vectorpair wb_dim_tags = wb_builder.DrawWellboreByShell();

            double lx = -3.0;
            double ly = -6.0;
            double lz = -2.0;
            double lx_e = 10.0;
            double ly_e = 10.0;
            double lz_e = 4.0;
            int box_tag = gmsh::model::occ::addBox(lx,ly,lz, lx_e,ly_e,lz_e);
            
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
            gmsh::model::occ::cut({{3, box_tag}},wellbore_volume, ov, ovv);
        }

    }
    
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
        
        
        int n_tags_0d = entities_0d.size();
        int n_tags_1d = entities_1d.size();
        int n_tags_2d = entities_2d.size();
        int n_tags_3d = entities_3d.size();
        std::vector<int> tags_0d(n_tags_0d);
        std::vector<int> tags_1d(n_tags_1d);
        std::vector<int> tags_2d(n_tags_2d);
        std::vector<int> tags_3d(n_tags_3d);
        
        int c=0;
        for (auto i: entities_0d) {
            tags_0d[c] = i.second;
            c++;
        }
        
        c=0;
        for (auto i: entities_1d) {
            tags_1d[c] = i.second;
            c++;
        }
        
        c=0;
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
        int physical_tag_0d, physical_tag_1d, physical_tag_2d,physical_tag_3d;
        physical_tag_3d = gmsh::model::addPhysicalGroup(3, tags_3d);
        physical_tag_2d = gmsh::model::addPhysicalGroup(2, tags_2d);
        physical_tag_1d = gmsh::model::addPhysicalGroup(1, tags_1d);
        physical_tag_0d = gmsh::model::addPhysicalGroup(0, tags_0d);
        
    }
    
    /// Meshing directives
    gmsh::model::mesh::generate(3);
//    gmsh::model::mesh::refine();
    //    gmsh::model::mesh::setOrder(2);
    //gmsh::model::mesh::partition(4)
    gmsh::write("wellbore_geo.msh");
    
    
    
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
    std::string name = "wellbore";
    TPZGeoMeshBluider::PrintGeometry(gmesh,name);
    
    gmsh::finalize();
    
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
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.4);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.4);
    
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
    //    gmsh::model::mesh::refine();
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




