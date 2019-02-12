
#include <iostream>
#include <fstream>

#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshBluider.h"


#include <cmath>
#include <thread>
#include "gmsh.h"




// Geometry construction from: http://en.wikipedia.org/wiki/Constructive_solid_geometry
int main()
{
    
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
    
    
    gmesh->BuildConnectivity();
    TPZGeoMeshBluider::PrintGeometry(gmesh);
    
    gmsh::finalize();
    
    return 0;
}





