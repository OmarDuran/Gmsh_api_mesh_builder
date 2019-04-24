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
#include <set>
#include <cmath>
#include <thread>
#include <algorithm>
#include "gmsh.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;

/// structure to reconstruc base fracture and associated microfractures
class EntityBinaryTree {
    
public:
    
    /// entity tag
    int m_entity_tag = -1;
    
    /// Left fracture
    EntityBinaryTree * m_left_fracture = nullptr;
    
    /// Right fracture
    EntityBinaryTree * m_right_fracture = nullptr;
    
    /// Default constructo
    EntityBinaryTree(){
        m_entity_tag = -1;
        m_left_fracture = nullptr;
        m_right_fracture = nullptr;
    }
    
    /// Constructor based on a new tag
    EntityBinaryTree(int entity_tag){
        m_entity_tag = entity_tag;
        m_left_fracture = nullptr;
        m_right_fracture = nullptr;
    }
    
    /// Copy constructor
    EntityBinaryTree(const EntityBinaryTree & other){
        m_entity_tag        = other.m_entity_tag;
        m_left_fracture     = other.m_left_fracture;
        m_right_fracture    = other.m_right_fracture;
    }
    
    /// Assignement constructor
    const EntityBinaryTree & operator=(const EntityBinaryTree & other){
        /// check for self-assignment
        if(&other == this){
            return *this;
        }
        m_entity_tag        = other.m_entity_tag;
        m_left_fracture     = other.m_left_fracture;
        m_right_fracture    = other.m_right_fracture;
        return *this;
    }
    
    /// Assignement constructor overload for pointer data
    const EntityBinaryTree & operator=(const EntityBinaryTree * other){
        /// check for self-assignment
        if(other == this){
            return *this;
        }
        m_entity_tag        = other->m_entity_tag;
        m_left_fracture     = other->m_left_fracture;
        m_right_fracture    = other->m_right_fracture;
        return *this;
    }
    
    /// Function to get the count of leaf nodes in a binary tree
    static unsigned int getLeafCount(EntityBinaryTree * tree)
    {
        if(tree == NULL){
            return 0;
        }
        if(tree->m_left_fracture == NULL && tree->m_right_fracture == NULL){
            return 1;
        }
        else{
            unsigned int c = getLeafCount(tree->m_left_fracture)+
            getLeafCount(tree->m_right_fracture);
            return c;
        }
    }
    
    /// Function to get the leaves in a binary tree
    static std::vector<int> getLeaves(EntityBinaryTree * tree)
    {
        std::vector<int> tags;
        if(tree == NULL){
            return tags;
        }
        if(tree->m_left_fracture == NULL && tree->m_right_fracture == NULL){
            tags.push_back(tree->m_entity_tag);
            return tags;
        }else if (tree->m_left_fracture != NULL && tree->m_right_fracture == NULL){
            std::vector<int> left_tags = getLeaves(tree->m_left_fracture);
            for (auto i : left_tags) {
                tags.push_back(i);
            }
            return tags;
        }
        else{
            std::vector<int> left_tags = getLeaves(tree->m_left_fracture);
            std::vector<int> right_tags = getLeaves(tree->m_right_fracture);
            for (auto i : left_tags) {
                tags.push_back(i);
            }
            for (auto i : right_tags) {
                tags.push_back(i);
            }
            return tags;
        }
    }
    
    /// Function to set the leaves in a binary tree
    static void setLeaves(EntityBinaryTree * tree, EntityBinaryTree left, EntityBinaryTree right)
    {
        if(tree == NULL){
            return;
        }
        if(tree->m_left_fracture == NULL && tree->m_right_fracture == NULL){
            tree->m_left_fracture = new EntityBinaryTree(left);
            tree->m_right_fracture = new EntityBinaryTree(right);
            return;
        }
        else if(tree->m_left_fracture != NULL && tree->m_right_fracture == NULL){
            tree->m_left_fracture = new EntityBinaryTree(left);
            return;
        }
        else{
            setLeaves(tree->m_left_fracture,left,right);
            setLeaves(tree->m_right_fracture,left,right);
            return;
        }
    }
    
};

// structure to reconstruc base fracture and associated microfractures
class EntityList {
    
public:
    
    /// entity tag
    int m_entity_tag = -1;
    
    /// Left fracture
    std::vector<EntityList> * m_micro_fractures = nullptr;
    
    /// Default constructor
    EntityList(){
        m_entity_tag = -1;
        m_micro_fractures = nullptr;
    }

    /// Constructor based on a new tag
    EntityList(int entity_tag){
        m_entity_tag = entity_tag;
        m_micro_fractures = nullptr;
    }

    /// Copy constructor
    EntityList(const EntityList & other){
        m_entity_tag        = other.m_entity_tag;
        m_micro_fractures   = other.m_micro_fractures;
    }

    /// Assignement constructor
    const EntityList & operator=(const EntityList & other){
        /// check for self-assignment
        if(&other == this){
            return *this;
        }
        m_entity_tag        = other.m_entity_tag;
        m_micro_fractures   = other.m_micro_fractures;
        return *this;
    }

    /// Assignement constructor overload for pointer data
    const EntityList & operator=(const EntityList * other){
        /// check for self-assignment
        if(other == this){
            return *this;
        }
        m_entity_tag        = other->m_entity_tag;
        m_micro_fractures   = other->m_micro_fractures;
        return *this;
    }

    /// Function to get the count of leaf nodes in a binary tree
    static unsigned int getLeafCount(EntityList * tree)
    {
        if(tree == NULL){
            return 0;
        }
        if(tree->m_micro_fractures == NULL){
            return 1;
        }
        else{
            unsigned int c = 0;
            for (auto i : *tree->m_micro_fractures) {
                c += getLeafCount(&i);
            }
            return c;
        }
    }

    /// Function to get the leaves in a binary tree
    static std::vector<int> getLeaves(EntityList * tree)
    {
        std::vector<int> tags;
        if(tree == NULL){
            return tags;
        }
        if(tree->m_micro_fractures == NULL){
            tags.push_back(tree->m_entity_tag);
            return tags;
        }
        else{
            for (auto i : *tree->m_micro_fractures) {
                std::vector<int> miroc_f_tags = getLeaves(&i);
                for (auto mf : miroc_f_tags) {
                    tags.push_back(mf);
                }
            }
            return tags;
        }
    }

    /// Function to set the leaves in a binary tree
    static void setLeaves(EntityList * tree, std::vector<EntityList> micro_fractures)
    {
        if(tree == NULL){
            return;
        }
        if(tree->m_micro_fractures == NULL){
            tree->m_micro_fractures = new std::vector<EntityList>(micro_fractures);
            return;
        }
        else{
            for (auto i : *tree->m_micro_fractures) {
                setLeaves(&i,micro_fractures);
            }
            return;
        }
    }
    
};

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
    
    /// Vector of external wire point tags
    std::vector<Point> m_internal_wire_cgal_points;
    
    /// Vector of internal wire curve tags
    std::vector<int> m_internal_wire_curve_tags;
    
    /// Internal wire loop tag
    int m_internal_wire_curve_loop;
    
    /// Vector of external wire point tags
    std::vector<int> m_external_wire_point_tags;
    
    /// Vector of external wire point tags
    std::vector<Point> m_external_wire_cgal_points;
    
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
    std::map<int, std::vector<int> > m_boundary_curve_tags;
    
    /// Map of base fracture curve tag to subfracture curve tags
    std::map<int, std::vector<int> > m_fracture_curve_tags;
    
    /// Map of internal boundary tree curve tags
    std::map<int, EntityBinaryTree > m_internal_boundary_tree_tags;
    
    /// Map of external boundary tree curve tags
    std::map<int, EntityBinaryTree > m_external_boundary_tree_tags;
    
    /// Map of fracture tree curve tags
    std::map<int, EntityBinaryTree > m_fracture_tree_tags;
    
    /// Map of internal boundary tree curve tags
    std::map<int, EntityList > m_internal_boundary_list_tags;
    
    /// Map of external boundary tree curve tags
    std::map<int, EntityList > m_external_boundary_list_tags;
    
    /// Map of fracture tree curve tags
    std::map<int, EntityList > m_fracture_list_tags;
    
    /// Entity tag for the wellbore region and boundaries
    std::vector<int>  m_wellbore_region_tags;
    
    /// Physical tag for the wellbore region and boundaries
    std::vector<int>  m_wellbore_region_physical_tags;
    
    /// Entity tag for the fractures that are outside of computational domain
    std::set<int>  m_fracture_tags_to_remove;
    
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
    
    /// Draw points and DFN's base fractures 
    void DrawBaseFractures();
    
    /// Computes the intersection between two fractures and return false when the intersection is empty.
    bool IntersectLines(int object_tag, int tool_tag, std::vector<gmsh::vectorpair> & map_dim_tags);
    
    bool ComputeFracturesIntersections(std::map<int,std::vector<int>> & objects, std::map<int,std::vector<int>>& tools, std::map<int,std::vector<int>>  & fractures);
    
    std::vector<int> ComputeAssociatedMicroFractures(std::pair<int,std::vector<int>> & fracture_data,std::map<int,std::vector<int>>  & fractures);
    
    bool ComputeFractureBCIntersections(std::map<int,std::vector<int>> & objects, std::map<int,std::vector<int>>& tools, std::map<int,std::vector<int>>  & fractures, std::map<int,std::vector<int>>  & boundaries);
    
    

    
//    /// Helper function that allocates a new node with the given data and NULL left and right pointers.
//    EntityBinaryTree node * newNode(int data)
//    {
//        struct node* node = (struct node*)
//        malloc(sizeof(struct node));
//        node->data = data;
//        node->left = NULL;
//        node->right = NULL;
//
//        return(node);
//    }
    
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
    
    /// Draws a computable DFN
    void DrawDFN();
    
    /// Compute DFN physical tags (fractures and end points)
    void ComputeDFNPhysicalTags();
    
    /// Embed the DFN inside reservoir
    void EmbedDFNInsideReservoir();
    
    /// Divide wellbore boundary at least in n_points
    void RefineWellboreElements(double omega, double size_ratio, int n_base_points);

    /// Divide DFN fractrures according to the weigth omega.
    void RefineDFN(double omega, double size_ratio, int n_base_points);
    
    const bool IsMemeberQ(Point pt, Polygon_2 & polygon) const
    {
        bool is_member_Q = false;
        switch(polygon.bounded_side(pt)) {
            case CGAL::ON_BOUNDED_SIDE :
                is_member_Q = true;
                break;
            case CGAL::ON_BOUNDARY:
                is_member_Q = true;
                break;
            case CGAL::ON_UNBOUNDED_SIDE:
                is_member_Q = false;
                break;
        }
        return is_member_Q;
    }
    
    const void PrintPolygon(Polygon_2 & polygon) const
    {
        int n_vertices = polygon.size();
        for(std::size_t i = 0; i < n_vertices; i++){
            Point p = polygon.vertex(i);
            std::cout << "Coordinate number " << i << std::endl;
            std::cout << "x = " <<  p.x() << std::endl;
            std::cout << "y = " <<  p.y() << std::endl;
        }
    }
    
};


#endif /* TGeometryBuilder_h */
