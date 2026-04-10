#ifndef BASEMESH_H
#define BASEMESH_H
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/IO/MeshIO.hh>
#include"OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include"vertex_remove_data.h"
#include<optional>
#include"Maps.h"
#include<Eigen/Dense>
#include<tuple>
#include<iostream>
#include<queue>
#include <set>

using namespace Eigen;
using namespace std;

typedef OpenMesh::TriMesh_ArrayKernelT<>::VertexHandle VertexHandle;
typedef OpenMesh::TriMesh_ArrayKernelT<>::FaceHandle FaceHandle;
typedef OpenMesh::TriMesh_ArrayKernelT<>::HalfedgeHandle  HalfedgeHandle;
typedef OpenMesh::TriMesh_ArrayKernelT<>::Point Point;
typedef OpenMesh::TriMesh_ArrayKernelT<>::Normal Normal;
typedef OpenMesh::Vec2d Vec2d;
typedef OpenMesh::Vec3d Vec3d;
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;


using Barycoor = std::array<std::pair<unsigned int,double>,3>;
struct mapstraits:public OpenMesh::DefaultTraits
{
public:
    FaceTraits {
        double area;
        std::vector<unsigned int> vertices;
    };

    VertexTraits{
        unsigned int ori_id;
        double curvature;
        bool canberemoved = false;
        bool isdeleted = false;
        double weight;
        double ringArea;
        //update times
        unsigned int c=1;
        std::optional<Barycoor>barycoor;

		//feature edge flag
        int feature_edge_count = 0;
        bool is_dart = false;
        bool is_corner = false;
    };

    EdgeTraits{
        bool is_feature = false; //feature edge;
	};
};

// template <typename Trait = OpenMesh::DefaultTraits>
// class maps:public OpenMesh::TriMesh_ArrayKernelT<Trait>
// {
// public:
//     maps(){};
//     ~maps(){};
// };

class BaseMesh:public maps<mapstraits>
{
public:
    struct CompareFun {
    private:
        BaseMesh* _mesh;

    public:
        explicit CompareFun(BaseMesh* mesh) : _mesh(mesh){};
        bool operator()(const VertexHandle& L, const VertexHandle& R) {
            return _mesh->data(L).weight < _mesh->data(R).weight;
        }
    };
public:
    BaseMesh():removal_P(CompareFun(this)) {};
    ~BaseMesh()override {

    };
    BaseMesh(const BaseMesh& mapMesh) = default;

    struct EdgeKey {
        unsigned int v0, v1;
        EdgeKey(unsigned int a, unsigned int b) {
            v0 = std::min(a, b);
            v1 = std::max(a, b);
        }
        bool operator<(const EdgeKey& other) const {
            return (v0 < other.v0) || (v0 == other.v0 && v1 < other.v1);
        }
    };





public:

    void Initialization();
	void detect_feature_edges(double threshold_degree);
    bool check_uv_flip(std::array<Vec2d, 3>uv_face);
    bool check_all_orifaces_uv_flip(std::vector<FaceHandle>ring_faces, std::vector<std::pair<VertexHandle, Vec2d>> coordinates);

    double cal_area(const Vec2d p1,const Vec2d p2,const Point P0,const Point P1,const Point P2);// a face area and 1-ring area sum
    double cal_Angles(Point &v0,Point &v1,Point& v2);
    bool In_2Dtriangle(const Vec2d& point2d,const std::array<Vec2d,3>&f2d);

    std::tuple<double,double,double>cal_bar_2dcoord(const Vec2d& point,const std::array<Vec2d,3>&f2d);
	void compute_curvature(const VertexHandle& v_h);
	void compute_ringArea(const VertexHandle& v_h);
    void initial_weights(double lamda);
    void map2plane(VertexHandle v_h,std::vector<std::pair<VertexHandle,Vec2d>>&coordinates);
    void CGAL_CDT(const std::vector<std::pair<VertexHandle, Vec2d>>& coordinates, std::vector<std::array<VertexHandle, 3>>& faces,
        std::pair<VertexHandle,VertexHandle> feature_constraint = {VertexHandle(-1),VertexHandle(-1)});
    void Retriangle(const VertexHandle &deletevertex,unsigned int remaining_nfs);
    void construt_bm(unsigned int remaining_nfs);
    void save_pts();
    void save_feature_lines(const std::string& filename);

    void MidSubdivision(unsigned int level);
    void remesh_(unsigned int level);

    bool areElementsUnique(const std::vector<std::array<unsigned int, 3>>& vec) {
        std::set<std::array<unsigned int, 3>> seen(vec.begin(), vec.end());
        return seen.size() == vec.size();
    }
    




private:
    double maxringArea;
    double maxCurvature;


    std::vector<Eigen::Vector3d> uvs_;
    unsigned int nv_;
	unsigned int nf_;
    std::vector<FaceHandle> originFaces;
    std::vector<VertexHandle>original_vhs;
    std::map<FaceHandle,std::vector<VertexHandle>>originFaceVertices;
    std::map < VertexHandle, std::vector<FaceHandle>>originVertexFaces;
    std::map<VertexHandle, int>ori_vids;
    std::map<int, Point>positions;
    std::vector<Barycoor> barycoordinates;
    std::vector<bool>is_deleted;
    std::priority_queue<VertexHandle, std::vector<VertexHandle>, CompareFun>removal_P;
    std::map < std::array<unsigned int, 3>, unsigned int> Face2IM;
    std::vector<vertex_remove_data>vremoveIM;
	std::vector<std::array<unsigned int, 3>> all_orifaces,test_faces;

public:
    unsigned int current_nfs;
};












#endif // BASEMESH_H
