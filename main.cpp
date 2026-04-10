#include "BaseMesh.h"
#include<iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace std;

int main()
{
    BaseMesh mesh;
    BaseMesh originalmesh;
    OpenMesh::IO::Options opt;
    if (!OpenMesh::IO::read_mesh(mesh, "D:/workplace/MAPS/data/39050__sf.obj", opt))
    {
        std::cerr<<"Cannot open file"<<std::endl;
        exit(1);
    }
    else
    {
        std::cout<<"Successfully read mesh"<<std::endl;
        std::cerr << "Vertices: " << mesh.n_vertices() << ", Faces: " << mesh.n_faces() << std::endl;
    }
    //remaining faces
    unsigned int target_faces = 400;
    int max_iter = 100;
    originalmesh = mesh;
    mesh.Initialization();

    mesh.request_face_status();
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.detect_feature_edges(75.0);
    mesh.save_feature_lines("original_features.obj");
    unsigned int num = 0;
    while(mesh.current_nfs > target_faces && num<max_iter)
    { 
        mesh.construt_bm(target_faces);
        num++;
    }
    mesh.garbage_collection();
    mesh.save_feature_lines("base_features.obj");
    int mapping_points = 0;
    for (BaseMesh::FaceIter f_iter = mesh.faces_begin(); f_iter != mesh.faces_end(); f_iter++)
    {
        size_t n = mesh.data(*f_iter).vertices.size();
        std::vector<unsigned int>vertices = mesh.data(*f_iter).vertices;
   

        mapping_points += n;
        
    }
    std::cout << "delete points number: "<<mapping_points << "  basemesh vertices:"<<mesh.n_vertices()<<"  original mesh: "<<originalmesh.n_vertices()<<std::endl;

    mesh.save_pts();
    mesh.remesh_(3);
    mesh.MidSubdivision(3);
    OpenMesh::IO::write_mesh(mesh, "../data/result.obj");
	std::cout << "final mesh vertices:" << mesh.n_vertices() << " faces: " << mesh.n_faces() << std::endl;

}
