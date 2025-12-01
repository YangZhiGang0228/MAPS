#include "BaseMesh.h"
#include<iostream>
#include <filesystem>
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace std;
namespace fs = std::filesystem;

int main(int argc,char* argv[])
{
    string mesh_path,output_dir;
	unsigned int target_faces,num_subdivs;
    if (argc == 5)
    {
		mesh_path = argv[1];
		target_faces = stoi(argv[2]);
        num_subdivs = stoi(argv[3]);
		output_dir = argv[4];
    }
    else
    {
		throw invalid_argument("Usage: ./program [mesh_path] [target_face] [num_subdivs] [output_dir]");
    }

    fs::path mesh_fs(mesh_path);
    std::string base_name = mesh_fs.stem().string();

    BaseMesh mesh;
    BaseMesh originalmesh;
    if(!OpenMesh::IO::read_mesh(mesh,mesh_path))
    {
        std::cerr<<"Cannot open file"<<std::endl;
        exit(1);
    }
    else
    {
        std::cout<<"Successfully read mesh"<<std::endl;
    }
    if (mesh.n_faces() >= 1000)
    {
        originalmesh = mesh;
        mesh.Initialization();
        mesh.request_face_status();
        mesh.request_vertex_status();
        mesh.request_edge_status();
        unsigned int max_iter = 100;
        unsigned int iter = 0;
        while (mesh.current_nfs > target_faces && iter < max_iter)
        {
            mesh.construt_bm(target_faces);
            iter++;
        }
        mesh.garbage_collection();

        int mapping_points = 0;
        for (BaseMesh::FaceIter f_iter = mesh.faces_begin(); f_iter != mesh.faces_end(); f_iter++)
        {
            size_t n = mesh.data(*f_iter).vertices.size();
            std::vector<unsigned int>vertices = mesh.data(*f_iter).vertices;


            mapping_points += n;

        }
        std::cout << "delete points number: " << mapping_points << "  basemesh vertices:" << mesh.n_vertices() << "  original mesh: " << originalmesh.n_vertices() << std::endl;

        //mesh.save_pts();
        if (mesh.n_faces() <= target_faces)
        {
            fs::create_directories(output_dir + "/simplified");
            for (unsigned int i = 1; i <= num_subdivs; i++)
            {
                fs::create_directories(output_dir + "/subd" + std::to_string(i));
                mesh.remesh_(i, output_dir, base_name);
            }

            //mesh.MidSubdivision(4);
            string simplified_path = output_dir + "/simplified/" + base_name + "_001.obj";
            OpenMesh::IO::write_mesh(mesh, simplified_path);
            std::cout << "Finished remeshing!" << std::endl;
        }
        else
        {
            fs::create_directories(output_dir + "/ignore_mesh");
            string output_name = output_dir + "/ignore_mesh/" + base_name + ".obj";
            std::cout << "Can't simplified mesh to target faces  !" << std::endl;
            OpenMesh::IO::write_mesh(originalmesh, output_name);
        }
    }
    else
    {
        fs::create_directories(output_dir + "/ignore_mesh");
        string output_name = output_dir + "/ignore_mesh/" + base_name + ".obj";
        std::cout << "The input mesh is too small to be processed!" << std::endl;
        OpenMesh::IO::write_mesh(mesh, output_name);
    }
}
