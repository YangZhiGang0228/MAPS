#ifndef MAPS_H
#define MAPS_H
#include <OpenMesh/Core/IO/MeshIO.hh>
#include"OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

template <typename Trait = OpenMesh::DefaultTraits>
class maps:public OpenMesh::TriMesh_ArrayKernelT<Trait>
{
public:
    maps(){};
    ~maps(){};
};
#endif // MAPS_H
