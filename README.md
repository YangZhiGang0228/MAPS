# MAPS
This project is the implementation of the MAPS: Multiresolution Adaptive Parameterization of Surfaces.

## 1 Mesh simplification and self-parameterization
### 1.1 Vertex removal
The order of vertex removal is determined based on the curvature magnitude of vertices, giving priority to removing vertices with smaller curvature.We can decide how many faces or vertices to retain in the final triangle mesh.

Some results as follow：
<img width="794" height="365" alt="image" src="https://github.com/user-attachments/assets/3b631c8a-775b-40c9-8a94-fcc8bad33506" />


