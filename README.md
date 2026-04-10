# MAPS
This project is an improved implementation of MAPS: Multiresolution Adaptive Parameterization of Surfaces.

## 1 Mesh simplification and self-parameterization
### 1.1 Simplification
The order of vertex removal is determined based on the curvature magnitude of vertices, giving priority to removing vertices with smaller curvature.We can decide how many faces or vertices to retain in the final triangle mesh.

Some results as follow：

<img width="758" height="455" alt="image" src="https://github.com/user-attachments/assets/651d3bf6-245c-4bde-a57b-39e1db0584a7" />

<img width="879" height="500" alt="image" src="https://github.com/user-attachments/assets/63e1cc36-f8d3-4a37-8de0-e31e96dddb27" />

### 1.2 self-parameterization
During the vertex removal process, through surface parameterization, we find which specific triangle the removed vertex lies in and compute its barycentric coordinates.Visualize result：

<img width="885" height="360" alt="image" src="https://github.com/user-attachments/assets/df4265c1-fba8-4b67-89fd-1801b54affa2" />

## 2 Mesh subdividion and remesh
### 2.1 Subdivision
We subdivide the simplified mesh using midpoint subdivision.

<img width="395" height="473" alt="image" src="https://github.com/user-attachments/assets/2c40523a-fb1a-4cfe-8448-4e2e782e91e5" />


### 2.2 Remesh

Because the positional relationships between the original vertices and the hierarchical levels are recorded during the simplification process, the corresponding location of a subdivided vertex in the original mesh can be found, thereby enabling remeshing.

<img width="901" height="351" alt="image" src="https://github.com/user-attachments/assets/d8df71e0-99fe-4c89-8aaa-a727ba7684cc" />


## 3 Feature Line

<img width="1075" height="768" alt="image" src="https://github.com/user-attachments/assets/5fe08f99-59b1-412b-8026-1cb77f708fa3" />




















