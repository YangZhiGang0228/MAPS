#ifndef VERTEX_REMOVE_DATA
#define VERTEX_REMOVE_DATA

#include <vector>
#include <Eigen/Core>

struct vertex_remove_data
{
	std::map<unsigned int, unsigned int>g2l;
	std::vector<std::array<double, 2>>uv_coords;                 //Include center vetrex(0,0)
	std::vector<unsigned int>subsetvids;                         //1-ring vertices include center point
	std::vector<std::array<unsigned int, 3>> FUV_pre, FUV_post;  //before vertex remove and after uv faces
};

#endif //VERTEX_REMOVE_DATA
