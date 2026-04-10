#include"BaseMesh.h"
#include <algorithm>

namespace MyCGALStuff {
    template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
    class My_vertex_base
        : public  Vb
    {
        typedef Vb                              Base;
    public:
        typedef typename Vb::Vertex_handle      Vertex_handle;
        typedef typename Vb::Face_handle        Face_handle;
        typedef typename Vb::Point              Point;
        template < typename TDS2 >
        struct Rebind_TDS {
            typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
            typedef My_vertex_base<Gt, Vb2>                           Other;
        };

    public:
        My_vertex_base() : Base(), id(-1) {}
        My_vertex_base(const Point& p) : Base(p), id(-1) {}
        My_vertex_base(const Point& p, Face_handle f) : Base(f, p), id(-1) {}
        My_vertex_base(Face_handle f) : Base(f), id(-1) {}

        int id;
    };

    struct FaceInfo2
    {
        FaceInfo2() {}
        int nesting_level;
        bool in_domain() {
            //return nesting_level % 2 == 1;
            return nesting_level > 0;
        }
    };
    typedef CGAL::Simple_cartesian<double> K;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
    typedef CGAL::Exact_predicates_tag                                           Itag;
    typedef My_vertex_base<K>                                                    Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>              Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>                  Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                         Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>             CDT;
    typedef CDT::Point              Point;
    typedef CDT::Vertex_handle      Vhandle;
    typedef CDT::Face_handle        Face_handle;

    // the following two functions are drawn from CGAL website:
    // https://doc.cgal.org/latest/Triangulation_2/index.html#Section_2D_Triangulations_Constrained
    void mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border) {
        if (start->info().nesting_level != -1) {
            return;
        }
        std::list<Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; i++) {
                    CDT::Edge e(fh, i);
                    Face_handle n = fh->neighbor(i);
                    if (n->info().nesting_level == -1) {
                        if (ct.is_constrained(e)) border.push_back(e);
                        else queue.push_back(n);
                    }
                }
            }
        }
    }

    void mark_domains(CDT& cdt) {
        for (CDT::Face_handle f : cdt.all_face_handles()) {
            f->info().nesting_level = -1;
        }
        std::list<CDT::Edge> border;
        mark_domains(cdt, cdt.infinite_face(), 0, border);
        while (!border.empty()) {
            CDT::Edge e = border.front();
            border.pop_front();
            Face_handle n = e.first->neighbor(e.second);
            if (n->info().nesting_level == -1) {
                mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
        }
    }
}

void BaseMesh::Initialization()
{
    originFaces = std::vector<FaceHandle>(faces_begin(), faces_end());
    for (const auto& face : originFaces)
    {
		std::array<unsigned int, 3> f_vids;
		unsigned int vidx = 0;
        originFaceVertices[face] = std::vector<VertexHandle>(fv_begin(face), fv_end(face));
        for(FaceVertexIter fv = fv_iter(face); fv.is_valid(); fv++)
        {
			f_vids[vidx] = fv->idx();
			vidx++;
		}
		all_orifaces.push_back(f_vids);
    }
    for (VertexIter vh = vertices_begin(); vh != vertices_end(); vh++)
    {
        for (VertexFaceIter vf = vf_iter(*vh); vf.is_valid(); vf++)
        {
            originVertexFaces[*vh].push_back(*vf);
        }
    }
    unsigned int i = 0;
    nv_ = n_vertices();
	nf_ = n_faces();
    uvs_.resize(nv_);
    barycoordinates.resize(nv_);
	is_deleted.resize(nv_, false);
	current_nfs = nf_;
    for (auto& v : vertices())
    {
        original_vhs.push_back(v);
        ori_vids[v] = v.idx();
        data(v).ori_id = v.idx();
        positions[i] = point(v);
        uvs_[i] = { point(v)[0],point(v)[1],point(v)[2] };
        i++;
    }
   
}

void BaseMesh::detect_feature_edges(double threshold_degree)
{
    if(!has_face_normals())
    {
        request_face_normals();
	}
	update_face_normals();

	double threshold_cos = std::cos(threshold_degree * M_PI / 180.0);

    for (EdgeIter e_it = edges_begin(); e_it != edges_end(); e_it++)
    {
        //boundary is a feature edge
        if (is_boundary(*e_it))
        {
            data(*e_it).is_feature = true;
            continue;
        }
		HalfedgeHandle heh = halfedge_handle(*e_it, 0);
        FaceHandle fh1 = face_handle(heh);
        FaceHandle fh2 = face_handle(opposite_halfedge_handle(heh));
        Normal normal1 = normal(fh1);
        Normal normal2 = normal(fh2);
        double cos_angle = std::max(-1.0, std::min(1.0, static_cast<double>(normal1.normalized().dot(normal2.normalized()))));
        if (cos_angle < threshold_cos)
        {
            data(*e_it).is_feature = true;
		}
        else {
            data(*e_it).is_feature = false;
        }
    }
    for (VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
    {
        int incident_features = 0;
        for (VertexEdgeIter ve_it = ve_iter(*v_it); ve_it.is_valid(); ++ve_it)
        {
            if (data(*ve_it).is_feature) {
                incident_features++;
            }
        }

        data(*v_it).feature_edge_count = incident_features;

        if (incident_features == 1) {
            data(*v_it).is_dart = true;
        }
        else if (incident_features > 2) {
            data(*v_it).is_corner = true;
        }
    }
}

bool BaseMesh::check_uv_flip(std::array<Vec2d, 3>uv_face)
{
    Vec2d v1 = uv_face[1] - uv_face[0];
    Vec2d v2 = uv_face[2] - uv_face[0];
    double signarea = v1[0] * v2[1] - v1[1] * v2[0];
    if (signarea < 1e-10 || isnan(signarea))
        return true;
    return false;
}

bool BaseMesh::check_all_orifaces_uv_flip(std::vector<FaceHandle>ring_faces, std::vector<std::pair<VertexHandle, Vec2d>> coordinates)
{
    unsigned int n_uvs = coordinates.size();
    std::vector<bool>vmask(nv_, false);
    std::vector<unsigned int>all_v_ids;
    for (unsigned int i = 0; i < ring_faces.size(); i++)
    {
        std::vector<unsigned int>v_ids = data(ring_faces[i]).vertices;
        for (unsigned int j = 0; j < v_ids.size(); j++)
        {
            vmask[v_ids[j]] = true;
            all_v_ids.push_back(v_ids[j]);
        }
    }
    std::vector<FaceHandle>all_in;
    for (unsigned int i = 0; i < all_v_ids.size(); i++)
    {
        VertexHandle vh = original_vhs[i];
        std::vector<FaceHandle>r_faces = originVertexFaces[vh];
        for (unsigned int j = 0; j < r_faces.size(); j++)
        {
            std::vector<VertexHandle>r_vhs = originFaceVertices[r_faces[j]];
            assert(r_vhs.size() == 3);
            bool flag = false;
            for (unsigned int k = 0; k < 3; k++)
            {
                if (!vmask[ori_vids[r_vhs[k]]])
                {
                    flag = true;
                }
            }
            if (!flag)
            {
                all_in.push_back(r_faces[j]);
            }
        }
    }
    for (unsigned int i = 0; i < all_in.size(); i++)
    {
        std::vector<VertexHandle>f_vs = originFaceVertices[all_in[i]];
        std::array<Vec2d, 3>face;
        assert(f_vs.size() == 3);
        for (int j = 0; j < 3; j++)
        {
            VertexHandle vh = f_vs[j];
            Vec2d uv = Vec2d(0, 0);
            for (unsigned int n = 0; n < 3; n++)
            {
                for (unsigned int k = 0; k < n_uvs; k++)
                {
                    if (data(vh).barycoor.has_value())
                    {
                        VertexHandle v_h = original_vhs[data(vh).barycoor.value()[n].first];
                        if (coordinates[k].first == v_h)
                        {
                            uv += coordinates[k].second * data(vh).barycoor.value()[n].second;
                        }
                    }
                }
            }
            face[j] = uv;
        }
        if (!check_uv_flip(face))
        {
            std::cout << "original face flip,cancel this step !" << std::endl;
            return false;
        }
    }
    return true;
}

double BaseMesh::cal_area(const Vec2d p1,const Vec2d p2,const Point P0,const Point P1,const Point P2)
{
    Point P03d = P1-P0;
    Point P13d = P2-P0;
    double z = P03d.cross(P13d)[2];
    double area1 = p1[0]*p2[1]-p2[0]*p1[1];
    double sign = z*area1;

    return sign;
}

double BaseMesh::cal_Angles(Point& v0, Point& v1, Point& v2)
{
    double eps = 1e-8;

    double a = (v1 - v0).norm();
    double b = (v2 - v0).norm();
    double c = (v2 - v1).norm();
    double theta = (a * a + b * b - c * c) / (2 * a * b + eps);
    if (theta > 1)
    {
        theta = 1;
    }
    if (theta < -1)
    {
        theta = -1;
    }


    return std::acos(theta);
}

//bool BaseMesh::In_2Dtriangle(const Vec2d& point2d,const std::array<Vec2d,3>&f2d)
//{
//    /*double esp = 1e-10;
//    Point point0 (point2d[0],point2d[1],0);
//    Point point1 (f2d[0][0],f2d[0][1],0);
//    Point point2 (f2d[1][0],f2d[1][1],0);
//    Point point3 (f2d[2][0],f2d[2][1],0);
//
//    double z1 = (point0-point1).cross(point2-point1)[2];
//    double z2 = (point0-point2).cross(point3-point2)[2];
//    double z3 = (point0-point3).cross(point1-point3)[2];
//
//    if((z1>=-esp && z2 >= -esp && z3 >= -esp) || (z1<=esp && z2 <= esp && z3 <= esp))
//    {
//        return true;
//    }
//    return false;*/
//
//	Point_2 p(point2d[0], point2d[1]);
//	Point_2 p1(f2d[0][0], f2d[0][1]);
//	Point_2 p2(f2d[1][0], f2d[1][1]);
//	Point_2 p3(f2d[2][0], f2d[2][1]);
//
//	CGAL::Orientation o1 = CGAL::orientation(p1, p2, p);
//	CGAL::Orientation o2 = CGAL::orientation(p2, p3, p);
//	CGAL::Orientation o3 = CGAL::orientation(p3, p1, p);
//
//    return ((o1 == o2 && o2 == o3) ||
//        (o1 == CGAL::COLLINEAR && o2 == o3) ||
//        (o2 == CGAL::COLLINEAR && o1 == o3) ||
//        (o3 == CGAL::COLLINEAR && o1 == o2));
//}

std::tuple<double,double,double> BaseMesh::cal_bar_2dcoord(const Vec2d& point,const std::array<Vec2d,3>&f2d)
{
    /*Point point0 = Point(point[0],point[1],0);
    Point point1 = Point(f2d[0][0],f2d[0][1],0);
    Point point2 = Point(f2d[1][0],f2d[1][1],0);
    Point point3 = Point(f2d[2][0],f2d[2][1],0);
    double eps = 1e-10;

    double area = (point2-point1).cross(point3-point1)[2];
    double alpha = (point0-point2).cross(point0-point3)[2]/area;
    double beta = (point0-point1).cross(point0-point3)[2]/area;
    double gamma = 1.0f-alpha-beta;


    return {alpha,beta,gamma};*/

    /*Point_2 p(point[0], point[1]);
    Point_2 p1(f2d[0][0], f2d[0][1]);
    Point_2 p2(f2d[1][0], f2d[1][1]);
    Point_2 p3(f2d[2][0], f2d[2][1]);
    std::array<typename Point_2::R::FT, 3> coords;
    CGAL::Barycentric_coordinates::triangle_coordinates_2(
        p1, p2, p3, p, coords.begin()
    );

    return { CGAL::to_double(coords[0]),
            CGAL::to_double(coords[1]),
            CGAL::to_double(coords[2]) };*/
    const Vec2d& a = f2d[0];
    const Vec2d& b = f2d[1];
    const Vec2d& c = f2d[2];

    auto signed_area = [](const Vec2d& u, const Vec2d& v, const Vec2d& w) -> double {
        return 0.5 * ((v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0]));
        };

    double area_total = signed_area(a, b, c);

    if (std::abs(area_total) < 1e-20)
        return { 0.0, 0.0, 0.0 };

    double area_pbc = signed_area(point, b, c);
    double area_pca = signed_area(point, c, a);
    double area_pab = signed_area(point, a, b);

    double alpha = area_pbc / area_total;
    double beta = area_pca / area_total;
    double gamma = area_pab / area_total;

    return { alpha, beta, gamma };

}

void BaseMesh::compute_curvature(const VertexHandle& v_h)
{
    double angleSum = 0.0;
    for (auto vf_it = vf_iter(v_h); vf_it.is_valid(); ++vf_it)
    {
        HalfedgeHandle heh = halfedge_handle(*vf_it);
        VertexHandle v0 = to_vertex_handle(heh);
        VertexHandle v1 = from_vertex_handle(heh);
        VertexHandle v2 = to_vertex_handle(next_halfedge_handle(heh));
        Point p0 = point(v0);
        Point p1 = point(v1);
        Point p2 = point(v2);
        Vec3d u(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
        Vec3d v(p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]);
        double cosTheta = std::max(-1.0, std::min(1.0, u.normalized().dot(v.normalized())));
        double theta = std::acos(cosTheta);
        angleSum += theta;
    }
    double k = 2 * M_PI - angleSum;
    data(v_h).curvature = k;
}

void BaseMesh::compute_ringArea(const VertexHandle& v_h)
{
    double ringArea = 0;
    for (auto faceIter = vf_begin(v_h); faceIter != vf_end(v_h); faceIter++)
    {
        ringArea += calc_face_area(*faceIter);
    }
    data(v_h).ringArea = ringArea;
}


void BaseMesh::initial_weights(double lamda)
{
    while(!removal_P.empty())
    {
        removal_P.pop();
    }
    maxCurvature = (std::numeric_limits<double>::min());
    maxringArea = (std::numeric_limits<double>::min());

    for(auto i = vertices_begin();i != vertices_end();i++)
    {
		compute_curvature(*i);
		compute_ringArea(*i);
        double curvature = data(*i).curvature;
		double ringArea = data(*i).ringArea;
        if(curvature >= maxCurvature)
        {
            maxCurvature = curvature;
		}
        if (ringArea >= maxringArea)
        {
            maxringArea = ringArea;
        }
    }
    for (auto i = vertices_begin(); i != vertices_end(); i++)
    {
		data(*i).weight = lamda * (data(*i).curvature / maxCurvature) + (1 - lamda) * (data(*i).ringArea / maxringArea);
		removal_P.push(*i);
    }
}

void BaseMesh::map2plane(VertexHandle v_h,std::vector<std::pair<VertexHandle,Vec2d>>&coordinates)
{
    coordinates.clear();
    auto vertex3D = point(v_h);
    double K_i = 0.0;
    std::vector<VertexHandle>one_ring_vs(vv_begin(v_h), vv_end(v_h));
    unsigned int n = one_ring_vs.size();
    std::vector<double>angles(n);
    for (unsigned int i = 0; i < n;i++) {
        Point v1 = point(one_ring_vs[(i + n - 1) % n]);
        Point v2 = point(one_ring_vs[i]);
        angles[i] = cal_Angles(vertex3D, v1, v2);
        K_i += angles[i];
    }
    double a = 2 * M_PI / K_i;
    double k_i = 0.0;
    for (unsigned int i = 0; i < n;i++) {
     
        double r_k;
        VertexHandle rv_h = one_ring_vs[i];
        r_k = (point(rv_h) - vertex3D).norm();
        k_i += angles[i];

        Vec2d point2D(std::pow(r_k, a) * std::cos(k_i * a), std::pow(r_k, a) * std::sin(k_i * a));
        coordinates.emplace_back(rv_h, point2D);
    }
}

void BaseMesh::CGAL_CDT(const std::vector<std::pair<VertexHandle, Vec2d>>& coordinates, std::vector<std::array<VertexHandle, 3>>& faces,
    std::pair<VertexHandle, VertexHandle> feature_constraint)
{
    MyCGALStuff::CDT cdt;
    std::vector<MyCGALStuff::Vhandle>cdt_vhandles;
    std::map<VertexHandle, MyCGALStuff::Vhandle> vh_to_cdt_vh;
    for (unsigned int i = 0; i < coordinates.size(); i++) {
        double x = coordinates[i].second[0];
        double y = coordinates[i].second[1];

        MyCGALStuff::Vhandle v = cdt.insert(MyCGALStuff::Point(x, y));
        if (v->id != -1)
        {
            printf("CDT insert point error\n");
            return;
        }
        else
        {
            v->id = i;
            cdt_vhandles.push_back(v);
            vh_to_cdt_vh[coordinates[i].first] = v;
        }
    }
    for (unsigned int i = 0; i < cdt_vhandles.size(); i++)
    {
        cdt.insert_constraint(cdt_vhandles[i], cdt_vhandles[(i + 1) % cdt_vhandles.size()]);
    }
	//insert feature constraint
    if (feature_constraint.first.is_valid() && feature_constraint.second.is_valid()) {
        if (vh_to_cdt_vh.count(feature_constraint.first) && vh_to_cdt_vh.count(feature_constraint.second)) {
            cdt.insert_constraint(vh_to_cdt_vh[feature_constraint.first], vh_to_cdt_vh[feature_constraint.second]);
        }
    }
    assert(cdt.is_valid());
    MyCGALStuff::mark_domains(cdt);
    std::array<VertexHandle, 3>cdt_face;
    for (const MyCGALStuff::CDT::Face_handle& f : cdt.finite_face_handles())
    {
        if (f->info().in_domain())
        {
            cdt_face[0]=coordinates[f->vertex(0)->id].first;
            cdt_face[1] = coordinates[f->vertex(2)->id].first;
            cdt_face[2] = coordinates[f->vertex(1)->id].first;
            faces.push_back(cdt_face);
        }
    }
}

void BaseMesh::Retriangle(const VertexHandle &deletevertex,unsigned int remaining_nfs)
{
    
    std::vector<FaceHandle> ringFaces(vf_begin(deletevertex), vf_end(deletevertex));
    std::vector < std::vector<VertexHandle>>one_ring_vhs;
    for (unsigned int i = 0; i < ringFaces.size(); i++)
    {
        FaceHandle f_h = ringFaces[i];
        std::vector<VertexHandle>v_hs(fv_begin(f_h), fv_end(f_h));
        one_ring_vhs.push_back(v_hs);
    }
    std::vector<std::pair<VertexHandle,Vec2d>> coordinates;
    map2plane(deletevertex, coordinates);

    bool flag = false;
    for (unsigned int i = 0; i < one_ring_vhs.size(); i++)
    {
        std::array<Vec2d, 3>UV_face_pre;
        for (unsigned int j = 0; j < 3; j++)
        {
            for (unsigned int k = 0; k < coordinates.size(); k++)
            {
                if (coordinates[k].first == one_ring_vhs[i][j])
                {
                    UV_face_pre[j] = coordinates[k].second;
                }
                else if (one_ring_vhs[i][j] == deletevertex)
                {
                    UV_face_pre[j] = Vec2d(0.0, 0.0);
                }
            }
        }
        if (!check_uv_flip(UV_face_pre))
        {
            flag = true;
            break;
        }
    }
    //check original face flip
    if (!check_all_orifaces_uv_flip(ringFaces, coordinates))
    {
        flag = true;
    }
    //insert new triangulation, if this vertex is on feature edge and has two feature neighbors, add feature constraint
    std::pair<VertexHandle, VertexHandle> feature_constraint = { VertexHandle(-1), VertexHandle(-1) };
    if (data(deletevertex).feature_edge_count == 2 && !data(deletevertex).is_dart && !data(deletevertex).is_corner) {
        std::vector<VertexHandle> feature_neighbors;
        for (VertexEdgeIter ve_it = ve_iter(deletevertex); ve_it.is_valid(); ++ve_it) {
            if (data(*ve_it).is_feature) {
                HalfedgeHandle he = halfedge_handle(*ve_it, 0);
                if (to_vertex_handle(he) == deletevertex) {

                    feature_neighbors.push_back(from_vertex_handle(he));
                }
                else {
                    feature_neighbors.push_back(to_vertex_handle(he));
                }
            }
        }
        if (feature_neighbors.size() == 2) {
            feature_constraint.first = feature_neighbors[0];
            feature_constraint.second = feature_neighbors[1];
        }
    }
    //delete one ring of this vertex
    for (unsigned int i = 0; i < ringFaces.size(); i++)
    {
        delete_face(ringFaces[i],false);
    }
    std::vector<std::array<VertexHandle, 3>> newFaces;
    CGAL_CDT(coordinates, newFaces, feature_constraint);
    //add new triangulation，check the new face isn't complex edge and uv flip
    std::vector<VertexHandle> newfhandles;
    std::vector<FaceHandle> testfhandles;
    unsigned int new_face_size = newFaces.size();
    for (unsigned int i = 0; i < new_face_size; i++)
    {
        std::array<Vec2d,3>UV_face_post;
        newfhandles.clear();
        for (unsigned int j = 0; j < 3; j++)
        {
            newfhandles.push_back(newFaces[i][j]);
            for (unsigned int k = 0; k < coordinates.size(); k++)
            {
                if (coordinates[k].first == newFaces[i][j])
                {
                    UV_face_post[j]=(coordinates[k].second);
                }
            }
        }
        if (!check_uv_flip(UV_face_post))
        {
            flag = true;
        }
        FaceHandle f = add_face(newfhandles);
        if (!f.is_valid())
        {
            flag = true;
        }
        else
        {
            testfhandles.push_back(f);
            if (feature_constraint.first.is_valid() && feature_constraint.second.is_valid()) {
                for (FaceEdgeIter fe_it = fe_iter(f); fe_it.is_valid(); ++fe_it) {
                    VertexHandle v0 = to_vertex_handle(halfedge_handle(*fe_it, 0));
                    VertexHandle v1 = from_vertex_handle(halfedge_handle(*fe_it, 0));

					//mark feature edge in the new triangulation
                    if ((v0 == feature_constraint.first && v1 == feature_constraint.second) ||
                        (v1 == feature_constraint.first && v0 == feature_constraint.second)) {
                        data(*fe_it).is_feature = true;
                    }
                }
            }
        }
    }
    bool is_find = false;

    unsigned int ncoords = coordinates.size();
    if (!flag)
    {
        std::vector<unsigned int>subvids;
        std::vector<std::array<double, 2>>vuv_coords;
        std::vector<std::array<unsigned int, 3>> fuv_pre, fuv_post;
        vertex_remove_data vrd;
        
		subvids.push_back(deletevertex.idx());
		vuv_coords.push_back({ 0.0,0.0 });
		vrd.g2l[deletevertex.idx()] = 0;

		is_deleted[deletevertex.idx()] = true;
        delete_vertex(deletevertex, true);
        data(deletevertex).isdeleted = true;
		unsigned int delta = ringFaces.size() - new_face_size;
		current_nfs -= delta;
        //collect vertex remove information
        for (unsigned int i = 0; i < coordinates.size(); i++)
        {
			subvids.push_back(coordinates[i].first.idx());
			vuv_coords.push_back({ coordinates[i].second[0],coordinates[i].second[1] });
			vrd.g2l[coordinates[i].first.idx()] = i + 1;
        }
        for (unsigned int i = 0; i < one_ring_vhs.size(); i++)
        {
            std::array<unsigned int, 3>fuv;
            std::vector <VertexHandle>fvhs = one_ring_vhs[i];
            fuv[0] = fvhs[0].idx();
            fuv[1] = fvhs[1].idx();
            fuv[2] = fvhs[2].idx();
            fuv_pre.push_back(fuv);
        }
        for (unsigned int i = 0; i < new_face_size; i++)
        {
            std::array<unsigned int, 3>fuv;
            std::array<VertexHandle, 3>fvs = newFaces[i];
            fuv[0] = fvs[0].idx();
            fuv[1] = fvs[1].idx();
            fuv[2] = fvs[2].idx();
            fuv_post.push_back(fuv);
        }
        vrd.subsetvids = subvids;
        vrd.uv_coords = vuv_coords;
        vrd.FUV_pre = fuv_pre;
        vrd.FUV_post = fuv_post;
        vremoveIM.push_back(vrd);
        for (unsigned int i = 0; i < fuv_post.size(); i++)
        {
            Face2IM[fuv_post[i]] = vremoveIM.size() - 1;
			test_faces.push_back(fuv_post[i]);
        }

        //update barycenter coordinates
        std::array<Vec2d, 3> triangle2D;
        double min_distance = 1.0;
        std::array<VertexHandle, 3>cloest_f;
        std::array<double, 3> closest_bary;
        unsigned int cloest_i = 0;
        for (unsigned int i = 0; i < new_face_size; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                for (unsigned int k = 0; k < ncoords; k++)
                {
                    if (coordinates[k].first == newFaces[i][j])
                    {
                        triangle2D[j] = coordinates[k].second;
                    }
                }
            }
            auto [alpha, beta, gamma] = cal_bar_2dcoord(Vec2d(0, 0), triangle2D);
            double min_coord = std::min(alpha, std::min(beta, gamma));
            double distance = -min_coord;
            if (distance < min_distance)
            {
                min_distance = distance;
                cloest_f = newFaces[i];
                closest_bary = { alpha, beta, gamma };
                cloest_i = i;
            }
        }
        double alpha = std::max(0.0, closest_bary[0]);
        double beta = std::max(0.0, closest_bary[1]);
        double gamma = std::max(0.0, closest_bary[2]);
        double sum = alpha + beta + gamma;
        Barycoor barycoord;
        barycoord[0] = std::make_pair(ori_vids[cloest_f[0]], alpha/sum);
        barycoord[1] = std::make_pair(ori_vids[cloest_f[1]], beta/sum);
        barycoord[2] = std::make_pair(ori_vids[cloest_f[2]], gamma/sum);
        data(deletevertex).barycoor = barycoord;
        barycoordinates[ori_vids[deletevertex]] = barycoord;
        data(testfhandles[cloest_i]).vertices.push_back(ori_vids[deletevertex]);
            
        //recompute 2d coordinates in the new 1-ring
        std::vector<unsigned int>local_vertices;
        for (unsigned int j = 0; j < ringFaces.size(); j++)
        {
            local_vertices = data(ringFaces[j]).vertices;

            for (unsigned int i = 0; i < local_vertices.size(); i++)
            {
                VertexHandle &inner_vertex = original_vhs[local_vertices[i]];
                Vec2d uv = Vec2d(0,0);
                for (unsigned int n = 0; n < 3; n++)
                {
                    for (unsigned int k = 0; k < ncoords; k++)
                    {
                        if (data(inner_vertex).barycoor.has_value())
                        {
                            VertexHandle v_h = original_vhs[data(inner_vertex).barycoor.value()[n].first];
                            if (coordinates[k].first == v_h)
                            {
                                uv += coordinates[k].second * data(inner_vertex).barycoor.value()[n].second;
                            }
                        }
                    }
                }
                double min_distance1 = 1.0;
                std::array<VertexHandle, 3>cloest_f1;
                std::array<double, 3> closest_bary1;
                unsigned int cloest_ii = 0;
                for (unsigned int s = 0; s < new_face_size; s++)
                {
                    std::array<Vec2d, 3>new_f;
                    std::array<VertexHandle, 3>new_fv = newFaces[s];
                    for (unsigned int t = 0; t < 3; t++)
                    {
                        for (unsigned int n = 0; n < ncoords; n++)
                        {
                            if (coordinates[n].first == new_fv[t])
                            {
                                new_f[t] = coordinates[n].second;
                            }
                        }
                    }
                    auto [alpha, beta, gamma] = cal_bar_2dcoord(uv, new_f);
                    double min_coord1 = std::min(alpha, std::min(beta, gamma));
                    double distance1 = -min_coord1;
                    if (distance1 < min_distance1)
                    {
                        min_distance1 = distance1;
                        cloest_f1 = newFaces[s];
                        closest_bary1 = { alpha, beta, gamma };
                        cloest_ii = s;
                    }
                }
				double alpha1 = std::max(0.0, closest_bary1[0]);
				double beta1 = std::max(0.0, closest_bary1[1]);
				double gamma1 = std::max(0.0, closest_bary1[2]);
				double sum1 = alpha1 + beta1 + gamma1;
                barycoord[0] = std::make_pair(ori_vids[cloest_f1[0]], alpha1/sum1);
                barycoord[1] = std::make_pair(ori_vids[cloest_f1[1]], beta1/sum1);
                barycoord[2] = std::make_pair(ori_vids[cloest_f1[2]], gamma1/sum1);
                data(inner_vertex).barycoor = barycoord;
                barycoordinates[ori_vids[inner_vertex]] = barycoord;
                data(testfhandles[cloest_ii]).vertices.push_back(ori_vids[inner_vertex]);
            }
        }

    }
    else
    {
        for (unsigned int j = 0; j < testfhandles.size(); j++)
        {
            delete_face(testfhandles[j], false);
        }
        data(deletevertex).isdeleted = true;
        for (unsigned int j = 0; j < one_ring_vhs.size(); j++)
        {
            newfhandles.clear();
            for (unsigned int k = 0; k < 3; k++)
            {
                newfhandles.push_back(one_ring_vhs[j][k]);
            }

            FaceHandle f = add_face(newfhandles);
            std::vector<unsigned int>vertices;
            vertices = data(ringFaces[j]).vertices;
            for (unsigned int k = 0; k < vertices.size(); k++)
            {
                data(f).vertices.push_back(vertices[k]);
            }
        }
    }

}



void BaseMesh::construt_bm(unsigned int remaining_nfs)
{
    for(auto vertexiter = vertices_begin();vertexiter != vertices_end();vertexiter++)
    {
        if (data(*vertexiter).is_dart || data(*vertexiter).is_corner) {
            data(*vertexiter).canberemoved = false;
        }
        else {
            data(*vertexiter).canberemoved = true;
        }
    }

    initial_weights(0.5);
    while(!removal_P.empty() && current_nfs>remaining_nfs)
    {
		VertexHandle p_temp;
        auto vertexhandle = removal_P.top();
        removal_P.pop();
		/*int nth_vertex = 1 + (rand() % 100);
		std::list<VertexHandle> temp_vertices;
        if(nth_vertex > removal_P.size())
        {
            nth_vertex = removal_P.size()-1;
		}
        
        for (int i = 0; i < nth_vertex; i++)
        {
			p_temp = removal_P.top();
			temp_vertices.push_back(p_temp);
			removal_P.pop();
        }

        if (removal_P.empty()) {
            break;  
        }
		auto vertexhandle = removal_P.top();
        removal_P.pop();

        while (temp_vertices.size() > 0)
        {
			removal_P.push(temp_vertices.front());
			temp_vertices.pop_front();
        }*/

        if(!data(vertexhandle).canberemoved or data(vertexhandle).isdeleted) continue;

        for(auto ringvertex = vv_begin(vertexhandle);ringvertex != vv_end(vertexhandle);ringvertex++)
        {
            data(*ringvertex).canberemoved = false;
        }
        Retriangle(vertexhandle,remaining_nfs);

    }


}

void BaseMesh::save_pts()
{
    unsigned int n = 0;
    unsigned int n1 = 0;
    for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); f_iter++)
    {
        std::vector<unsigned int> vertices = data(*f_iter).vertices;
        for (unsigned int i = 0; i < vertices.size(); i++)
        {
            unsigned int v_id = vertices[i];
            Barycoor v_idbarycoord = barycoordinates[v_id];
            Point v_idpoint(0, 0, 0);
            for (unsigned int j = 0; j < 3; j++)
            {
                v_idpoint += positions[v_idbarycoord[j].first] * v_idbarycoord[j].second;
            }
            uvs_[v_id] = { v_idpoint[0],v_idpoint[1],v_idpoint[2] };
        }
    }
    std::ofstream out("../data/mapping_pts.obj");
    for (unsigned int i = 0; i < nv_; i++) {
        out << "v " << uvs_[i](0) << " " << uvs_[i](1) << " " << uvs_[i](2) << std::endl;
    }
    out.close();
}

void BaseMesh::save_feature_lines(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) return;

    // 使用 ori_id -> OBJ局部索引 的映射，完美规避 garbage_collection 带来的句柄变化
    std::map<unsigned int, unsigned int> vertex_map;
    unsigned int current_v_idx = 1;

    std::vector<std::pair<unsigned int, unsigned int>> lines_to_write;

    // 遍历所有边
    for (EdgeIter e = edges_begin(); e != edges_end(); ++e) {
        // 【关键】：必须跳过已删除的边，并检查是否为特征边
        if (!status(*e).deleted() && data(*e).is_feature) {
            HalfedgeHandle he = halfedge_handle(*e, 0);
            VertexHandle v1 = from_vertex_handle(he);
            VertexHandle v2 = to_vertex_handle(he);

            // 【核心修复】：使用持久化属性 ori_id，而不是可能会改变的 v.idx()
            unsigned int id1 = data(v1).ori_id;
            unsigned int id2 = data(v2).ori_id;

            if (vertex_map.find(id1) == vertex_map.end()) {
                vertex_map[id1] = current_v_idx++;
                out << "v " << point(v1)[0] << " " << point(v1)[1] << " " << point(v1)[2] << "\n";
            }
            if (vertex_map.find(id2) == vertex_map.end()) {
                vertex_map[id2] = current_v_idx++;
                out << "v " << point(v2)[0] << " " << point(v2)[1] << " " << point(v2)[2] << "\n";
            }

            // 记录拓扑关系
            lines_to_write.push_back({ vertex_map[id1], vertex_map[id2] });
        }
    }

    // 写入线段
    for (const auto& line : lines_to_write) {
        out << "l " << line.first << " " << line.second << "\n";
    }

    out.close();
    std::cout << "Feature lines saved to " << filename
        << " (Vertices: " << (current_v_idx - 1)
        << ", Lines: " << lines_to_write.size() << ")" << std::endl;
}

void BaseMesh::MidSubdivision(unsigned int level)
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<unsigned int, 3>> faces;

    for (VertexIter v = vertices_begin(); v != vertices_end(); ++v) {
        vertices.emplace_back(
            point(*v)[0],
            point(*v)[1],
            point(*v)[2]
        );
    }

    for (FaceIter f = faces_begin(); f != faces_end(); ++f) {
        std::array<unsigned int, 3> fv;
        int k = 0;
        for (FaceVertexIter fv_it = fv_iter(*f); fv_it.is_valid(); ++fv_it)
            fv[k++] = fv_it->idx();
        faces.push_back(fv);
    }

    for (unsigned int lv = 0; lv < level; ++lv) {

        unsigned int V = vertices.size();
        unsigned int F = faces.size();

        std::map<EdgeKey, unsigned int> edge2mid;
        std::vector<Eigen::Vector3d> new_vertices = vertices;

        for (const auto& f : faces) {
            EdgeKey e01(f[0], f[1]);
            EdgeKey e12(f[1], f[2]);
            EdgeKey e20(f[2], f[0]);

            edge2mid[e01] = 0;
            edge2mid[e12] = 0;
            edge2mid[e20] = 0;
        }

        for (auto& it : edge2mid) {
            unsigned int v0 = it.first.v0;
            unsigned int v1 = it.first.v1;
            unsigned int mid_id = new_vertices.size();

            it.second = mid_id;
            new_vertices.push_back(
                0.5 * (vertices[v0] + vertices[v1])
            );
        }

        std::vector<std::array<unsigned int, 3>> new_faces;
        new_faces.reserve(F * 4);
        std::vector<std::array<unsigned int, 3>> corner0_faces, corner1_faces, corner2_faces, center_faces;

        corner0_faces.reserve(faces.size());
        corner1_faces.reserve(faces.size());
        corner2_faces.reserve(faces.size());
        center_faces.reserve(faces.size());

        for (const auto& f : faces) {
            unsigned int v0 = f[0];
            unsigned int v1 = f[1];
            unsigned int v2 = f[2];

            unsigned int e01 = edge2mid[EdgeKey(v0, v1)];
            unsigned int e12 = edge2mid[EdgeKey(v1, v2)];
            unsigned int e20 = edge2mid[EdgeKey(v2, v0)];

            corner0_faces.push_back({ v0, e01, e20 });
            corner1_faces.push_back({ v1, e12, e01 });
            corner2_faces.push_back({ v2, e20, e12 });
            center_faces.push_back({ e12, e20, e01 });
        }
        new_faces.insert(new_faces.end(), corner0_faces.begin(), corner0_faces.end());
        new_faces.insert(new_faces.end(), corner1_faces.begin(), corner1_faces.end());
        new_faces.insert(new_faces.end(), corner2_faces.begin(), corner2_faces.end());
        new_faces.insert(new_faces.end(), center_faces.begin(), center_faces.end());
        vertices.swap(new_vertices);
        faces.swap(new_faces);
    }

    std::string out_path ="../data/up_sample.obj";
    std::ofstream out(out_path);

    for (const auto& v : vertices) {
        out << "v " << v(0) << " " << v(1) << " " << v(2) << "\n";
    }

    for (const auto& f : faces) {
        out << "f "
            << f[0] + 1 << " "
            << f[1] + 1 << " "
            << f[2] + 1 << "\n";
    }

    out.close();
}

//void BaseMesh::remesh_(unsigned int level)
//{
    //valid test
   // if(areElementsUnique(test_faces)==false)
    //{
     //   std::cout<<"The test faces are not unique!"<<std::endl;
	//}



    //int nv = n_vertices();
    //int nf = n_faces();
    //int ne = n_edges();

    //int phi = std::pow(2, level);
    //int re_nv = nv + nf * ((phi - 2) * (phi - 1) / 2) + ne * (phi - 1);
    //int n_inner_f = nf * ((phi - 2) * (phi - 1) / 2);

    //std::vector<Eigen::Vector3d> points(re_nv);
    //std::vector<unsigned int> basefvids(nf * 3);
    //std::vector<bool> is_marked_semiv(re_nv, false);

    //for (VertexIter v = vertices_begin(); v != vertices_end(); v++) {
    //    is_marked_semiv[v->idx()] = true;
    //    for (unsigned int k = 0; k < 3; k++) points[v->idx()](k) = point(*v)[k];//remain points
    //}

    //for (FaceIter f = faces_begin(); f != faces_end(); f++) {
    //    unsigned int k = 0;
    //    for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); fv++) {
    //        basefvids[3 * f->idx() + k] = fv->idx();
    //        k++;
    //    }
    //}
    //int num = 0;
    //for (FaceIter fit = faces_begin(); fit != faces_end(); fit++)
    //{
    //    //face map vertex remove information
    //    std::array<unsigned int, 3>fuv;
    //    
    //    // get subdivied mesh of vertices from each face of base domain in a row-by-row way
    //    std::vector<unsigned int> ids;
    //    std::vector<std::vector<unsigned int> > eids(3);
    //    HalfedgeHandle he;
    //    EdgeHandle e;
    //    for (FaceHalfedgeIter fh = fh_iter(*fit); fh.is_valid(); fh++)
    //    {
    //        VertexHandle v1, v2;
    //        v1 = from_vertex_handle(*fh);
    //        v2 = to_vertex_handle(*fh);
    //        if (v1.idx() == (int)basefvids[3 * fit->idx()])
    //        {
    //            he = *fh;
    //            break;
    //        }
    //    }
    //    for (unsigned int k = 0; k < 3; k++) {
    //        e = edge_handle(he);
    //        if (he.idx() == halfedge_handle(e, 0).idx()) {
    //            for (int i = 0; i < phi - 1; i++) {
    //                eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
    //            }
    //        }
    //        else {
    //            for (int i = phi - 2; i >= 0; i--) {
    //                eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
    //            }
    //        }
    //        he = next_halfedge_handle(he);
    //    }
    //    std::vector<std::vector<unsigned int> > subdivide_ids(phi + 1);
    //    subdivide_ids[0].push_back(basefvids[3 * fit->idx()]);
    //    unsigned int pos = 0;
    //    for (int k = 1; k < phi + 1; k++) {
    //        if (k == phi) {
    //            subdivide_ids[k].push_back(basefvids[3 * fit->idx() + 1]);
    //            for (int i = 0; i < (int)eids[1].size(); i++) {
    //                subdivide_ids[k].push_back(eids[1][i]);
    //            }
    //            subdivide_ids[k].push_back(basefvids[3 * fit->idx() + 2]);
    //        }
    //        else {
    //            subdivide_ids[k].push_back(eids[0][k - 1]);
    //            for (int j = 0; j < k - 1; j++) {
    //                subdivide_ids[k].push_back(nv + fit->idx() * ((phi - 2) * (phi - 1) / 2) + pos);
    //                pos++;
    //            }
    //            subdivide_ids[k].push_back(eids[2][phi - 1 - k]);
    //        }
    //    }
    //    // subdivide points location
    //    unsigned int removeIM_id ;
    //    vertex_remove_data dataIM ;
    //    std::array<double, 3> barycoord_values;
    //    for (unsigned int i = 0; i < subdivide_ids.size(); i++) {
    //        for (unsigned int j = 0; j < subdivide_ids[i].size(); j++) {
    //            unsigned int k = subdivide_ids[i][j];
    //            if (!is_marked_semiv[k]) {
    //                Point query(0, 0, 0);
    //                barycoord_values = { (phi - i) * 1.0 / phi ,(i - j) * 1.0 / phi ,j * 1.0 / phi };
    //                int vsize = 0;
    //                for (FaceVertexIter fv = fv_iter(*fit); fv.is_valid(); fv++)
    //                {
    //                    fuv[vsize] = data(*fv).ori_id;
    //                    vsize++;
    //                }
    //                while (std::find(all_orifaces.begin(), all_orifaces.end(), fuv) == all_orifaces.end())
    //                {
    //                    bool found_triangle = false;
    //                    if (Face2IM.find(fuv) == Face2IM.end()) {
    //                        std::cerr << "Error: Face not found in Face2IM" << std::endl;
    //                        break;
    //                    }
    //                    removeIM_id = Face2IM[fuv];
    //                    if (removeIM_id >= vremoveIM.size()) {
    //                        std::cerr << "Error: Invalid removeIM_id" << std::endl;
    //                        break;
    //                    }
    //                    dataIM = vremoveIM[removeIM_id];

    //                    for (unsigned int s = 0; s < 2; s++)
    //                    {
    //                        query[s] = barycoord_values[0] * dataIM.uv_coords[dataIM.g2l[fuv[0]]][s] +
    //                            barycoord_values[1] * dataIM.uv_coords[dataIM.g2l[fuv[1]]][s] +
    //                            barycoord_values[2] * dataIM.uv_coords[dataIM.g2l[fuv[2]]][s];
    //                    }
    //                    //find cloest triangle
    //                    double min_distance = 1.0;
    //                    std::array<unsigned int, 3> closest_fuv = fuv;
    //                    std::array<double, 3> closest_bary = barycoord_values;

    //                    for (unsigned int m = 0; m < dataIM.FUV_pre.size(); m++)
    //                    {
    //                        std::array<Vec2d, 3> pre_fuv;
    //                        for (unsigned int n = 0; n < 3; n++)
    //                        {
    //                            unsigned int vid = dataIM.FUV_pre[m][n];
    //                            pre_fuv[n] = Vec2d(dataIM.uv_coords[dataIM.g2l[vid]][0], dataIM.uv_coords[dataIM.g2l[vid]][1]);
    //                        }

    //                        auto [alpha, beta, gamma] = cal_bar_2dcoord(Vec2d(query[0], query[1]), pre_fuv);
    //                        double min_coord = std::min(alpha,std::min( beta, gamma ));
    //                        double distance = -min_coord;
    //                        if (distance < min_distance)
    //                        {
    //                            min_distance = distance;
    //                            closest_fuv = dataIM.FUV_pre[m];
    //                            closest_bary = { alpha, beta, gamma };
    //                        }
    //                        
    //                    }
    //                    double alpha = std::max(0.0, closest_bary[0]);
    //                    double beta = std::max(0.0, closest_bary[1]);
    //                    double gamma = std::max(0.0, closest_bary[2]);
    //                    double sum = alpha + beta + gamma;
    //                    barycoord_values[0] = alpha / sum;
    //                    barycoord_values[1] = beta / sum;
    //                    barycoord_values[2] = gamma / sum;
    //                    fuv = closest_fuv;
    //                    //std::cout << "barycoord_values[0]: " << barycoord_values[0] << "barycoord_values[1]: " << barycoord_values[1] << "barycoord_values[2]: " << barycoord_values[2] << std::endl;
				//		//std::cout << "sum: " << sum << std::endl;
    //                    
    //                    /*std::ofstream out("../data/before_error_fuv.obj");
    //                    for (unsigned int s = 0; s < dataIM.uv_coords.size(); s++)
    //                    {
    //                        out << "v " << dataIM.uv_coords[s][0] << " " << dataIM.uv_coords[s][1] << " 0" << std::endl;
    //                    }
    //                    unsigned int nvs = dataIM.uv_coords.size() - 1;
    //                    for (unsigned int m = 0; m < nvs; m++)
    //                    {
    //                        out << "f 1 " << (m + 1) % nvs + 2 << " " << m + 2 << std::endl;
    //                    }
    //                    out.close();

    //                    std::ofstream out1("../data/before_error_query.obj");
    //                    out1 << "v " << query[0] << " " << query[1] << " 0" << std::endl;
    //                    out1.close();*/
    //                }
    //                //std::cout << "face id: " <<fit->idx()<< std::endl;
    //                points[k] = Eigen::Vector3d::Zero(0);
				//	Point p =  barycoord_values[0] * positions[fuv[0]]+
    //                    barycoord_values[1]*positions[fuv[1]] + 
    //                    barycoord_values[2] * positions[fuv[2]];
				//	for (unsigned int s = 0; s < 3; s++) points[k](s) = p[s];
    //                is_marked_semiv[k] = true;
    //            }

    //        }

    //    }
    //}
    //int count = 0;
    //for (int i = 0; i < re_nv; i++) {
    //    if (!is_marked_semiv[i]) {
    //        count++;
    //        printf("missing id: %d\n", i);
    //    }
    //}
    //if (count > 0) {
    //    printf("# %d points location is failed\n", count);
    //    return;
    //}
    //std::ofstream out("../data/subdivision.obj");
    //for (int i = 0; i < re_nv; i++) {
    //    out << "v " << points[i](0) << " " << points[i](1) << " " << points[i](2) << std::endl;
    //}
    //for (FaceIter f = faces_begin(); f != faces_end(); f++) {
    //    // get subdivied mesh of vertices from each face of base domain in a row-by-row way.
    //    std::vector<unsigned int> ids;
    //    std::vector<std::vector<unsigned int> > eids(3);
    //    HalfedgeHandle he;
    //    EdgeHandle e;
    //    for (FaceHalfedgeIter fh = fh_iter(*f); fh.is_valid(); fh++) {
    //        VertexHandle v1, v2;
    //        v1 = from_vertex_handle(*fh);
    //        v2 = to_vertex_handle(*fh);
    //        if (v1.idx() == (int)basefvids[3 * f->idx()]) {
    //            he = *fh;
    //            break;
    //        }
    //    }
    //    for (unsigned int k = 0; k < 3; k++) {
    //        e = edge_handle(he);
    //        if (he.idx() == halfedge_handle(e, 0).idx()) {
    //            for (int i = 0; i < phi - 1; i++) {
    //                eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
    //            }
    //        }
    //        else {
    //            for (int i = phi - 2; i >= 0; i--) {
    //                eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
    //            }
    //        }
    //        he = next_halfedge_handle(he);
    //    }
    //    std::vector<std::vector<unsigned int> > subdivide_ids(phi + 1);
    //    subdivide_ids[0].push_back(basefvids[3 * f->idx()]);
    //    unsigned int pos = 0;
    //    for (int k = 1; k < phi + 1; k++) {
    //        if (k == phi) {
    //            subdivide_ids[k].push_back(basefvids[3 * f->idx() + 1]);
    //            for (int i = 0; i < (int)eids[1].size(); i++) {
    //                subdivide_ids[k].push_back(eids[1][i]);
    //            }
    //            subdivide_ids[k].push_back(basefvids[3 * f->idx() + 2]);
    //        }
    //        else {
    //            subdivide_ids[k].push_back(eids[0][k - 1]);
    //            for (int j = 0; j < k - 1; j++) {
    //                subdivide_ids[k].push_back(nv + f->idx() * ((phi - 2) * (phi - 1) / 2) + pos);
    //                pos++;
    //            }
    //            subdivide_ids[k].push_back(eids[2][phi - 1 - k]);
    //        }
    //    }
    //    for (int i = 0; i < phi; i++) {
    //        for (int j = 0; j < phi - 1 - i; j++) {
    //            out << "f " << subdivide_ids[phi - i - 1][j] + 1 << " " << subdivide_ids[phi - i][j + 1] + 1 << " " << subdivide_ids[phi - i - 1][j + 1] + 1 << std::endl;
    //        }
    //        for (int j = 0; j < phi - i; j++) {
    //            out << "f " << subdivide_ids[phi - i][j] + 1 << " " << subdivide_ids[phi - i][j + 1] + 1 << " " << subdivide_ids[phi - i - 1][j] + 1 << std::endl;
    //        }
    //    }
    //}
    //out.close();

    
//}
//void BaseMesh::remesh_(unsigned int level)
//{
//    // Valid test
//    if (areElementsUnique(test_faces) == false) {
//        std::cout << "The test faces are not unique!" << std::endl;
//    }
//
//    // 1. Store Base Mesh Data (Reference for Projection)
//    std::vector<std::array<unsigned int, 3>> base_faces_ori_ids;
//    std::map<unsigned int, std::vector<unsigned int>> ori_id_to_base_faces_idx;
//
//    unsigned int f_idx = 0;
//    for (FaceIter f = faces_begin(); f != faces_end(); ++f) {
//        std::array<unsigned int, 3> fuv;
//        int k = 0;
//        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); ++fv) {
//            unsigned int oid = data(*fv).ori_id;
//            fuv[k++] = oid;
//            ori_id_to_base_faces_idx[oid].push_back(f_idx);
//        }
//        base_faces_ori_ids.push_back(fuv);
//        f_idx++;
//    }
//
//    // 2. Initialize Dynamic Mesh Data (Topology & Geometry)
//    std::vector<Eigen::Vector3d> vertices;
//    std::vector<std::array<unsigned int, 3>> current_faces; // Indices into 'vertices'
//    std::vector<std::map<unsigned int, double>> vert_tracking; // lineage: ori_id -> weight
//
//    // Init vertices and tracking
//    for (VertexIter v = vertices_begin(); v != vertices_end(); ++v) {
//        vertices.emplace_back(point(*v)[0], point(*v)[1], point(*v)[2]);
//        std::map<unsigned int, double> tracker;
//        tracker[data(*v).ori_id] = 1.0;
//        vert_tracking.push_back(tracker);
//    }
//
//    // Init current faces (using indices 0..N)
//    for (FaceIter f = faces_begin(); f != faces_end(); ++f) {
//        std::array<unsigned int, 3> f_vidx;
//        int k = 0;
//        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); ++fv) {
//            f_vidx[k++] = fv->idx(); // Use internal index corresponding to 'vertices' vector order
//        }
//        current_faces.push_back(f_vidx);
//    }
//
//
//    // 3. Subdivision Loop
//    for (unsigned int lv = 0; lv < level; ++lv) {
//
//        std::vector<Eigen::Vector3d> new_vertices = vertices;
//        std::vector<std::array<unsigned int, 3>> new_faces;
//        std::vector<std::map<unsigned int, double>> new_vert_tracking = vert_tracking;
//        new_faces.reserve(current_faces.size() * 4);
//
//        std::map<EdgeKey, unsigned int> edge2mid;
//
//        for (const auto& f : current_faces) {
//            for (int e = 0; e < 3; ++e) {
//                unsigned int v0 = f[e];
//                unsigned int v1 = f[(e + 1) % 3];
//
//                EdgeKey key(v0, v1);
//                if (edge2mid.count(key)) continue;
//
//                unsigned int mid_id = new_vertices.size();
//                edge2mid[key] = mid_id;
//
//                // --- A. Compute Provenance (Weights relative to Base Mesh) ---
//                std::map<unsigned int, double> mid_track;
//                for (auto const& [oid, w] : vert_tracking[v0]) mid_track[oid] += 0.5 * w;
//                for (auto const& [oid, w] : vert_tracking[v1]) mid_track[oid] += 0.5 * w;
//                new_vert_tracking.push_back(mid_track);
//
//                // --- B. Find Base Face ---
//                // Find a Base Face that contains the Base Vertices contributing to this point.
//                std::vector<unsigned int> candidates;
//                bool first = true;
//                for (auto const& [oid, w] : mid_track) {
//                    if (w < 1e-4) continue;
//                    if (first) {
//                        if (ori_id_to_base_faces_idx.find(oid) != ori_id_to_base_faces_idx.end())
//                            candidates = ori_id_to_base_faces_idx[oid];
//                        first = false;
//                    }
//                    else {
//                        if (ori_id_to_base_faces_idx.find(oid) != ori_id_to_base_faces_idx.end()) {
//                            std::vector<unsigned int> next_set = ori_id_to_base_faces_idx[oid];
//                            std::vector<unsigned int> intersection;
//                            std::sort(candidates.begin(), candidates.end());
//                            std::sort(next_set.begin(), next_set.end());
//                            std::set_intersection(candidates.begin(), candidates.end(),
//                                next_set.begin(), next_set.end(),
//                                std::back_inserter(intersection));
//                            candidates = intersection;
//                        }
//                    }
//                    if (candidates.empty()) break;
//                }
//
//                // --- C. Inverse Map Projection ---
//                Point final_p(0, 0, 0);
//                bool projection_success = false;
//
//                if (!candidates.empty()) {
//                    unsigned int base_face_idx = candidates[0];
//                    std::array<unsigned int, 3> base_fuv = base_faces_ori_ids[base_face_idx];
//
//                    // Reconstruct full barycentrics for this specific Base Face
//                    std::array<double, 3> bary = { 0,0,0 };
//                    for (int k = 0; k < 3; ++k) {
//                        if (mid_track.find(base_fuv[k]) != mid_track.end()) {
//                            bary[k] = mid_track[base_fuv[k]];
//                        }
//                    }
//                    // Normalize (just in case)
//                    double sum = bary[0] + bary[1] + bary[2];
//                    if (sum > 1e-9) { bary[0] /= sum; bary[1] /= sum; bary[2] /= sum; }
//
//                    // Execute IM Logic
//                    std::array<unsigned int, 3> curr_fuv = base_fuv;
//                    std::array<double, 3> curr_bary = bary;
//
//                    // Logic from original code: walk back the history
//                    while (std::find(all_orifaces.begin(), all_orifaces.end(), curr_fuv) == all_orifaces.end())
//                    {
//                        if (Face2IM.find(curr_fuv) == Face2IM.end()) {
//                            // Fallback if lookup fails
//                            break;
//                        }
//                        unsigned int removeIM_id = Face2IM[curr_fuv];
//                        if (removeIM_id >= vremoveIM.size()) break;
//
//                        vertex_remove_data& dataIM = vremoveIM[removeIM_id];
//
//                        // 1. Compute 2D position in the flattened plane
//                        Point query(0, 0, 0);
//                        for (unsigned int s = 0; s < 2; s++) {
//                            query[s] = curr_bary[0] * dataIM.uv_coords[dataIM.g2l[curr_fuv[0]]][s] +
//                                curr_bary[1] * dataIM.uv_coords[dataIM.g2l[curr_fuv[1]]][s] +
//                                curr_bary[2] * dataIM.uv_coords[dataIM.g2l[curr_fuv[2]]][s];
//                        }
//
//                        // 2. Find closest triangle in the "previous" mesh (FUV_pre)
//                        double min_dist = 1.0;
//                        std::array<unsigned int, 3> best_fuv = curr_fuv;
//                        std::array<double, 3> best_bary = curr_bary;
//                        bool found = false;
//
//                        for (unsigned int m = 0; m < dataIM.FUV_pre.size(); m++) {
//                            std::array<Vec2d, 3> tri;
//                            for (int n = 0; n < 3; n++) {
//                                unsigned int vid = dataIM.FUV_pre[m][n];
//                                tri[n] = Vec2d(dataIM.uv_coords[dataIM.g2l[vid]][0],
//                                    dataIM.uv_coords[dataIM.g2l[vid]][1]);
//                            }
//
//                            auto [a0, a1, a2] = cal_bar_2dcoord(Vec2d(query[0], query[1]), tri);
//
//                            // Check if inside (all positive) or closest
//                            double d = -std::min({ a0, a1, a2 });
//                            if (d < min_dist) {
//                                min_dist = d;
//                                best_fuv = dataIM.FUV_pre[m];
//                                best_bary = { a0, a1, a2 };
//                                found = true;
//                            }
//                        }
//
//                        if (found) {
//                            double s_sum = std::max(1e-12, best_bary[0] + best_bary[1] + best_bary[2]);
//                            curr_bary[0] = std::max(0.0, best_bary[0]) / s_sum;
//                            curr_bary[1] = std::max(0.0, best_bary[1]) / s_sum;
//                            curr_bary[2] = std::max(0.0, best_bary[2]) / s_sum;
//                            curr_fuv = best_fuv;
//                        }
//                        else {
//                            break;
//                        }
//                    }
//
//                    // Compute final position using original vertex positions
//                    // curr_fuv should now be in all_orifaces (Original Mesh Faces)
//                    final_p = curr_bary[0] * positions[curr_fuv[0]] +
//                        curr_bary[1] * positions[curr_fuv[1]] +
//                        curr_bary[2] * positions[curr_fuv[2]];
//
//                    new_vertices.emplace_back(final_p[0], final_p[1], final_p[2]);
//                    projection_success = true;
//                }
//
//                if (!projection_success) {
//                    // Fallback to Linear Subdivision if IM fails
//                    Eigen::Vector3d p = 0.5 * (vertices[v0] + vertices[v1]);
//                    new_vertices.emplace_back(p[0], p[1], p[2]);
//                }
//            }
//        }
//
//        std::vector<std::array<unsigned int, 3>> corner0_faces, corner1_faces, corner2_faces, center_faces;
//
//        corner0_faces.reserve(current_faces.size());
//        corner1_faces.reserve(current_faces.size());
//        corner2_faces.reserve(current_faces.size());
//        center_faces.reserve(current_faces.size());
//
//        // --- D. Topology Reassembly (Standard Loop) ---
//        for (const auto& f : current_faces) {
//            unsigned int v0 = f[0];
//            unsigned int v1 = f[1];
//            unsigned int v2 = f[2];
//
//            unsigned int e01 = edge2mid[EdgeKey(v0, v1)];
//            unsigned int e12 = edge2mid[EdgeKey(v1, v2)];
//            unsigned int e20 = edge2mid[EdgeKey(v2, v0)];
//
//            corner0_faces.push_back({ v0, e01, e20 });
//            corner1_faces.push_back({ v1, e12, e01 });
//            corner2_faces.push_back({ v2, e20, e12 });
//            center_faces.push_back({ e12, e20, e01 });
//        }
//
//        new_faces.insert(new_faces.end(), corner0_faces.begin(), corner0_faces.end());
//        new_faces.insert(new_faces.end(), corner1_faces.begin(), corner1_faces.end());
//        new_faces.insert(new_faces.end(), corner2_faces.begin(), corner2_faces.end());
//        new_faces.insert(new_faces.end(), center_faces.begin(), center_faces.end());
//        vertices.swap(new_vertices);
//        current_faces.swap(new_faces);
//        vert_tracking.swap(new_vert_tracking);
//    }
//
//    // Output
//    std::ofstream out("loop_IM_" + std::to_string(level) + ".obj");
//    std::cout << "face size: " << current_faces.size() << std::endl;
//
//    for (const auto& v : vertices)
//        out << "v " << v(0) << " " << v(1) << " " << v(2) << "\n";
//
//    for (const auto& f : current_faces)
//        out << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << "\n";
//
//    out.close();
//}

void BaseMesh::remesh_(unsigned int level)
{
    if (areElementsUnique(test_faces) == false) {
        std::cout << "The test faces are not unique!" << std::endl;
    }

    std::vector<std::array<unsigned int, 3>> base_faces_ori_ids;
    std::map<unsigned int, std::vector<unsigned int>> ori_id_to_base_faces_idx;

    unsigned int f_idx = 0;
    for (FaceIter f = faces_begin(); f != faces_end(); ++f) {
        std::array<unsigned int, 3> fuv;
        int k = 0;
        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); ++fv) {
            unsigned int oid = data(*fv).ori_id;
            fuv[k++] = oid;
            ori_id_to_base_faces_idx[oid].push_back(f_idx);
        }
        base_faces_ori_ids.push_back(fuv);
        f_idx++;
    }

    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<unsigned int, 3>> current_faces;
    std::vector<std::map<unsigned int, double>> vert_tracking;

    for (VertexIter v = vertices_begin(); v != vertices_end(); ++v) {
        vertices.emplace_back(point(*v)[0], point(*v)[1], point(*v)[2]);
        std::map<unsigned int, double> tracker;
        tracker[data(*v).ori_id] = 1.0;
        vert_tracking.push_back(tracker);
    }

    for (FaceIter f = faces_begin(); f != faces_end(); ++f) {
        std::array<unsigned int, 3> f_vidx;
        int k = 0;
        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); ++fv) {
            f_vidx[k++] = fv->idx();
        }
        current_faces.push_back(f_vidx);
    }

    for (unsigned int lv = 0; lv < level; ++lv) {
        std::vector<Eigen::Vector3d> new_vertices = vertices;
        std::vector<std::array<unsigned int, 3>> new_faces;
        std::vector<std::map<unsigned int, double>> new_vert_tracking = vert_tracking;

        std::set<EdgeKey> unique_edges;
        for (const auto& f : current_faces) {
            unique_edges.insert(EdgeKey(f[0], f[1]));
            unique_edges.insert(EdgeKey(f[1], f[2]));
            unique_edges.insert(EdgeKey(f[2], f[0]));
        }

        std::map<EdgeKey, unsigned int> edge2mid;
        for (const auto& key : unique_edges) {
            unsigned int mid_id = new_vertices.size();
            edge2mid[key] = mid_id;

            std::map<unsigned int, double> mid_track;
            for (auto const& [oid, w] : vert_tracking[key.v0]) mid_track[oid] += 0.5 * w;
            for (auto const& [oid, w] : vert_tracking[key.v1]) mid_track[oid] += 0.5 * w;
            new_vert_tracking.push_back(mid_track);

            std::vector<unsigned int> candidates;
            bool first = true;
            for (auto const& [oid, w] : mid_track) {
                if (w < 1e-4) continue;
                if (first) {
                    if (ori_id_to_base_faces_idx.count(oid))
                        candidates = ori_id_to_base_faces_idx[oid];
                    first = false;
                }
                else {
                    if (ori_id_to_base_faces_idx.count(oid)) {
                        std::vector<unsigned int> next_set = ori_id_to_base_faces_idx[oid];
                        std::vector<unsigned int> intersection;
                        std::sort(candidates.begin(), candidates.end());
                        std::sort(next_set.begin(), next_set.end());
                        std::set_intersection(candidates.begin(), candidates.end(),
                            next_set.begin(), next_set.end(), std::back_inserter(intersection));
                        candidates = intersection;
                    }
                }
                if (candidates.empty()) break;
            }

            Point final_p(0, 0, 0);
            bool projection_success = false;

            if (!candidates.empty()) {
                unsigned int base_face_idx = candidates[0];
                std::array<unsigned int, 3> curr_fuv = base_faces_ori_ids[base_face_idx];

                std::array<double, 3> curr_bary = { 0,0,0 };
                for (int k = 0; k < 3; ++k) {
                    if (mid_track.count(curr_fuv[k])) curr_bary[k] = mid_track[curr_fuv[k]];
                }
                double sum = curr_bary[0] + curr_bary[1] + curr_bary[2];
                if (sum > 1e-9) { curr_bary[0] /= sum; curr_bary[1] /= sum; curr_bary[2] /= sum; }

                while (std::find(all_orifaces.begin(), all_orifaces.end(), curr_fuv) == all_orifaces.end()) {
                    if (Face2IM.find(curr_fuv) == Face2IM.end()) break;
                    unsigned int removeIM_id = Face2IM[curr_fuv];
                    if (removeIM_id >= vremoveIM.size()) break;

                    vertex_remove_data& dataIM = vremoveIM[removeIM_id];
                    Point query(0, 0, 0);
                    for (unsigned int s = 0; s < 2; s++) {
                        query[s] = curr_bary[0] * dataIM.uv_coords[dataIM.g2l[curr_fuv[0]]][s] +
                            curr_bary[1] * dataIM.uv_coords[dataIM.g2l[curr_fuv[1]]][s] +
                            curr_bary[2] * dataIM.uv_coords[dataIM.g2l[curr_fuv[2]]][s];
                    }

                    double min_dist = 1e10;
                    std::array<unsigned int, 3> best_fuv = curr_fuv;
                    std::array<double, 3> best_bary = curr_bary;
                    bool found = false;

                    for (unsigned int m = 0; m < dataIM.FUV_pre.size(); m++) {
                        std::array<Vec2d, 3> tri;
                        for (int n = 0; n < 3; n++) {
                            unsigned int vid = dataIM.FUV_pre[m][n];
                            tri[n] = Vec2d(dataIM.uv_coords[dataIM.g2l[vid]][0], dataIM.uv_coords[dataIM.g2l[vid]][1]);
                        }
                        auto [a0, a1, a2] = cal_bar_2dcoord(Vec2d(query[0], query[1]), tri);
                        double d = -std::min({ a0, a1, a2 });
                        if (d < min_dist) {
                            min_dist = d;
                            best_fuv = dataIM.FUV_pre[m];
                            best_bary = { a0, a1, a2 };
                            found = true;
                        }
                    }

                    if (found) {
                        double s_sum = std::max(1e-12, best_bary[0] + best_bary[1] + best_bary[2]);
                        curr_bary = { std::max(0.0, best_bary[0]) / s_sum, std::max(0.0, best_bary[1]) / s_sum, std::max(0.0, best_bary[2]) / s_sum };
                        curr_fuv = best_fuv;
                    }
                    else break;
                }
                final_p = curr_bary[0] * positions[curr_fuv[0]] + curr_bary[1] * positions[curr_fuv[1]] + curr_bary[2] * positions[curr_fuv[2]];
                new_vertices.emplace_back(final_p[0], final_p[1], final_p[2]);
                projection_success = true;
            }

            if (!projection_success) {
                Eigen::Vector3d p = 0.5 * (vertices[key.v0] + vertices[key.v1]);
                new_vertices.emplace_back(p[0], p[1], p[2]);
            }
        }

        std::vector<std::array<unsigned int, 3>> corner0_faces, corner1_faces, corner2_faces, center_faces;

        corner0_faces.reserve(current_faces.size());
        corner1_faces.reserve(current_faces.size());
        corner2_faces.reserve(current_faces.size());
        center_faces.reserve(current_faces.size());

        for (const auto& f : current_faces) {
            unsigned int v0 = f[0], v1 = f[1], v2 = f[2];
            unsigned int e01 = edge2mid[EdgeKey(v0, v1)]; 
            unsigned int e12 = edge2mid[EdgeKey(v1, v2)]; 
            unsigned int e20 = edge2mid[EdgeKey(v2, v0)];

            corner0_faces.push_back({ v0, e01, e20 });
            corner1_faces.push_back({ v1, e12, e01 });
            corner2_faces.push_back({ v2, e20, e12 });
            center_faces.push_back({ e12, e20, e01 });
        }

        new_faces.insert(new_faces.end(), corner0_faces.begin(), corner0_faces.end());
        new_faces.insert(new_faces.end(), corner1_faces.begin(), corner1_faces.end());
        new_faces.insert(new_faces.end(), corner2_faces.begin(), corner2_faces.end());
        new_faces.insert(new_faces.end(), center_faces.begin(), center_faces.end());
        vertices.swap(new_vertices);
        current_faces.swap(new_faces);
        vert_tracking.swap(new_vert_tracking);
    }

    std::ofstream out("../data/loop_IM_sync_" + std::to_string(level) + ".obj");
    for (const auto& v : vertices) out << "v " << v(0) << " " << v(1) << " " << v(2) << "\n";
    for (const auto& f : current_faces) out << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << "\n";
    out.close();
}

