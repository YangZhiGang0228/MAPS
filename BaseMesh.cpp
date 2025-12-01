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
            return nesting_level % 2 == 1;
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

bool BaseMesh::In_2Dtriangle(const Vec2d& point2d,const std::array<Vec2d,3>&f2d)
{
    /*double esp = 1e-10;
    Point point0 (point2d[0],point2d[1],0);
    Point point1 (f2d[0][0],f2d[0][1],0);
    Point point2 (f2d[1][0],f2d[1][1],0);
    Point point3 (f2d[2][0],f2d[2][1],0);

    double z1 = (point0-point1).cross(point2-point1)[2];
    double z2 = (point0-point2).cross(point3-point2)[2];
    double z3 = (point0-point3).cross(point1-point3)[2];

    if((z1>=-esp && z2 >= -esp && z3 >= -esp) || (z1<=esp && z2 <= esp && z3 <= esp))
    {
        return true;
    }
    return false;*/

	Point_2 p(point2d[0], point2d[1]);
	Point_2 p1(f2d[0][0], f2d[0][1]);
	Point_2 p2(f2d[1][0], f2d[1][1]);
	Point_2 p3(f2d[2][0], f2d[2][1]);

	CGAL::Orientation o1 = CGAL::orientation(p1, p2, p);
	CGAL::Orientation o2 = CGAL::orientation(p2, p3, p);
	CGAL::Orientation o3 = CGAL::orientation(p3, p1, p);

    return ((o1 == o2 && o2 == o3) ||
        (o1 == CGAL::COLLINEAR && o2 == o3) ||
        (o2 == CGAL::COLLINEAR && o1 == o3) ||
        (o3 == CGAL::COLLINEAR && o1 == o2));
}

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

void BaseMesh::CGAL_CDT(const std::vector<std::pair<VertexHandle, Vec2d>>& coordinates, std::vector<std::array<VertexHandle, 3>>& faces)
{
    MyCGALStuff::CDT cdt;
    std::vector<MyCGALStuff::Vhandle>cdt_vhandles;
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
        }
    }
    for (unsigned int i = 0; i < cdt_vhandles.size(); i++)
    {
        cdt.insert_constraint(cdt_vhandles[i], cdt_vhandles[(i + 1) % cdt_vhandles.size()]);
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
    std::vector<std::pair<VertexHandle, Vec2d>> coordinates;
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
    //delete one ring of this vertex
    for (unsigned int i = 0; i < ringFaces.size(); i++)
    {
        delete_face(ringFaces[i], false);
    }

    std::vector<std::array<VertexHandle, 3>> newFaces;
    CGAL_CDT(coordinates, newFaces);
    //add new triangulation£¬check the new face isn't complex edge and uv flip
    std::vector<VertexHandle> newfhandles;
    std::vector<FaceHandle> testfhandles;
    unsigned int new_face_size = newFaces.size();
    for (unsigned int i = 0; i < new_face_size; i++)
    {
        std::array<Vec2d, 3>UV_face_post;
        newfhandles.clear();
        for (unsigned int j = 0; j < 3; j++)
        {
            newfhandles.push_back(newFaces[i][j]);
            for (unsigned int k = 0; k < coordinates.size(); k++)
            {
                if (coordinates[k].first == newFaces[i][j])
                {
                    UV_face_post[j] = (coordinates[k].second);
                }
            }
        }
        if (!check_uv_flip(UV_face_post))
        {
            flag = true;
            std::cout << "hello" << std::endl;
        }
        FaceHandle f = add_face(newfhandles);
        if (!f.is_valid())
        {
            flag = true;
        }
        else
        {
            testfhandles.push_back(f);
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
        barycoord[0] = std::make_pair(ori_vids[cloest_f[0]], alpha / sum);
        barycoord[1] = std::make_pair(ori_vids[cloest_f[1]], beta / sum);
        barycoord[2] = std::make_pair(ori_vids[cloest_f[2]], gamma / sum);
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
                VertexHandle& inner_vertex = original_vhs[local_vertices[i]];
                Vec2d uv = Vec2d(0, 0);
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
                barycoord[0] = std::make_pair(ori_vids[cloest_f1[0]], alpha1 / sum1);
                barycoord[1] = std::make_pair(ori_vids[cloest_f1[1]], beta1 / sum1);
                barycoord[2] = std::make_pair(ori_vids[cloest_f1[2]], gamma1 / sum1);
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
        data(*vertexiter).canberemoved = true;

    }

    initial_weights(0.5);
    while(!removal_P.empty() && current_nfs>remaining_nfs)
    {
        auto vertexhandle = removal_P.top();
        removal_P.pop();

        if (!data(vertexhandle).canberemoved || data(vertexhandle).isdeleted)
        {
            continue;
        }

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

void BaseMesh::MidSubdivision(unsigned int level)
{
    BaseMesh remesh_;
    int nv = n_vertices();
    int nf = n_faces();
    int ne = n_edges();
    int phi = std::pow(2, level);
    int re_nv = nv + nf * ((phi - 2) * (phi - 1) / 2) + ne * (phi - 1);

    std::vector<Eigen::Vector3d> points(re_nv);
    std::vector<std::vector<unsigned int> > fvids(nf);
    // nv
    for (VertexIter v = vertices_begin(); v != vertices_end(); v++) {
        for (unsigned int k = 0; k < 3; k++) points[v->idx()](k) = point(*v)[k];
    }
    unsigned int pos = 0;
    for (FaceIter f = faces_begin(); f != faces_end(); f++) {
        fvids[f->idx()].resize(3);
        pos = 0;
        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); fv++) {
            fvids[f->idx()][pos] = fv->idx();
            pos++;
        }
    }
    pos = 0;
    // nf * ((phi-2)*(phi-1)/2)
    Eigen::Vector3d p;
    for (FaceIter f = faces_begin(); f != faces_end(); f++) {
        for (int i = 0; i < phi - 2; i++) {
            for (int j = 0; j < i + 1; j++) {
                p = (phi - 2 - i) * 1.0 / phi * points[fvids[f->idx()][0]] +
                    (i + 1 - j) * 1.0 / phi * points[fvids[f->idx()][1]] +
                    (j + 1) * 1.0 / phi * points[fvids[f->idx()][2]];
                points[nv + pos] = p;
                pos++;
            }
        }
    }
    pos += nv;
    // ne * (phi-1)
    for (EdgeIter e = edges_begin(); e != edges_end(); e++) {
        HalfedgeHandle he = halfedge_handle(*e, 0);
        VertexHandle v1 = from_vertex_handle(he);
        VertexHandle v2 = to_vertex_handle(he);
        for (int i = 0; i < phi - 1; i++) {
            p = (phi - 1 - i) * 1.0 / phi * points[v1.idx()] +
                (i + 1) * 1.0 / phi * points[v2.idx()];
            points[pos] = p;
            pos++;
        }
    }

    // 2. obtain subdivided mesh connectivity in each face of base domain//
    unsigned int n_inner_f = nf * ((phi - 2) * (phi - 1) / 2);
    std::vector<VertexHandle> vhandles(points.size());
    for (unsigned int i = 0; i < points.size(); i++) {
        VertexHandle v = remesh_.add_vertex(Point(points[i](0), points[i](1), points[i](2)));
        vhandles[i] = v;
    }
    for (FaceIter f = faces_begin(); f != faces_end(); f++) {
        // get subdivied mesh of vertices from each face of base domain in a row-by-row way.
        std::vector<unsigned int> ids;
        std::vector<std::vector<unsigned int> > eids(3);
        HalfedgeHandle he;
        EdgeHandle e;
        for (FaceHalfedgeIter fh = fh_iter(*f); fh.is_valid(); fh++) {
           VertexHandle v1, v2;
            v1 = from_vertex_handle(*fh);
            v2 = to_vertex_handle(*fh);
            if (v1.idx() == (int)fvids[f->idx()][0]) {
                he = *fh;
                break;
            }
        }
        for (unsigned int k = 0; k < 3; k++) {
            e = edge_handle(he);
            if (he.idx() == halfedge_handle(e, 0).idx()) {
                for (int i = 0; i < phi - 1; i++) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            else {
                for (int i = phi - 2; i >= 0; i--) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            he = next_halfedge_handle(he);
        }
        std::vector<std::vector<unsigned int> > subdivide_ids(phi + 1);
        subdivide_ids[0].push_back(fvids[f->idx()][0]);
        pos = 0;
        for (int k = 1; k < phi + 1; k++) {
            if (k == phi) {
                subdivide_ids[k].push_back(fvids[f->idx()][1]);
                for (int i = 0; i < (int)eids[1].size(); i++) {
                    subdivide_ids[k].push_back(eids[1][i]);
                }
                subdivide_ids[k].push_back(fvids[f->idx()][2]);
            }
            else {
                subdivide_ids[k].push_back(eids[0][k - 1]);
                for (int j = 0; j < k - 1; j++) {
                    subdivide_ids[k].push_back(nv + f->idx() * ((phi - 2) * (phi - 1) / 2) + pos);
                    pos++;
                }
                subdivide_ids[k].push_back(eids[2][phi - 1 - k]);
            }
        }

        std::vector<VertexHandle> fhandles;
        for (int i = 0; i < phi; i++) {
            for (int j = 0; j < phi - 1 - i; j++) {
                fhandles.clear();
                fhandles.push_back(vhandles[subdivide_ids[phi - i - 1][j]]);
                fhandles.push_back(vhandles[subdivide_ids[phi - i][j + 1]]);
                fhandles.push_back(vhandles[subdivide_ids[phi - i - 1][j + 1]]);
                remesh_.add_face(fhandles);
            }
            for (int j = 0; j < phi - i; j++) {
                fhandles.clear();
                fhandles.push_back(vhandles[subdivide_ids[phi - i][j]]);
                fhandles.push_back(vhandles[subdivide_ids[phi - i][j + 1]]);
                fhandles.push_back(vhandles[subdivide_ids[phi - i - 1][j]]);
                remesh_.add_face(fhandles);
            }
        }
    }
    OpenMesh::IO::write_mesh(remesh_, "../data/semi_regular_basedomain.obj");
}

void BaseMesh::remesh_(unsigned int level,string output_dir,string basename)
{
    //valid test
    if(areElementsUnique(test_faces)==false)
    {
        std::cout<<"The test faces are not unique!"<<std::endl;
	}

    int nv = n_vertices();
    int nf = n_faces();
    int ne = n_edges();

    int phi = std::pow(2, level);
    int re_nv = nv + nf * ((phi - 2) * (phi - 1) / 2) + ne * (phi - 1);
    int n_inner_f = nf * ((phi - 2) * (phi - 1) / 2);

    std::vector<Eigen::Vector3d> points(re_nv);
    std::vector<unsigned int> basefvids(nf * 3);
    std::vector<bool> is_marked_semiv(re_nv, false);

    for (VertexIter v = vertices_begin(); v != vertices_end(); v++) {
        is_marked_semiv[v->idx()] = true;
        for (unsigned int k = 0; k < 3; k++) points[v->idx()](k) = point(*v)[k];//remain points
    }

    for (FaceIter f = faces_begin(); f != faces_end(); f++) {
        unsigned int k = 0;
        for (FaceVertexIter fv = fv_iter(*f); fv.is_valid(); fv++) {
            basefvids[3 * f->idx() + k] = fv->idx();
            k++;
        }
    }
    int num = 0;
    for (FaceIter fit = faces_begin(); fit != faces_end(); fit++)
    {
        //face map vertex remove information
        std::array<unsigned int, 3>fuv;
        
        // get subdivied mesh of vertices from each face of base domain in a row-by-row way
        std::vector<unsigned int> ids;
        std::vector<std::vector<unsigned int> > eids(3);
        HalfedgeHandle he;
        EdgeHandle e;
        for (FaceHalfedgeIter fh = fh_iter(*fit); fh.is_valid(); fh++)
        {
            VertexHandle v1, v2;
            v1 = from_vertex_handle(*fh);
            v2 = to_vertex_handle(*fh);
            if (v1.idx() == (int)basefvids[3 * fit->idx()])
            {
                he = *fh;
                break;
            }
        }
        for (unsigned int k = 0; k < 3; k++) {
            e = edge_handle(he);
            if (he.idx() == halfedge_handle(e, 0).idx()) {
                for (int i = 0; i < phi - 1; i++) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            else {
                for (int i = phi - 2; i >= 0; i--) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            he = next_halfedge_handle(he);
        }
        std::vector<std::vector<unsigned int> > subdivide_ids(phi + 1);
        subdivide_ids[0].push_back(basefvids[3 * fit->idx()]);
        unsigned int pos = 0;
        for (int k = 1; k < phi + 1; k++) {
            if (k == phi) {
                subdivide_ids[k].push_back(basefvids[3 * fit->idx() + 1]);
                for (int i = 0; i < (int)eids[1].size(); i++) {
                    subdivide_ids[k].push_back(eids[1][i]);
                }
                subdivide_ids[k].push_back(basefvids[3 * fit->idx() + 2]);
            }
            else {
                subdivide_ids[k].push_back(eids[0][k - 1]);
                for (int j = 0; j < k - 1; j++) {
                    subdivide_ids[k].push_back(nv + fit->idx() * ((phi - 2) * (phi - 1) / 2) + pos);
                    pos++;
                }
                subdivide_ids[k].push_back(eids[2][phi - 1 - k]);
            }
        }
        // subdivide points location
        unsigned int removeIM_id ;
        vertex_remove_data dataIM ;
        std::array<double, 3> barycoord_values;
        for (unsigned int i = 0; i < subdivide_ids.size(); i++) {
            for (unsigned int j = 0; j < subdivide_ids[i].size(); j++) {
                unsigned int k = subdivide_ids[i][j];
                if (!is_marked_semiv[k]) {
                    Point query(0, 0, 0);
                    barycoord_values = { (phi - i) * 1.0 / phi ,(i - j) * 1.0 / phi ,j * 1.0 / phi };
                    int vsize = 0;
                    for (FaceVertexIter fv = fv_iter(*fit); fv.is_valid(); fv++)
                    {
                        fuv[vsize] = data(*fv).ori_id;
                        vsize++;
                    }
                    while (std::find(all_orifaces.begin(), all_orifaces.end(), fuv) == all_orifaces.end())
                    {
                        bool found_triangle = false;
                        if (Face2IM.find(fuv) == Face2IM.end()) {
                            std::cerr << "Error: Face not found in Face2IM" << std::endl;
                            break;
                        }
                        removeIM_id = Face2IM[fuv];
                        if (removeIM_id >= vremoveIM.size()) {
                            std::cerr << "Error: Invalid removeIM_id" << std::endl;
                            break;
                        }
                        dataIM = vremoveIM[removeIM_id];

                        for (unsigned int s = 0; s < 2; s++)
                        {
                            query[s] = barycoord_values[0] * dataIM.uv_coords[dataIM.g2l[fuv[0]]][s] +
                                barycoord_values[1] * dataIM.uv_coords[dataIM.g2l[fuv[1]]][s] +
                                barycoord_values[2] * dataIM.uv_coords[dataIM.g2l[fuv[2]]][s];
                        }
                        //find cloest triangle
                        double min_distance = 1.0;
                        std::array<unsigned int, 3> closest_fuv = fuv;
                        std::array<double, 3> closest_bary = barycoord_values;

                        for (unsigned int m = 0; m < dataIM.FUV_pre.size(); m++)
                        {
                            std::array<Vec2d, 3> pre_fuv;
                            for (unsigned int n = 0; n < 3; n++)
                            {
                                unsigned int vid = dataIM.FUV_pre[m][n];
                                pre_fuv[n] = Vec2d(dataIM.uv_coords[dataIM.g2l[vid]][0], dataIM.uv_coords[dataIM.g2l[vid]][1]);
                            }

                            auto [alpha, beta, gamma] = cal_bar_2dcoord(Vec2d(query[0], query[1]), pre_fuv);
                            double min_coord = std::min(alpha,std::min( beta, gamma ));
                            double distance = -min_coord;
                            if (distance < min_distance)
                            {
                                min_distance = distance;
                                closest_fuv = dataIM.FUV_pre[m];
                                closest_bary = { alpha, beta, gamma };
                            }
                            
                        }
                        double alpha = std::max(0.0, closest_bary[0]);
                        double beta = std::max(0.0, closest_bary[1]);
                        double gamma = std::max(0.0, closest_bary[2]);
                        double sum = alpha + beta + gamma;
                        barycoord_values[0] = alpha / sum;
                        barycoord_values[1] = beta / sum;
                        barycoord_values[2] = gamma / sum;
                        fuv = closest_fuv;
                        //std::cout << "barycoord_values[0]: " << barycoord_values[0] << "barycoord_values[1]: " << barycoord_values[1] << "barycoord_values[2]: " << barycoord_values[2] << std::endl;
						//std::cout << "sum: " << sum << std::endl;
                        
                        /*std::ofstream out("../data/before_error_fuv.obj");
                        for (unsigned int s = 0; s < dataIM.uv_coords.size(); s++)
                        {
                            out << "v " << dataIM.uv_coords[s][0] << " " << dataIM.uv_coords[s][1] << " 0" << std::endl;
                        }
                        unsigned int nvs = dataIM.uv_coords.size() - 1;
                        for (unsigned int m = 0; m < nvs; m++)
                        {
                            out << "f 1 " << (m + 1) % nvs + 2 << " " << m + 2 << std::endl;
                        }
                        out.close();

                        std::ofstream out1("../data/before_error_query.obj");
                        out1 << "v " << query[0] << " " << query[1] << " 0" << std::endl;
                        out1.close();*/
                    }
                    //std::cout << "face id: " <<fit->idx()<< std::endl;
                    points[k] = Eigen::Vector3d::Zero(0);
					Point p =  barycoord_values[0] * positions[fuv[0]]+
                        barycoord_values[1]*positions[fuv[1]] + 
                        barycoord_values[2] * positions[fuv[2]];
					for (unsigned int s = 0; s < 3; s++) points[k](s) = p[s];
                    is_marked_semiv[k] = true;
                }

            }

        }
    }
    int count = 0;
    for (int i = 0; i < re_nv; i++) {
        if (!is_marked_semiv[i]) {
            count++;
            printf("missing id: %d\n", i);
        }
    }
    if (count > 0) {
        printf("# %d points location is failed\n", count);
        return;
    }
	string sub_path = output_dir + "/subd" + std::to_string(level) + "/" + basename + "_001.obj";
    std::ofstream out(sub_path);
    for (int i = 0; i < re_nv; i++) {
        out << "v " << points[i](0) << " " << points[i](1) << " " << points[i](2) << std::endl;
    }
    for (FaceIter f = faces_begin(); f != faces_end(); f++) {
        // get subdivied mesh of vertices from each face of base domain in a row-by-row way.
        std::vector<unsigned int> ids;
        std::vector<std::vector<unsigned int> > eids(3);
        HalfedgeHandle he;
        EdgeHandle e;
        for (FaceHalfedgeIter fh = fh_iter(*f); fh.is_valid(); fh++) {
            VertexHandle v1, v2;
            v1 = from_vertex_handle(*fh);
            v2 = to_vertex_handle(*fh);
            if (v1.idx() == (int)basefvids[3 * f->idx()]) {
                he = *fh;
                break;
            }
        }
        for (unsigned int k = 0; k < 3; k++) {
            e = edge_handle(he);
            if (he.idx() == halfedge_handle(e, 0).idx()) {
                for (int i = 0; i < phi - 1; i++) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            else {
                for (int i = phi - 2; i >= 0; i--) {
                    eids[k].push_back(nv + n_inner_f + e.idx() * (phi - 1) + i);
                }
            }
            he = next_halfedge_handle(he);
        }
        std::vector<std::vector<unsigned int> > subdivide_ids(phi + 1);
        subdivide_ids[0].push_back(basefvids[3 * f->idx()]);
        unsigned int pos = 0;
        for (int k = 1; k < phi + 1; k++) {
            if (k == phi) {
                subdivide_ids[k].push_back(basefvids[3 * f->idx() + 1]);
                for (int i = 0; i < (int)eids[1].size(); i++) {
                    subdivide_ids[k].push_back(eids[1][i]);
                }
                subdivide_ids[k].push_back(basefvids[3 * f->idx() + 2]);
            }
            else {
                subdivide_ids[k].push_back(eids[0][k - 1]);
                for (int j = 0; j < k - 1; j++) {
                    subdivide_ids[k].push_back(nv + f->idx() * ((phi - 2) * (phi - 1) / 2) + pos);
                    pos++;
                }
                subdivide_ids[k].push_back(eids[2][phi - 1 - k]);
            }
        }
        for (int i = 0; i < phi; i++) {
            for (int j = 0; j < phi - 1 - i; j++) {
                out << "f " << subdivide_ids[phi - i - 1][j] + 1 << " " << subdivide_ids[phi - i][j + 1] + 1 << " " << subdivide_ids[phi - i - 1][j + 1] + 1 << std::endl;
            }
            for (int j = 0; j < phi - i; j++) {
                out << "f " << subdivide_ids[phi - i][j] + 1 << " " << subdivide_ids[phi - i][j + 1] + 1 << " " << subdivide_ids[phi - i - 1][j] + 1 << std::endl;
            }
        }
    }
    out.close();
}
