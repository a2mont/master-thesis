#ifndef OPENFLIPPER_MAINLOOP_HH
#define OPENFLIPPER_MAINLOOP_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
#include "Smoothing.hh"
#include "eigen3/Eigen/Dense"

class MainLoop
{

public:
    MainLoop(TriMesh& _mesh, double _q_min) :
        mesh_(_mesh),
        q_min_(_q_min)
    {
        if(!mesh_.get_property_handle(face_visited_, "face visited"))
            mesh_.add_property(face_visited_, "face visited");
        for(auto fh: mesh_.faces()){
            mesh_.property(face_visited_, fh) = false;
        }
        if(!mesh_.get_property_handle(cavity_edge_, "cavity edge"))
            mesh_.add_property(cavity_edge_, "cavity edge");
        for(auto vh: mesh_.vertices()){
            mesh_.property(cavity_edge_, vh) = false;
        }
    }
    ~MainLoop(){
        mesh_.remove_property(face_visited_);
        mesh_.remove_property(cavity_edge_);
    }

    struct Triangle{
        OpenMesh::SmartFaceHandle face_handle_;
        double quality_;
        Triangle(OpenMesh::SmartFaceHandle _face_handle, double _quality): face_handle_(_face_handle), quality_(_quality){}
        std::string toString() const {return "Face id: " + std::to_string(face_handle_.idx()) + " Quality: " + std::to_string(quality_);}
    };
    struct CompareQuality{
        bool operator()(Triangle const& t1,Triangle const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    using PriorityQueue = std::priority_queue<Triangle, std::vector<Triangle>, CompareQuality>;
    using Point = ACG::Vec3d;


public:
    void loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose= false, int _max_iter=1);
private:
    void reset_queue(PriorityQueue& _queue);

    bool topologial_pass(PriorityQueue* _A);
    void edge_contraction_pass(PriorityQueue* _A);
    void insertion_pass(PriorityQueue* _A);
    void smoothing_pass(PriorityQueue* _A);

    void improve_mesh(PriorityQueue _badTriangles);
    void improve_triangle(Triangle _t);

    void find_faces_with_p(std::vector<OpenMesh::SmartFaceHandle> &_list, OpenMesh::SmartFaceHandle _fh, const Point _p);
    bool contains_p(OpenMesh::SmartFaceHandle _fh, const Point _p);

private:
    TriMesh& mesh_;

    // used in insertion pass

    // true means the circumcircle contains the new vertex p
    OpenMesh::FPropHandleT<bool> face_visited_;
    // true means the vertex is at the border of the cavity
    OpenMesh::VPropHandleT<bool> cavity_edge_;

    PriorityQueue quality_queue_;


    const double q_min_;
};

#endif // OPENFLIPPER_MAINLOOP_HH
