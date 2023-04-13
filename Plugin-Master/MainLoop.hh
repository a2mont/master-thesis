#ifndef OPENFLIPPER_MAINLOOP_HH
#define OPENFLIPPER_MAINLOOP_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
#include "Smoothing.hh"

class MainLoop
{

public:
    MainLoop(TriMesh& _mesh, double _q_min) :
        mesh_(_mesh),
        vd_(_mesh),
        smoother_(_mesh),
        q_min_(_q_min){}
    ~MainLoop(){}

    struct Triangle{
        int face_id_;
        double quality_;
        Triangle(int _face_id, double _quality): face_id_(_face_id), quality_(_quality){}
        std::string toString() const {return "Face id: " + std::to_string(face_id_) + " Quality: " + std::to_string(quality_);}
    };
    struct CompareQuality{
        bool operator()(Triangle const& t1,Triangle const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    using PriorityQueue = std::priority_queue<Triangle, std::vector<Triangle>, CompareQuality>;


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

private:
    TriMesh& mesh_;
    VertexDisplacement vd_;
    Smoothing smoother_;

    PriorityQueue quality_queue_;


    const double q_min_;
};

#endif // OPENFLIPPER_MAINLOOP_HH
