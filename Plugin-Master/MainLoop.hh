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
        qe_(_mesh),
        smoother_(_mesh),
        q_min_(_q_min){}
    ~MainLoop(){}


public:
    void loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose= false, int _max_iter=1);
    QualityEvaluation& get_quality_evaluation() {return qe_;}
private:
    bool topologial_pass(QualityEvaluation::PriorityQueue* _A);
    void edge_contraction_pass(QualityEvaluation::PriorityQueue* _A);
    void insertion_pass(QualityEvaluation::PriorityQueue* _A);
    void smoothing_pass(QualityEvaluation::PriorityQueue* _A);

    void improve_mesh(QualityEvaluation::PriorityQueue _badTriangles);
    void improve_triangle(QualityEvaluation::Triangle _t);

private:
    TriMesh& mesh_;
    VertexDisplacement vd_;
    QualityEvaluation qe_;
    Smoothing smoother_;


    const double q_min_;
};

#endif // OPENFLIPPER_MAINLOOP_HH
