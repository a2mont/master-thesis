#ifndef OPENFLIPPER_MAINLOOP_HH
#define OPENFLIPPER_MAINLOOP_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
class MainLoop
{

public:
    MainLoop(TriMesh &_mesh) :
        mesh_(_mesh),
        vd_(_mesh),
        qe_(_mesh){}
    ~MainLoop(){}


public:
    void loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose= false, int _max_iter=1);
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

    const double q_min = 0.15;
};

#endif // OPENFLIPPER_MAINLOOP_HH
