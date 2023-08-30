#ifndef OPENFLIPPER_QUALITYEVALUATION_HH
#define OPENFLIPPER_QUALITYEVALUATION_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <queue>
#include "eigen3/Eigen/Dense"

class QualityEvaluation
{
private:
    QualityEvaluation();

public:
    static double evaluate(const OpenMesh::SmartFaceHandle _face, TriMesh& _mesh, const bool _verbose = false);
    static double evaluate(const OpenVolumeMesh::CellHandle _cell, TetrahedralMesh& _mesh, const bool _verbose = false);
private:
    static double calculate_area(TriMesh& _mesh, const OpenMesh::SmartFaceHandle _face);
    static double calculate_volume(TetrahedralMesh& _mesh, OpenVolumeMesh::CellHandle _cellHandle);
    static double calculate_l_harm(std::vector<double> _edge_lengths);
    static double calculate_l_rms(std::vector<double> _edge_lengths);

    static double calculate_edge_length(TetrahedralMesh& _mesh, OpenVolumeMesh::EdgeHandle _edge);
};


#endif // OPENFLIPPER_QUALITYEVALUATION_HH
