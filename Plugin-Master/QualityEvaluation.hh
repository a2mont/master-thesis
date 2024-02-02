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
    // 3D
    static double evaluate(const OpenVolumeMesh::CellHandle _cell, TetrahedralMesh& _mesh, const bool _verbose = false);
    static double evaluate(
            OpenVolumeMesh::VertexHandle _p0,
            OpenVolumeMesh::VertexHandle _p1,
            OpenVolumeMesh::VertexHandle _p2,
            OpenVolumeMesh::VertexHandle _p3,
            TetrahedralMesh &_mesh,
            const bool _verbose = false);
    static void scaleMesh(TetrahedralMesh &_mesh);
    // 2D
    static double evaluate(const OpenMesh::SmartFaceHandle _face, TriMesh& _mesh, const bool _verbose = false);
    // 2d & 3D
    static double getMaxQuality();

private:
    // 3D
    static double symmetric_dirichlet_energy(const TetrahedralMesh &mesh, const std::vector<OpenVolumeMesh::VertexHandle> &_c_verts);
    static double calculate_volume(TetrahedralMesh& _mesh, OpenVolumeMesh::CellHandle _cellHandle);
    static double calculate_volume(TetrahedralMesh& _mesh,
                                   OpenVolumeMesh::VertexHandle _v0,
                                   OpenVolumeMesh::VertexHandle _v1,
                                   OpenVolumeMesh::VertexHandle _v2,
                                   OpenVolumeMesh::VertexHandle _v3);
    static double calculate_edge_length(TetrahedralMesh& _mesh, OpenVolumeMesh::EdgeHandle _edge);
    static double calculate_edge_length(TetrahedralMesh &_mesh,
                                        OpenVolumeMesh::VertexHandle _v0,
                                        OpenVolumeMesh::VertexHandle _v1);
    static double computeQuality(TetrahedralMesh &_mesh,
                                 OpenVolumeMesh::VertexHandle _v0,
                                 OpenVolumeMesh::VertexHandle _v1,
                                 OpenVolumeMesh::VertexHandle _v2,
                                 OpenVolumeMesh::VertexHandle _v3,
                                 std::vector<double> _edge_lengths,
                                 bool _verbose = false);
    // 2D
    static double calculate_area(TriMesh& _mesh, const OpenMesh::SmartFaceHandle _face);
    // 2D & 3D
    static double calculate_l_harm(std::vector<double> _edge_lengths);
    static double calculate_l_rms(std::vector<double> _edge_lengths);

};


#endif // OPENFLIPPER_QUALITYEVALUATION_HH
