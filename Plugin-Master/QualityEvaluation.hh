#ifndef OPENFLIPPER_QUALITYEVALUATION_HH
#define OPENFLIPPER_QUALITYEVALUATION_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <queue>

class QualityEvaluation
{
private:
    QualityEvaluation();

public:
    static double evaluate(const OpenMesh::SmartFaceHandle _face, TriMesh& _mesh, const bool _verbose = false);
private:
    static double calculate_area(TriMesh& _mesh, const OpenMesh::SmartFaceHandle _face);
    static double calculate_l_harm(std::vector<double> _edge_lengths);
    static double calculate_l_rms(std::vector<double> _edge_lengths);

};


#endif // OPENFLIPPER_QUALITYEVALUATION_HH
