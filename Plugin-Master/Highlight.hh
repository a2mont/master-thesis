#ifndef OPENFLIPPER_HIGHLIGHT_HH
#define OPENFLIPPER_HIGHLIGHT_HH


#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "QualityEvaluation.hh"

class Highlight
{
public:
    Highlight();

public :
    static void highlight_worst_triangles(TriMesh& _mesh, double _treshold);

};

#endif // OPENFLIPPER_HIGHLIGHT_HH
