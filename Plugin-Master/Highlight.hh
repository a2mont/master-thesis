#ifndef OPENFLIPPER_HIGHLIGHT_HH
#define OPENFLIPPER_HIGHLIGHT_HH


#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "QualityEvaluation.hh"
#include "MainLoop.hh"

class Highlight
{
public:
    Highlight();

public :
    static void highlight_triangles(TriMesh& _mesh, const bool _verbose = false);
};

#endif // OPENFLIPPER_HIGHLIGHT_HH
