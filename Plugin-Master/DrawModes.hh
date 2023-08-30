#pragma once
//
// Created by Martin Heistermann on 14.05.18.
//
// Access drawmodes by method, so typos are impossible.
//
// This should go into OpenFlipper sometime
//

#include "OpenFlipper/common/GlobalOptions.hh"
#include "ACG/Scenegraph/DrawModes.hh"

enum class BlendMode {
      Opaque,
      Transparent,
  };
using DrawMode = ACG::SceneGraph::DrawModes::DrawMode;

struct DrawModes {
    static DrawMode Faces_flat()       { return getDrawMode("Faces (flat shaded)"); }
    static DrawMode Faces_colored()    { return getDrawMode("Faces (colored per face)"); }
    static DrawMode Faces_colored_flat() { return getDrawMode("Faces (colored per face, flat shaded)"); }
    static DrawMode Halffaces_flat()   { return getDrawMode("Halffaces (flat shaded)"); }
    static DrawMode Edges_colored()    { return getDrawMode("Edges (colored per edge)"); }
    static DrawMode Edges_wireframe()  { return getDrawMode("Edges (wireframe)"); }
    static DrawMode Edges_hiddenline() { return getDrawMode("Edges (hidden line)"); }
    static DrawMode Vertices_colored() { return getDrawMode("Vertices (colored)"); }
    static DrawMode Vertices()         { return getDrawMode("Vertices"); }
    static DrawMode Wireframe()        { return getDrawMode("Wireframe"); }
    static DrawMode Solid_perVertex()  { return getDrawMode("Solid (colored per-vertex)"); }
    static DrawMode Cells_transparent(){ return getDrawMode("Cells (transparent)"); }
    static DrawMode Cells_flat(){ return getDrawMode("Cells (flat shaded)"); }
private:
    static DrawMode getDrawMode (const std::string &name) {
        if (!OpenFlipper::Options::gui()) {
            return DrawMode();
        }
        auto mode = ACG::SceneGraph::DrawModes::getDrawMode(name);
        return mode;
    }
};
