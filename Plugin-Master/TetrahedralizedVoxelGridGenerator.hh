#pragma once


template <class TET_MESH_>
class TetrahedralizedVoxelGridGenerator
{
public:

    static TET_MESH_ generate_mesh(int n,
                                   TET_MESH_& mesh);

    static TET_MESH_ generate_mesh(int width,
                                   int height,
                                   int depth,
                                   TET_MESH_& mesh);


private:
    TetrahedralizedVoxelGridGenerator(int width,
                                      int height,
                                      int depth,
                                      TET_MESH_& mesh);

    void generate_five_tet_voxel();

    void generate_vertices();

    /* (x,y,z) = voxel coordinates */
    void generate_cube(int x, int y, int z,
                       bool orientation);

    void dig_hole(int start_x, int start_y, int start_z,
                  int end_x,   int end_y,   int end_z);


    bool is_digged(int x, int y, int z);

    int coordinates_to_corner_vertex_idx(int x, int y, int z);

    int coordinates_to_voxel_idx(int x, int y, int z);

    TET_MESH_& mesh_;
    const int width_;
    const int height_;
    const int depth_;
};

#if !defined(TET_VOXEL_GRID_GEN_C)
#define TET_VOXEL_GRID_GEN_TEMPLATES
#include "TetrahedralizedVoxelGridGenerator.cc"
#endif
