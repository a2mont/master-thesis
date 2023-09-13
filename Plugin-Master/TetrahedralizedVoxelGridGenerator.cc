#define TET_VOXEL_GRID_GEN_C


#include "TetrahedralizedVoxelGridGenerator.hh"
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <iostream>

using namespace OpenVolumeMesh;


template<class TET_MESH_>
TET_MESH_ TetrahedralizedVoxelGridGenerator<TET_MESH_>::generate_mesh(int n, TET_MESH_& mesh){
    return generate_mesh(n,n,n, mesh);
}

template<class TET_MESH_>
TET_MESH_ TetrahedralizedVoxelGridGenerator<TET_MESH_>::generate_mesh(int width,
                                                                      int height,
                                                                      int depth,
                                                                      TET_MESH_& mesh){

    TetrahedralizedVoxelGridGenerator<TET_MESH_> generator(width, height, depth, mesh);



    std::cout<<" generating block..."<<std::endl;
    generator.generate_five_tet_voxel();

    std::cout<<" generated voxel grid."<<std::endl;
    std::cout<<" -> vertices count = "<<generator.mesh_.n_vertices()<<std::endl;
    std::cout<<" ->    edges count = "<<generator.mesh_.n_edges()<<std::endl;
    std::cout<<" ->    faces count = "<<generator.mesh_.n_faces()<<std::endl;
    std::cout<<" ->    cells count = "<<generator.mesh_.n_cells()<<std::endl;

    generator.mesh_.collect_garbage();

    return generator.mesh_;
}



template<class TET_MESH_>
TetrahedralizedVoxelGridGenerator<TET_MESH_>::TetrahedralizedVoxelGridGenerator(int width,
                                                                                int height,
                                                                                int depth,
                                                                                TET_MESH_& mesh)
    : mesh_(mesh), width_(width), height_(height), depth_(depth){}




template<class TET_MESH_>
void TetrahedralizedVoxelGridGenerator<TET_MESH_>::generate_five_tet_voxel(){

    std::cout<<" - generating vertices..."<<std::endl;
    generate_vertices();

    std::cout<<" - generating cubes..."<<std::endl;
    for(auto x(0); x < width_; x++){
        for(auto y(0); y < height_; y++){
            for(auto z(0); z < depth_; z++){
                generate_cube(x,y,z, (x+y+z) % 2);
            }
        }
    }
}


template<class TET_MESH_>
void TetrahedralizedVoxelGridGenerator<TET_MESH_>::generate_vertices(){

    for(auto z(0); z <= depth_; z++){
        for(auto y(0); y <= height_; y++){
            for(auto x(0); x <= width_; x++){
                auto v = mesh_.add_vertex({x,y,-z});
                //std::cout<<" -- added vertex "<<v<<" at "<<mesh_.vertex(v)<<std::endl;
            }
        }
    }
}


template<class TET_MESH_>
void TetrahedralizedVoxelGridGenerator<TET_MESH_>::generate_cube(int x, int y, int z,
                                                      bool orientation){

    //std::cout<<" -- generating cube at "<<x<<", "<<y<<", "<<z<<" with orientation "<<orientation<<std::endl;

    const int i = coordinates_to_corner_vertex_idx(x,y,z);
    //std::cout<<" -- corner vertex idx = "<<i<<std::endl;

    const int w(width_ + 1), h(height_ + 1), wh(w * h);

    if(orientation){
        mesh_.add_cell(VertexHandle(     w     + i),
                       VertexHandle(wh + w + 1 + i),
                       VertexHandle(wh         + i),
                       VertexHandle(wh + w     + i));

        mesh_.add_cell(VertexHandle(     w     + i),
                       VertexHandle(         1 + i),
                       VertexHandle(wh + w + 1 + i),
                       VertexHandle(     w + 1 + i));

        mesh_.add_cell(VertexHandle(             i),
                       VertexHandle(         1 + i),
                       VertexHandle(wh         + i),
                       VertexHandle(     w     + i));

        mesh_.add_cell(VertexHandle(     w     + i),
                       VertexHandle(         1 + i),
                       VertexHandle(wh         + i),
                       VertexHandle(wh + w + 1 + i));

        mesh_.add_cell(VertexHandle(         1 + i),
                       VertexHandle(wh     + 1 + i),
                       VertexHandle(wh         + i),
                       VertexHandle(wh + w + 1 + i));
    }else{

        mesh_.add_cell(VertexHandle(           + i),
                       VertexHandle(     w + 1 + i),
                       VertexHandle(wh + w     + i),
                       VertexHandle(     w     + i));

        mesh_.add_cell(VertexHandle(             i),
                       VertexHandle(         1 + i),
                       VertexHandle(wh     + 1 + i),
                       VertexHandle(     w + 1 + i));

        mesh_.add_cell(VertexHandle(             i),
                       VertexHandle(wh     + 1 + i),
                       VertexHandle(wh         + i),
                       VertexHandle(wh + w     + i));

        mesh_.add_cell(VertexHandle(     w + 1 + i),
                       VertexHandle(wh     + 1 + i),
                       VertexHandle(wh + w     + i),
                       VertexHandle(wh + w + 1 + i));

        mesh_.add_cell(VertexHandle(             i),
                       VertexHandle(     w + 1 + i),
                       VertexHandle(wh     + 1 + i),
                       VertexHandle(wh + w     + i));
    }

}

template<class TET_MESH_>
int TetrahedralizedVoxelGridGenerator<TET_MESH_>::coordinates_to_corner_vertex_idx(int x, int y, int z){

    return x + y * (width_ + 1) + z * (width_ + 1) * (height_ + 1);
}


template<class TET_MESH_>
int TetrahedralizedVoxelGridGenerator<TET_MESH_>::coordinates_to_voxel_idx(int x, int y, int z){
    return x + y * width_ + z * width_ * height_;
}



