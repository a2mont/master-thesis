#include "Experiments3D.hh"

void Experiment3D::generate_torsion_mesh(double torsion_turns_count){

    TetrahedralMesh mesh;
    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::generate_mesh(10, 10, 10, mesh);

    double target_angle = torsion_turns_count * 2 * M_PI;

    Eigen::Vector3d min_pos(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Eigen::Vector3d max_pos(0.0, 0.0, 0.0);

    for(auto v: mesh_.vertices()){
        auto init_pos = mesh_.vertex(v);

        for(int i(0); i<3; i++){
            min_pos[i] = std::min(min_pos[i], init_pos[i]);
            max_pos[i] = std::max(max_pos[i], init_pos[i]);
        }
    }

    auto min_z = min_pos[2];
    auto max_z = max_pos[2];

    std::cout<<" min/max z = "<<min_z<<" / "<<max_z<<std::endl;


    auto center = 0.5 * (max_pos + min_pos);
    std::cout<<" - center at "<<center.transpose()<<std::endl;

    for(auto v: mesh.vertices()){
        mesh.set_vertex(v, mesh.vertex(v) - TetrahedralMesh::PointT(center[0], center[1], 0.0));
    }

    for(auto v: mesh.vertices()){
        auto init_pos = mesh.vertex(v);
        Eigen::Vector2f xy;
        xy <<init_pos[0], init_pos[1];
        auto z = init_pos[2];

        auto z_prime = (z - min_z) / (max_z - min_z);

        auto theta = z_prime * target_angle;

        Eigen::Matrix2f rot;
        rot << cos(theta), -sin(theta),
               sin(theta), cos(theta);

        auto target_xy = rot * xy;

        mesh.set_vertex(v, {target_xy[0], target_xy[1], z});

    }
    mesh_ = mesh;
}
