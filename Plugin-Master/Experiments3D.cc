#include "Experiments3D.hh"

void Experiment3D::generate_torsion_mesh(double torsion_turns_count, bool _withRemesh){

//    TetrahedralMesh mesh;
//    TetrahedralizedVoxelGridGenerator<TetrahedralMesh>::generate_mesh(10, 10, 10, mesh);
    double target_angle = torsion_turns_count * 2 * M_PI;
    double boundaryOnly(false);

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

//    std::cout<<" min/max z = "<<min_z<<" / "<<max_z<<std::endl;


    auto center = 0.5 * (max_pos + min_pos);
//    std::cout<<" - center at "<<center.transpose()<<std::endl;

    for(auto v: mesh_.vertices()){
        mesh_.set_vertex(v, mesh_.vertex(v) - TetrahedralMesh::PointT(center[0], center[1], 0.0));
    }

    for(auto v: mesh_.vertices()){
        if(boundaryOnly){
            if(!mesh_.is_boundary(v))continue;
        }
        auto init_pos = mesh_.vertex(v);
        Eigen::Vector2f xy;
        xy <<init_pos[0], init_pos[1];
        auto z = init_pos[2];

        auto z_prime = (z - min_z) / (max_z - min_z);

        auto theta = z_prime * target_angle;

        Eigen::Matrix2f rot;
        rot << cos(theta), -sin(theta),
               sin(theta), cos(theta);

        auto target_xy = rot * xy;

        mesh_.set_vertex(v, {target_xy[0], target_xy[1], z});

    }
    if(_withRemesh){
        loop_.loop();
    }
}

void Experiment3D::generate_stretch_mesh(double stretch_factor,
                                         bool _withRemesh){
    bool printDebug(false);
    Eigen::Vector3d min_pos(std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max());
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

    if(printDebug){
        std::cout<<" min/max z = "<<min_z<<" / "<<max_z<<std::endl;
    }


    auto center = 0.5 * (max_pos + min_pos);
    if(printDebug){
        std::cout<<" - center at "<<center.transpose()<<std::endl;
    }

    double d = stretch_factor;
    double A = 1. * (1.- 1./(1 * d));

    for(auto v: mesh_.vertices()){
        ACG::Vec3d new_pos = mesh_.vertex(v) - TetrahedralMesh::PointT(center[0], center[1], 0.0);
        mesh_.set_vertex(v, new_pos);
    }
    if(printDebug){
        std::cout<<"amplitude: "<<A<<std::endl;
    }
        for(auto v: mesh_.vertices()){
            auto init_pos = mesh_.vertex(v);
            Eigen::Vector2f xy;
            xy <<init_pos[0], init_pos[1];
            auto z = init_pos[2];

            double s = 1 + (4.0 * A/(min_z * min_z)) * (z - min_z) * (z - max_z);
            if(printDebug){
                std::cout<<" - z = "<<z<<", s = "<<s<<std::endl;
            }

            auto target_xy = s * xy;

            mesh_.set_vertex(v, {target_xy[0], target_xy[1], z * d});

        }
    if(_withRemesh){
        loop_.loop();
    }

}






void Experiment3D::generate_compress_mesh(double stretch_factor,
                                         bool _withRemesh){
    bool printDebug(false);
    Eigen::Vector3d min_pos(std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max());
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

    if(printDebug){
        std::cout<<" min/max z = "<<min_z<<" / "<<max_z<<std::endl;
    }


    auto center = 0.5 * (max_pos + min_pos);
    if(printDebug){
        std::cout<<" - center at "<<center.transpose()<<std::endl;
    }

    double d = -stretch_factor;
    double A = 1. * (1.- 1./(1 * d));

    for(auto v: mesh_.vertices()){
        ACG::Vec3d new_pos = mesh_.vertex(v) - TetrahedralMesh::PointT(center[0], center[1], 0.0);
        mesh_.set_vertex(v, new_pos);
    }
    if(printDebug){
        std::cout<<"amplitude: "<<A<<std::endl;
    }
        for(auto v: mesh_.vertices()){
            auto init_pos = mesh_.vertex(v);
            Eigen::Vector2f xy;
            xy <<init_pos[0], init_pos[1];
            auto z = init_pos[2];

            double s = 1 + (4.0 * A/(min_z * min_z)) * (z - min_z) * (z - max_z);
            if(printDebug){
                std::cout<<" - z = "<<z<<", s = "<<s<<std::endl;
            }

            auto target_xy = s * xy;

            mesh_.set_vertex(v, {target_xy[0], target_xy[1], z * d});

        }
    if(_withRemesh){
//        loop_.loop();
    }

}
