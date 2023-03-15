#include "MainLoop.hh"

void MainLoop::loop(ACG::Vec3d _displacement, int _constraint_vh, bool _verbose,int _max_iter ){
    for(int i=0; i < _max_iter; i++){
        vd_.displace(_displacement, _constraint_vh, _verbose);
        for(auto fh: mesh_.faces()){
            qe_.evaluate(fh, _verbose);
        }
    }

    std::cout << "Worst Triangle\t" << qe_.get_face_quality_queue().top().toString() << std::endl;
    qe_.reset_queue();

}
