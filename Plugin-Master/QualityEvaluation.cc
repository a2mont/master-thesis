#include "QualityEvaluation.hh"

// ------------------------------------ 3D ----------------------------------------------


double QualityEvaluation::evaluate(const OpenVolumeMesh::CellHandle _cell, TetrahedralMesh& _mesh, const bool _verbose){
    double quality = 0.;

    std::vector<double> edge_lengths;
    std::vector<OpenVolumeMesh::VertexHandle> vertices;

    for(auto ce_it = _mesh.ce_iter(_cell); ce_it.valid(); ++ce_it){
        double edgeLength = calculate_edge_length(_mesh, *ce_it);
        edge_lengths.push_back(edgeLength);
    }
    for(auto vh: _mesh.get_cell_vertices(_cell)){
        vertices.push_back(vh);
    }

    quality = computeQuality(_mesh, vertices[0], vertices[1], vertices[2], vertices[3], edge_lengths,_verbose);

    return quality;
}

double QualityEvaluation::evaluate(OpenVolumeMesh::VertexHandle _v0,
                                   OpenVolumeMesh::VertexHandle _v1,
                                   OpenVolumeMesh::VertexHandle _v2,
                                   OpenVolumeMesh::VertexHandle _v3,
                                   TetrahedralMesh& _mesh,
                                   const bool _verbose){
    double quality = 0.;
    double edgeLength = 0.;
    std::vector<double> edge_lengths;

    edgeLength = calculate_edge_length(_mesh, _v0,_v1);
    edge_lengths.push_back(edgeLength);
    edgeLength = calculate_edge_length(_mesh, _v0,_v2);
    edge_lengths.push_back(edgeLength);
    edgeLength = calculate_edge_length(_mesh, _v0,_v3);
    edge_lengths.push_back(edgeLength);
    edgeLength = calculate_edge_length(_mesh, _v1,_v2);
    edge_lengths.push_back(edgeLength);
    edgeLength = calculate_edge_length(_mesh, _v1,_v3);
    edge_lengths.push_back(edgeLength);
    edgeLength = calculate_edge_length(_mesh, _v2,_v3);
    edge_lengths.push_back(edgeLength);

    quality = computeQuality(_mesh, _v0, _v1, _v2, _v3, edge_lengths, _verbose);

    return quality;
}

double QualityEvaluation::computeQuality(TetrahedralMesh& _mesh,
                                         OpenVolumeMesh::VertexHandle _v0,
                                         OpenVolumeMesh::VertexHandle _v1,
                                         OpenVolumeMesh::VertexHandle _v2,
                                         OpenVolumeMesh::VertexHandle _v3,
                                         std::vector<double> _edge_lengths,
                                         bool _verbose){
    double quality = 0.;
    double volume = 0.;
    double l_harm = 0.;
    double l_rms = 1.;

    volume = calculate_volume(_mesh, _v0, _v1, _v2, _v3);

    l_harm = calculate_l_harm(_edge_lengths);
    l_rms = calculate_l_rms(_edge_lengths);

    if(!std::isnan(volume) && !std::isnan(l_harm) && !std::isnan(l_rms)){
      quality = 6 * sqrt(2) * volume * l_harm / pow(l_rms, 4);
    }

    if(_verbose){
        std::cout <<
                     "Tet with vertices: " << _v0 << ", " << _v1 << ", " << _v2 << ", " << _v3 <<
                     "\nEdge lengths: ";
        for(auto e: _edge_lengths){
                     std::cout << e << ", ";
        }
        std::cout <<
                     "\nVolume: "<< volume <<
                     "\nHarmonic length: " << l_harm <<
                     "\nRoot-mean-squared edge length: " << l_harm <<
                     "\nQuality: " << quality <<
                     "\n------------------------------"
        << std::endl;
    }

    return quality;
}

double QualityEvaluation::calculate_edge_length(TetrahedralMesh& _mesh, OpenVolumeMesh::EdgeHandle _edge){
    double length = 0;
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    for(auto vh: _mesh.edge_vertices(_edge)){
        vertices.push_back(vh);
    }
    length = calculate_edge_length(_mesh, vertices[0], vertices[1]);
    return length;
}

double QualityEvaluation::calculate_edge_length(TetrahedralMesh& _mesh,
                                                OpenVolumeMesh::VertexHandle _v0,
                                                OpenVolumeMesh::VertexHandle _v1){
    double length = 0;
    auto vector = _mesh.vertex(_v1) - _mesh.vertex(_v0);

    length = sqrt(pow(vector[0],2)+pow(vector[1],2)+pow(vector[2],2));
//    std::cout << "Edge length: " <<length <<", Vector: " << vector <<std::endl;
    return length;
}

double QualityEvaluation::calculate_volume(TetrahedralMesh& _mesh, OpenVolumeMesh::CellHandle _cellHandle){
    double volume = 0;
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    for(auto vh : _mesh.get_cell_vertices(_cellHandle)){
        vertices.push_back(vh);
    }
    volume = calculate_volume(_mesh, vertices[0], vertices[1], vertices[2], vertices[3]);

    return volume;
}

double QualityEvaluation::calculate_volume(TetrahedralMesh& _mesh,
                                           OpenVolumeMesh::VertexHandle _v0,
                                           OpenVolumeMesh::VertexHandle _v1,
                                           OpenVolumeMesh::VertexHandle _v2,
                                           OpenVolumeMesh::VertexHandle _v3){
    double volume = 0;
    ACG::Vec3d e0, e1, e2;
    e0 =_mesh.vertex(_v1)-_mesh.vertex(_v0);
    e1 = _mesh.vertex(_v2)-_mesh.vertex(_v0);
    e2 = _mesh.vertex(_v3)-_mesh.vertex(_v0);
    Eigen::Vector3d AB, AC, AD;
    AB = Eigen::Vector3d(e0[0], e0[1], e0[2]);
    AC = Eigen::Vector3d(e1[0], e1[1], e1[2]);
    AD = Eigen::Vector3d(e2[0], e2[1], e2[2]);

    Eigen::Matrix3d matrix;
    matrix.col(0) = AB;
    matrix.col(1) = AC;
    matrix.col(2) = AD;


    volume = matrix.determinant()/6;

//    if(volume < 10e-7){
//        std::cout << "\033[1;31m--------- Inverted tet with vertices: \033[0m "
//                  <<  _v0 << ", " << _v1 << ", " << _v2 << ", " << _v3
//                  << ", volume= "
//                  << volume << "\033[1;31m ---------\033[0m"
//                  << std::endl;
//    }

    return volume;
}

// ------------------------------------ 2D ----------------------------------------------

double QualityEvaluation::evaluate(const OpenMesh::SmartFaceHandle _face, TriMesh& _mesh, const bool _verbose){
    double quality = 0.;
    double area = 0.;

    std::vector<double> edge_lengths;

    for(auto feh : _face.edges()){
        edge_lengths.push_back(_mesh.calc_edge_length(feh));
    }

    area = calculate_area(_mesh, _face);
    double squaredEdgeLength = 0.;
    for(auto edge: edge_lengths){
        squaredEdgeLength += pow(edge,2);
    }

    if(!std::isnan(area) && squaredEdgeLength > 0.){
        quality = 4 * sqrt(3) * area / squaredEdgeLength;
    }

    if(_verbose)
        std::cout <<
                     "Triangle " << _face.idx() <<
                     "\nEdge lengths: " << edge_lengths[0] << " " << edge_lengths[1] << " " << edge_lengths[2] <<
                     "\nArea: "<< area <<
                     "\nTotal squared edges length: " << squaredEdgeLength <<
                     "\nQuality: " << quality <<
                     "\n------------------------------"
        << std::endl;


    return quality;
}

double QualityEvaluation::calculate_area(TriMesh& _mesh, OpenMesh::SmartFaceHandle _faceHandle){
    double area = 0.;
    std::vector<ACG::Vec3d> points;
    ACG::Vec3d v,w;
    for(auto vh: _faceHandle.vertices()){
        points.emplace_back(_mesh.point(vh));
    }
    v = points[1] - points[0];
    w = points[2] - points[0];
    ACG::Vec3d xproduct = w.cross(v);
    area = xproduct[2]/2;
    return area;
}

// ------------------------------------ 2D&3D -------------------------------------------

double QualityEvaluation::calculate_l_harm(std::vector<double> _edge_lengths){
    //https://en.wikipedia.org/wiki/Harmonic_mean
    double l_harm = 0.;
    double div = 0;
    for(auto edge: _edge_lengths){
        div += 1/edge;
    }
    l_harm = _edge_lengths.size()/div;

    return l_harm;
}

double QualityEvaluation::calculate_l_rms(std::vector<double> _edge_lengths){
    double l_rms = 1.;
    double mean_squared = 0;
    for(auto edge: _edge_lengths){
        mean_squared += pow(edge,2);
    }
    l_rms = sqrt(mean_squared / _edge_lengths.size());

    return l_rms;
}

// Can be modified if metric changes
double QualityEvaluation::getMaxQuality(){
    return 1.;
}
