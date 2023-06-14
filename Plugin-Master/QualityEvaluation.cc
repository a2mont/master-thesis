#include "QualityEvaluation.hh"

double QualityEvaluation::evaluate(const OpenMesh::SmartFaceHandle _face, TriMesh& _mesh, const bool _verbose){
    double quality = 0.;
    double area = 0.;
    double l_harm = 0.;
    double l_rms = 1.;

    std::vector<double> edge_lengths;

    for(auto feh : _face.edges()){
        edge_lengths.push_back(_mesh.calc_edge_length(feh));
    }

    area = calculate_area(edge_lengths);
    double squaredEdgeLength = 0.;
    for(auto edge: edge_lengths){
        squaredEdgeLength += pow(edge,2);
    }
//    l_harm = calculate_l_harm(edge_lengths); used for tets
//    l_rms = calculate_l_rms(edge_lengths); used for tets

    if(!std::isnan(area) && squaredEdgeLength > 0.){ //&& !std::isnan(l_harm) && !std::isnan(l_rms))
//      quality = 6 * sqrt(2) * area * l_harm / pow(l_rms, 4); used for tets
        quality = 4 * sqrt(3) * area / squaredEdgeLength;
    }

    if(_verbose)
        std::cout <<
                     "Triangle " << _face.idx() <<
                     "\nEdge lengths: " << edge_lengths[0] << " " << edge_lengths[1] << " " << edge_lengths[2] <<
                     "\nArea: "<< area <<
//                     "\nHarmonic length: " << l_harm <<
//                     "\nRoot-mean-squared edge length: " << l_harm <<
                     "\nTotal squared edges length: " << squaredEdgeLength <<
                     "\nQuality: " << quality <<
                     "\n------------------------------"
        << std::endl;


    return quality;
}


double QualityEvaluation::calculate_area(std::vector<double> _edge_lengths){
    double area = 0.;
    double x1,x2,x3,semi_perimeter,total;
    x1 = _edge_lengths[0];
    x2 = _edge_lengths[1];
    x3 = _edge_lengths[2];
    semi_perimeter = (x1+x2+x3)/2;

    total = semi_perimeter*(semi_perimeter-x1)*(semi_perimeter-x2)*(semi_perimeter-x3);
    total = total < 0. ? 0. : total;

    area = sqrt(total);
    return area;
}

double QualityEvaluation::calculate_l_harm(std::vector<double> _edge_lengths){
    //https://en.wikipedia.org/wiki/Harmonic_mean
    double l_harm = 0.;
    double x1,x2,x3;
    x1 = _edge_lengths[0];
    x2 = _edge_lengths[1];
    x3 = _edge_lengths[2];

    l_harm = 3*x1*x2*x3 / (x1*x2 + x1*x3 + x2*x3);
    return l_harm;
}

double QualityEvaluation::calculate_l_rms(std::vector<double> _edge_lengths){
    double l_rms = 1.;
    double x1,x2,x3;
    x1 = _edge_lengths[0];
    x2 = _edge_lengths[1];
    x3 = _edge_lengths[2];

    l_rms = sqrt(pow(x1,2)+pow(x2,2)+pow(x3,2));

    return l_rms;
}
