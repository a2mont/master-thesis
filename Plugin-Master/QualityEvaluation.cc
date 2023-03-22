#include "QualityEvaluation.hh"

double QualityEvaluation::evaluate(const OpenMesh::SmartFaceHandle _face, const bool _verbose){
    double quality = 0.;
    double area = 0.;
    double l_harm = 0.;
    double l_rms = 1.;


    std::vector<double> edge_legnths;
    for(auto fvh : _face.vertices()){
        auto edge = fvh.outgoing_halfedges().to_vector().front();
        edge_legnths.push_back(mesh_.calc_edge_length(edge));
    }
    area = calculate_area(edge_legnths);
    l_harm = calculate_l_harm(edge_legnths);
    l_rms = calculate_l_rms(edge_legnths);
    if(!std::isnan(area) && !std::isnan(area) && !std::isnan(area))
        quality = 3 * sqrt(2) * area * l_harm / pow(l_rms, 4);

    if(_verbose)
        std::cout <<
                     "Edge lengths: " << edge_legnths[0] << " " << edge_legnths[1] << " " << edge_legnths[2] <<
                     "\nArea: "<< area <<
                     "\nHarmonic length: " << l_harm <<
                     "\nRoot-mean-squared edge length: " << l_harm <<
                     "\nQuality: " << quality <<
                     "\n------------------------------"
        << std::endl;

    face_quality_queue_.emplace(Triangle(_face.idx(), quality));

    return quality;
}

void QualityEvaluation::reset_queue(){
    PriorityQueue empty;
    std::swap(face_quality_queue_, empty);
}

double QualityEvaluation::calculate_area(std::vector<double> _edge_lengths){
    double area = 0.;
    double x1,x2,x3,semi_perimeter;
    x1 = _edge_lengths[0];
    x2 = _edge_lengths[1];
    x3 = _edge_lengths[2];
    semi_perimeter = (x1+x2+x3)/2;

    area = sqrt(semi_perimeter*(semi_perimeter-x1)*(semi_perimeter-x2)*(semi_perimeter-x3));
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
    double l_rms = sqrt(
                std::accumulate(_edge_lengths.begin(),_edge_lengths.end(),decltype(_edge_lengths)::value_type(0))
                /_edge_lengths.size()
    );

    return l_rms;
}
