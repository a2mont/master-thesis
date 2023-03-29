#ifndef OPENFLIPPER_QUALITYEVALUATION_HH
#define OPENFLIPPER_QUALITYEVALUATION_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <queue>

class QualityEvaluation
{
public:
    QualityEvaluation(TriMesh& _mesh) : mesh_(_mesh){}
    ~QualityEvaluation(){}


    struct Triangle{
        int face_id_;
        double quality_;
        Triangle(int _face_id, double _quality): face_id_(_face_id), quality_(_quality){}
        std::string toString() const {return "Face id: " + std::to_string(face_id_) + " Quality: " + std::to_string(quality_);}
    };
    struct CompareQuality{
        bool operator()(Triangle const& t1,Triangle const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    using PriorityQueue = std::priority_queue<Triangle, std::vector<Triangle>, CompareQuality>;

public:
    PriorityQueue& get_face_quality_queue() {return face_quality_queue_;}
    double evaluate(const OpenMesh::SmartFaceHandle _face, const bool _verbose = false);
    void reset_queue();
private:
    double calculate_area(std::vector<double> _edge_lengths);
    double calculate_l_harm(std::vector<double> _edge_lengths);
    double calculate_l_rms(std::vector<double> _edge_lengths);

private:
    TriMesh& mesh_;
    std::priority_queue<Triangle, std::vector<Triangle>, CompareQuality> face_quality_queue_;

};


#endif // OPENFLIPPER_QUALITYEVALUATION_HH
