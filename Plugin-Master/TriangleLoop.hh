#ifndef OPENFLIPPER_MAINLOOP_HH
#define OPENFLIPPER_MAINLOOP_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
#include "Smoothing.hh"
#include "Logger.hh"
#include "eigen3/Eigen/Dense"

class TriangleLoop
{

public:
    TriangleLoop(TriMesh& _mesh, TriMesh& _worldMesh, double _q_min, std::map<int,int>& _constraint_vhs, const bool _logs = false) :
        mesh_(_mesh),
        worldSpaceMesh_(_worldMesh),
        constraint_vhs_(_constraint_vhs),
        includeLogs_(_logs),
        q_min_(_q_min)
    {
        if(!mesh_.get_property_handle(face_visited_, "face visited"))
            mesh_.add_property(face_visited_, "face visited");
        for(auto fh: mesh_.faces()){
            mesh_.property(face_visited_, fh) = false;
        }
        if(!mesh_.get_property_handle(cavity_edge_, "cavity edge"))
            mesh_.add_property(cavity_edge_, "cavity edge");
        for(auto heh: mesh_.halfedges()){
            mesh_.property(cavity_edge_, heh) = false;
        }
        if(!mesh_.get_property_handle(deleted_border_, "deleted border"))
            mesh_.add_property(deleted_border_, "deleted border");
        for(auto vh: mesh_.vertices()){
            mesh_.property(deleted_border_, vh) = false;
        }
        if(_logs){
            logsAddress_ = LOGS_BASE + std::to_string(q_min_).substr(0,4) + LOGS_EXTENSION;
            std::string logsTimesteps = LOGS_STEP + std::to_string(q_min_).substr(0,4) + LOGS_EXTENSION;
//            logger_= new Logger(logsAddress_, q_min_);
//            timeStepLogger_ = new Logger(logsTimesteps, q_min_);
        }

    }
    ~TriangleLoop(){
        mesh_.remove_property(face_visited_);
        mesh_.remove_property(cavity_edge_);
    }

    struct Triangle{
        OpenMesh::SmartFaceHandle face_handle_;
        double quality_;
        Triangle(OpenMesh::SmartFaceHandle _face_handle, double _quality): face_handle_(_face_handle), quality_(_quality){}
        std::string toString() const {return "Face id: " + std::to_string(face_handle_.idx()) + " Quality: " + std::to_string(quality_);}
    };
    struct CompareQuality{
        bool operator()(Triangle const& t1,Triangle const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    using PriorityQueue = std::priority_queue<Triangle, std::vector<Triangle>, CompareQuality>;
    using Point = ACG::Vec3d;


public:
    void loop(int _max_iter=1);
private:
    void static reset_queue(PriorityQueue& _queue);

    bool topologial_pass(PriorityQueue* _A);
    void edge_contraction_pass(PriorityQueue* _A);
    void insertion_pass(PriorityQueue* _A);
    void smoothing_pass(PriorityQueue* _A, int iterations);

    void improve_mesh(PriorityQueue _badTriangles);
    void improve_triangle(Triangle _t);

    void find_faces_with_p(std::vector<OpenMesh::SmartFaceHandle> &_list, OpenMesh::SmartFaceHandle _fh, const Point _p);
    bool contains_p(OpenMesh::SmartFaceHandle _fh, const Point _p);
    void log( Logger* _logger, bool _endOfLine = false);
    void computeQuality();

private:
    TriMesh& mesh_;
    TriMesh& worldSpaceMesh_;

    // true means the circumcircle contains the new vertex p
    OpenMesh::FPropHandleT<bool> face_visited_;
    // true means the halfedge is at the border of the cavity
    OpenMesh::HPropHandleT<bool> cavity_edge_;
    // true means the vertex from a deleted face is at the border of the mesh
    OpenMesh::VPropHandleT<bool> deleted_border_;

    PriorityQueue quality_queue_;

    Logger *logger_;
    Logger *timeStepLogger_;

    std::map<int,int>& constraint_vhs_;
    std::string logsAddress_;
    const bool includeLogs_;
    const double q_min_;
    const std::string LOGS_BASE = "../../../../Plugin-Master/logs/2D/quality_logs_";
    const std::string LOGS_EXTENSION = ".csv";
    const std::string LOGS_STEP = "../../../../Plugin-Master/logs/2D/quality_timesteps_";

};

#endif // OPENFLIPPER_MAINLOOP_HH
