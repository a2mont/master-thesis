#ifndef TETLOOP_HH
#define TETLOOP_HH

#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>

#include "TopologicalLink.hh"
#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
#include "Smoothing.hh"
#include "Logger.hh"
#include "eigen3/Eigen/Dense"

using namespace OpenVolumeMesh;

class TetLoop
{
public:
    TetLoop(TetrahedralMesh& _mesh, double _q_min, std::map<int,int>& _constraint_vhs, const bool _logs = false) :
        mesh_(_mesh),
        constraint_vhs_(_constraint_vhs),
        includeLogs_(_logs),
        q_min_(_q_min)
    {
        if(_logs){
            logsAddress_ = LOGS_BASE + std::to_string(q_min_).substr(0,4) + LOGS_EXTENSION;
            std::string logsTimesteps = LOGS_STEP + std::to_string(q_min_).substr(0,4) + LOGS_EXTENSION;
            logger_= new Logger(logsAddress_, q_min_);
            timeStepLogger_ = new Logger(logsTimesteps, q_min_);
        }

    }
    ~TetLoop(){}

    struct Tet{
        CellHandle cell_handle_;
        double quality_;
        Tet(CellHandle _cell_handle, double _quality): cell_handle_(_cell_handle), quality_(_quality){}
        std::string toString() const {return "Cell id: " + std::to_string(cell_handle_.idx()) + " Quality: " + std::to_string(quality_);}
    };
    struct CompareQuality{
        bool operator()(Tet const& t1,Tet const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    using PriorityQueue = std::priority_queue<Tet, std::vector<Tet>, CompareQuality>;
    using Point = ACG::Vec3d;

private:
    struct FaceWithChildren
    {
        FaceHandle fh_;
        std::set<FaceHandle> children_;
        FaceWithChildren(){}
        FaceWithChildren(FaceHandle _fh): fh_(_fh) {}
        FaceWithChildren(FaceHandle _fh, std::set<FaceHandle> _children): fh_(_fh), children_(_children) {}
    };
    struct TestNeighborResult {
        double o_;
        double n_;
        std::vector<FaceWithChildren> h_;
        TestNeighborResult(double _o, double _n, std::vector<FaceWithChildren> _h):
            o_(_o),
            n_(_n),
            h_(_h)
        {}
    };

public:
    void loop(int _max_iter=1);
    void static reset_queue(PriorityQueue& _queue);
    bool edgeRemoval(EdgeHandle _eh);
    VertexHandle contractEdge(EdgeHandle _eh, std::vector<CellHandle>& _tetsAltered);
private:

    bool topologial_pass(PriorityQueue* _A);
    void edge_contraction_pass(PriorityQueue* _A);
    void insertion_pass(PriorityQueue* _A);
    void smoothing_pass(PriorityQueue* _A, int _iterations);

    void improve_mesh(PriorityQueue _badTets);
    void improve_tet(Tet _t);

    void log( Logger* _logger, bool _endOfLine = false);
    void computeQuality();
    bool faceRemoval(FaceHandle _fh);
    void multiFace(FaceHandle _fh);
    TestNeighborResult testNeighbor(VertexHandle a,
                      VertexHandle b,
                      VertexHandle u,
                      VertexHandle w);
    double orient3D(VertexHandle _a, VertexHandle _b, VertexHandle _c, VertexHandle _d);

private:
    TetrahedralMesh& mesh_;

    PriorityQueue quality_queue_;

    Logger *logger_;
    Logger *timeStepLogger_;

    std::map<int,int>& constraint_vhs_;

    std::string logsAddress_;
    const bool includeLogs_;
    const double q_min_;
    const std::string LOGS_BASE = "../../../../Plugin-Master/logs/3D/quality_logs_";
    const std::string LOGS_EXTENSION = ".csv";
    const std::string LOGS_STEP = "../../../../Plugin-Master/logs/3D/quality_timesteps_";

};

#endif // TETLOOP_HH

