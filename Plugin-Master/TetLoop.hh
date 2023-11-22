#ifndef TETLOOP_HH
#define TETLOOP_HH

#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

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
        cavityEdge_(mesh_.request_halfface_property<bool>("Cavity Edge")),
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
    struct Star{
        ACG::Vec3d center_;
        std::set<HalfFaceHandle> bounds_;
        std::set<CellHandle> tets_;
        Star(std::set<CellHandle> _tets) :
            tets_(_tets){}
        Star(std::set<CellHandle> _tets, std::set<HalfFaceHandle> _bounds) :
            bounds_(_bounds),
            tets_(_tets){}
        Star(std::set<CellHandle> _tets, std::set<HalfFaceHandle> _bounds, ACG::Vec3d _center) :
            center_(_center),
            bounds_(_bounds),
            tets_(_tets){}
        Star(){}
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
    bool edgeRemoval(EdgeHandle _eh, std::vector<CellHandle>* const _cellsAdded, bool _verbose = true);
    bool edgeRemoval(EdgeHandle _eh, bool _verbose = true);
    bool faceRemoval(FaceHandle _fh, bool _verbose = true);
    VertexHandle contractEdge(EdgeHandle _eh, std::vector<CellHandle>& _tetsAltered);
    bool find_chebyshev_center(const std::vector<HalfFaceHandle> &constraint_hfs, const double radius_lower_bound, ACG::Vec3d &new_position);

    static void reset_queue(PriorityQueue& _queue);
    static void computeQuality(TetLoop::PriorityQueue &_queue, TetrahedralMesh &_mesh);
private:

    bool topologial_pass(PriorityQueue* _A);
    void edge_contraction_pass(PriorityQueue* _A);
    void insertion_pass(PriorityQueue* _A);
    void smoothing_pass(PriorityQueue* _A, int _iterations);

    void improve_mesh(PriorityQueue _badTets);
    void improve_tet(Tet _t);

    void log( Logger* _logger, bool _endOfLine = false);
    void computeQuality();
    void multiFace(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded);
    TestNeighborResult testNeighbor(VertexHandle a,
                      VertexHandle b,
                      VertexHandle u,
                      VertexHandle w);
    double orient3D(VertexHandle _a, VertexHandle _b, VertexHandle _c, VertexHandle _d);
    void flip32Recurse(TetLoop::FaceWithChildren _g, TetLoop::FaceWithChildren _parent, std::vector<CellHandle>* const _cellsAdded);
    void flip32(EdgeHandle _eh, std::vector<CellHandle>* const _cellsAdded);
    void flip23(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded);
    void flip22(FaceHandle _fh, std::vector<CellHandle>* const _cellsAdded);

    CellHandle findNextCell(Star& _startStar);
    void findCavityBoundary(std::set<HalfFaceHandle> &_cavityBoundary);
    void findCavityBoundary(std::set<HalfFaceHandle> &_cavityBoundary, Star _constraint);
    bool checkStarConditions(Star _star);
    bool checkCenter(Star _star);
    void createBoundaryMesh(std::set<HalfFaceHandle> &_cavityBoundary);
    void triangle_normal_and_centroid(const HalfFaceHandle &hf, ACG::Vec3d &normal, ACG::Vec3d &centroid);

private:
    TetrahedralMesh& mesh_;
    HalfFacePropertyT<bool> cavityEdge_;
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
    const std::string LOGS_MESH = "../../../../Plugin-Master/logs/meshes/";

    static constexpr double CHEBY_THRESHOLD = 10e-2;
};

#endif // TETLOOP_HH

