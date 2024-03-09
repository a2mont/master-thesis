#ifndef TETLOOP_HH
#define TETLOOP_HH

#include <ObjectTypes/TetrahedralMesh/TetrahedralMesh.hh>
#include "TopologicalLink.hh"
#include "VertexDisplacement.hh"
#include "QualityEvaluation.hh"
#include "Smoothing.hh"
#include "Logger.hh"
#include "eigen3/Eigen/Dense"
#include "OpenVolumeMesh/FileManager/FileManager.hh"
#include "ortools/linear_solver/linear_solver.h"
#include <chrono>

using namespace OpenVolumeMesh;

class TetLoop
{
public:
    TetLoop(TetrahedralMesh& _mesh, double _q_min, std::map<int,int>& _constraint_vhs, const bool _logs = false) :
        mesh_(_mesh),
        world_mesh_(_mesh),
        cavityEdge_(mesh_.request_halfface_property<bool>("Cavity Edge")),
        cavityEdge_wm_(world_mesh_.request_halfface_property<bool>("Cavity Edge")),
        constraint_vhs_(_constraint_vhs),
        includeLogs_(_logs),
        q_min_(_q_min)
    {
        if(_logs){
            logsAddress_ = LOGS_BASE + std::to_string(q_min_).substr(0,6) + LOGS_EXTENSION;
            logger_= new Logger(logsAddress_, q_min_, stats_.names_list());
            computeQuality();
            addToStats(Stats::INIT, quality_queue_.top().quality_);
            double avgQuality(0.);
            while(!quality_queue_.empty()){
                Tet top = quality_queue_.top();
                quality_queue_.pop();
                avgQuality += top.quality_;
            }
            avgQuality /= mesh_.n_logical_cells();
            addToStats(Stats::INIT_AVG, avgQuality);

        }
    }
    TetLoop(TetrahedralMesh& _mesh,
            double _q_min,
            std::map<int,int>& _constraint_vhs,
            const bool _useWorldSpace,
            const bool _logs = false) :
        mesh_(_mesh),
        world_mesh_(_mesh),
        cavityEdge_(mesh_.request_halfface_property<bool>("Cavity Edge")),
        cavityEdge_wm_(world_mesh_.request_halfface_property<bool>("Cavity Edge")),
        constraint_vhs_(_constraint_vhs),
        includeLogs_(_logs),
        q_min_(_q_min),
        useWorldSpace_(_useWorldSpace)
    {
        if(_logs){
            logsAddress_ = LOGS_BASE + std::to_string(q_min_).substr(0,6) + LOGS_EXTENSION;
            logger_= new Logger(logsAddress_, q_min_, stats_.names_list());
        }
        // record initial mesh quality to stats
        computeQuality();
        addToStats(Stats::INIT, quality_queue_.top().quality_);
        double avgQuality(0.);
        while(!quality_queue_.empty()){
            Tet top = quality_queue_.top();
            quality_queue_.pop();
            avgQuality += top.quality_;
        }
        avgQuality /= mesh_.n_logical_cells();
        addToStats(Stats::INIT_AVG, avgQuality);
        std::cout << "Starting tet loop with parameters:\n\t-Q min: "<<
                     q_min_ << " -World space: "<< useWorldSpace_ << " -Use logs: "<< includeLogs_
                  << std::endl;

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
        std::vector<std::vector<VertexHandle>> reconstructionVectors_;
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
        void describe() const {    std::cout << "-Star\n\t-Center = " <<
                                                center_ << "\n\t-Bounds = "<<
                                                iterableToString<std::set<HalfFaceHandle>>(bounds_) <<
                                                   " (total = "<< bounds_.size()<< ")"
                                                << "\n\t-Tets = " <<
                                                iterableToString<std::set<CellHandle>>(tets_)
                                                << " (total= "<< tets_.size() << ")"
                                         << std::endl;
                                }
    };

    struct CompareQuality{
        bool operator()(Tet const& t1,Tet const& t2){
            return t1.quality_ > t2.quality_;
        }
    };

    struct Stats{
        enum StatType {
            INIT,
            INIT_AVG,
            BEFORE,
            BEFORE_AVG,
            AFTER,
            AFTER_AVG,
            TOPOLOGY,
            TOPOLOGY_REJECT,
            CONTRACTION,
            CONTRACTION_REJECT,
            INSERTION,
            INSERTION_REJECT,
            SMOOTHING,
            SMOOTHING_REJECT,
            LAST /* used to get number of entries (HAS to be last !)*/};
        std::map<StatType, const char*> types_titles{
            {INIT, "Initial mesh quality"},
            {INIT_AVG, "Initial average mesh quality"},
            {BEFORE, "Quality before remeshing"},
            {BEFORE_AVG, "Average quality before remeshing"},
            {AFTER, "Quality after remeshing"},
            {AFTER_AVG, "Average quality after remeshing"},
            {TOPOLOGY, "Quality improvements for topological pass"},
            {TOPOLOGY_REJECT, "Quality decrease rejected for topological pass"},
            {CONTRACTION, "Quality improvements for contraction pass"},
            {CONTRACTION_REJECT, "Quality decrease rejected for contraction pass"},
            {INSERTION, "Quality improvements for insertion pass"},
            {INSERTION_REJECT, "Quality decrease rejected for insertion pass"},
            {SMOOTHING, "Quality improvements for smoothing pass"},
            {SMOOTHING_REJECT, "Quality decrease rejected for smoothing pass"},
            {LAST, "Used for enumeration, do NOT fill"}
        };
        std::map<StatType, const char*> types_names{
            {INIT, "Initial"},
            {INIT_AVG, "Initial_avg"},
            {BEFORE, "Before"},
            {BEFORE_AVG, "Before_avg"},
            {AFTER, "After"},
            {AFTER_AVG, "After_avg"},
            {TOPOLOGY, "Topological"},
            {TOPOLOGY_REJECT, "Topological_reject"},
            {CONTRACTION, "Contraction"},
            {CONTRACTION_REJECT, "Contraction_reject"},
            {INSERTION, "Insertion"},
            {INSERTION_REJECT, "Insertion_reject"},
            {SMOOTHING, "Smoothing"},
            {SMOOTHING_REJECT, "Smoothing_reject"},
            {LAST, "Ignore"}
        };
        std::map<StatType, std::vector<double>> stat_data_;



        Stats(){
            for(StatType i = INIT; i != LAST; i = StatType(i+1)){
                StatType type = static_cast<StatType>(i);
                stat_data_[type] = std::vector<double>();
            }
        }
        std::vector<std::string> names_list(){
            std::vector<std::string> names;
            for(auto kv: types_names){
                names.push_back(kv.second);
            }
            return names;
        }
        void flush(){
            // empty content
            for(auto& kv_pair: stat_data_){
                if(kv_pair.first == INIT || kv_pair.first == INIT_AVG) continue;
                kv_pair.second.clear();
            }
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
    bool edgeRemoval(EdgeHandle _eh,
                     std::vector<CellHandle>& _cellsAdded,
                     double& _qualityDelta,
                     double& _rejectedTotal,
                     bool _verbose = true);
    bool edgeRemoval(EdgeHandle _eh,
                     double& _qualityDelta,
                     double& _rejectedTotal,
                     bool _verbose = true);
    bool faceRemoval(FaceHandle _fh,
                     std::vector<CellHandle>& _cellsAdded,
                     std::vector<int>& _counter,
                     double& _qualityDelta,
                     double& _rejectedTotal,
                     bool _verbose = true);
    bool faceRemoval(FaceHandle _fh,
                     double& _qualityDelta,
                     double& _rejectedTotal,
                     bool _verbose = true);
    VertexHandle contractEdge(EdgeHandle _eh,
                              std::vector<CellHandle>& _tetsAltered,
                              double& _qualityDelta,
                              double& _rejectTotal);
    bool find_chebyshev_center(const TetrahedralMesh &_mesh,
                               const std::set<HalfFaceHandle> &constraint_hfs,
                               const double radius_lower_bound,
                               ACG::Vec3d &new_position,
                               bool printDebug = false);

    void flip32(TetrahedralMesh& _mesh, EdgeHandle _eh, std::vector<CellHandle>& _cellsAdded);
    void flip23(TetrahedralMesh& _mesh, FaceHandle _fh, std::vector<CellHandle>& _cellsAdded);
    void flip22(TetrahedralMesh& _mesh, FaceHandle _fh, std::vector<CellHandle>& _cellsAdded);
    void multiFace(TetrahedralMesh& _mesh, FaceHandle _fh, std::vector<CellHandle>& _cellsAdded);

    static void reset_queue(PriorityQueue& _queue);
    static void computeQuality(TetLoop::PriorityQueue &_queue, TetrahedralMesh &_mesh);
    template<typename T>
    static void computeQuality(TetLoop::PriorityQueue &_queue,
                        TetrahedralMesh &_mesh,
                        const T &_iterableToEvaluate);
    template<typename T>
    static void printIterable(T _toPrint, bool _eol = false);
    template<typename T>
    static std::string iterableToString(T _toPrint);
    template<typename T>
    static std::string vectorToString(std::vector<T> _toPrint);
    template<typename T>
    static std::string setToString(std::set<T> _toPrint);
    template<typename T>
    static void clearDuplicates(std::vector<T> &_vector);

    static constexpr double CHEBY_THRESHOLD = 0.01;
    static void cleanMesh(TetrahedralMesh& _mesh, bool _keepBoundary = false);
    static void cleanQualityQueue(PriorityQueue& _queue, TetrahedralMesh& _mesh);
    static void displayIterationTime(std::chrono::system_clock::time_point& _begin, std::string _name);
    static void logStats(Stats& _stats, Logger& _logger);
private:

    void subdivision_pass(PriorityQueue &_A, double &_qualityDelta);
    bool topological_pass(PriorityQueue& _A, double& _qualityDelta, double& _rejectedTotal);
    void edge_contraction_pass(PriorityQueue& _A, double& _qualityDelta);
    void insertion_pass(PriorityQueue& _A, double& _qualityDelta);
    void smoothing_pass(PriorityQueue& _A, double& _qualityDelta, int _iterations = 1);

    void improve_mesh(PriorityQueue& _badTets);
    void improve_tet(Tet _t);

    void log(Logger* _logger, bool _endOfLine = false);
    void computeQuality();
    TestNeighborResult testNeighbor(TetrahedralMesh& _mesh,
                                    VertexHandle a,
                                    VertexHandle b,
                                    VertexHandle u,
                                    VertexHandle w);
    double orient3D(TetrahedralMesh& _mesh,
                    VertexHandle _a,
                    VertexHandle _b,
                    VertexHandle _c,
                    VertexHandle _d);
    void flip32Recurse(TetrahedralMesh& _mesh,
                       TetLoop::FaceWithChildren _g,
                       FaceHandle _parent,
                       std::vector<CellHandle>& _cellsAdded);

    CellHandle findNextCell(Star& _startStar, std::vector<Star>& _galaxy, bool& _resultsLeft);
    void findCavityBoundary(std::set<HalfFaceHandle> &_cavityBoundary);
    void findCavityBoundary(Star& _constraint, bool _isWorldMesh = false);
    bool checkStarConditions(Star& _star, CellHandle _lastAdded);
    bool checkCenter(Star _star);
    void createBoundaryMesh(std::set<HalfFaceHandle> &_cavityBoundary);
    void cavityMesh3D(Star& _star);
    void triangle_normal_and_centroid(const HalfFaceHandle &hf,
                                      ACG::Vec3d &normal,
                                      ACG::Vec3d &centroid,
                                      bool printDebug = false);
    void triangle_normal_and_centroid(const TetrahedralMesh &_mesh,
                                      const HalfFaceHandle &hf,
                                      ACG::Vec3d &normal,
                                      ACG::Vec3d &centroid,
                                      bool printDebug = false);
    void recoverBadStar(Star &_star, std::string _meshName);
    void compareStarWithMesh(Star& _star, TetrahedralMesh& _mesh);
    bool validateWorldMesh(TetrahedralMesh &_mesh);
    double computeQualityDelta(double _before);
    void addToStats(Stats::StatType _statName, double _quality);
    void saveMesh(TetrahedralMesh &_mesh, std::set<CellHandle> &_cellsToKeep, std::string _name);
    void saveMesh(TetrahedralMesh &_mesh, std::string _name);
    void compute_volumes_and_min_heights(TetrahedralMesh &mesh, double &min_height, double &min_vol);

private:
    TetrahedralMesh& mesh_;
    TetrahedralMesh world_mesh_;
    HalfFacePropertyT<bool> cavityEdge_;
    HalfFacePropertyT<bool> cavityEdge_wm_;
    PriorityQueue quality_queue_;

    Logger *logger_;
    Stats stats_;


    std::map<int,int>& constraint_vhs_;

    std::string logsAddress_;
    const bool includeLogs_;
    const double q_min_;
    const bool useWorldSpace_ = false;


    const std::string LOGS_BASE = "../../../../Plugin-Master/logs/3D/logs_";
    const std::string LOGS_EXTENSION = ".csv";
    const std::string LOGS_MESH = "../../../../Plugin-Master/logs/meshes/";

};

#endif // TETLOOP_HH

