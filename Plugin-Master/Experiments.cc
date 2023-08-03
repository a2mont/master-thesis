#include "Experiments.hh"

void Experiment::stretch2D(const int _timesteps, const int _currentT){
    if(_currentT > _timesteps)
        return;

    // Physical properties
    // from https://www.bu.edu/moss/mechanics-of-materials-strain/
    double poissonRatio = 0.5;
    double deltaLength = 1;
    double length = 2;
    double axialStrain = deltaLength/length;
    double lateralStrain = - poissonRatio * axialStrain;

    //  calculate target position, executed only once
    if(!initialized_){
        updateBoundaries();
        // Top
        auto topStartPoint = mesh_.point(topBoundary_[0]);
        auto topEndPoint = mesh_.point(topBoundary_.back());
        topEndPoint[0] += deltaLength;
        auto topMidPoint = Point((topEndPoint[0] - topStartPoint[0])/2, topStartPoint[1] + lateralStrain, 0);
        // Bot
        auto botStartPoint = mesh_.point(bottomBoundary_[0]);
        auto botEndPoint = mesh_.point(bottomBoundary_.back());
        botEndPoint[0] += deltaLength;
        auto botMidPoint = Point((botEndPoint[0] - botStartPoint[0])/2, botStartPoint[1] - lateralStrain, 0);

        for(auto vh: rightBoundary_){
            mesh_.set_color(vh, ACG::Vec4f(1,1,0,1));
            for(int i = 0; i < _timesteps; ++i){
                auto step = basePoints_[vh.idx()] + Point(deltaLength/_timesteps * i,0,0);
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(int j =1; j < bottomBoundary_.size() - 1; ++j){
            auto vh = bottomBoundary_[j];
            double lateralTarget = (length+deltaLength)/bottomBoundary_.size();
            for(int i = 0; i < _timesteps; ++i){
                Point pt = basePoints_[vh.idx()];
                Point mid = botMidPoint;
                mid[1] = mid[1]/_timesteps * i;
                double deltaLat = lateralTarget * j - pt[0];
                pt[0] += deltaLat/_timesteps * i;
                auto parabola = generateParabola(botStartPoint, mid, botEndPoint);
                pt[1] = parabola[0] * pow(pt[0], 2) + parabola[1] * pt[0] + parabola[2];
                Point step = pt;
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(int j =1; j < topBoundary_.size() - 1; ++j){
            auto vh = topBoundary_[j];
            double lateralTarget = (length+deltaLength)/topBoundary_.size();
            for(int i = 0; i < _timesteps; ++i){
                Point pt = basePoints_[vh.idx()];
                Point mid = botMidPoint;
                mid[1] = -mid[1]/_timesteps * i;
                double deltaLat = lateralTarget * j - pt[0];
                pt[0] += deltaLat/_timesteps * i;
                auto parabola = generateParabola(botStartPoint, mid, botEndPoint);
                parabola[2] = 2;
                pt[1] = parabola[0] * pow(pt[0], 2) + parabola[1] * pt[0] + parabola[2];
                Point step = pt;
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(auto vh: leftBoundary_){
            for(int i = 0; i < _timesteps; ++i){
                // left border does not move
                auto step = basePoints_[vh.idx()] + Point(0,0,0);
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        initialized_ = true;
    }

    // select vertices to move
    for(auto nextPoint: targetPoints_){
        Point target = nextPoint.second[_currentT];
        auto vh = mesh_.vertex_handle(nextPoint.first);
        mesh_.point(vh) = target;
    }
    // remesh
    loop2D_.loop();
}

void Experiment::compress2D(const int _timesteps, const int _currentT){
    if(_currentT > _timesteps)
        return;
    // Physical properties
    // from https://www.bu.edu/moss/mechanics-of-materials-strain/
    double poissonRatio = 0.5;
    double deltaLength = -0.65;
    double length = 2;
    double axialStrain = deltaLength/length;
    double lateralStrain = 0.5;

    //  calculate target position, executed only once
    if(!initialized_){
        updateBoundaries();
        // Top
        auto topStartPoint = mesh_.point(topBoundary_[0]);
        auto topEndPoint = mesh_.point(topBoundary_.back());
        topEndPoint[0] += deltaLength;
        auto topMidPoint = Point((topEndPoint[0] - topStartPoint[0])/2, topStartPoint[1] - lateralStrain, 0);
        // Bot
        auto botStartPoint = mesh_.point(bottomBoundary_[0]);
        auto botEndPoint = mesh_.point(bottomBoundary_.back());
        botEndPoint[0] += deltaLength;
        auto botMidPoint = Point((botEndPoint[0] - botStartPoint[0])/2, botStartPoint[1] + lateralStrain, 0);

        for(auto vh: rightBoundary_){
            mesh_.set_color(vh, ACG::Vec4f(1,1,0,1));
            for(int i = 0; i < _timesteps; ++i){
                auto step = basePoints_[vh.idx()] + Point(deltaLength/_timesteps * i,0,0);
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(int j =1; j < bottomBoundary_.size() - 1; ++j){
            auto vh = bottomBoundary_[j];
            double lateralTarget = (length+deltaLength)/bottomBoundary_.size();
            for(int i = 0; i < _timesteps; ++i){
                Point pt = basePoints_[vh.idx()];
                Point mid = botMidPoint;
                mid[1] = -mid[1]/_timesteps * i;
                double deltaLat = lateralTarget * j - pt[0];
                pt[0] += deltaLat/_timesteps * i;
                auto parabola = generateParabola(botStartPoint, mid, botEndPoint);
                pt[1] = parabola[0] * pow(pt[0], 2) + parabola[1] * pt[0] + parabola[2];
                Point step = pt;
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(int j =1; j < topBoundary_.size() - 1; ++j){
            auto vh = topBoundary_[j];
            double lateralTarget = (length+deltaLength)/topBoundary_.size();
            for(int i = 0; i < _timesteps; ++i){
                Point pt = basePoints_[vh.idx()];
                Point mid = botMidPoint;
                mid[1] = mid[1]/_timesteps * i;
                double deltaLat = lateralTarget * j - pt[0];
                pt[0] += deltaLat/_timesteps * i;
                auto parabola = generateParabola(botStartPoint, mid, botEndPoint);
                parabola[2] = 2;
                pt[1] = parabola[0] * pow(pt[0], 2) + parabola[1] * pt[0] + parabola[2];
                Point step = pt;
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        for(auto vh: leftBoundary_){
            for(int i = 0; i < _timesteps; ++i){
                // left border does not move
                auto step = basePoints_[vh.idx()] + Point(0,0,0);
                targetPoints_[vh.idx()].push_back(step);
            }
        }
        initialized_ = true;
    }

    // select vertices to move
    for(auto nextPoint: targetPoints_){
        Point target = nextPoint.second[_currentT];
        auto vh = mesh_.vertex_handle(nextPoint.first);
        mesh_.point(vh) = target;
    }
    // remesh
    loop2D_.loop();
}

Eigen::Vector3f Experiment::generateParabola(Point _a, Point _b, Point _c){
    Eigen::Matrix3f A;
    Eigen::Vector3f b,x;
    A <<
        pow(_a[0],2),_a[0], 1,
        pow(_b[0],2),_b[0], 1,
        pow(_c[0],2),_c[0], 1;
    b << _a[1], _b[1], _c[1];
    x = A.colPivHouseholderQr().solve(b);
    return x;
}

Experiment::Point Experiment::findMidPoint(Point _position, std::vector<OpenMesh::SmartVertexHandle> _boundary){
    Point mid(0,0,0);
    double minDistance = 100;
    for(auto vh: _boundary){
        auto pt = mesh_.point(vh);
        auto distance = std::abs(pt[0]- _position[0]);
        if(distance < minDistance){
            minDistance = distance;
            mid = pt;
        }
    }
    return mid;
}

void Experiment::updateBoundaries(){
    leftBoundary_.clear();
    rightBoundary_.clear();
    topBoundary_.clear();
    bottomBoundary_.clear();

    for(auto vh: mesh_.vertices()){
        if(!vh.is_boundary()){
            continue;
        }
        auto pt = mesh_.point(vh);
        targetPoints_[vh.idx()].push_back(pt);
        if(pt[0] <= 0)
            leftBoundary_.push_back(vh);
        else if (pt[0] >= 2)
            rightBoundary_.push_back(vh);
        if(pt[1] <= 0)
            bottomBoundary_.push_back(vh);
        else if (pt[1] >= 2)
            topBoundary_.push_back(vh);
    }

}
