#include "Logger.hh"

void Logger::logQuality(double _quality){
    if(!file_.is_open())
        file_.open(filename_, std::ofstream::out | std::ios::app);
    file_ << _quality << ",";
    file_.close();
}

void Logger::nextLine(){
    if(!file_.is_open())
        file_.open(filename_, std::ofstream::out | std::ios::app);
    file_ << "\n";
    file_.close();
}

void Logger::logLine(const std::vector<double> _line, const bool _keepFile){
    if(!file_.is_open())
        file_.open(filename_, std::ofstream::out | std::ios::app);
    std::string data = q_min_string_ + ",";
    size_t nb(0);
    for(auto l: _line){
        data += std::to_string(l);
        if(++nb < _line.size()){
            data+=",";
        }
    }

    file_ << data << "\n";
    if(!_keepFile){
        file_.close();
    }
}

void Logger::close(){
    if(file_.is_open()){
        file_.close();
    }
}
