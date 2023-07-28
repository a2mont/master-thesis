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

