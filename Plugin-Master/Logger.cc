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

void Logger::logLine(std::vector<double> _line, bool _keepFile){
    if(!file_.is_open())
        file_.open(filename_, std::ofstream::out | std::ios::app);
    if(_line.size() == 5){
        std::string data = q_min_string_ + ","
                + std::to_string(_line[0]) + ","
                + std::to_string(_line[1]) + ","
                + std::to_string(_line[2]) + ","
                + std::to_string(_line[3]) + ","
                + std::to_string(_line[4]);
        file_ << data << "\n";
    }else{
        std::cout << "Line contains too many elements !" << std::endl;
    }
    if(!_keepFile){
        file_.close();
    }
}

void Logger::close(){
    if(file_.is_open()){
        file_.close();
    }
}
