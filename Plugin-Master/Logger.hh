#ifndef OPENFLIPPER_LOGGER_HH
#define OPENFLIPPER_LOGGER_HH

#include <fstream>
#include <iostream>
#include <vector>


class Logger
{
public:
    Logger(std::string _filename, const double _q_min): filename_(_filename), q_min_(_q_min) {
        file_.open(_filename, std::ofstream::out | std::ios::trunc);
        if(file_){
            file_ << "Quality Min,Quality,Topological,Contraction,Insertion,Smoothing\n";
            file_.close();
        }else
            std::cout << "\033[1;4;33mFailed to open file "<<_filename<<"!!\033[0m" << std::endl;
    }
    ~Logger(){
        if(file_.is_open())
            file_.close();
    }


    void logQuality(double _quality);
    void nextLine();
    void logAverage();
    void logLine(std::vector<double> _line, bool _keepFile = false);
    void close();
private:
    std::fstream file_;
    std::string filename_;

    double q_min_;
    std::string q_min_string_ = std::to_string(q_min_);
};

#endif // OPENFLIPPER_LOGGER_HH
