#ifndef OPENFLIPPER_LOGGER_HH
#define OPENFLIPPER_LOGGER_HH

#include <fstream>
#include <iostream>

class Logger
{
public:
    Logger(std::string _filename, const double _q_min): filename_(_filename) {
        file_.open(_filename, std::ofstream::out | std::ios::trunc);
        if(file_){
            file_ << "Quality Min\n";
            file_ << _q_min << "\n";
            file_.close();
        }else
            std::cout << "\033[1;4;33mFailed to open file !!\033[0m" << std::endl;
    }
    ~Logger(){
        if(file_.is_open())
            file_.close();
    }


    void logQuality(double _quality);
    void nextLine();
private:
    std::fstream file_;
    std::string filename_;
};

#endif // OPENFLIPPER_LOGGER_HH
