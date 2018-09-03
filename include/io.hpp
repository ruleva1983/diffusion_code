#ifndef IO_H_
#define IO_H_

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <fstream>
#include <string>


using namespace boost::posix_time;
using namespace boost::gregorian;

class KP{
public:
    KP(std::string filename){
        
        std::ifstream infile(filename);
        std::string line;
        while (std::getline(infile, line)){
            std::istringstream iss(line);
            std::string date, time, datetime;
            float value;
            if (!(iss >> date >> time >> value)) { break; } // error
            datetime = date + " " + time;
            timestamp.push_back(boost::posix_time::time_from_string("2002-01-20 23:59:59.000"));

            }
    }
    
private:
    std::vector<float> values;
    std::vector<boost::posix_time::ptime> timestamp;
};


void read_kp(){
    
}


#endif
