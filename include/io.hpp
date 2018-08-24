#ifndef IO_H_
#define IO_H_

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <fstream>


using namespace boost::posix_time;
using namespace boost::gregorian;

class KP{
public:
    KP();
    
private:
    std::vector<float> values;
    std::vector<ptime> timestamp;
};


void read_kp(){
    
}


#endif
