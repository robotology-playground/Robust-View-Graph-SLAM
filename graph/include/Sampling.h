

#include <string>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <yarp/os/LogStream.h>
#include <vector>

class Sampling {
    struct Sample {
        double value;
        double time;
        Sample() { value = 0.0; time = -1.0;}
        Sample(double v, double t) { value = v; time = t; }
    };

private:
    std::vector<Sample> samples;
    unsigned int logSize;
    bool relativeTime;
    double firstTime;
public:
    Sampling(bool relativeTime=false, unsigned int logSize=0) {
        Sampling::firstTime = 0.0;
        Sampling::relativeTime = relativeTime;
        Sampling::logSize = logSize;
        if(logSize > 0)
            samples.resize(logSize);
    }

    void add(double value) { add(value, yarp::os::Time::now()); }

    void add(double value, double time) {
        if(samples.size()==0 && relativeTime) firstTime = time;
        if(logSize == 0) {
            samples.push_back(Sample(value, time - firstTime));
            return;
        }
        try {
            samples[samples.size()] = Sample(value, time - firstTime);
        }
        catch (const std::out_of_range& e) {
            yWarning() << "Samples count is out of maximum range. error: " << e.what();
        }
    }

    const std::vector<Sample>& data() { return samples; }

    bool save(std::string filename) {
        std::ofstream file(filename.c_str());
        if(!file.is_open()) {
            yError() << "cannot open file"<<filename;
            return false;
        }
        std::vector<Sample>::iterator itr = samples.begin();
        file << std::fixed << std::showpoint;
        while(itr++ != samples.end()-1)
            file << std::setprecision(6) << (*itr).time <<" "<< std::setprecision(6) << (*itr).value<<std::endl;
        file.close();
    }

    ~Sampling() { }
};


