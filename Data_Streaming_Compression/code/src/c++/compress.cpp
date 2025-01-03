#include "piecewise-approximation/constant.h"
#include "piecewise-approximation/linear.h"
#include "piecewise-approximation/polynomial.h"
#include "piecewise-approximation/model-seletion.h"

using namespace std;

long Monitor::page_size;
bool Monitor::flag = false;


TimeSeries loadTimeseries(string input) {
    TimeSeries timeseries;
    CSVObj* obj = BatchIO::readCSV(input);
    CSVObj* d_obj = nullptr;

    while (obj != nullptr) {
        time_t time = (time_t) stol(obj->getData(0));
        float value = stof(obj->getData(1));        
        timeseries.push(new Univariate<float>(time, value));

        d_obj = obj;
        obj = (CSVObj*) obj->getNext();
        delete d_obj;
    }

    return timeseries;
}


int main(int argc, char** argv) {
    const string INPUT = argv[1];
    const string OUTPUT = argv[2];
    const string OUT_MONITOR = argv[3];
    const float ERROR = atof(argv[4]);
    const string ALGO = argv[5];

    TimeSeries timeseries = loadTimeseries(INPUT);
    std::thread monitor(&Monitor::monitor, OUT_MONITOR);
    Monitor::start();

    if (ALGO == "pmc") {
        PMC::compress(timeseries, argv[6], ERROR, OUTPUT);
    }
    else if (ALGO == "hybrid-pca") {
        HybridPCA::compress(timeseries, atoi(argv[6]), atoi(argv[7]), ERROR, OUTPUT);
    }
    else if (ALGO == "swing-filter") {
        SwingFilter::compress(timeseries, ERROR, OUTPUT);
    }
    else if (ALGO == "slide-filter") {
        SlideFilter::compress(timeseries, ERROR, OUTPUT);
    }
    else if (ALGO == "optimal-pla") {
        OptimalPLA::compress(timeseries, ERROR, OUTPUT);
    }
    else if (ALGO == "normal-equation") {
        NormalEquation::compress(timeseries, argv[6], atoi(argv[7]), ERROR, OUTPUT);
    }
    else if (ALGO == "mix-piece") {
        MixPiece::compress(timeseries, atoi(argv[6]), ERROR, OUTPUT);
    }

    timeseries.finalize();
    Monitor::stop();
    monitor.join();

    return 0;
}