#include "piecewise-approximation/constant.hpp"
#include "piecewise-approximation/linear.hpp"
#include "piecewise-approximation/polynomial.hpp"
#include "model-selection/polynomial.hpp"

using namespace std;

Monitor Monitor::instance;

TimeSeries loadTimeseries(string input) {
    TimeSeries timeseries;
    CSVObj* obj = BatchIO::readCSV(input);
    CSVObj* curr_obj = obj;

    while (curr_obj != nullptr) {
        time_t time = (time_t) stol(curr_obj->getData(0));
        float value = stof(curr_obj->getData(1));        
        timeseries.push(new Univariate(time, value));

        curr_obj = (CSVObj*) curr_obj->getNext();
    }

    IOObj::clear(obj);
    return timeseries;
}


int main(int argc, char** argv) {
    const string INPUT = argv[1];
    const string OUTPUT = argv[2];
    const string ALGO = argv[3];

    TimeSeries timeseries = loadTimeseries(INPUT);
    Monitor::instance.start(OUTPUT+".mon");

    if (ALGO == "pmc") {
        PMC::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "hybrid-pca") {
        HybridPCA::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "swing-filter") {
        SwingFilter::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "slide-filter") {
        SlideFilter::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "optimal-pla") {
        OptimalPLA::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "semi-optimal-pla") {
        SemiOptimalPLA::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "semi-mixed-pla") {
        SemiMixedPLA::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "normal-equation") {
        NormalEquation::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "mix-piece") {
        MixPiece::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    else if (ALGO == "smart-grid-compression") {
        SmartGridCompression::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }
    // else if (ALGO == "unbounded") {
    //     Unbounded::compress(timeseries, ERROR, OUTPUT);
    // }
    else if (ALGO == "bounded") {
        Bounded::Compression compressor(OUTPUT);
        compressor.initialize(argc - 4, &argv[4]);
        compressor.run(timeseries);
    }

    timeseries.finalize();
    Monitor::instance.stop();

    return 0;
}