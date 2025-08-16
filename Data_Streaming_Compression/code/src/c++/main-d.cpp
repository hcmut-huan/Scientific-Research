#include "piecewise-approximation/constant.hpp"
#include "piecewise-approximation/linear.hpp"
#include "piecewise-approximation/polynomial.hpp"
#include "model-selection/polynomial.hpp"


using namespace std;

Monitor Monitor::instance;

int main(int argc, char** argv) {
    const string INPUT = argv[1];
    const string OUTPUT = argv[2];
    const int INTERVAL = atoi(argv[3]);
    const string ALGO = argv[4];

    Monitor::instance.start(OUTPUT+".mon"); 

    if (ALGO == "pmc") {
        PMC::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "hybrid-pca") {
        HybridPCA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "cov-pla") {
        CovariancePLA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "swing-filter") {
        SwingFilter::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "slide-filter") {
        SlideFilter::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "optimal-pla") {
        OptimalPLA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "conn-I-pla") {
        ConnIPLA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "semi-optimal-pla") {
        SemiOptimalPLA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "semi-mixed-pla") {
        SemiMixedPLA::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "normal-equation") {
        NormalEquation::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "mix-piece") {
        MixPiece::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    else if (ALGO == "smart-grid-compression") {
        SmartGridCompression::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }
    // else if (ALGO == "unbounded") {
    //     Unbounded::decompress(INPUT, OUTPUT, INTERVAL);
    // }
    else if (ALGO == "bounded") {
        Bounded::Decompression decompressor(INPUT, OUTPUT, INTERVAL);
        decompressor.initialize();
        decompressor.run();
    }

    Monitor::instance.stop();

    return 0;
}