#include "piecewise-approximation/constant.hpp"
#include "piecewise-approximation/linear.hpp"
#include "piecewise-approximation/polynomial.hpp"
#include "model-selection/polynomial.hpp"

using namespace std;
using namespace std::chrono;


Monitor Monitor::instance;
std::queue<high_resolution_clock::time_point> time_stream;

TimeSeries loadTimeseries(string input) {
    TimeSeries timeseries;
    CSVObj* head_obj = BatchIO::readCSV(input);
    CSVObj* curr_obj = head_obj;

    while (curr_obj != nullptr) {
        time_t time = (time_t) stol(curr_obj->getData(0));
        float value = stof(curr_obj->getData(1));        
        timeseries.push(new Univariate(time, value));

        curr_obj = (CSVObj*) curr_obj->getNext();
    }

    IOObj::clear(head_obj);
    return timeseries;
}


int main(int argc, char** argv) {
    const string INPUT = argv[1];
    const string COM_OUTPUT = argv[2];
    const string DECOM_OUTPUT = argv[3];
    const int INTERVAL = stoi(argv[4]);
    const string ALGO = argv[5];

    TimeSeries data_stream = loadTimeseries(INPUT);
    Monitor::instance.start(".mon");
    
    BaseDecompression* decompressor = nullptr;
    BaseCompression* compressor = nullptr;

    if (ALGO == "pmc") {
        compressor = new PMC::Compression(COM_OUTPUT);       
        decompressor = new PMC::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "hybrid-pca") {
        compressor = new HybridPCA::Compression(COM_OUTPUT);
        decompressor = new HybridPCA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "cov-pla") {
        compressor = new CovariancePLA::Compression(COM_OUTPUT);
        decompressor = new CovariancePLA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "swing-filter") {
        compressor = new SwingFilter::Compression(COM_OUTPUT);
        decompressor = new SwingFilter::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "slide-filter") {
        compressor = new SlideFilter::Compression(COM_OUTPUT);
        decompressor = new SlideFilter::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "optimal-pla") {
        compressor = new OptimalPLA::Compression(COM_OUTPUT);
        decompressor = new OptimalPLA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "conn-I-pla") {
        compressor = new ConnIPLA::Compression(COM_OUTPUT);
        decompressor = new ConnIPLA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "semi-optimal-pla") {
        compressor = new SemiOptimalPLA::Compression(COM_OUTPUT);
        decompressor = new SemiOptimalPLA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "semi-mixed-pla") {
        compressor = new SemiMixedPLA::Compression(COM_OUTPUT);
        decompressor = new SemiMixedPLA::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "normal-equation") {
        compressor = new NormalEquation::Compression(COM_OUTPUT);
        decompressor = new NormalEquation::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "mix-piece") {
        compressor = new MixPiece::Compression(COM_OUTPUT);
        decompressor = new MixPiece::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "smart-grid-compression") {
        compressor = new SmartGridCompression::Compression(COM_OUTPUT);
        decompressor = new SmartGridCompression::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }
    else if (ALGO == "bounded") {
        compressor = new Bounded::Compression(COM_OUTPUT);
        decompressor = new Bounded::Decompression(DECOM_OUTPUT, INTERVAL, ((Univariate*) data_stream.get(0))->get_time()); 
    }

    compressor->initialize(argc - 6, &argv[6]);
    decompressor->initialize(argc - 6, &argv[6]);

    long count = 0;
    long max_latency = -1e18;
    double sum_latency = 0;

    Clock com_clock;
    Clock decom_clock;
    while (data_stream.hasNext()) {
        Univariate* data = (Univariate*) data_stream.next();
        
        com_clock.tick();
        time_stream.push(high_resolution_clock::now());
        BinObj* obj = compressor->process(data);
        com_clock.tick();

        while (obj != nullptr) {
            decom_clock.tick();
            high_resolution_clock::time_point curr_time = high_resolution_clock::now();
            long length = decompressor->process(obj);

            if (length > 0 && time_stream.size() > 0) {
                for (int i=0; i<length; i++) {
                    long latency = duration_cast<nanoseconds>(curr_time - time_stream.front()).count();
                    sum_latency += latency;
                    max_latency = max_latency < latency ? latency : max_latency;
                    
                    time_stream.pop();
                }
            }

            count += length;
            decom_clock.tick();
            obj = (BinObj*) obj->getNext();
        }
        compressor->clear_buffer();
    }

    // Finalizing
    BinObj* obj = compressor->complete();
    while (obj != nullptr) {
        decom_clock.tick();
        high_resolution_clock::time_point curr_time = high_resolution_clock::now();
        long length = decompressor->process(obj);

        if (length > 0) {
            for (int i=0; i<length; i++) {
                long latency = duration_cast<nanoseconds>(curr_time - time_stream.front()).count();
                sum_latency += latency;
                max_latency = max_latency < latency ? latency : max_latency;
                time_stream.pop();
            }
        }
        
        count += length;
        decom_clock.tick();
        obj = (BinObj*) obj->getNext();
    }
    decompressor->complete();
    data_stream.finalize();

    // Profiling time
    std::cout << std::fixed << "Average compress time (ns): " << com_clock.getAvgDuration() << "\n"; 
    std::cout << std::fixed << "Average latency (ns): " << (sum_latency / count) << "\n"; 
    std::cout << std::fixed << "Max latency (ns): " << max_latency << "\n";
    std::cout << std::fixed << "Average decompress segment time (ns): " << decom_clock.getAvgDuration() << "\n"; 
    std::cout << std::fixed << "Max decompress segment time (ns): " << decom_clock.getMaxDuration() << "\n"; 

    IterIO timeFile(".time", false);
    timeFile.write("Average compress time (ns): " + std::to_string(com_clock.getAvgDuration()));
    timeFile.write("Average latency (ns): " + std::to_string(sum_latency / count));
    timeFile.write("Max latency (ns): " + std::to_string(max_latency));
    timeFile.write("Average decompress segment time (ns): " + std::to_string(decom_clock.getAvgDuration()));
    timeFile.write("Max decompress segment time (ns): " + std::to_string(decom_clock.getMaxDuration()));
    timeFile.close();
    
    delete compressor;
    delete decompressor;
    Monitor::instance.stop();

    return 0;
}