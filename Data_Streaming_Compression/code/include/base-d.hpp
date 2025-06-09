#ifndef BASE_DECOMPRESSION_HPP
#define BASE_DECOMPRESSION_HPP

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <eigen/Eigen>

#include "system/io.hpp"
#include "system/monitor.hpp"
#include "algebraic/sdlp.hpp"
#include "algebraic/convex.hpp"
#include "algebraic/matrix.hpp"
#include "algebraic/function.hpp"
#include "timeseries.hpp"

class BaseDecompression {
    private:
        Clock clock;
        IterIO* timeFile = nullptr;
        IterIO* outputFile = nullptr;

    protected:
        int interval = 0;
        time_t basetime = 0;
        BinObj* compress_data = nullptr;

        virtual CSVObj* decompress() = 0;

    public:
        virtual void initialize() = 0;
        virtual void finalize() = 0;

        BaseDecompression(std::string input, std::string output, int interval) {
            this->interval = interval;
            this->outputFile = new IterIO(output, false);
            this->timeFile = new IterIO(output + ".time", false);
            
            IterIO inputFile(input, true, true);
            this->compress_data = inputFile.readBin();
            inputFile.close();
        }

        ~BaseDecompression() {
            IOObj::clear(compress_data);
            this->outputFile->close();
            delete this->outputFile;
            this->timeFile->close();
            delete this->timeFile;
        }

        void run() {
            this->basetime = this->compress_data->getLong();
            this->clock.start();
            
            while (this->compress_data->getSize() != 0) {
                CSVObj* data = this->decompress();
                if (data != nullptr) {
                    this->outputFile->write(data);
                    IOObj::clear(data);
                    this->clock.tick();
                }
            }
            this->finalize();

            // Profile average latency
            std::cout << std::fixed << "Average time taken for each segment (ns): " << this->clock.getAvgDuration() << "\n";
            std::cout << std::fixed << "Max time taken for each segment (ns): " << this->clock.getMaxDuration() << "\n";

            this->timeFile->write("Time taken for each segment (ns): " + std::to_string(this->clock.getAvgDuration()));
            this->timeFile->write("Max time taken for each segment (ns): " + std::to_string(this->clock.getMaxDuration()));
        }
};

#endif