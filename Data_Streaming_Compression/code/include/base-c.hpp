#ifndef BASE_COMPRESSION_HPP
#define BASE_COMPRESSION_HPP

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <eigen/Eigen>

#include "system/io.hpp"
#include "system/monitor.hpp"
#include "algebraic/function.hpp"
#include "algebraic/convex.hpp"
#include "algebraic/sdlp.hpp"
#include "timeseries.hpp"


class BaseCompression {
    private:
        Clock clock;
        BinObj* currObj = nullptr;
        IterIO* timeFile = nullptr;
        IterIO* outputFile = nullptr;
        
    protected:
        BinObj* headObj = nullptr;

        virtual void compress(Univariate* data) = 0;
        virtual BinObj* serialize() = 0;

        void yield() {
            this->currObj->setNext(this->serialize());
            this->currObj = (BinObj*) this->currObj->getNext();
            this->clock.tick();
        }

    public:
        virtual void initialize(int count, char** params) = 0;
        virtual void finalize() = 0;

        BaseCompression(std::string output) {
            this->headObj = new BinObj;
            this->currObj = this->headObj;
            this->outputFile = new IterIO(output, false);
            this->timeFile = new IterIO(output + ".time", false);
        }

        ~BaseCompression() {
            IOObj::clear(this->headObj);
            this->outputFile->close();
            delete this->outputFile;
            this->timeFile->close();
            delete this->timeFile;
        }

        void run(TimeSeries& timeseries) {
            this->clock.start();
            this->headObj->put(((Univariate*) timeseries.get(0))->get_time());
            
            while (timeseries.hasNext()) {
                Univariate* data = (Univariate*) timeseries.next();
                this->compress(data);
            }
            this->finalize();

            long total_time = this->clock.stop();
            double avg_time = (double) total_time / (double) timeseries.size();

            // Profile average latency
            std::cout << std::fixed << "Time taken for each data point (ns): " << avg_time << "\n";
            std::cout << std::fixed << "Average latency (ns): " << this->clock.getAvgDuration() << "\n";
            std::cout << std::fixed << "Max latency (ns): " << this->clock.getMaxDuration() << "\n";
            
            this->timeFile->write("Time taken for each data point (ns): " + std::to_string(avg_time));
            this->timeFile->write("Average latency (ns): " + std::to_string(this->clock.getAvgDuration()));
            this->timeFile->write("Max latency (ns): " + std::to_string(this->clock.getMaxDuration()));

            // Write compression data
            this->outputFile->writeBin(this->headObj);
        }
};

#endif