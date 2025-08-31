#ifndef BASE_DECOMPRESSION_HPP
#define BASE_DECOMPRESSION_HPP

#include <map>
#include <queue>
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

using namespace std::chrono;

class BaseDecompression {
    private:
        IterIO* decomFile = nullptr;

    protected:
        int interval = 0;
        time_t basetime = 0;

        virtual CSVObj* decompress(BinObj* compress_data) = 0;

    public:
        virtual void initialize(int count, char** params) = 0;
        virtual void finalize() = 0;

        BaseDecompression(std::string output, int interval, std::time_t basetime) {
            this->interval = interval;
            this->basetime = basetime;
            this->decomFile = new IterIO(output, false);
        }

        ~BaseDecompression() {            
            this->decomFile->close();
        }

        long process(BinObj* obj) {
            long length = 0;
            while (obj->getSize() != 0) {
                // std::cout << obj->getSize() << "\n";
                CSVObj* data = this->decompress(obj);
                if (data != nullptr) {
                    length += this->decomFile->write(data);
                }
            }

            return length;
        }

        void complete() {
            this->finalize();
        }
};

#endif