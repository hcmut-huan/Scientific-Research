#include "base-c.hpp"
#include "base-d.hpp"

namespace PMC {
    // Source paper: Capturing Sensor-Generated Time Series with Quality Guarantees
    // Source path: src/piecewise-approximation/constant/pmc.cpp

    class Compression : public BaseCompression {
        private:
            double error = 0;
            std::string mode = "";

            double min = INFINITY;
            double max = -INFINITY;
            double value = 0;
            int length = 0;

            long index = 0;
            long prev_end = -1;
            long curr_end = 0;

            void __mean(Univariate* data);
            void __midrange(Univariate* data);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            long index = -1;

        protected:
            CSVObj* decompress() override;

        public:
            Decompression(std::string input, std::string output, int interval) : BaseDecompression(input, output, interval) {} 
            void initialize() override;
            void finalize() override;
    };
};


namespace HybridPCA {
    // Source paper: Improved Piecewise Constant Approximation Method for Compressing Data Streams
    // Source path: src/piecewise-approximation/constant/hybrid-pmc

    struct Window {
        double min = INFINITY;
        double max = -INFINITY;
        std::vector<Point2D*> data;

        ~Window();
        int size();
        void append(Point2D& p);
    };

    struct Buffer {
        double min = INFINITY;
        double max = -INFINITY;
        std::vector<Window*> windows;

        int size();
        void pop();
        void clear();
        double value();
        void append(Window* window);
        bool is_appendable(Window* window, float bound);
    };

    class Compression : public BaseCompression {
        private:
            double error = 0;
            int w_size = 0;
            int n_window = 0;

            int length = 0;
            double value = 0;
            long index = 0;
            Buffer buffer;
            Window* window = nullptr;
            
            void __pmc();

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        protected:
            CSVObj* decompress() override;

        public:
            Decompression(std::string input, std::string output, int interval) : BaseDecompression(input, output, interval) {} 
            void initialize() override;
            void finalize() override;
    };
};