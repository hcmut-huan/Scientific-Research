#include "base-c.hpp"
#include "base-d.hpp"

namespace SwingFilter {
    // Source paper: Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees
    // Source path: src/piecewise-approximation/linear/swing-filter.cpp
    
    class Compression : public BaseCompression {
        private:
            double error = 0;

            double A_num = 0;
            double A_den = 0;
            Point2D* pivot = nullptr;
            Line* u_line = nullptr;
            Line* l_line = nullptr;

            long index = 0;
            Point2D* prev_end = nullptr;
            Point2D* curr_end = nullptr;

            void __fit(Point2D& p);

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
            long index = 0;
            Point2D* prev_end = nullptr;

        protected:
            CSVObj* decompress() override;

        public:
            Decompression(std::string input, std::string output, int interval) : BaseDecompression(input, output, interval) {}
            void initialize() override;
            void finalize() override;
    };
};

// namespace SlideFilter {
//     // Source paper: Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees
//     // Source path: src/piecewise-approximation/linear/slide-filter.cpp
//     void compress(TimeSeries& timeseries, float bound, std::string output);
//     void decompress(std::string input, std::string output, int interval);
// };

namespace OptimalPLA {
    // Source paper: Maximum error-bounded Piecewise Linear Representation for online stream approximation
    // Source path: src/piecewise-approximation/linear/optimal-pla.cpp
    class Compression : public BaseCompression {
        private:
            double error = 0;

            Point2D* pivot = nullptr;
            Line* u_line = nullptr;
            Line* l_line = nullptr;
            UpperHull u_cvx;
            LowerHull l_cvx;

            int length = 0;

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


namespace MixPiece {
    // Source paper: Flexible grouping of linear segments for highly accurate lossy compression of time series data
    // Source path: src/piecewise-approximation/linear/mix-piece.cpp
    struct Segment {
        long t;
        double aMin;
        double aMax;
        double b;
        double a;

        Segment(long t, double a, double b);
        Segment(long t, double aMin, double aMax, double b);
    };

    struct BBlock {

        struct Block {
            double a_u;
            double a_l;
            std::vector<long> t;

            Block();
            Block(double a_u, double a_l, long t);
        };

        double b;
        std::vector<Block> blocks;

        BBlock(double b);
        BBlock(double b, double a_u, double a_l, long t);

        void sort();
        Segment pop_back();
        void push_back(double a_u, double a_l, long t);
        bool is_intersect(double a_u, double a_l);
        void intersect(double a_u, double a_l, long t);
    };

    struct ABlock {

        struct Block {
            double b;
            long t;

            Block(double b, long t);
        };

        double a_u;
        double a_l;
        std::vector<Block> blocks;

        ABlock();
        ABlock(double b, double a_u, double a_l, long t);

        void sort();
        bool is_intersect(double a_u, double a_l);
        void intersect(double b, double a_u, double a_l, long t);
    };

    class Compression : public BaseCompression {
        private:
            double error = 0;
            int n_segment = 0;

            int flag = 0;
            bool floor_flag = true;
            bool ceil_flag = true;
            double slp_u_1 = INFINITY;
            double slp_u_2 = INFINITY;
            double slp_l_1 = -INFINITY;
            double slp_l_2 = -INFINITY;
            long t_s = 0;
            double b_1 = 0;
            double b_2 = 0;

            std::vector<BBlock> BBlocks;
            std::vector<ABlock> ABlocks;
            std::vector<Segment> r_blocks;
            std::map<double, std::vector<Segment>> intervals;

            long index = 0;
            int curr_segment = 0;
            
            void __group();
            void __serializeBblock(BinObj* obj);
            void __serializeAblock(BinObj* obj);
            void __serializeRblock(BinObj* obj);

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
            long start = 0;
            std::pair<CSVObj*, CSVObj*> __decompress_segment(Segment& segment);

        protected:
            CSVObj* decompress() override;

        public:
            Decompression(std::string input, std::string output, int interval) : BaseDecompression(input, output, interval) {}
            void initialize() override;
            void finalize() override;
    };
};