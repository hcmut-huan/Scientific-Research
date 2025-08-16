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

namespace SlideFilter {
    // Source paper: Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees
    // Source path: src/piecewise-approximation/linear/slide-filter.cpp
    class Compression : public BaseCompression {
        private:
            double error = 0;

            double A_num = 0;
            double A_den = 0;
            int segment_pos = 0;
            
            Point2D* pivot = nullptr;
            Line* u_line = nullptr;
            Line* l_line = nullptr;
            ConvexHull cvx;
            std::vector<Point2D> segments;

            long index = 0;
            Line* prev_g = nullptr;
            Line* prev_u = nullptr;
            Line* prev_l = nullptr;
            Point2D* prev_end = nullptr;

            Line __fit();
            double __checkConnected();

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
            const double EPSILON = 1e-6;
            
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

namespace CovariancePLA {
    // Source paper: Maximum error-bounded Piecewise Linear Representation for online stream approximation
    // Source path: src/piecewise-approximation/linear/optimal-pla.cpp
    class Compression : public BaseCompression {
        private:
            double error = 0;

            double accumulate = 0;
            double average_x = 0;
            double average_y = 0;
            long accumulate_square = 0;
            Line* line = nullptr;
            
            UpperHull u_cvx;
            LowerHull l_cvx;

            int length = 0;

            bool __error_bound_verify(Line& line);

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

            bool __yield_10_bytes(BinObj* obj);
            bool __yield_8_bytes(BinObj* obj);
            bool __yield_6_bytes(BinObj* obj);

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

namespace ConnIPLA {
    // Source paper: Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees
    // Source path: src/piecewise-approximation/linear/swing-filter.cpp

    class LinearSegment {
        private:
            bool isOptimal;
            UpperHull u_cvx;
            LowerHull l_cvx;
            Point2D* pivot = nullptr;
            
            void __optimal_approximate(Point2D& p, double error);
            void __fsw_approximate(Point2D& p, double error);

        public:
            long t_start;
            bool is_complete = false;

            Line* u_line = nullptr;
            Line* l_line = nullptr;
            
            LinearSegment(Point2D p);
            LinearSegment(Point2D u_p, Point2D l_p);
            ~LinearSegment();

            Line getLine();
            void approximate(Point2D& p, double error);
    };

    
    class Compression : public BaseCompression {
        private:
            double error = 0;

            int phase = 0;
            LinearSegment* s_1 = nullptr;
            LinearSegment* s_2 = nullptr;

            long index = 0;
            Point2D* prev_end = nullptr;
            Point2D* curr_end = nullptr;

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

namespace SemiOptimalPLA {
    // Source paper: An Optimal Online Semi-Connected PLA Algorithm With Maximum Error Bound
    // Source path: src/piecewise-approximation/linear/semi-optimal-pla.cpp
    class OptimalPLA {
        private:
            UpperHull u_cvx;
            LowerHull l_cvx;
            Point2D* pivot = nullptr;

            std::vector<Point2D> u_points;
            std::vector<Point2D> l_points;
            std::vector<Point2D> window;

            std::vector<std::pair<long, Line>> u_lines;
            std::vector<std::pair<long, Line>> l_lines;

        public:
            long t_start = -1;
            long t_end = -1;
            int status = 0; // 1 is up and -1 is down
            bool is_complete = false;

            Line* extrm = nullptr;
            Line* u_line = nullptr;
            Line* l_line = nullptr;

            ~OptimalPLA();

            void reconstructCvx();
            Point2D shrink(OptimalPLA* prev_seg, double error);
            void extendBackward(OptimalPLA* prev_seg, double error);
            void updateExtrm();
            Point2D popLast(bool flag);
            bool isSemiConnected(Line* line, double bound);            
            void approximate(Point2D& p, double error);
    };

    class Compression : public BaseCompression {
        private:
            double error = 0;

            OptimalPLA* seg_1 = nullptr;
            OptimalPLA* seg_2 = nullptr;
            std::deque<Point2D> buffer;

            long index = 0;
            Point2D* prev_end = nullptr;
            Point2D* curr_end = nullptr;

            int __check(); // check if two segments belong to which case

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

namespace SemiMixedPLA {
    // Source paper: An online PLA algorithm with maximum error bound for generating optimal mixed‐segments
    // Source path: src/piecewise-approximation/linear/semi-mixed-pla.cpp
    class OptimalPLA {
        private:
            UpperHull u_cvx;
            LowerHull l_cvx;
            Point2D* pivot = nullptr;

            std::vector<Point2D> u_points;
            std::vector<Point2D> l_points;

            std::vector<std::pair<long, Line>> u_lines;
            std::vector<std::pair<long, Line>> l_lines;
        
        public:
            long t_start = -1;
            long t_end = -1;
            int status = 0; // 1 is up and -1 is down
            bool is_complete = false;
            std::vector<Point2D> window;

            Line* extrm = nullptr;
            Line* u_line = nullptr;
            Line* l_line = nullptr;

            ~OptimalPLA();
            
            void reconstructCvx();
            Point2D shrink(OptimalPLA* prev_seg, double error);
            void extendBackward(OptimalPLA* prev_seg, double error);
            void updateExtrm();
            Point2D popLast(bool flag);
            bool isSemiConnected(Line* line, double bound);            
            void approximate(Point2D& p, double error);
    };

    class SemiOptimalPLA {
        private:
            int __check(double error); // check if two segments belong to which case

        public:
            int seg_count = 0;
            OptimalPLA* seg_1 = nullptr;
            OptimalPLA* seg_2 = nullptr;
            std::deque<Point2D> buffer;

            BinObj* inter = nullptr;
            Point2D* prev_end = nullptr;
            
            SemiOptimalPLA();
            SemiOptimalPLA(OptimalPLA* p_seg, Point2D* p_point);
            ~SemiOptimalPLA();

            void finalize();
            long getLastEnd();
            BinObj* getInter();
            OptimalPLA* getPendSeg();
            Point2D* getPendPoint();
            void makeSemiConnect(int semi_case, double error);
            // return 0: segments are not completed
            // return 1: case 1 of semi connected
            // return 2: case 2 of semi connected
            // return 3: case 3 of semi connected
            int approximate(Point2D& p, double error);
    };

    class Compression : public BaseCompression {
        private:
            double error = 0;

            SemiOptimalPLA* semi_segs_1 = nullptr;
            SemiOptimalPLA* semi_segs_2 = nullptr;

            // state = 0: fit semi 
            // state = 1: fit semi and mixed
            // state = 2: semi is selected
            // state = 3: mixed is selected
            int state = 0;
            long index = 0;

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
            const double EPSILON = 1e-6;
            
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