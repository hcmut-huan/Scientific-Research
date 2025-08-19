#include "base-c.hpp"
#include "base-d.hpp"

namespace SmartGridCompression {
    // Source paper: A time-series compression technique and its application to the smart grid
    // Source path: src/model-selection/polynomial/smart-grid-compression.cpp

    class Model {
        public:
            int degree;
            unsigned short length;

            virtual void clear() = 0;
            virtual float getCompressionRatio() = 0;
            virtual bool approximate(float bound, std::vector<Point2D>& segment) = 0;
    };

    // Constant approximate
    class ConstantModel : public Model {
        private:
            float value;
            float min;
            float max;

        public:
            ConstantModel();
            float getValue();

            void clear() override;
            float getCompressionRatio() override;
            bool approximate(float bound, std::vector<Point2D>& segment) override;
    };

    // Linear approximate    
    class LinearModel : public Model {
        private:
            Line* line;
            ConvexHull cvx;

            float __cal_error(Line& line, std::vector<Point2D>& segment);
            float __distance(Line line, Point2D p);
            int __x_external(int x, int x1, int x2);
            int __search(int m, int p, int j);
            Line __approx(const std::vector<Point2D>& segment);

        public:
            LinearModel();
            Line* getLine();

            void clear() override;
            float getCompressionRatio() override;
            bool approximate(float bound, std::vector<Point2D>& segment) override;
    };

    // Polynomial approximate
    class PolynomialModel : public Model {
        private:
            Polynomial* polynomial;

            float __cal_error(std::vector<Point2D>& segment);

        public:
            PolynomialModel(int degree);
            Polynomial* getPolynomial();

            void clear() override;
            float getCompressionRatio() override;
            bool approximate(float bound, std::vector<Point2D>& segment) override;
    };

    class Compression : public BaseCompression {
        private:
            int max_degree = 0;
            double error = 0;

            int curr_degree = 0;
            int chosen_model = -1;
            std::vector<Point2D> window;
            std::vector<Model*> models;

            bool __approximate(int degree);
            void __choose_best_model();

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
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

// namespace Unbounded {
//     // Source path: src/model-selection/polynomial/unbounded.cpp
//     void compress(TimeSeries& timeseries, float bound, std::string output);
//     void decompress(std::string input, std::string output, int interval);
// };

namespace Bounded {
    // Source path: src/model-selection/polynomial/bounded.cpp
    class Model {
        public:
            virtual Polynomial* get() = 0;
            virtual int getDegree () = 0;        
            virtual bool fit(float bound, std::vector<Point2D>& data, int start) = 0;
            virtual bool test(float bound, std::vector<Point2D>& data, int start) = 0;
    };

    class LinearModel : public Model {
        private:
            UpperHull u_cvx; 
            LowerHull l_cvx;
            Line* u_line;
            Line* l_line;
        
        public:        
            LinearModel();
            ~LinearModel();

            Polynomial* get() override;
            int getDegree() override;
            bool test(float bound, std::vector<Point2D>& data, int start) override;
            bool fit(float bound, std::vector<Point2D>& data, int start) override;
    };

    class PolyModel : public Model {
        private:
            int degree;
            Polynomial* func;

        public:
            PolyModel(int degree);
            ~PolyModel();

            Polynomial* get() override;
            int getDegree() override;
            bool test(float bound, std::vector<Point2D>& data, int start) override;
            bool fit(float bound, std::vector<Point2D>& data, int start) override;
    };

    struct Segment {
        int start;
        short length;
        int triedDegree;
        int extreme;

        int childLength;
        int childCoeffByte;

        bool isComplete;
        Model* model;

        Segment(int start, Model* model);
        ~Segment();

        bool trigger();
        void updateChild(Segment* child);
        void updateModel(Model* model, short length);
    };
    
    class Compression : public BaseCompression {
        private:
            int max_degree = 0;
            double error = 0;

            int model_index = 0;
            int trigger_index = -1;
            std::vector<Point2D> window;
            std::vector<Segment*> segments;

            Line __translate(float slp, float intercept, int step);
            int __check(Point2D& p1, Point2D& p2, Point2D& p3, Point2D& p4, Point2D& p5, Line& line_1, Line& line_2);

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
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};