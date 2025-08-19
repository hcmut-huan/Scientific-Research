#include "base-c.hpp"
#include "base-d.hpp"


namespace NormalEquation {
    // Source paper: Fast Piecewise Polynomial Fitting of Time-Series Data for Streaming Computing
    // Source path: src/piecewise-approximation/polynomial/normal-equation.cpp
    
    struct Model {
        int length = 0;
        float error = -1;
        Polynomial* function = nullptr;

        Model(Polynomial* function);
        ~Model();
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            double error = 0;
            std::string mode = "";

            Model* model = nullptr;
            std::vector<Point2D> window;

            // std::map<int, Matrix<double>*> cache;    // Use our matrix implementation
            std::map<int, Eigen::MatrixXd> cache;       // Use eigen library
            
            Polynomial* __calPolynomial();
            bool __approxSuccess(Model* model);

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
            int degree = 0;

        protected:
            CSVObj* decompress(BinObj* compress_data) override;
        
        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};