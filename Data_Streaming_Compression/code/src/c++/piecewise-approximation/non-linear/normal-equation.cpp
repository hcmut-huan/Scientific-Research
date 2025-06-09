#include "piecewise-approximation/polynomial.hpp"

namespace NormalEquation {

    Model::Model(Polynomial* function) {
        this->function = function;
    }

    Model::~Model() {
        delete this->function;
        this->error = -1;
        this->length = 0;
    }

    // Begin: compression
    // Using eigen library
    Polynomial* Compression::__calPolynomial() {
            int N = this->window.size();
            int n = this->degree;
            Eigen::VectorXd theta;

            if (this->cache.find(this->window.size()) == this->cache.end()) {
                Eigen::MatrixXd X(N, n + 1);
                Eigen::VectorXd y(N);

                for (int i = 0; i < N; i++) {
                    for (int k=0; k<n+1; k++) {
                        X(i, k) = pow(i, k);
                    } 
                    y(i) = this->window[i].y;
                }

                Eigen::MatrixXd X_T_X_inv_X_T = (X.transpose() * X).inverse() * X.transpose();
                this->cache.insert({window.size(), X_T_X_inv_X_T});

                theta = X_T_X_inv_X_T * y;
            }
            else {
                Eigen::VectorXd y(N);
                for (int i=0; i<N; i++) {
                    y(i) = this->window[i].y;
                }

                Eigen::MatrixXd X_T_X_inv_X_T = this->cache.find(this->window.size())->second;
                theta = X_T_X_inv_X_T * y;
            }
            
            float* coeffs = new float[n+1];
            for (int i = 0; i < n+1; ++i) {
                coeffs[i] = theta[i];
            }

            Polynomial* model = new Polynomial(degree, coeffs);
            
            delete coeffs;
            return model;
    }

    // // Using our own matrix implementation
    // Polynomial* Compression::__calPolynomial() {
    //     Matrix<double>* theta = nullptr;
    //     if (cache.find(this->window.size()) == cache.end()) {
    //         Matrix<double> *X = new Matrix<double>(this->window.size(), this->degree+1);
    //         Matrix<double> *y = new Matrix<double>(this->window.size(), 1);

    //         for (int i = 0; i < this->window.size(); i++) {
    //             for (int k=0; k<degree+1; k++) {
    //                 X->cell[i][k] = pow(i, k);
    //             } 
    //             y->cell[i][0] = this->window[i].y;
    //         }

    //         Matrix<double>* X_T = X->transpose();
    //         Matrix<double>* X_T_X = Matrix<double>::matrix_outter_product(X_T, X);
    //         Matrix<double>* X_T_X_inv = X_T_X->inverse();
    //         Matrix<double>* X_T_X_inv_X_T = Matrix<double>::matrix_outter_product(X_T_X_inv, X_T);
    //         theta = Matrix<double>::matrix_outter_product(X_T_X_inv_X_T, y);

    //         cache.insert({this->window.size(), X_T_X_inv_X_T});

    //         delete X; delete y; delete X_T; 
    //         delete X_T_X; delete X_T_X_inv;
    //     }
    //     else {
    //         Matrix<double> *y = new Matrix<double>(this->window.size(), 1);
    //         for (int i=0; i<this->window.size(); i++) {
    //             y->cell[i][0] = this->window[i].y;
    //         }
    //         Matrix<double>* X_T_X_inv_X_T = cache.find(this->window.size())->second;
    //         theta = Matrix<double>::matrix_outter_product(X_T_X_inv_X_T, y);
    //     }
        
    //     double* coeffs = theta->toVec();
    //     Polynomial* model = new Polynomial(degree, coeffs);
    //     delete coeffs; delete theta;
    //     return model;
    // }

    bool Compression::__approxSuccess(Model* model) {
        if (this->mode == "individual") {
            if (model->error == -1) {
                model->error = -INFINITY;
                for (int i = 0; i < this->window.size(); i++) {
                    float error = std::abs(model->function->subs(i) - this->window[i].y);
                    model->error = model->error < error ? error : model->error;
                }
            }
            else {
                float error = std::abs(model->function->subs(this->window.size()-1) - this->window[this->window.size()-1].y);
                model->error = model->error < error ? error : model->error;
            }
            
            return model->error < this->error;
        }
        else if (mode == "accumulate") {
            if (model->error == -1) {
                model->error = 0;
                for (int i = 0; i < this->window.size(); i++) {
                    model->error += std::abs(model->function->subs(i) - this->window[i].y);
                }
            }
            else {
                model->error += std::abs(model->function->subs(this->window.size()-1) - this->window[this->window.size()-1].y);
            }
            
            return model->error < this->error * this->window.size();
        }

        return false;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->degree = atoi(params[2]);

        this->headObj->put((unsigned char) this->degree);
    }

    void Compression::finalize() {
        if (this->window.size() >= this->degree + 1) {
            this->yield();
            delete this->model;
        }
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(Point2D(this->window.size(), data->get_value()));

        if (this->window.size() == this->degree + 1) {
            model = new Model(__calPolynomial());
            this->model->length = this->window.size();
        }
        if (this->window.size() > this->degree + 1) {
            if (!__approxSuccess(this->model)) {
                Model* n_model = new Model(__calPolynomial());
                if (!__approxSuccess(n_model)) {
                    this->yield();
                    window = {Point2D(0, data->get_value())};

                    delete this->model;
                    delete n_model;
                }
                else {
                    delete this->model;
                    this->model = n_model;
                    this->model->length = this->window.size();
                }
            }
            else {
                this->model->length = this->window.size();
            }
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        
        obj->put((short) this->model->length);
        for (int i = 0; i <= this->degree; i++) {
            obj->put((float) this->model->function->coefficients[i]);
        }

        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize() {
        this->degree = (int) this->compress_data->getByte();        
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress() {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        unsigned short length = this->compress_data->getShort();
        float* coefficients = new float[this->degree+1];

        for (int i = 0; i <= this->degree; i++) {
            coefficients[i] = this->compress_data->getFloat();
        }
        Polynomial model(this->degree, coefficients);

        for (int i=0; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(model.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(model.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        delete[] coefficients;
        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};