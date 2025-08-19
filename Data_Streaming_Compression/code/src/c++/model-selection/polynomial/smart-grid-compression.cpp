#include "model-selection/polynomial.hpp"

namespace SmartGridCompression {
    // Begin: material
    ConstantModel::ConstantModel() {
        this->degree = 0;
        this->length = 0;
        
        this->value = 0;
        this->min = INFINITY;
        this->max = -INFINITY;
    }

    float ConstantModel::getValue() {
        return this->value;
    }

    float ConstantModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else return (float) this->length / (2 + 4);
    }

    void ConstantModel::clear() {
        this->degree = 0;
        this->length = 0;
        
        this->min = INFINITY;
        this->max = -INFINITY;
    }

    bool ConstantModel::approximate(float bound, std::vector<Point2D>& segment) {
        if (this->length == 0) {
            for (Point2D& p : segment) {
                this->min = this->min < p.y ? this->min : p.y;
                this->max = this->max > p.y ? this->max : p.y;
                
                if (this->max - this->min <= 2 * bound) {
                    this->value = (this->max + this->min) / 2;
                    this->length++;
                }
                else {
                    return false;
                }
            }
            
            return true;
        }
        else {
            Point2D& p = segment.back();
            this->min = this->min < p.y ? min : p.y;
            this->max = this->max > p.y ? max : p.y;
            
            if (max - min <= 2 * bound) {
                this->value = (max + min) / 2;
                this->length++;
                return true;
            }
            else {
                return false;
            }
        }
    }

    float LinearModel::__cal_error(Line& line, std::vector<Point2D>& segment) {
        float max_error = 0;
        for (Point2D& p : segment) {
            float error = std::abs(line.subs(p.x) - p.y);
            max_error = error > max_error ? error : max_error;
        }
        return max_error;
    }

    float LinearModel::__distance(Line line, Point2D p) {
        return std::abs(p.y - line.subs(p.x)) 
            / sqrt(line.get_slope() * line.get_slope() + 1);
    }

    int LinearModel::__x_external(int x, int x1, int x2) {
        // 0:  x_external
        // 1:  x_external to the right
        // -1: x_external to the left

        if (x >= x1 && x >= x2) return 1;
        else if (x <= x1 && x <= x2) return -1;
        else return 0;
    }

    int LinearModel::__search(int m, int p, int j) {
        Line line = Line::line(this->cvx.at_ccw(p), this->cvx.at_ccw(j));

        // search in ccw order
        // if pj(i+1) < pj(i) return i
        int v = -1; float prev_dis = -INFINITY;
        for (int i=0; i<this->cvx.size(); i++) {
            float dis = this->__distance(line, this->cvx.at_ccw((m+i)%this->cvx.size()));
            if (prev_dis < dis) {
                prev_dis = dis;
                v = (m+i)%this->cvx.size();
            }
            else break;
        }

        return v;
    }

    Line LinearModel::__approx(const std::vector<Point2D>& segment) {
        if (this->cvx.size() <= 2) return Line::line(segment[0], segment[1]);
        else if (this->cvx.size() > 2) {
            // find the A, B and C
            Point2D A(0, 0), B(0,0), C(0,0);
            int prev_v_l = this->cvx.rightmost_index_ccw();

            // find the first v(l1)
            for (int i=0; i<this->cvx.size(); i++) {
                int v_l = this->__search(prev_v_l, i, (i+1)%this->cvx.size());                        
                int external = this->__x_external(
                    this->cvx.at_ccw(v_l).x, this->cvx.at_ccw(i).x, this->cvx.at_ccw((i+1)%this->cvx.size()).x
                );

                if (external == 1) {
                    prev_v_l = v_l;
                }
                else if (external == 0) {
                    A = this->cvx.at_ccw(i);
                    B = this->cvx.at_ccw((i+1)%this->cvx.size());
                    C = this->cvx.at_ccw(v_l);
                    break;
                }
                else {
                    int step = v_l > prev_v_l ? v_l - prev_v_l : v_l - prev_v_l + this->cvx.size();
                    for (int j=0; j<step; j++) {
                        int external = this->__x_external(
                            this->cvx.at_ccw(i).x,
                            this->cvx.at_ccw((prev_v_l+j)%this->cvx.size()).x,
                            this->cvx.at_ccw((prev_v_l+j+1)%this->cvx.size()).x
                        );

                        if (external == 0) {
                            A = this->cvx.at_ccw((prev_v_l+j)%this->cvx.size());
                            B = this->cvx.at_ccw((prev_v_l+j+1)%this->cvx.size());
                            C = this->cvx.at_ccw(i);
                            break;
                        }
                    }
                    break;
                }
            }
            
            // find the optimal approximation line
            Line line = Line::line(A, B);
            float a = line.get_slope();
            float b = line.get_intercept();
            float c = C.y > A.x ? b + this->__distance(line, C)/2*sqrt(a*a+1) 
                                : b - this->__distance(line, C)/2*sqrt(a*a+1);
            
            return Line(a, c);
        }

        return Line(this->line->get_slope(), this->line->get_intercept());
    }

    LinearModel::LinearModel() {
        this->degree = 1;
        this->length = 0;
        this->line = nullptr;
    }

    Line* LinearModel::getLine() {
        return this->line;
    }

    float LinearModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else return (float) this->length / (2 + 4 + 4);
    }

    void LinearModel::clear() {
        this->length = 0;
        this->cvx.clear();

        if (this->line != nullptr) {
            delete this->line;
            this->line = nullptr;
        }
    }

    bool LinearModel::approximate(float bound, std::vector<Point2D>& segment) {
        if (this->length == 0) {
            for (Point2D& p : segment) {
                this->cvx.append(p);
                
                if (this->cvx.size() > 1) {
                    Line line = this->__approx(segment);
                    if (this->__cal_error(line, segment) > bound) return false;
                    else {
                        if (this->line != nullptr) delete this->line;
                        this->line = new Line(line.get_slope(), line.get_intercept());
                    }
                }
                this->length++;
            }
            
            return true;
        }
        else {
            Point2D& p = segment.back();
            this->cvx.append(p);

            if (this->cvx.size() > 1) {
                Line line = this->__approx(segment);
                if (this->__cal_error(line, segment) > bound) return false;
                else {
                    if (this->line != nullptr) delete this->line;
                    this->line = new Line(line.get_slope(), line.get_intercept());
                } 
            }

            this->length++;

            return true;
        }
    }

    float PolynomialModel::__cal_error(std::vector<Point2D>& segment) {
        float max_error = 0;
        for (Point2D& p : segment) {
            float error = std::abs(this->polynomial->subs(p.x) - p.y);
            max_error = error > max_error ? error : max_error;
        }
        return max_error;
    }

    PolynomialModel::PolynomialModel(int degree) {
        this->degree = degree;
        this->length = 0;   
        this->polynomial = nullptr;
    }

    Polynomial* PolynomialModel::getPolynomial() {
        return this->polynomial;
    }

    float PolynomialModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else return (float) this->length / (2 + 4 * (this->degree + 1));
    }

    void PolynomialModel::clear() {
        this->length = 0;
        if (this->polynomial != nullptr) {
            delete this->polynomial;
            this->polynomial = nullptr;
        }
    }

    bool PolynomialModel::approximate(float bound, std::vector<Point2D>& segment) {
        Eigen::VectorXd x(this->degree + 2);                    // decision variables
        Eigen::VectorXd c(this->degree + 2);                    // objective coefficients
        Eigen::MatrixXd A(2*segment.size(), this->degree + 2);  // constraint matrix
        Eigen::VectorXd b(2*segment.size());                    // constraint bound

        for (int i=0; i<this->degree+2; i++) {
            if (i != this->degree + 1) c(i) = 0.0;
            else c(i) = 1.0;
        }

        for (int i=0; i<segment.size(); i++) {
            for (int j=this->degree; j>=0; j--) {
                A(2*i, this->degree-j) = -pow(segment[i].x, j);
                A(2*i+1, this->degree-j) = pow(segment[i].x, j);
            }
            A(2*i, this->degree+1) = -1.0;
            A(2*i+1, this->degree+1) = -1.0;

            b(2*i) = -segment[i].y;
            b(2*i+1) = segment[i].y;
        }

        double minobj = sdlp::linprog(c, A, b, x);
        if (minobj == INFINITY || minobj == -INFINITY || minobj > bound) {
            return false;
        }
        else {
            float* coefficients = new float[this->degree+1];
            for (int i = 0; i <= this->degree; i++) {
                coefficients[this->degree-i] = x(i);
            }

            if (this->polynomial != nullptr) delete this->polynomial;
            this->polynomial = new Polynomial(this->degree, coefficients);

            this->length = segment.size();
            return true;
        }
    }

    // End: material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->max_degree = atoi(params[1]);

        for (int i=0; i<=this->max_degree; i++) {
            if (i == 0) this->models.push_back(new ConstantModel());
            else if (i == 1) this->models.push_back(new LinearModel());
            else this->models.push_back(new PolynomialModel(i));
        }
    }

    void Compression::finalize() {
        for (Point2D& p : this->window) {
            bool success = false;
            for (int k=this->curr_degree; k<=this->max_degree; k++) {
                if (this->__approximate(k)) {
                    this->curr_degree = k;
                    success = true;
                }
            }

            if (!success) {
                this->__choose_best_model();
                this->yield();

                this->curr_degree = 0;
                this->chosen_model = -1;
                this->window.erase(this->window.begin(), this->window.begin() + this->models[this->chosen_model]->length);
                for (int i=0; i<this->window.size(); i++) {
                    this->window[i].x = i;
                }
                for (int i=0; i<=this->max_degree; i++) {
                    this->models[i]->clear();
                    delete this->models[i];
                }
            }
        }

        if (this->window.size() > 0) {
            this->__choose_best_model();
            this->yield();
            for (int i=0; i<=this->max_degree; i++) {
                this->models[i]->clear();
                delete this->models[i];
            }
        }
    }

    bool Compression::__approximate(int degree) {
        bool success = true;

        if (degree == 0) {
            ConstantModel* constantModel = (ConstantModel*) this->models[degree];
            success = constantModel->approximate(this->error, this->window);
        }
        else if (degree == 1) {
            LinearModel* linearModel = (LinearModel*) this->models[degree];
            success = linearModel->approximate(this->error, this->window);
        }
        else {
            if (this->window.size() > degree) {
                PolynomialModel* polynomialModel = (PolynomialModel*) this->models[degree];
                success = polynomialModel->approximate(this->error, this->window);
            }
        }

        return success;
    }

    void Compression::__choose_best_model() {
        int index = 0;
        float max_gain = -INFINITY;

        for (int i=0; i<this->models.size(); i++) {
            float gain = this->models[i]->getCompressionRatio();
            if (gain > max_gain) {
                index = i;
                max_gain = gain;
            }
        }
        this->chosen_model = index;
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Model* model = this->models[this->chosen_model];

        if (model->degree == 0) {
            short degree_length = (unsigned short) model->length | (0 << 14);

            obj->put(degree_length);
            obj->put((float) ((ConstantModel*) model)->getValue());
        }
        else if (model->degree == 1) {
            Line* line = ((LinearModel*) model)->getLine();
            short degree_length = (unsigned short) model->length | (1 << 14);

            obj->put(degree_length);
            obj->put((float) line->get_intercept());
            obj->put((float) line->get_slope());
        }
        else {
            Polynomial* polynomial = ((PolynomialModel*) model)->getPolynomial();
            short degree_length = (unsigned short) model->length | (polynomial->degree << 14);
            obj->put(degree_length);
            for (int i = 0; i <= polynomial->degree; i++) {
                obj->put((float) polynomial->coefficients[i]);
            }
        }

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->window.size(), data->get_value());
        this->window.push_back(p);

        bool success = false;
        for (int k=this->curr_degree; k<=this->max_degree; k++) {
            if (this->__approximate(k)) {
                this->curr_degree = k;
                success = true;
                break;
            }
        }

        if (!success) {
            this->__choose_best_model();
            this->yield();

            this->window.erase(this->window.begin(), this->window.begin() + this->models[this->chosen_model]->length);
            for (int i=0; i<this->window.size(); i++) {
                this->window[i].x = i;
            }
            for (int i=0; i<=this->max_degree; i++) {
                this->models[i]->clear();
            }

            this->curr_degree = 0;
            this->chosen_model = -1;
        }
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        // Do nothing
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        unsigned short degree_length = compress_data->getShort();
        int degree = degree_length >> 14;
        unsigned short length = degree_length & (0xffff >> 2);

        float* coefficients = new float[degree+1];
        for (int i = 0; i <= degree; i++) {
            coefficients[i] = compress_data->getFloat();
        }

        Polynomial polynomial(degree, coefficients);
        for (int i = 0; i < length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(polynomial.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(polynomial.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        this->basetime += length * this->interval;

        delete[] coefficients;
        return base_obj;
    }
    // End: decompression
}