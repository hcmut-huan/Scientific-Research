#include "model-selection/polynomial.hpp"

namespace Bounded {

    int count = 0;
    int total = 0;
    int e_count = 0;

    void permutate(int level, std::vector<std::vector<double>>& result,
        const std::vector<std::vector<double>>& arrays, std::vector<double>& current) {
        
        int n = arrays.size();
        if (level == n) {
            // Base case: one permutation is complete
            std::vector<double> permutation;
            for (double val : current)
                permutation.push_back(val);
            result.push_back(permutation);
            return;
        }
    
        // Iterate over each element in the current array
        for (int i = 0; i < arrays[level].size(); ++i) {
            current[level] = arrays[level][i];
            permutate(level + 1, result, arrays, current);
        }
    }

    // Begin: material
    LinearModel::LinearModel() {
        this->u_line = nullptr;
        this->l_line = nullptr;
    }

    LinearModel::~LinearModel() {
        if (this->l_line != nullptr) delete this->u_line;
        if (this->u_line != nullptr) delete this->l_line;
        this->u_cvx.clear();
        this->l_cvx.clear();
    }

    Polynomial* LinearModel::get() {
        if (this->u_line != nullptr && this->l_line != nullptr) {
            float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
            float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;
            
            float coefficients[2] = {intercept, slope};
            return new Polynomial(1, coefficients);
        }
        
        return nullptr;
    }

    int LinearModel::getDegree() {
        return 1;
    }

    bool LinearModel::test(float bound, std::vector<Point2D>& data, int start) {
        return true;
    }

    bool LinearModel::fit(float bound, std::vector<Point2D>& data, int start) {
        Point2D p(data.back().x - start - 1, data.back().y);

        if (this->u_cvx.size() == 0 && this->u_cvx.size() == 0) {
            this->u_cvx.append(Point2D(p.x, p.y-bound));
            this->l_cvx.append(Point2D(p.x, p.y+bound));
        }
        else if (this->u_line == nullptr && this->l_line == nullptr) {
            Line u = Line::line(this->u_cvx.at(0), Point2D(p.x, p.y+bound));
            Line l = Line::line(this->l_cvx.at(0), Point2D(p.x, p.y-bound));

            u_cvx.append(Point2D(p.x, p.y-bound));
            l_cvx.append(Point2D(p.x, p.y+bound));
            this->u_line = new Line(u.get_slope(), u.get_intercept());
            this->l_line = new Line(l.get_slope(), l.get_intercept());
        }
        else {
            if (this->l_line->subs(p.x) > p.y + bound || p.y - bound > this->u_line->subs(p.x)) {
                return false;
            }
            else {
                bool update_u = p.y + bound < this->u_line->subs(p.x);
                bool update_l = p.y - bound > this->l_line->subs(p.x);

                if (update_u) {
                    int index = 0;
                    float min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y+bound));
                        if (line.get_slope() < min_slp) {
                            min_slp = line.get_slope();
                            index = i;
                            
                            delete this->u_line;
                            this->u_line = new Line(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->u_cvx.erase_from_begin(index);
                }
                if (update_l) {
                    int index = 0;
                    float max_slp = -INFINITY;

                    for (int i=0; i<this->l_cvx.size(); i++) {
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y-bound));
                        if (line.get_slope() > max_slp) {
                            max_slp = line.get_slope();
                            index = i;

                            delete this->l_line;
                            this->l_line = new Line(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y+bound));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y-bound));
            }
        }

        return true;
    }

    PolyModel::PolyModel(int degree) {
        this->degree = degree;
        this->func = nullptr;
    }

    PolyModel::~PolyModel() {
        if (this->func != nullptr) delete this->func;
    }

    Polynomial* PolyModel::get() {
        return this->func;
    }

    int PolyModel::getDegree() {
        return this->degree;
    }

    bool PolyModel::test(float bound, std::vector<Point2D>& data, int start) {
        if (std::abs(this->func->subs(data.back().x - start - 1) - data.back().y) <= bound) {
            return true;
        }
        else {
            return false;
        }
    }

    bool PolyModel::fit(float bound, std::vector<Point2D>& data, int start) {
        if (this->func != nullptr) {
            if (std::abs(this->func->subs(data.back().x - start - 1) - data.back().y) <= bound) {
                return true;
            }  
        }

        Eigen::VectorXd x(degree + 2);                    
        Eigen::VectorXd c(degree + 2);                    
        Eigen::MatrixXd A(2*(data.size()-start), degree + 2);  
        Eigen::VectorXd b(2*(data.size()-start));              

        for (int i=0; i<degree+2; i++) {
            if (i != degree + 1) c(i) = 0.0;
            else c(i) = 1.0;
        }
        
        for (int i=0; i<data.size()-start; i++) {                                                              
            for (int j=degree; j>=0; j--) {
                A(2*i, degree-j) = -pow(i, j);
                A(2*i+1, degree-j) = pow(i, j);
            }
            A(2*i, degree+1) = -1.0;
            A(2*i+1, degree+1) = -1.0;

            b(2*i) = -data[i+start].y;
            b(2*i+1) = data[i+start].y;
        }

        double minobj = sdlp::linprog(c, A, b, x);
        
        if (minobj <= bound) {
            float* coefficients = new float[degree+1];
            for (int i = 0; i <= degree; i++) {
                coefficients[degree-i] = x(i);
            }

            std::vector<std::vector<double>> vecs(degree+1, std::vector<double>(3));
            for (int i=0; i<=degree; i++) {
                vecs[i][0] = std::ceil(coefficients[degree-i]) - coefficients[degree-i];
                vecs[i][1] = coefficients[degree-i] - std::floor(coefficients[degree-i]);
                vecs[i][2] = 0;
            }

            std::vector<std::vector<double>> permutations;
            std::vector<double> current(degree + 1);
            permutate(0, permutations, vecs, current);
            
            bool flag = false;
            double remain = bound - minobj;
            Polynomial* base_poly = new Polynomial(degree, coefficients);
            Polynomial* poly = nullptr;

            for (int i=0; i<permutations.size(); i++) {
                for (int j=0; j<=degree; j++) {
                    coefficients[degree-j] = permutations[i][j];
                }

                if (i == permutations.size() - 1) break;

                poly = new Polynomial(degree, coefficients);
                double x1 = std::abs(poly->subs(0));
                double x2 = std::abs(poly->subs(1));
                double x3 = std::abs(poly->subs(data.size()-start-1));

                if (x1 <= remain && x2 <= remain && x3 <= remain) {
                    flag = true;
                    break;
                }
                delete poly;
            }

            total++;
            if (flag) count++;

            if (this->func != nullptr) delete this->func;
            this->func = new Polynomial(degree, coefficients);
            delete coefficients;
            
            return true;
        }
        
        return false;
    }

    Segment::Segment(int start, Model* model) {
        this->start = start;
        this->length = 0;
        this->triedDegree = 1;
        this->extreme = 0;

        this->childLength = 0;
        this->childCoeffByte = 0;

        this->model = model; 
        this->isComplete = false;
    }

    Segment::~Segment() {
        if (this->model != nullptr) delete this->model;
    }

    void Segment::updateChild(Segment* child) {
        this->childLength += child->length;
        this->childCoeffByte += 6 + child->model->getDegree() * 4;
    }

    bool Segment::trigger() {
        return (this->triedDegree * 4 + 10) <
            (this->childCoeffByte + 6 + this->model->getDegree() * 4);
    }

    void Segment::updateModel(Model* model, short length) {
        this->length = length;
        this->triedDegree = model->getDegree();
        this->childLength = 0;
        this->childCoeffByte = 0;
        this->isComplete = false;

        delete this->model;
        this->model = model;
    }
    // End: material

    // Begin: compression
    Line Compression::__translate(float slp, float intercept, int step) {
        return Line(slp, intercept-slp*step);
    }

    int Compression::__check(Point2D& p1, Point2D& p2, Point2D& p3, Point2D& p4, Point2D& p5, Line& line_1, Line& line_2) {
        if (p2.x <= p1.x || p2.x >= p3.x) return false;
        else {
            Eigen::MatrixXd A(3, 3);
            A << p1.x*p1.x, p1.x, 1,
                p2.x*p2.x, p2.x, 1,
                p3.x*p3.x, p3.x, 1;
        
            Eigen::VectorXd b(3);
            b << p1.y, p2.y, p3.y;
        
            Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

            float extreme_x = (-x(1)) / (2*x(0));
            float extreme_y = x(0)*extreme_x*extreme_x + x(1)*extreme_x + x(2);

            if (line_1.get_slope() * line_2.get_slope() > 0) {
                if (extreme_x > p1.x && extreme_x < p3.x) return 0;
                else if (std::abs(x(0)*p4.x*p4.x + x(1)*p4.x + x(2) - p4.y) > this->error
                    || std::abs(x(0)*p5.x*p5.x + x(1)*p5.x + x(2) - p5.y) > this->error) return 0;
            }
            else {
                if (extreme_x >= p2.x && std::abs(extreme_y - line_2.subs(extreme_x)) > this->error) return 1; 
                else if (extreme_x < p2.x && std::abs(extreme_y - line_1.subs(extreme_x)) > this->error) return 1;
                else if (std::abs(x(0)*p4.x*p4.x + x(1)*p4.x + x(2) - p4.y) > this->error
                    || std::abs(x(0)*p5.x*p5.x + x(1)*p5.x + x(2) - p5.y) > this->error) return 1;
            }
            
            return (line_1.get_slope() * line_2.get_slope() > 0) ? 2 : 3;
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->max_degree = atoi(params[1]);

        this->segments = {new Segment(0, new LinearModel())};
    }

    void Compression::finalize() {
        for (int i=0; i<this->segments.size(); i++) {
            Segment* seg = this->segments[i];
            this->model_index = i;
            this->yield();

            delete seg;
        }

        std::cout << "Total: " << total << "\n";
        std::cout << "Count: " << count << "\n";
        std::cout << "Error: " << e_count << "\n";
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;

        Segment* seg = this->segments[this->model_index];
        Polynomial* polynomial = seg->model->get();

        short degree_length = seg->length | (seg->model->getDegree() << 14);
        obj->put((short) degree_length);
        for (int i = 0; i <= seg->model->getDegree(); i++) {
            obj->put((float) polynomial->coefficients[i]);
        }

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->window.size() + 1, data->get_value());
        this->window.push_back(p);

        for (int i=0; i<this->segments.size(); i++) {
            Segment* seg = this->segments[i];
            if(!seg->isComplete) {
                bool flag = (seg->model->getDegree() > 1)
                    ? seg->model->test(this->error, this->window, seg->start)
                    : seg->model->fit(this->error, this->window, seg->start);
                
                if (!flag) {
                    seg->isComplete = true;

                    Segment* n_seg = new Segment(this->window.size()-1, new LinearModel());
                    this->segments.push_back(n_seg);

                    this->trigger_index = i;
                    for (int j=0; j<i; j++) this->segments[j]->updateChild(seg);
                }
                else seg->length++;
            }
            
            if ((this->window.size() - seg->start > 16000) || 
                (seg->triedDegree == this->max_degree && seg->isComplete)) {
                this->model_index = i;
                this->yield();

                for (int j=i+1; j<this->segments.size(); j++) 
                    this->segments[j]->start -= seg->length;
                for (int j=seg->length; j<this->window.size(); j++)
                    this->window[j].x -= seg->length;
                this->segments.erase(segments.begin() + i--);
                this->window.erase(this->window.begin(), this->window.begin() + seg->length);

                delete seg;
            }
        }

        for (int i=0; i<this->trigger_index; i++) {
            Segment* seg = this->segments[i];
            
            int flag = -1;
            if (seg->model->getDegree() == 1 && this->segments[i+1]->model->getDegree() == 1) {
                Segment* s1 = this->segments[i];
                Segment* s2 = this->segments[i+1];
                
                Line line_1(s1->model->get()->coefficients[1], s1->model->get()->coefficients[0]);
                Line line_2 = this->__translate(s2->model->get()->coefficients[1], s2->model->get()->coefficients[0], s1->length);

                Point2D p1(0, line_1.subs(0));
                Point2D p2 = Line::intersection(line_1, line_2);
                Point2D p3(s1->length + s2->length - 1, line_2.subs(s1->length + s2->length - 1));
                Point2D p4(s1->length-1, line_1.subs(s1->length-1));
                Point2D p5(s1->length, line_2.subs(s2->start));
                
                flag = this->__check(p1, p2, p3, p4, p5, line_1, line_2);
                if (flag == 0) seg->triedDegree = this->max_degree;
                else if (flag == 1) seg->triedDegree = 2;
            }

            while (seg->trigger() && seg->triedDegree < this->max_degree) {
                PolyModel* n_model = new PolyModel(seg->triedDegree+1);
                if (n_model->getDegree() > 2) {
                    int d = (seg->model->getDegree() % 2 == 0) ? -1 : 1;
                    int a = (seg->model->get()->coefficients[seg->model->getDegree()] > 0) ? -1 : 1;
                    int e = (seg->extreme % 2 == 0) ? -1 : 1;

                    if (d * a * e == -1) {
                        // slope decrease
                        float p1 = seg->model->get()->subs(seg->length);
                        float p2 = this->segments[i+1]->model->get()->subs(0);
                        if (this->window[seg->start+seg->length-1].y <= seg->model->get()->subs(seg->length-1) && p2 > p1) {
                            seg->triedDegree++;
                            continue;
                        }
                    } 
                    else {
                        // slope increase
                        float p1 = seg->model->get()->subs(seg->length);
                        float p2 = this->segments[i+1]->model->get()->subs(0);
                        if (this->window[seg->start+seg->length-1].y >= seg->model->get()->subs(seg->length-1) && p1 > p2) {
                            seg->triedDegree++;
                            continue;
                        }   

                    }
                }

                if (n_model->fit(this->error, this->window, seg->start)) {
                    if (flag == 3) seg->extreme++;
                    seg->updateModel(n_model, this->window.size() - seg->start);
                    this->segments.erase(this->segments.begin() + i + 1, this->segments.end());
                    break;
                }
                else {
                    seg->triedDegree++;
                }
            }

            // Fit success -> break
            if (!seg->isComplete) break;
        }
        this->trigger_index = -1;
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
};