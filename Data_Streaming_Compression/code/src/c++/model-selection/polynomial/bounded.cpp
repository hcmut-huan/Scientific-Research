// #include "model-selection/polynomial.hpp"

// namespace Bounded {

//     int total = 0;
//     int count = 0;
//     Clock clock;

//     class Model {
//         public:
//             virtual Polynomial* get() = 0;
//             virtual int getDegree () = 0;        
//             virtual bool fit(float bound, std::vector<float>& data, int start) = 0;
//             virtual bool test(float bound, std::vector<float>& data, int start) = 0;
//     };

//     class LinearModel : public Model {
//         private:
//             UpperHull u_cvx; 
//             LowerHull l_cvx;
//             Line* u_line;
//             Line* l_line;
        
//         public:        
//             LinearModel() {
//                 this->u_line = nullptr;
//                 this->l_line = nullptr;
//             }

//             ~LinearModel() {
//                 if (this->l_line != nullptr) delete this->u_line;
//                 if (this->u_line != nullptr) delete this->l_line;
//                 this->u_cvx.clear();
//                 this->l_cvx.clear();
//             }

//             Polynomial* get() override {
//                 if (this->u_line != nullptr && this->l_line != nullptr) {
//                     float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
//                     float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;
                    
//                     float coefficients[2] = {intercept, slope};
//                     return new Polynomial(1, coefficients);
//                 }
                
//                 return nullptr;
//             }

//             int getDegree() override {
//                 return 1;
//             }

//             bool test(float bound, std::vector<float>& data, int start) override {
//                 return true;
//             }

//             bool fit(float bound, std::vector<float>& data, int start) override {
//                 Point2D p(data.size()-start-1, data.back());

//                 if (this->u_cvx.size() == 0 && this->u_cvx.size() == 0) {
//                     this->u_cvx.append(Point2D(p.x, p.y-bound));
//                     this->l_cvx.append(Point2D(p.x, p.y+bound));
//                 }
//                 else if (this->u_line == nullptr && this->l_line == nullptr) {
//                     Line u = Line::line(this->u_cvx.at(0), Point2D(p.x, p.y+bound));
//                     Line l = Line::line(this->l_cvx.at(0), Point2D(p.x, p.y-bound));

//                     u_cvx.append(Point2D(p.x, p.y-bound));
//                     l_cvx.append(Point2D(p.x, p.y+bound));
//                     this->u_line = new Line(u.get_slope(), u.get_intercept());
//                     this->l_line = new Line(l.get_slope(), l.get_intercept());
//                 }
//                 else {
//                     if (this->l_line->subs(p.x) > p.y + bound || p.y - bound > this->u_line->subs(p.x)) {
//                         return false;
//                     }
//                     else {
//                         bool update_u = p.y + bound < this->u_line->subs(p.x);
//                         bool update_l = p.y - bound > this->l_line->subs(p.x);
    
//                         if (update_u) {
//                             int index = 0;
//                             float min_slp = INFINITY;
    
//                             for (int i=0; i<this->u_cvx.size(); i++) {
//                                 Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y+bound));
//                                 if (line.get_slope() < min_slp) {
//                                     min_slp = line.get_slope();
//                                     index = i;
                                    
//                                     delete this->u_line;
//                                     this->u_line = new Line(line.get_slope(), line.get_intercept());
//                                 }
//                             }
//                             this->u_cvx.erase_from_begin(index);
//                         }
//                         if (update_l) {
//                             int index = 0;
//                             float max_slp = -INFINITY;
    
//                             for (int i=0; i<this->l_cvx.size(); i++) {
//                                 Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y-bound));
//                                 if (line.get_slope() > max_slp) {
//                                     max_slp = line.get_slope();
//                                     index = i;

//                                     delete this->l_line;
//                                     this->l_line = new Line(line.get_slope(), line.get_intercept());
//                                 }
//                             }
//                             this->l_cvx.erase_from_begin(index);
//                         }
    
//                         if (update_u) this->l_cvx.append(Point2D(p.x, p.y+bound));
//                         if (update_l) this->u_cvx.append(Point2D(p.x, p.y-bound));
//                     }
//                 }

//                 return true;
//             }
//     };

//     class PolyModel : public Model {
//         private:
//             int degree;
//             Polynomial* func;

//         public:
//             PolyModel(int degree) {
//                 this->degree = degree;
//                 this->func = nullptr;
//             }

//             ~PolyModel() {
//                 if (this->func != nullptr) delete this->func;
//             }

//             Polynomial* get() override {
//                 return this->func;
//             }

//             int getDegree() override {
//                 return this->degree;
//             }

//             bool test(float bound, std::vector<float>& data, int start) override {
//                 if (std::abs(this->func->subs(data.size()-start-1) - data.back()) <= bound) {
//                     return true;
//                 }
//                 else {
//                     return false;
//                 }
//             }

//             bool fit(float bound, std::vector<float>& data, int start) override {
//                 if (this->func != nullptr) {
//                     if (std::abs(this->func->subs(data.size()-start-1) - data.back()) <= bound) {
//                         return true;
//                     }  
//                 }

//                 Eigen::VectorXd x(degree + 2);                    
//                 Eigen::VectorXd c(degree + 2);                    
//                 Eigen::MatrixXd A(2*(data.size()-start), degree + 2);  
//                 Eigen::VectorXd b(2*(data.size()-start));              
        
//                 for (int i=0; i<degree+2; i++) {
//                     if (i != degree + 1) c(i) = 0.0;
//                     else c(i) = 1.0;
//                 }
                
//                 for (int i=0; i<data.size()-start; i++) {                                                              
//                     for (int j=degree; j>=0; j--) {
//                         A(2*i, degree-j) = -pow(i, j);
//                         A(2*i+1, degree-j) = pow(i, j);
//                     }
//                     A(2*i, degree+1) = -1.0;
//                     A(2*i+1, degree+1) = -1.0;
        
//                     b(2*i) = -data[i+start];
//                     b(2*i+1) = data[i+start];
//                 }
        
//                 double minobj = sdlp::linprog(c, A, b, x);
                
//                 if (minobj <= bound) {
//                     float* coefficients = new float[degree+1];
//                     for (int i = 0; i <= degree; i++) {
//                         coefficients[degree-i] = x(i);
//                     }
                    
//                     if (this->func != nullptr) delete this->func;
//                     this->func = new Polynomial(degree, coefficients);
//                     delete coefficients;
                    
//                     return true;
//                 }
                
//                 return false;
//             }
//     };

//     struct Segment {
//         int start;
//         short length;
//         int triedDegree;
//         int extreme;

//         int childLength;
//         int childCoeffByte;

//         bool isComplete;
//         Model* model;

//         Segment(int start, Model* model) {
//             this->start = start;
//             this->length = 0;
//             this->triedDegree = 1;
//             this->extreme = 0;

//             this->childLength = 0;
//             this->childCoeffByte = 0;

//             this->model = model; 
//             this->isComplete = false;
//         }

//         ~Segment() {
//             if (this->model != nullptr) delete this->model;
//         }

//         void updateChild(Segment* child) {
//             this->childLength += child->length;
//             this->childCoeffByte += 6 + child->model->getDegree() * 4;
//         }

//         bool trigger() {
//             return (this->triedDegree * 4 + 10) <
//                 (this->childCoeffByte + 6 + this->model->getDegree() * 4);
//         }

//         void updateModel(Model* model, short length) {
//             this->length = length;
//             this->triedDegree = model->getDegree();
//             this->childLength = 0;
//             this->childCoeffByte = 0;
//             this->isComplete = false;

//             delete this->model;
//             this->model = model;
//         }

//     };

//     void __yield(BinObj* obj, Model* model, short length) {
//         Polynomial* polynomial = model->get();
//         if (polynomial == nullptr) return;

//         short degree_length = length | ((model->getDegree() - 1) << 14) ;
//         obj->put(degree_length);
//         for (int i = 0; i <= model->getDegree(); i++) {
//             obj->put(polynomial->coefficients[i]);
//         }
//     }

//     Line translate(float slp, float intercept, int step) {
//         return Line(slp, intercept-slp*step);
//     }

//     int __check(Point2D& p1, Point2D& p2, Point2D& p3, Point2D& p4, Point2D& p5, Line& line_1, Line& line_2, float bound) {
//         if (p2.x <= p1.x || p2.x >= p3.x) return false;
//         else {
//             Eigen::MatrixXd A(3, 3);
//             A << p1.x*p1.x, p1.x, 1,
//                 p2.x*p2.x, p2.x, 1,
//                 p3.x*p3.x, p3.x, 1;
        
//             Eigen::VectorXd b(3);
//             b << p1.y, p2.y, p3.y;
        
//             Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

//             float extreme_x = (-x(1)) / (2*x(0));
//             float extreme_y = x(0)*extreme_x*extreme_x + x(1)*extreme_x + x(2);

//             if (line_1.get_slope() * line_2.get_slope() > 0) {
//                 if (extreme_x > p1.x && extreme_x < p3.x) return 0;
//                 else if (std::abs(x(0)*p4.x*p4.x + x(1)*p4.x + x(2) - p4.y) > bound
//                     || std::abs(x(0)*p5.x*p5.x + x(1)*p5.x + x(2) - p5.y) > bound) return 0;
//             }
//             else {
//                 if (extreme_x >= p2.x && std::abs(extreme_y - line_2.subs(extreme_x)) > bound) return 1; 
//                 else if (extreme_x < p2.x && std::abs(extreme_y - line_1.subs(extreme_x)) > bound) return 1;
//                 else if (std::abs(x(0)*p4.x*p4.x + x(1)*p4.x + x(2) - p4.y) > bound
//                     || std::abs(x(0)*p5.x*p5.x + x(1)*p5.x + x(2) - p5.y) > bound) return 1;
//             }
            
//             return (line_1.get_slope() * line_2.get_slope() > 0) ? 2 : 3;
//         }
//     }

//     void compress(TimeSeries& timeseries, int max_degree, float bound, std::string output) {
//         clock.start();
//         IterIO outputFile(output, false);
//         BinObj* compress_data = new BinObj;

//         Univariate* data = (Univariate*) timeseries.next();
//         time_t basetime = data->get_time();
//         compress_data->put(basetime);

//         std::vector<Segment*> segments = {new Segment(0, new LinearModel())};
//         std::vector<float> buffer = {data->get_value()};
//         segments[0]->model->fit(bound, buffer, 0);
//         segments[0]->length++;

//         int trigger_index = -1;
//         while (timeseries.hasNext()) {
//             buffer.push_back(((Univariate*) timeseries.next())->get_value());

//             for (int i=0; i<segments.size(); i++) {
//                 Segment* seg = segments[i];
//                 if(!seg->isComplete) {
//                     bool flag = (seg->model->getDegree() > 1)
//                         ? seg->model->test(bound, buffer, seg->start)
//                         : seg->model->fit(bound, buffer, seg->start);
                    
//                     if (!flag) {
//                         seg->isComplete = true;

//                         Segment* n_seg = new Segment(buffer.size()-1, new LinearModel());
//                         segments.push_back(n_seg);

//                         trigger_index = i;
//                         for (int j=0; j<i; j++) segments[j]->updateChild(seg);
//                     }
//                     else seg->length++;
//                 }
                
//                 if ((buffer.size() - seg->start > 16000) || 
//                     (seg->triedDegree == max_degree && seg->isComplete)) {
                    
//                     __yield(compress_data, seg->model, seg->length);
//                     for (int j=i+1; j<segments.size(); j++) 
//                         segments[j]->start -= seg->length;
//                     segments.erase(segments.begin() + i--);
//                     buffer.erase(buffer.begin(), buffer.begin() + seg->length);
//                 }
//             }

//             for (int i=0; i<trigger_index; i++) {
//                 Segment* seg = segments[i];
                
//                 int flag = -1;
//                 if (seg->model->getDegree() == 1 && segments[i+1]->model->getDegree() == 1) {
//                     Segment* s1 = segments[i];
//                     Segment* s2 = segments[i+1];
                    
//                     Line line_1(s1->model->get()->coefficients[1], s1->model->get()->coefficients[0]);
//                     Line line_2 = translate(s2->model->get()->coefficients[1], s2->model->get()->coefficients[0], s1->length);

//                     Point2D p1(0, line_1.subs(0));
//                     Point2D p2 = Line::intersection(line_1, line_2);
//                     Point2D p3(s1->length + s2->length - 1, line_2.subs(s1->length + s2->length - 1));
//                     Point2D p4(s1->length-1, line_1.subs(s1->length-1));
//                     Point2D p5(s1->length, line_2.subs(s2->start));
                    
//                     flag = __check(p1, p2, p3, p4, p5, line_1, line_2, bound);
//                     if (flag == 0) seg->triedDegree = max_degree;
//                     else if (flag == 1) seg->triedDegree = 2;
//                 }

//                 // Fit because of condition 1
//                 while (seg->trigger() && seg->triedDegree < max_degree) {
//                     PolyModel* n_model = new PolyModel(seg->triedDegree+1);
//                     if (n_model->getDegree() > 2) {
//                         int d = (seg->model->getDegree() % 2 == 0) ? -1 : 1;
//                         int a = (seg->model->get()->coefficients[seg->model->getDegree()] > 0) ? -1 : 1;
//                         int e = (seg->extreme % 2 == 0) ? -1 : 1;

//                         if (d * a * e == -1) {
//                             // slope decrease
//                             float p1 = seg->model->get()->subs(seg->length);
//                             float p2 = segments[i+1]->model->get()->subs(0);
//                             if (buffer[seg->start+seg->length-1] <= seg->model->get()->subs(seg->length-1) && p2 > p1) {
//                                 seg->triedDegree++;
//                                 continue;
//                             }
//                         } 
//                         else {
//                             // slope increase
//                             float p1 = seg->model->get()->subs(seg->length);
//                             float p2 = segments[i+1]->model->get()->subs(0);
//                             if (buffer[seg->start+seg->length-1] >= seg->model->get()->subs(seg->length-1) && p1 > p2) {
//                                 seg->triedDegree++;
//                                 continue;
//                             }   

//                         }
//                     }

//                     if (n_model->fit(bound, buffer, seg->start)) {
//                         if (flag == 3) seg->extreme++;
//                         seg->updateModel(n_model, buffer.size() - seg->start);
//                         segments.erase(segments.begin() + i + 1, segments.end());
//                         break;
//                     }
//                     else {
//                         seg->triedDegree++;
//                     }
//                 }

//                 // Fit success -> break
//                 if (!seg->isComplete) break;
//             }
//             trigger_index = -1;
//         }

//         for (Segment* seg : segments) {
//             __yield(compress_data, seg->model, seg->length);
//             delete seg;
//         }

//         outputFile.writeBin(compress_data);
//         outputFile.close();
//         delete compress_data;
        
//         clock.tick();
//         double avg_time = clock.getAvgDuration() / timeseries.size();

//         std::cout << total << " " << count << "\n";

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each data point (ns): " << avg_time << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each data point (ns): " + std::to_string(avg_time));
//         timeFile.close();
//     }

//     void __decompress_segment(IterIO& file, int interval, time_t basetime, int length, Polynomial& polynomial) {
//         for (int i = 0; i < length; i++) {
//             CSVObj obj;
//             obj.pushData(std::to_string(basetime + i * interval));
//             obj.pushData(std::to_string(polynomial.subs(i)));
//             file.write(&obj);
//         }
//     }
    
//     void decompress(std::string input, std::string output, int interval) {
//         IterIO inputFile(input, true, true);
//         IterIO outputFile(output, false);
//         BinObj* compress_data = inputFile.readBin();

//         time_t basetime = compress_data->getLong();
//         clock.start();
//         while (compress_data->getSize() != 0) {
//             unsigned short degree_length = compress_data->getShort();
//             int degree = (degree_length >> 14) + 1;
//             unsigned short length = degree_length & (0xffff >> 2);

//             float* coefficients = new float[degree+1];
//             for (int i = 0; i <= degree; i++) {
//                 coefficients[i] = compress_data->getFloat();
//             }

//             Polynomial polynomial(degree, coefficients);
//             __decompress_segment(outputFile, interval, basetime, length, polynomial);
//             delete[] coefficients;
            
//             basetime += length * interval;
//             clock.tick();
//         }

//         delete compress_data;
//         inputFile.close();
//         outputFile.close();

//         std::cout << total << " " << count << "\n";

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each segment (ns): " << clock.getAvgDuration() << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each segment (ns): " + std::to_string(clock.getAvgDuration()));
//         timeFile.close();
//     }
    
// };