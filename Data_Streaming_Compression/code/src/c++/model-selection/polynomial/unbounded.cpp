// #include "model-selection/polynomial.hpp"

// namespace Unbounded {

//     Clock clock;

//     // class LinearModel {
//     //     private:
//     //         Line* u_line;
//     //         Line* l_line;
//     //         Line* line;
//     //         Point2D* p1;

//     //         float A_num = 0;
//     //         float A_den = 0;

//     //     public:
//     //         bool isComplete;
//     //         short length;
            
//     //         LinearModel() {
//     //             this->length = 0;
//     //             this->line = nullptr;
//     //             this->u_line = nullptr;
//     //             this->l_line = nullptr;
//     //             this->isComplete = false;
//     //         }

//     //         ~LinearModel() {
//     //             if (this->p1 != nullptr) delete this->p1;
//     //             if (this->line != nullptr) delete this->line;
//     //             if (this->l_line != nullptr) delete this->u_line;
//     //             if (this->u_line != nullptr) delete this->l_line;
//     //         }

//     //         void translation(Point2D& p) {
//     //             Line* n_line = new Line(this->line->get_slope(), this->line->subs(p.x));
//     //             if (this->line != nullptr) delete this->line;

//     //             this->line = n_line; 
//     //             this->p1->x = 0;
//     //         }

//     //         Line* getLine() {
//     //             if (this->line == nullptr) {
//     //                 float A_ig = this->A_num / this->A_den;
//     //                 float temp = A_ig > this->u_line->get_slope() ? this->u_line->get_slope() : A_ig;
//     //                 double a_ig = temp > this->l_line->get_slope() ? temp : this->l_line->get_slope();
//     //                 double b_ig = this->p1->y - a_ig * this->p1->x;

//     //                 this->line = new Line(a_ig, b_ig);
//     //             }
            
//     //             return this->line;
//     //         }
            
//     //         void fit(float bound, Point2D& p)  {
//     //             if (this->length == 0) {
//     //                 this->p1 = new Point2D(p.x, p.y);
//     //             }
//     //             else if (this->length == 1) {
//     //                 Line l = Line::line(*p1, Point2D(p.x, p.y-bound));
//     //                 Line u = Line::line(*p1, Point2D(p.x, p.y+bound));

//     //                 this->u_line = new Line(u.get_slope(), u.get_intercept());
//     //                 this->l_line = new Line(l.get_slope(), l.get_intercept());
//     //             }
//     //             else {
//     //                 if (this->l_line->subs(p.x) > p.y + bound || p.y - bound > this->u_line->subs(p.x)) {
//     //                     this->isComplete = true;
//     //                     return;
//     //                 }
//     //                 else {
//     //                     if (p.y + bound < this->u_line->subs(p.x)) {
//     //                         Line l = Line::line(*p1, Point2D(p.x, p.y+bound));

//     //                         delete this->u_line;
//     //                         this->u_line = new Line(l.get_slope(), l.get_intercept());
//     //                     }
                        
//     //                     if (p.y - bound > this->l_line->subs(p.x)) {
//     //                         Line l = Line::line(*p1, Point2D(p.x, p.y-bound));

//     //                         delete this->l_line;
//     //                         this->l_line = new Line(l.get_slope(), l.get_intercept());
//     //                     }        
//     //                 }
//     //             }
                
//     //             this->length++;
//     //             this->A_num += (p.y - this->p1->y) * (p.x - this->p1->x);
//     //             this->A_den += (p.x - this->p1->x) * (p.x - this->p1->x);
//     //         }
//     // };

//     // Linear approximate    
//     class LinearModel {
//         private:
//             Line* line;
//             UpperHull u_cvx; 
//             LowerHull l_cvx;
//             Line* u_line;
//             Line* l_line;
        
//         public:
//             bool isComplete;
//             short length;
            
//             LinearModel() {
//                 this->length = 0;
//                 this->u_line = nullptr;
//                 this->l_line = nullptr;
//                 this->line = nullptr;
//                 this->isComplete = false;
//             }

//             ~LinearModel() {
//                 if (this->l_line != nullptr) delete this->u_line;
//                 if (this->u_line != nullptr) delete this->l_line;
//                 if (this->line != nullptr) delete this->line;
//                 this->u_cvx.clear();
//                 this->l_cvx.clear();
//             }

//             void translation(Point2D& p) {
//                 Line* n_line = new Line(this->line->get_slope(), this->line->subs(p.x));                
//                 if (this->line != nullptr) delete this->line;
//                 this->line = n_line;
//             }

//             Line* getLine() {
//                 if (this->line == nullptr) { 
//                     this->line = new Line(
//                         (this->u_line->get_slope() + this->l_line->get_slope()) / 2,
//                         (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2
//                     );
//                 }
                
//                 return this->line;
//             }

//             void fit(float bound, Point2D& p)  {
//                 if (this->length == 0) {
//                     this->u_cvx.append(Point2D(p.x, p.y-bound));
//                     this->l_cvx.append(Point2D(p.x, p.y+bound));
//                 }
//                 else if (this->length == 1) {
//                     Line u = Line::line(this->u_cvx.at(0), Point2D(p.x, p.y+bound));
//                     Line l = Line::line(this->l_cvx.at(0), Point2D(p.x, p.y-bound));

//                     u_cvx.append(Point2D(p.x, p.y-bound));
//                     l_cvx.append(Point2D(p.x, p.y+bound));
//                     this->u_line = new Line(u.get_slope(), u.get_intercept());
//                     this->l_line = new Line(l.get_slope(), l.get_intercept());
//                 }
//                 else {
//                     if (this->l_line->subs(p.x) > p.y + bound || p.y - bound > this->u_line->subs(p.x)) {
//                         this->isComplete = true;
//                         return;
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
                
//                 this->length++;
//             }
//     };

//     // Polynomial approximate
//     class PolynomialModel {
//         private:
//             Polynomial* polynomial;

//         public:
//             int degree;
//             short length;

//             PolynomialModel(int degree) {
//                 this->degree = degree;
//                 this->length = 0;   
//                 this->polynomial = nullptr;
//             }

//             ~PolynomialModel() {
//                 this->length = 0;
//                 delete this->polynomial;
//             }

//             Polynomial* getPolynomial() {
//                 return this->polynomial;
//             }

//             float getCompressionRatio() {
//                 if (this->length == 0) return -1;       // Ignore this model
//                 else return (float) this->length / (1 + 2 + 4 * (this->degree + 1));
//             }

//             void fit(float bound, std::vector<Point2D>& segment, int pivot = 0)  {                
//                 if (segment.size() - pivot < this->degree + 1) {
//                     return;
//                 }

//                 int n = segment.size() - pivot;
//                 Eigen::MatrixXd A(n, this->degree + 1);
//                 Eigen::VectorXd b(n);

//                 for (int i = 0; i < n; ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
//                     b(i) = segment[i].y;
//                     for (int j = 0; j <= this->degree; ++j) {
//                         A(i, j) = std::pow(segment[i].x, j);
//                     }
//                 }

//                 Eigen::VectorXd coefficients = (((A.transpose() * A).inverse()) * A.transpose()) * b;
//                 if (this->polynomial != nullptr) delete this->polynomial;

//                 float* coeffs = new float[this->degree+1];
//                 for (int i=0; i<=this->degree; i++) coeffs[i] = coefficients(i);
//                 this->polynomial = new Polynomial(this->degree, coeffs);
//                 this->length = n;

//                 delete coeffs;
//             }
//     };

//     void __yield(BinObj* obj, LinearModel* model) {
//         short degree_length = model->length;
//         Line* line = model->getLine();

//         obj->put(degree_length);
//         obj->put(line->get_slope());
//         obj->put(line->get_intercept());
//     }

//     void __yield(BinObj* obj, PolynomialModel* model) {
//         short degree_length = model->length | ((model->degree - 1) << 14);
//         Polynomial* polynomial = model->getPolynomial();

//         obj->put(degree_length);
//         for (int i = 0; i <= model->degree; i++) {
//             obj->put(polynomial->coefficients[i]);
//         }
//     }

//     bool __check(Line* line_1, Line* line_2, float bound, Point2D& p1, Point2D& p2, Point2D& p3, Point2D& p4, Point2D& p5) {
//         Eigen::MatrixXd A(3, 3);
//         A << p1.x*p1.x, p1.x, 1,
//              p2.x*p2.x, p2.x, 1,
//              p3.x*p3.x, p3.x, 1;
    
//         Eigen::VectorXd b(3);
//         b << p1.y, p2.y, p3.y;
    
//         Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

//         float extreme_x = (-x(1)) / (2*x(0));
//         float extreme_y = x(0)*extreme_x*extreme_x + x(1)*extreme_x + x(2);

//         if (line_1->get_slope() * line_2->get_slope() > 0) {
//             if (extreme_x > p1.x && extreme_x < p3.x) return true;
//         }
//         else {  
//             if (extreme_x >= p2.x && std::abs(extreme_y - p2.y) > bound) return true; 
//             else if (extreme_x < p2.x && std::abs(extreme_y - line_1->subs(extreme_x)) > bound) return true;
//         }
        
//         return false;
//     }
    
//     void compress(TimeSeries& timeseries, float bound, std::string output) {
//         clock.start();
//         IterIO outputFile(output, false);
//         BinObj* compress_data = new BinObj;

//         Univariate* d = (Univariate*) timeseries.next();
//         compress_data->put(d->get_time());
//         std::vector<Point2D> segment = {Point2D(0, d->get_value())};
//         Point2D p1(-1, -1); Point2D p2(-1, -1); Point2D p3(-1, -1);
        
//         int phase = 1;
//         LinearModel* l_1 = new LinearModel(); l_1->fit(bound, p1);
//         LinearModel* l_2 = new LinearModel();

//         int index = 1;
//         bool flag = false;
//         int degree = 1;     // current degree of polynomial
//         int direction = 0;  // -1 decrease, 1 increase
//         while (timeseries.hasNext()) {
//             Point2D p(segment.size(), ((Univariate*) timeseries.next())->get_value());

//             if (phase == 1) {
//                 if (!l_1->isComplete) l_1->fit(bound, p);

//                 if (l_1->isComplete) {
//                     phase = 2;
//                     p1 = Point2D(0, l_1->getLine()->subs(0));
//                     p2 = Point2D(p.x, p.y);
//                     l_2->fit(bound, p);
//                 }
//             }
//             else if (phase == 2) {
//                 if (!l_2->isComplete) l_2->fit(bound, p);
//                 if (l_2->isComplete) {
//                     p3 = Point2D(
//                         p2.x + l_2->length - 1,
//                         l_2->getLine()->subs(p2.x + l_2->length - 1)
//                     );
//                 }
//             }

//             if (l_1->isComplete && l_2->isComplete) {
//                 Line* line_1 = l_1->getLine();
//                 Line* line_2 = l_2->getLine();

//                 int n_direction = line_1->get_slope() > line_2->get_slope() ? -1 : 1;
//                 Point2D intersection = Line::intersection(*line_1, *line_2);
//                 Point2D e1(l_1->length-1, line_1->subs(l_1->length-1));
//                 Point2D s2(l_1->length, line_2->subs(l_1->length));

//                 if (segment.size() > 16000 || intersection.x <= p1.x || intersection.x >= p3.x) flag = true;
//                 else {
//                     flag = __check(
//                         line_1, line_2, bound, p1, intersection, p3, e1, s2
//                     );
//                 }
//                 if (!flag) {
//                     if (direction != n_direction) {
//                         direction = n_direction;
//                         degree += 1;
//                     }
//                     if (degree > 4) {
//                         flag = true;
//                         degree = 4;
//                     }
//                     else {
//                         delete l_1; l_1 = l_2;
//                         l_2 = new LinearModel();
//                         l_2->fit(bound, p);

//                         p1 = Point2D(p2.x, l_1->getLine()->subs(p2.x));
//                         p2 = Point2D(p.x, p.y);
//                     }
//                 }
                 
//                 if (flag) {
//                     if (degree == 1) {
//                         __yield(compress_data, l_1);
//                     }
//                     else {
//                         PolynomialModel* polynomialModel = new PolynomialModel(degree);
//                         polynomialModel->fit(bound, segment, l_2->length);
                        
//                         __yield(compress_data, polynomialModel);
//                     }

//                     segment = { segment.end() - l_2->length, segment.end() };
//                     for (int i=0; i<segment.size(); i++) segment[i].x = i; 
//                     delete l_1; l_1 = l_2; l_1->translation(p2);

//                     phase = 2;
//                     p1 = Point2D(0, l_1->getLine()->subs(0));
//                     p2 = Point2D(segment.size(), p.y);
//                     p = Point2D(segment.size(), p.y);

//                     l_2 = new LinearModel(); l_2->fit(bound, p);

//                     degree = 1; 
//                     direction = 0;
//                     flag = false;
//                 }
//             }

//             index++;
//             segment.push_back(p);
//         }

//         delete l_1; delete l_2;

//         outputFile.writeBin(compress_data);
//         outputFile.close();
//         delete compress_data;
        
//         clock.tick();
//         double avg_time = clock.getAvgDuration() / timeseries.size();

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each data point (ns): " << avg_time << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each data point (ns): " + std::to_string(avg_time));
//         timeFile.close();
//     }

//     void __decompress_segment(IterIO& file, int interval, time_t basetime, int length, Line& line) {
//         for (int i = 0; i < length; i++) {
//             CSVObj obj;
//             obj.pushData(std::to_string(basetime + i * interval));
//             obj.pushData(std::to_string(line.subs(i)));
//             file.write(&obj);
//         }
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

//             if (degree == 1) {
//                 float slope = compress_data->getFloat();
//                 float intercept = compress_data->getFloat();
                
//                 Line line(slope, intercept);
//                 __decompress_segment(outputFile, interval, basetime, length, line);
//             }
//             else {
//                 float* coefficients = new float[degree+1];
//                 for (int i = 0; i <= degree; i++) {
//                     coefficients[i] = compress_data->getFloat();
//                 }

//                 Polynomial polynomial(degree, coefficients);
//                 __decompress_segment(outputFile, interval, basetime, length, polynomial);
//                 delete[] coefficients;
//             }
            
//             basetime += length * interval;
//             clock.tick();
//         }

//         delete compress_data;
//         inputFile.close();
//         outputFile.close();

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each segment (ns): " << clock.getAvgDuration() << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each segment (ns): " + std::to_string(clock.getAvgDuration()));
//         timeFile.close();
//     }
    
// };