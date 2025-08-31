#include "piecewise-approximation/linear.hpp"

namespace SlideFilter {

    // Begin: compression
    Line Compression::__fit() {
        double slope = 0;
        double intercept = 0;

        if (this->u_line->get_slope() == this->l_line->get_slope()) {
            slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;    
            intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;
        }
        else {
            double A_num = 0;
            double A_den = 0;
            Point2D p0 = Line::intersection(*this->u_line, *this->l_line);
            
            for (Point2D& p : this->segments) {
                A_num += (p.y-p0.y)*(p.x-p0.x);
                A_den += (p.x-p0.x)*(p.x-p0.x);
            }

            double A_ig = this->A_num / this->A_den;
            double temp = A_ig > this->u_line->get_slope() ? this->u_line->get_slope() : A_ig;
            
            slope = temp > this->l_line->get_slope() ? temp : this->l_line->get_slope();
            intercept = p0.y - slope * p0.x;
        }
        
        return Line(slope, intercept);
    }
    
    double Compression::__checkConnected() {
        double t_j_k = this->pivot->x;
        Point2D prev = Line::intersection(*this->prev_l, *this->prev_u);
        if (std::abs(this->u_line->get_slope() - this->l_line->get_slope()) < 0.0000001) {
            Point2D p = Line::intersection(*this->u_line, *this->prev_g);
            
            if (p.x < t_j_k && p.x > prev.x 
                && this->prev_l->subs(p.x) > p.y 
                && p.y > this->prev_u->subs(p.x)) {
                    
                return p.x;
            }
        }
        else { 
            Point2D p = Line::intersection(*this->u_line, *this->l_line);
            if (this->prev_g->subs(p.x)>p.y) {
                // point below g^{k-1}
                double f_i_k = Line::intersection(*this->l_line, *this->prev_g).x;
                if (f_i_k<t_j_k && this->l_line->subs(t_j_k-1)>prev_l->subs(t_j_k-1)) {
                    // find hyperplane s
                    Line s = Line::line(Point2D(t_j_k-1, this->prev_l->subs(t_j_k-1)), p);

                    double c_i = Line::intersection(*this->prev_g, *this->u_line).x;     // c_i^k
                    double d_i = Line::intersection(*this->prev_g, s).x;     // d_i^k
                    double max_c_vs_d = (c_i > d_i) ? c_i : d_i;

                    if ((c_i - f_i_k > 0.0000001 && d_i - f_i_k > 0.0000001) || max_c_vs_d - p.x > 0.0000001 || f_i_k - p.x > 0.0000001) return -1;
                    else if (c_i - f_i_k > 0.0000001 && d_i - prev.x > 0.0000001) return (d_i + f_i_k) / 2;
                    else if (d_i - f_i_k > 0.0000001 && c_i - prev.x > 0.0000001) return (c_i + f_i_k) / 2;
                    else return (max_c_vs_d + f_i_k) / 2;            
                }
            }
            else if (p.y>this->prev_g->subs(p.x)) {
                // point above g^{k-1}
                double f_i_k = Line::intersection(*this->u_line, *this->prev_g).x;                          
                if (f_i_k<t_j_k && this->prev_u->subs(t_j_k-1)>this->u_line->subs(t_j_k-1)) {
                    // find hyperplane q
                    Line q = Line::line(Point2D(t_j_k-1, this->prev_u->subs(t_j_k-1)), p);

                    double c_i = Line::intersection(*this->prev_g, *this->l_line).x;   // c_i^k'
                    double d_i = Line::intersection(*this->prev_g, q).x;   // d_i^k'
                    double max_c_vs_d = (c_i > d_i) ? c_i : d_i;

                    if ((c_i - f_i_k > 0.0000001 && d_i - f_i_k > 0.0000001) || max_c_vs_d - p.x > 0.0000001 || f_i_k - p.x > 0.0000001) return -1;
                    else if (c_i - f_i_k > 0.0000001 && d_i - prev.x > 0.0000001) return (d_i + f_i_k) / 2;
                    else if (d_i - f_i_k > 0.0000001 && c_i - prev.x > 0.0000001) return (c_i + f_i_k) / 2;
                    else return (max_c_vs_d + f_i_k) / 2;         
                }
            }
        }

        return -1;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        this->yield();
        this->segment_pos = 2;
        this->yield();

        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        if (this->prev_end != nullptr) delete this->prev_end;
        if (this->prev_g != nullptr) delete this->prev_g;
        if (this->prev_u != nullptr) delete this->prev_u;
        if (this->prev_l != nullptr) delete this->prev_l;
    }

    BinObj* Compression::serialize() {
        if (this->u_line == nullptr) return nullptr;
        BinObj* obj = new BinObj;

        if (this->segment_pos == 0) {           
            Line g_k = this->__fit(); 
            this->prev_g = new Line(g_k.get_slope(), g_k.get_intercept());
            this->prev_end = new Point2D(0, this->prev_g->subs(0));
            
            obj->put((short) this->prev_end->x);
            obj->put((float) this->prev_end->y);
        }
        else if (this->segment_pos == 2) {
            obj->put((float) (this->index - 1 - this->prev_end->x));
            obj->put((float) this->prev_g->subs(this->index-1));
        }
        else {
            double connected_endpoint = this->__checkConnected();
            if (connected_endpoint == -1) {
                Line g_k = this->__fit(); 
            
                obj->put((float) (this->prev_end->x - this->pivot->x));
                obj->put((float) this->prev_g->subs(this->pivot->x));
                obj->put((float) g_k.subs(this->pivot->x));

                this->prev_end->x = this->pivot->x;
                this->prev_end->y = g_k.subs(this->pivot->x);

                delete this->prev_g;
                this->prev_g = new Line(g_k.get_slope(), g_k.get_intercept());
            }
            else {
                Point2D p = Line::intersection(*this->u_line, *this->l_line);
                Point2D curr_end(connected_endpoint, this->prev_g->subs(connected_endpoint));
                Line g_k = Line::line(p, curr_end);

                obj->put((float) (connected_endpoint - this->prev_end->x));
                obj->put((float) curr_end.y);

                this->prev_end->x = curr_end.x;
                this->prev_end->y = curr_end.y;

                delete this->prev_g;
                this->prev_g = new Line(g_k.get_slope(), g_k.get_intercept());
            }
        }

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->index++, data->get_value());

        if (this->pivot == nullptr) {
            this->pivot = new Point2D(p.x, p.y);
            this->cvx.append(p);
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y-this->error), 
                Point2D(p.x, p.y+this->error)
            );
            Line l_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y+this->error), 
                Point2D(p.x, p.y-this->error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->cvx.append(p);
        }
        else {
            if (p.y-this->u_line->subs(p.x)>this->error || this->l_line->subs(p.x)-p.y>this->error) {
                this->yield();
                if (this->segment_pos == 0) this->segment_pos = 1;
                else {
                    delete this->prev_u;
                    delete this->prev_l;
                }

                this->prev_u = this->u_line; this->u_line = nullptr;
                this->prev_l = this->l_line; this->l_line = nullptr;
                this->cvx.clear();

                this->pivot->x = p.x; this->pivot->y = p.y;
                this->cvx.append(p);
                this->segments.clear();
                this->segments.push_back(p);
            }
            else {
                cvx.append(p);
                if (p.y-this->l_line->subs(p.x)>this->error) {
                    for (int i=0; i<cvx.size(); i++) {
                        Point2D cvx_p = cvx.at(i);
                        if (cvx_p.x == p.x) continue;

                        Line l = Line::line(Point2D(cvx_p.x, cvx_p.y+this->error), Point2D(p.x, p.y-this->error));
                        if (l.get_slope() > this->l_line->get_slope()) {
                            delete this->l_line;
                            this->l_line = new Line(l.get_slope(), l.get_intercept());
                        }
                    }
                }
                
                if (this->u_line->subs(p.x)-p.y>this->error) {
                    for (int i=0; i<cvx.size(); i++) {
                        Point2D cvx_p = cvx.at(i);
                        if (cvx_p.x == p.x) continue;

                        Line l = Line::line(Point2D(cvx_p.x, cvx_p.y-this->error), Point2D(p.x, p.y+this->error));
                        if (l.get_slope() < this->u_line->get_slope()) {
                            delete this->u_line;
                            this->u_line = new Line(l.get_slope(), l.get_intercept());
                        }
                    }
                }
            }
        }

        this->segments.push_back(p);
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        // Do nothing
    }

    void Decompression::finalize() {
        if (this->prev_end != nullptr) delete this->prev_end;
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        if (this->prev_end == nullptr) {
            unsigned short start = compress_data->getShort();
            float value = compress_data->getFloat();
            this->prev_end = new Point2D(start, value);

            return nullptr;
        }

        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        float delta = compress_data->getFloat();
        if (delta > 0) {
            float value = compress_data->getFloat();

            Point2D* curr_end = new Point2D(this->prev_end->x + delta, value);
            Line line = Line::line(*curr_end, *this->prev_end);

            delete this->prev_end;
            this->prev_end = curr_end;

            for (long i=this->index; curr_end->x - i > 0.0000001 || std::abs(curr_end->x - i) < 0.0000001; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * interval));
                    base_obj->pushData(std::to_string(line.subs(i)));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * interval));
                    obj->pushData(std::to_string(line.subs(i)));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }

            this->index = ((long) curr_end->x) + 1;
        }
        else {
            float val_1 = compress_data->getFloat();
            float val_2 = compress_data->getFloat();
            
            Point2D* curr_end = new Point2D(this->prev_end->x - delta, val_1);
            Line line = Line::line(*curr_end, *this->prev_end);

            delete this->prev_end;
            this->prev_end = curr_end;
            this->prev_end->y = val_2;

            for (long i=this->index; i<std::lround(curr_end->x); i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * interval));
                    base_obj->pushData(std::to_string(line.subs(i)));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * interval));
                    obj->pushData(std::to_string(line.subs(i)));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }

            this->index = std::lround(curr_end->x);
        }

        
        return base_obj;
    }
    // End: decompression
};

