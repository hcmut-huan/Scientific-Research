#include "piecewise-approximation/linear.hpp"

namespace ConnIPLA {
    // Begin: material
    void LinearSegment::__optimal_approximate(Point2D& p, double error) {
        if (this->u_line == nullptr) {
            Point2D u_p = this->u_cvx.at(0);
            Line u_line = Line::line(
                Point2D(u_p.x, u_p.y), 
                Point2D(p.x, p.y + error)
            );
            Point2D l_p = this->l_cvx.at(0);
            Line l_line = Line::line(
                Point2D(l_p.x, l_p.y), 
                Point2D(p.x, p.y - error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
        }
        else {
            if (this->l_line->subs(p.x) > p.y + error || p.y - error > this->u_line->subs(p.x)) {
                this->is_complete = true;
            }
            else {
                bool update_u = p.y + error < this->u_line->subs(p.x);
                bool update_l = p.y - error > this->l_line->subs(p.x);
                if (update_u) {
                    int index = 0;
                    double min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + error));
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
                    double max_slp = -INFINITY;

                    for (int i=0; i<this->l_cvx.size(); i++) {
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - error));
                        if (line.get_slope() > max_slp) {
                            max_slp = line.get_slope();
                            index = i;

                            delete this->l_line;
                            this->l_line = new Line(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y + error));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y - error));
            }
        }
    }

    void LinearSegment::__fsw_approximate(Point2D& p, double error) {
        if (this->u_line == nullptr) {
            Line u_line = Line::line(*this->pivot, Point2D(p.x, p.y + error));
            Line l_line = Line::line(*this->pivot, Point2D(p.x, p.y - error));

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
        }
        else {
            if (this->l_line->subs(p.x) > p.y + error || p.y - error > this->u_line->subs(p.x)) {
                this->is_complete = true;
            }
            else {
                if (p.y + error < this->u_line->subs(p.x)) {
                    Line u_line = Line::line(*this->pivot, Point2D(p.x, p.y + error));
                    this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
                }
                
                if (p.y - error > this->l_line->subs(p.x)) {
                    Line l_line = Line::line(*this->pivot, Point2D(p.x, p.y - error));
                    this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
                }
            }
        }
    }

    LinearSegment::LinearSegment(Point2D p) {
        this->t_start = p.x;
        this->isOptimal = false;
        this->pivot = new Point2D(p.x, p.y);
    }

    LinearSegment::LinearSegment(Point2D u_p, Point2D l_p) {
        this->t_start = u_p.x;
        this->isOptimal = true;
        this->u_cvx.append(Point2D(u_p.x, u_p.y));
        this->l_cvx.append(Point2D(l_p.x, l_p.y));
    }

    LinearSegment::~LinearSegment() {
        this->u_cvx.clear();
        this->l_cvx.clear();
        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
    }

    Line LinearSegment::getLine() {
        double slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
        double intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;

        return Line(slope, intercept);
    }

    void LinearSegment::approximate(Point2D& p, double error) {
        if (this->isOptimal) this->__optimal_approximate(p, error);
        else this->__fsw_approximate(p, error);
    }
    // End: material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        if (!this->s_1->is_complete) {
            Line line = this->s_1->getLine();
            this->curr_end = new Point2D(this->index - 1, line.subs(this->index - 1));
            this->yield();

            delete this->s_1;
            delete this->curr_end;
            delete this->prev_end;
        }
        else if (!this->s_2->is_complete) {
            Line line = this->s_2->getLine();
            Point2D p2 = Point2D(this->index - 1, line.subs(this->index - 1));
            Point2D p1 = Point2D(this->s_2->t_start, line.subs(this->s_2->t_start));
            
            this->curr_end = new Point2D(p1.x, p1.y);
            this->yield();
            delete this->prev_end;
            this->prev_end = this->curr_end;
            this->curr_end = new Point2D(p2.x, p2.y);
            this->yield();

            delete this->s_1;
            delete this->s_2;
            delete this->curr_end;
            delete this->prev_end;
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        long delta = this->curr_end->x - this->prev_end->x;
        
        obj->put(VariableByteEncoding::encode(delta));
        obj->put((float) this->curr_end->y);

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->index++, data->get_value());

        if (this->phase == 0) {
            this->s_1 = new LinearSegment(Point2D(p.x, p.y-this->error), Point2D(p.x, p.y+this->error));
            this->phase = 1;
        }
        else if (this->phase == 1) {
            if (!this->s_1->is_complete) {
                this->s_1->approximate(p, this->error);
                if (this->s_1->is_complete) {
                    this->s_2 = new LinearSegment(
                        Point2D(p.x - 1, this->s_1->l_line->subs(p.x - 1)),
                        Point2D(p.x - 1, this->s_1->u_line->subs(p.x - 1))
                    );
                    this->s_2->approximate(p, this->error);
                }
            }
            else if (!this->s_2->is_complete) {
                
                this->s_2->approximate(p, this->error);
                if (this->s_2->is_complete) {
                    Line line_2 = this->s_2->getLine();
                    Point2D p3 = Point2D(p.x - 1, line_2.subs(p.x - 1));
                    Point2D p2 = Point2D(this->s_2->t_start, line_2.subs(this->s_2->t_start));
                    Line line_1 = Line::line(p2, Line::intersection(*this->s_1->u_line, *this->s_1->l_line));
                    Point2D p1 = Point2D(this->s_1->t_start, line_1.subs(this->s_1->t_start));

                    this->prev_end = new Point2D(p1.x, p1.y);
                    this->curr_end = prev_end;
                    this->yield();
                    this->curr_end = new Point2D(p2.x, p2.y);
                    this->yield();
                    delete this->prev_end;
                    this->prev_end = this->curr_end;
                    this->curr_end = new Point2D(p3.x, p3.y);
                    this->yield();
                    delete this->prev_end;
                    this->prev_end = this->curr_end;

                    delete this->s_1;
                    delete this->s_2;
                    this->s_1 = new LinearSegment(p3);
                    this->s_1->approximate(p, this->error);
                    this->phase = 2;
                }
            }
        }
        else if (this->phase == 2) {
            if (!this->s_1->is_complete) {
                this->s_1->approximate(p, this->error);
                if (this->s_1->is_complete) {
                    this->s_2 = new LinearSegment(
                        Point2D(p.x - 1, this->s_1->l_line->subs(p.x - 1)),
                        Point2D(p.x - 1, this->s_1->u_line->subs(p.x - 1))
                    );
                    this->s_2->approximate(p, this->error);
                }
            }
            else if (!this->s_2->is_complete) {
                this->s_2->approximate(p, this->error);
                if (this->s_2->is_complete) {
                    Line line = this->s_2->getLine();
                    Point2D p2 = Point2D(p.x - 1, line.subs(p.x - 1));
                    Point2D p1 = Point2D(this->s_2->t_start, line.subs(this->s_2->t_start));
                    
                    this->curr_end = new Point2D(p1.x, p1.y);
                    this->yield();
                    delete this->prev_end;
                    this->prev_end = this->curr_end;
                    this->curr_end = new Point2D(p2.x, p2.y);
                    this->yield();
                    delete this->prev_end;
                    this->prev_end = this->curr_end;

                    delete this->s_1;
                    delete this->s_2;
                    this->s_1 = new LinearSegment(p2);
                    this->s_1->approximate(p, this->error);
                }
            }
        }
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
            long start = VariableByteEncoding::decode(compress_data);;
            float value = compress_data->getFloat();
            this->prev_end = new Point2D(start, value);

            return nullptr;
        }

        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        long length = VariableByteEncoding::decode(compress_data);
        float value = compress_data->getFloat();
        Point2D* curr_end = new Point2D(this->prev_end->x + length, value);
        Line line = Line::line(*curr_end, *this->prev_end);

        delete this->prev_end;
        this->prev_end = curr_end;

        for (long i=this->index; i<=curr_end->x; i++) {
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

        this->index = curr_end->x + 1;
        return base_obj;
    }
    // End: decompression
};