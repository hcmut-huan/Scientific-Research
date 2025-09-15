#include "piecewise-approximation/linear.hpp"

namespace ConnOptimalPLA {
    // Begin: material
    OptimalPLA::OptimalPLA(Point2D u_p, Point2D l_p) {
        this->t_start = u_p.x;
        this->u_cvx.append(l_p);
        this->l_cvx.append(u_p);
    }

    OptimalPLA::~OptimalPLA() {
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear(); this->l_cvx.clear();
    }

    Point2D* OptimalPLA::backCheck(OptimalPLA* prev) {
        long pivot = prev->t_end;

        if (this->l_line->subs(pivot) <= prev->u_line->subs(pivot) && 
            this->u_line->subs(pivot) >= prev->l_line->subs(pivot)) {

            return new Point2D(
                pivot, (this->l_line->subs(pivot) + this->u_line->subs(pivot)) / 2
            );
        }

        if (this->l_line->subs(pivot) <= prev->u_line->subs(pivot) && 
            this->l_line->subs(pivot) >= prev->l_line->subs(pivot)) {

            return new Point2D(
                pivot, (this->l_line->subs(pivot) + prev->l_line->subs(pivot)) / 2
            );
        }

        if (this->u_line->subs(pivot) <= prev->u_line->subs(pivot) && 
            this->u_line->subs(pivot) >= prev->l_line->subs(pivot)) {

            return new Point2D(
                pivot, (this->u_line->subs(pivot) + prev->u_line->subs(pivot)) / 2
            );
        }

        return nullptr;
    }

    void OptimalPLA::approximate(Point2D& p, double error) { 
        if (this->u_line == nullptr) {
            Line u_line = Line::line(this->u_cvx.at(0), Point2D(p.x, p.y + error));
            Line l_line = Line::line(this->l_cvx.at(0), Point2D(p.x, p.y - error));
            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
        }
        else {
            if ((this->l_line->subs(p.x) - p.y - error) > 0.0000001) {
                this->t_end = p.x - 1;
                this->is_complete = true;
                return;
            }
            else if (p.y - error - this->u_line->subs(p.x) > 0.0000001) {
                this->t_end = p.x - 1;
                this->is_complete = true;
                return;
            }
            else {
                bool update_u = p.y + error < this->u_line->subs(p.x);
                bool update_l = p.y - error > this->l_line->subs(p.x);

                if (update_u) {
                    int index = -1;
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

                    if (index >= 0) this->u_cvx.erase_from_begin(index);
                }
                if (update_l) {
                    int index = -1;
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

                    if (index >= 0) this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y + error));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y - error));
            }
        }
    }
    // End: material 

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->prev_end = new Point2D(0, -1);
    }

    void Compression::finalize() {
        if (this->seg_1 != nullptr) delete this->seg_1;
        if (this->seg_2 != nullptr) delete this->seg_2;
        if (this->seg_3 != nullptr) delete this->seg_3;
        if (this->seg_4 != nullptr) delete this->seg_4;
        if (this->prev_end != nullptr) delete this->prev_end;
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

        std::cout << p.x << " " << p.y << "\n";
        
        if (this->state == 0) {
            if (this->seg_1 == nullptr) {
                this->seg_1 = new OptimalPLA(
                    Point2D(p.x, p.y + this->error),
                    Point2D(p.x, p.y - this->error)
                );
            }
            else if (!this->seg_1->is_complete) {
                this->seg_1->approximate(p, this->error);
                if (this->seg_1->is_complete) {
                    this->state = 1;
                    long pivot = this->seg_1->t_start;
                    this->curr_end = new Point2D(
                        pivot, (this->seg_1->u_line->subs(pivot) + this->seg_1->l_line->subs(pivot)) / 2
                    );
                    this->yield();
                    
                    buffer.push_back(p);
                    this->seg_2 = new OptimalPLA(
                        Point2D(p.x, p.y + this->error),
                        Point2D(p.x, p.y - this->error)
                    );
                }
            }
        }
        else if (this->state == 1) {
            this->buffer.push_back(p);
            this->seg_2->approximate(p, this->error);
            if (this->seg_2->is_complete) {
                Point2D* back_point = this->seg_2->backCheck(this->seg_1);
                if (back_point != nullptr) {
                    // TODO: backward checking successfully
                    this->curr_end = back_point;
                    this->yield();

                    delete this->seg_1;
                    delete this->prev_end;
                    this->prev_end = this->curr_end;
                    this->seg_1 = this->seg_2;
                    this->seg_2 = new OptimalPLA(
                        Point2D(p.x, p.y + this->error),
                        Point2D(p.x, p.y - this->error)
                    );

                    this->buffer.clear();
                    this->buffer.push_back(p);
                }
                else {
                    long pivot = this->seg_1->t_end;
                    this->seg_3 = new OptimalPLA(
                        Point2D(pivot, this->seg_1->u_line->subs(pivot)),
                        Point2D(pivot, this->seg_1->l_line->subs(pivot))
                    );        

                    for (Point2D& p : this->buffer) {
                        if (!this->seg_3->is_complete) {
                            this->seg_3->approximate(p, this->error);
                            if (this->seg_3->is_complete) {
                                long pivot = this->seg_3->t_end;
                                this->seg_4 = new OptimalPLA(
                                    Point2D(pivot, this->seg_3->u_line->subs(pivot)),
                                    Point2D(pivot, this->seg_3->l_line->subs(pivot))
                                );
                                this->seg_4->approximate(p, this->error);
                            }
                        }
                        else {
                            this->seg_4->approximate(p, this->error);
                            if (this->seg_3->is_complete && this->seg_4->is_complete) {
                                this->state = 2;
                                break;
                            }
                        }
                    }

                    this->state = this->state == 1 ? 3 : this->state;
                    std::cout << "select: " << state << "\n";
                }
            }
        }
    
        if (this->state == 2) {
            // TODO: length check failed
            long pivot_1 = this->seg_1->t_end;
            this->curr_end = new Point2D(
                pivot_1, (this->seg_1->u_line->subs(pivot_1) + this->seg_1->l_line->subs(pivot_1)) / 2
            );
            this->yield();

            delete this->prev_end;
            this->prev_end = this->curr_end;

            long pivot_2 = this->seg_2->t_start;
            this->curr_end = new Point2D(
                pivot_2, (this->seg_2->u_line->subs(pivot_2) + this->seg_2->l_line->subs(pivot_2)) / 2
            );
            this->yield();

            delete this->prev_end;
            this->prev_end = this->curr_end;
            
            delete this->seg_1; this->seg_1 = nullptr;
            delete this->seg_3; this->seg_3 = nullptr;
            delete this->seg_4; this->seg_4 = nullptr;
            this->buffer.clear();
            this->buffer.push_back(p);
            
            this->state = 1;
            this->seg_1 = this->seg_2;
            this->seg_2 = new OptimalPLA(
                Point2D(p.x, p.y + this->error),
                Point2D(p.x, p.y - this->error)
            );
        }
        else if (this->state == 3) {
            // TODO: length check successfully
            long pivot = this->seg_1->t_end;
            this->curr_end = new Point2D(
                pivot, (this->seg_3->u_line->subs(pivot) + this->seg_3->l_line->subs(pivot)) / 2
            );
            this->yield();

            delete this->prev_end;
            this->prev_end = this->curr_end;

            delete this->seg_1; this->seg_1 = nullptr;
            delete this->seg_2; this->seg_2 = nullptr;
            delete this->seg_3; this->seg_3 = nullptr;

            this->state = 0;
            this->seg_1 = this->seg_4; this->seg_4 = nullptr;
        }

        std::cout << "end " << this->state << "\n";
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