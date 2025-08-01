#include "piecewise-approximation/linear.hpp"

namespace SemiMixedPLA {
    // Begin: material
    OptimalPLA::~OptimalPLA() {
        if (this->pivot != nullptr) delete this->pivot;         
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear(); this->l_cvx.clear();
        this->l_lines.clear(); this->u_lines.clear();
    }

    void OptimalPLA::updateExtrm() {
        if (this->status == -1) this->extrm = this->l_line;
        else if (this->status == 1) this->extrm = this->u_line;
    }

    bool OptimalPLA::isSemiConnected(Line* line, double bound) {
        Point2D inter = Line::intersection(*line, *this->extrm);
        
        if (inter.x > this->t_end + 1 || inter.x < this->t_start) return false;
        else if (inter.x > this->t_end) return true;

        int start = std::ceil(inter.x) - this->t_start;
        for (int i=start; i<this->window.size(); i++) {
            Point2D p = this->window.at(i);
            if (p.x > inter.x && std::abs(line->subs(p.x) - p.y) > bound) {
                return false;
            }
        }

        return true;
    }

    void OptimalPLA::reconstructCvx() {
        this->u_cvx.clear();
        for (Point2D& p : this->u_points) {
            if (p.x > this->t_end) break;
            else this->u_cvx.append(p);
        }

        this->l_cvx.clear();
        for (Point2D& p : this->l_points) {
            if (p.x > this->t_end) break;
            else this->l_cvx.append(p);
        }
    }

    Point2D OptimalPLA::shrink(OptimalPLA* prev_seg, double error) {
        Point2D last = this->popLast(true);
        this->t_end = last.x - 1;

        if (this->window.size() > 1) {
            if (last.x == this->l_lines.back().first) {
                delete this->l_line;
                this->l_lines.pop_back();
                this->l_line = new Line(
                    this->l_lines.back().second.get_slope(),
                    this->l_lines.back().second.get_intercept()
                );
            }
            if (last.x == this->u_lines.back().first) {
                delete this->u_line;
                this->u_lines.pop_back();
                this->u_line = new Line(
                    this->u_lines.back().second.get_slope(),
                    this->u_lines.back().second.get_intercept()
                );
            }
            return last;
        }
        else {
            delete this->u_line; this->u_line = new Line(INFINITY, 0);
            delete this->l_line; this->l_line = new Line(-INFINITY, 0);
            this->u_cvx.clear(); this->u_cvx.append(this->u_points.at(0));
            this->l_cvx.clear(); this->l_cvx.append(this->l_points.at(0));

            this->extendBackward(prev_seg, error);

            return last;
        }
    }

    void OptimalPLA::extendBackward(OptimalPLA* prev_seg, double error) {
        Point2D p = prev_seg->popLast(false);
        this->t_start = p.x;
        this->window.insert(this->window.begin(), p);

        bool update_l = this->window.size() == 2 ? true : p.y + error < this->l_line->subs(p.x);
        bool update_u = this->window.size() == 2 ? true : p.y - error > this->u_line->subs(p.x);

        if (update_u) {
            int index = 0;
            for (int i=0; i<this->l_cvx.size(); i++) {
                Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - error));
                if (line.get_slope() < this->u_line->get_slope()) {
                    index = i;
                    delete this->u_line;
                    this->u_line = new Line(line.get_slope(), line.get_intercept());
                }
            }
            if (index != this->l_cvx.size() - 1) this->l_cvx.erase_from_end(index + 1);
        }
        
        if (update_l) {
            int index = 0;
            for (int i=0; i<this->u_cvx.size(); i++) {
                Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + error));
                if (line.get_slope() > this->l_line->get_slope()) {
                    index = i;
                    delete this->l_line;
                    this->l_line = new Line(line.get_slope(), line.get_intercept());
                }
            }
            
            if (index != this->u_cvx.size() - 1) this->u_cvx.erase_from_end(index + 1);
            
        }

        if (update_u) {
            this->u_cvx.append_backward(Point2D(p.x, p.y - error));
        }
        if (update_l) {
            this->l_cvx.append_backward(Point2D(p.x, p.y + error));
        }
    }

    Point2D OptimalPLA::popLast(bool flag) {
        Point2D p = this->window.back();
        this->window.pop_back();

        if (flag) {
            if (p.x == this->u_points.back().x) this->u_points.pop_back();
            if (p.x == this->l_points.back().x) this->l_points.pop_back();
        }
        
        return p;
    }

    void OptimalPLA::approximate(Point2D& p, double error) {        
        if (this->pivot == nullptr) {
            this->t_start = p.x;
            this->pivot = new Point2D(p.x, p.y);
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
            this->u_points.push_back(Point2D(p.x, p.y - error));
            this->l_points.push_back(Point2D(p.x, p.y + error));
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y - error), 
                Point2D(p.x, p.y + error)
            );
            Line l_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y + error), 
                Point2D(p.x, p.y - error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
            this->u_points.push_back(Point2D(p.x, p.y - error));
            this->l_points.push_back(Point2D(p.x, p.y + error));
            this->u_lines.push_back(std::make_pair(p.x, Line(this->u_line->get_slope(), this->u_line->get_intercept())));
            this->l_lines.push_back(std::make_pair(p.x, Line(this->l_line->get_slope(), this->l_line->get_intercept())));
        }
        else {
            if (this->l_line->subs(p.x) > p.y + error) {
                this->t_end = p.x - 1;
                this->status = -1;
                this->is_complete = true;
                return;
            }
            else if (p.y - error > this->u_line->subs(p.x)) {
                this->t_end = p.x - 1;
                this->status = 1;
                this->is_complete = true;
                return;
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

                if (update_u) {
                    this->l_cvx.append(Point2D(p.x, p.y + error));
                    this->l_points.push_back(Point2D(p.x, p.y + error));
                    this->u_lines.push_back(std::make_pair(p.x, Line(this->u_line->get_slope(), this->u_line->get_intercept())));
                }
                if (update_l) {
                    this->u_cvx.append(Point2D(p.x, p.y - error));
                    this->u_points.push_back(Point2D(p.x, p.y - error));
                    this->l_lines.push_back(std::make_pair(p.x, Line(this->l_line->get_slope(), this->l_line->get_intercept())));
                }
            }
        }

        this->window.push_back(p);
    }

    int SemiOptimalPLA::__check(double error) {
        if (this->seg_2->status == -1) {
            if (this->seg_1->isSemiConnected(this->seg_2->l_line, error)) return 1;
            else if (this->seg_1->isSemiConnected(this->seg_2->u_line, error)) return 2;
            else return 3;
        }
        else {
            if (this->seg_1->isSemiConnected(this->seg_2->u_line, error)) return 1;
            else if (this->seg_1->isSemiConnected(this->seg_2->l_line, error)) return 2;
            else return 3;
        }
    }

    SemiOptimalPLA::SemiOptimalPLA() {
        this->seg_1 = new OptimalPLA;
        this->seg_2 = new OptimalPLA;
        this->inter = new BinObj;
    }

    SemiOptimalPLA::SemiOptimalPLA(OptimalPLA* p_seg, Point2D* p_point) {
        this->seg_1 = p_seg;
        this->seg_2 = new OptimalPLA;

        this->inter = new BinObj;
        this->inter->put((float) p_point->x);
        this->inter->put((float) p_point->y);

        double value = this->seg_1->extrm->subs(this->seg_1->t_start);
        this->inter->put((float) value);

        this->prev_end = new Point2D(this->seg_1->t_start, value);
        this->seg_count = 1;
    }

    SemiOptimalPLA::~SemiOptimalPLA() {
        if (this->seg_1 != nullptr) delete this->seg_1;
        if (this->seg_2 != nullptr) delete this->seg_2;
        if (this->prev_end != nullptr) delete this->prev_end;
    }

    void SemiOptimalPLA::finalize() {
        Point2D last_point = this->seg_1->popLast(false);
        
        this->inter->put((float) (last_point.x - this->prev_end->x));
        this->inter->put((float) (this->seg_1->extrm->subs(last_point.x)));
    }

    long SemiOptimalPLA::getLastEnd() {
        if (this->seg_2->is_complete) return this->seg_2->t_end;
        else return this->seg_1->t_end;
    }

    BinObj* SemiOptimalPLA::getInter() {
        BinObj* obj = this->inter;
        this->inter = new BinObj;
        this->seg_count = 1;

        return obj;
    }

    OptimalPLA* SemiOptimalPLA::getPendSeg() {
        OptimalPLA* p_seg = new OptimalPLA;
        p_seg->window = this->seg_2->window;
        p_seg->t_end = this->seg_2->t_end;
        p_seg->t_start = this->seg_2->t_start;
        p_seg->status = this->seg_2->status;
        p_seg->is_complete = true;
        p_seg->u_line = new Line(
            this->seg_2->u_line->get_slope(),
            this->seg_2->u_line->get_intercept()
        );
        p_seg->l_line = new Line(
            this->seg_2->l_line->get_slope(),
            this->seg_2->l_line->get_intercept()
        );
        p_seg->updateExtrm();
        
        return p_seg;
    }

    Point2D* SemiOptimalPLA::getPendPoint() {
        return new Point2D(
            this->prev_end->x - this->seg_2->t_start,
            this->seg_1->extrm->subs(this->seg_2->t_start)
        );
    }

    void SemiOptimalPLA::makeSemiConnect(int semi_case, double error) {
        switch (semi_case) {
            case 3:
                this->seg_2->status = -this->seg_1->status;
                do {
                    Point2D shrink_point = this->seg_2->shrink(this->seg_1, error);
                    this->buffer.push_front(shrink_point);
                } 
                while(this->__check(error) == 3);  
                this->seg_2->reconstructCvx();

            case 2:
                while (this->__check(error) == 2) {
                    this->seg_2->extendBackward(this->seg_1, error);
                }
            case 1:
                this->seg_2->updateExtrm();
                Point2D inter_p = Line::intersection(
                    *this->seg_1->extrm, *this->seg_2->extrm
                );

                this->inter->put((float) (inter_p.x - this->prev_end->x));
                this->inter->put((float) inter_p.y);

                delete this->prev_end;
                this->prev_end = new Point2D(inter_p.x, inter_p.y);
                break;
        }

        delete this->seg_1;
        this->seg_1 = this->seg_2;
        this->seg_2 = new OptimalPLA();

        while (!this->buffer.empty()) {
            Point2D p = this->buffer.front();
            this->buffer.pop_front();
            this->seg_2->approximate(p, error);
        }
    }

    int SemiOptimalPLA::approximate(Point2D& point, double error) {
        this->buffer.push_back(point);
        Point2D p = this->buffer.front();
        this->buffer.pop_front();

        if (!this->seg_1->is_complete) { 
            this->seg_1->approximate(p, error);
            if (this->seg_1->is_complete) {
                this->seg_1->updateExtrm();
                this->prev_end = new Point2D(
                    this->seg_1->t_start,
                    this->seg_1->extrm->subs(this->seg_1->t_start)
                );
                
                this->inter->put((float) prev_end->x);
                this->inter->put((float) prev_end->y);

                this->seg_count++;
            }
        }
        if (this->seg_1->is_complete) {
            this->seg_2->approximate(p, error);
            if (this->seg_2->is_complete) this->seg_2->updateExtrm();
        }

        if (this->seg_1->is_complete && this->seg_2->is_complete) {
            this->seg_count++;
            return this->__check(error);
        }
        
        return 0;
    }
    // End: material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->semi_segs_1 = new SemiOptimalPLA();
    }

    void Compression::finalize() {
        if (this->state == 0) {
            this->semi_segs_1->finalize();
        }
        if (this->state == 1) {
            if (this->semi_segs_1->seg_count > this->semi_segs_2->seg_count) {
                delete this->semi_segs_1;
                this->semi_segs_1 = this->semi_segs_2;
            }

            this->semi_segs_1->finalize();
        }
        
        this->state = 1;
        this->yield();

        delete this->semi_segs_1;
    }

    BinObj* Compression::serialize() {
        if (this->state == 1) return this->semi_segs_1->getInter();
        else if (this->state == 2) return this->semi_segs_1->getInter();
        else return this->semi_segs_2->getInter();

        return nullptr;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->index++, data->get_value());
        
        if (this->state == 0) {            
            int flag = this->semi_segs_1->approximate(p, this->error);
            if (flag == 1 || flag == 2) {
                this->semi_segs_1->makeSemiConnect(flag, this->error);
                this->semi_segs_1->approximate(p, this->error);
            }
            else if (flag == 3) {
                this->state = 1;
                this->yield();
                this->semi_segs_2 = new SemiOptimalPLA(
                    this->semi_segs_1->getPendSeg(),
                    this->semi_segs_1->getPendPoint()
                );
                this->semi_segs_1->makeSemiConnect(flag, this->error);
            }
        }
        if (this->state == 1) {
            int flag_1 = this->semi_segs_1->approximate(p, this->error);
            int flag_2 = this->semi_segs_2->approximate(p, this->error);
            
            if (flag_1 != 0) this->semi_segs_1->makeSemiConnect(flag_1, this->error);

            if (flag_1 == 3) {
                if (this->semi_segs_1->seg_count > this->semi_segs_2->seg_count + 1) {
                    this->state = 3;
                }
                else if (this->semi_segs_1->seg_count > this->semi_segs_2->seg_count) {
                    this->state = this->semi_segs_2->getLastEnd() >= this->semi_segs_1->getLastEnd() ? 3 : 2;
                }
                else if (this->semi_segs_1->seg_count == this->semi_segs_2->seg_count) {
                    this->state = this->semi_segs_2->seg_1->t_end >= this->semi_segs_1->getLastEnd() ? 3 : 2;
                }
            }

            if (flag_2 != 0) this->semi_segs_2->makeSemiConnect(flag_2, this->error);

            if (this->state == 2) {
                this->yield();
                this->state = 0;
                delete this->semi_segs_2;

                int flag = this->semi_segs_1->approximate(p, this->error);
                if (flag == 1 || flag == 2) {
                    this->semi_segs_1->makeSemiConnect(flag, this->error);
                    this->semi_segs_1->approximate(p, this->error);
                }
                else if (flag == 3) {
                    this->state = 1;
                    this->semi_segs_2 = new SemiOptimalPLA(
                        this->semi_segs_1->getPendSeg(),
                        this->semi_segs_1->getPendPoint()
                    );
                    this->semi_segs_1->makeSemiConnect(flag, this->error);
                    this->semi_segs_1->approximate(p, this->error);
                    this->semi_segs_2->approximate(p, this->error);
                }
            }
            else if (this->state == 3) {
                this->yield();
                this->state = 0;
                delete this->semi_segs_1;
                this->semi_segs_1 = this->semi_segs_2;

                if (flag_2 != 0) {
                    int flag = this->semi_segs_1->approximate(p, this->error);
                    if (flag == 1 || flag == 2) {
                        this->semi_segs_1->makeSemiConnect(flag, this->error);
                        this->semi_segs_1->approximate(p, this->error);
                    }
                    else if (flag == 3) {
                        this->state = 1;
                        this->semi_segs_2 = new SemiOptimalPLA(
                            this->semi_segs_1->getPendSeg(),
                            this->semi_segs_1->getPendPoint()
                        );
                        this->semi_segs_1->makeSemiConnect(flag, this->error);
                        this->semi_segs_1->approximate(p, this->error);
                        this->semi_segs_2->approximate(p, this->error);
                    }
                }
            }
            else {
                if (flag_1 != 0) {
                    this->semi_segs_1->approximate(p, this->error);
                }
                if (flag_2 != 0) {
                    int flag = this->semi_segs_2->approximate(p, this->error);
                    if (flag != 0) {
                        this->semi_segs_2->makeSemiConnect(flag, this->error);
                        this->semi_segs_2->approximate(p, this->error);
                    }
                }
            }
        }
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize() {
        // Do nothing
    }

    void Decompression::finalize() {
        if (this->prev_end != nullptr) delete this->prev_end;
    }

    CSVObj* Decompression::decompress() {
        if (this->prev_end == nullptr) {
            float start = this->compress_data->getFloat();
            float value = compress_data->getFloat();
            this->prev_end = new Point2D(start, value);

            return nullptr;
        }

        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        float delta = this->compress_data->getFloat();
        if (delta > 0) {
            float value = compress_data->getFloat();

            Point2D* curr_end = new Point2D(this->prev_end->x + delta, value);
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