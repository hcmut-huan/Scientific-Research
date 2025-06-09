#include "piecewise-approximation/linear.hpp"

namespace OptimalPLA {

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        if (this->length >= 2) this->yield();
        
        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear();
        this->l_cvx.clear();
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Point2D p = Line::intersection(*this->u_line, *this->l_line);

        float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
        float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;
        Line line = Line::line(slope, p);
        // short embedded = this->length | (0 << 15);

        obj->put((short) this->length);
        obj->put((float) slope);
        obj->put((float) intercept);
        
        // float x1 = - this->u_line->get_intercept() / this->u_line->get_slope();
        // float x2 = - this->l_line->get_intercept() / this->l_line->get_slope();

        // int index = 0;

        // if (this->u_line->get_slope() * this->l_line->get_slope() < 0) {
        //     if (p.y >= 0) index = std::floor(x1);
        //     else index = std::ceil(x1);
        // }
        // else {
        //     index = std::floor(x1 < x2 ? x2 : x1);
        // }

        // if (std::abs(x1) > 8000 || std::abs(x2) > 8000 || index == (this->length-1) || std::abs(x2-x1) < 1) {
        //     float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
        //     Line line = Line::line(slope, p);
        //     short embedded = this->length | (0 << 15);

        //     obj->put((short) embedded);
        //     obj->put((float) line.get_slope());
        //     obj->put((float) line.get_intercept());
        // }
        // else {
        //     Point2D p2(index, 0);
        //     Line line = Line::line(p2, p);
        //     short embedded = (this->length - 1) | (1 << 15);

        //     obj->put((short) embedded);
        //     obj->put((short) index);
        //     obj->put((float) line.subs(this->length - 1));
        // }

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->length++, data->get_value());

        if (this->pivot == nullptr) {
            this->pivot = new Point2D(p.x, p.y);
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
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
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
        }
        else {
            if (this->l_line->subs(p.x) > p.y + this->error || p.y - this->error > this->u_line->subs(p.x)) {
                this->length--;
                this->yield();
                
                delete this->pivot;
                delete this->u_line; this->u_line = nullptr;
                delete this->l_line; this->l_line = nullptr;
                this->u_cvx.clear(); this->l_cvx.clear();

                this->pivot = new Point2D(0, p.y);
                u_cvx.append(Point2D(0, p.y - this->error));
                l_cvx.append(Point2D(0, p.y + this->error));
                this->length = 1;
            }
            else {
                bool update_u = p.y + this->error < this->u_line->subs(p.x);
                bool update_l = p.y - this->error > this->l_line->subs(p.x);

                if (update_u) {
                    int index = 0;
                    double min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + this->error));
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
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - this->error));
                        if (line.get_slope() > max_slp) {
                            max_slp = line.get_slope();
                            index = i;

                            delete this->l_line;
                            this->l_line = new Line(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y + this->error));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y - this->error));
            }
        }
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize() {
        // Do nothing
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress() {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        short embedded = this->compress_data->getShort();
        bool flag = embedded >> 15;

        int length = 0;
        float slp = 0;
        float intercept = 0;

        if (!flag) {
            length = embedded & (0xffff >> 1);
            slp = this->compress_data->getFloat();
            intercept = this->compress_data->getFloat();
        }
        else {
            length = (embedded & (0xffff >> 1)) + 1;
            short index = this->compress_data->getShort();
            float value = this->compress_data->getFloat();

            Point2D p1(index, 0);
            Point2D p2(length - 1, value);
            Line line = Line::line(p1, p2);

            slp = line.get_slope();
            intercept = line.get_intercept();
        }

        for (int i=0; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(slp * i + intercept));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(slp * i + intercept));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};