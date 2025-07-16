#include "piecewise-approximation/linear.hpp"

namespace SwingFilter {

    // Begin: compression
    void Compression::__fit(Point2D& p) {
        this->A_num += (p.y-this->pivot->y)*(p.x-this->pivot->x);
        this->A_den += (p.x-this->pivot->x)*(p.x-this->pivot->x);
    }
    
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        double A_ig = this->A_num / this->A_den;
        double temp = A_ig > this->u_line->get_slope() ? this->u_line->get_slope() : A_ig;
        double a_ig = temp > this->l_line->get_slope() ? temp : this->l_line->get_slope();
        double b_ig = this->pivot->y - a_ig * this->pivot->x;

        Line line(a_ig, b_ig);
        this->curr_end = new Point2D(this->index-1, line.subs(this->index-1));
        this->yield();

        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        if (this->prev_end != nullptr) delete this->prev_end;
        if (this->curr_end != nullptr) delete this->curr_end;
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

        if (this->pivot == nullptr) {
            this->pivot = new Point2D(p.x, p.y);
            
            this->prev_end = new Point2D(0, p.y);
            this->curr_end = prev_end;
            this->yield();
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(*this->pivot, Point2D(p.x, p.y + this->error));
            Line l_line = Line::line(*this->pivot, Point2D(p.x, p.y - this->error));

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
        }
        else {
            if (this->l_line->subs(p.x) > p.y + this->error || p.y - this->error > this->u_line->subs(p.x)) {
                double A_ig = this->A_num / this->A_den;
                double temp = A_ig > this->u_line->get_slope() ? this->u_line->get_slope() : A_ig;
                double a_ig = temp > this->l_line->get_slope() ? temp : this->l_line->get_slope();
                double b_ig = this->pivot->y - a_ig * this->pivot->x;

                Line line(a_ig, b_ig);
                this->curr_end = new Point2D(p.x-1, line.subs(p.x-1));
                this->yield();
                delete this->prev_end;
                this->prev_end = this->curr_end; 
                this->curr_end = nullptr;

                delete this->pivot;
                delete this->u_line;
                delete this->l_line;

                this->pivot = new Point2D(p.x-1, line.subs(p.x-1));
                Line u_line = Line::line(*this->pivot, Point2D(p.x, p.y+this->error));
                Line l_line = Line::line(*this->pivot, Point2D(p.x, p.y-this->error));
                this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
                this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());

                this->A_den = 0;
                this->A_num = 0;
                this->__fit(*this->pivot);
            }
            else {
                if (p.y + this->error < this->u_line->subs(p.x)) {
                    Line u_line = Line::line(*this->pivot, Point2D(p.x, p.y+this->error));
                    this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
                }
                
                if (p.y - this->error > this->l_line->subs(p.x)) {
                    Line l_line = Line::line(*this->pivot, Point2D(p.x, p.y-this->error));
                    this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
                }
            }
        }

        this->__fit(p);
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
            long start = VariableByteEncoding::decode(this->compress_data);;
            float value = this->compress_data->getFloat();
            this->prev_end = new Point2D(start, value);

            return nullptr;
        }

        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        long length = VariableByteEncoding::decode(this->compress_data);
        float value = this->compress_data->getFloat();
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