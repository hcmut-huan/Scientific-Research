#include "piecewise-approximation/linear.hpp"

namespace CovariancePLA {
    // Begin: compression
    bool Compression::__error_bound_verify(Line& line) {
        for (int i=0; i<this->u_cvx.size(); i++) {
            Point2D p = this->u_cvx.at(i);
            if (line.subs(p.x) < p.y) 
                return false;
        }

        for (int i=0; i<this->l_cvx.size(); i++) {
            Point2D p = this->l_cvx.at(i);
            if (line.subs(p.x) > p.y) 
                return false;
        }

        return true;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        if (this->length >= 2) this->yield();
        if (this->line != nullptr) delete this->line;

        this->u_cvx.clear();
        this->l_cvx.clear();
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;

        obj->put(VariableByteEncoding::encode(this->length - 1));
        obj->put((float) this->line->get_slope());
        obj->put((float) this->line->get_intercept());

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->length++, data->get_value());
        this->u_cvx.append(Point2D(p.x, p.y - this->error));
        this->l_cvx.append(Point2D(p.x, p.y + this->error));

        this->accumulate += p.x*p.y;
        this->accumulate_square += p.x*p.x;
        this->average_x = (this->average_x * (this->length - 1) + p.x) / this->length;
        this->average_y = (this->average_y * (this->length - 1) + p.y) / this->length;

        if (this->length > 1) {
            double variance = (double) this->accumulate_square / this->length - this->average_x * this->average_x;
            double covariance = (double) this->accumulate / this->length - this->average_x * this->average_y;
            double slope = covariance / variance;
            double intercept = this->average_y - slope * this->average_x;

            Line n_line(slope, intercept);
            if (this->__error_bound_verify(n_line)) {
                if (this->line != nullptr) delete this->line;
                this->line = new Line(n_line.get_slope(), n_line.get_intercept());
            }
            else {
                this->yield();
                this->length = 1;
                this->average_x = 0;
                this->average_y = p.y;
                this->accumulate = 0;
                this->accumulate_square = 0;

                this->u_cvx.clear(); this->u_cvx.append(Point2D(0, p.y - this->error));
                this->l_cvx.clear(); this->l_cvx.append(Point2D(0, p.y + this->error));
            }
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

        int length = VariableByteEncoding::decode(compress_data);
        float slp = compress_data->getFloat();
        float intercept = compress_data->getFloat();
        
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