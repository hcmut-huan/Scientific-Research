#include "piecewise-approximation/constant.hpp"

namespace PMC {
    // Begin: compression
    void Compression::__mean(Univariate* data) {
        Point2D p(this->index++, data->get_value());   
        float n_value = (value * length + p.y) / (length + 1);
        this->min = this->min < p.y ? this->min : p.y;
        this->max = this->max > p.y ? this->max : p.y;
    
        if (this->max - n_value > this->error || n_value - this->min > this->error) {
            this->yield();
            this->prev_end = this->curr_end;

            this->min = p.y;
            this->max = p.y;
            this->value = p.y;
            length = 1;
        }
        else {
            value = n_value;
            length++;
        }

        this->curr_end = p.x;
    }

    void Compression::__midrange(Univariate* data) {
        Point2D p(this->index++, data->get_value());        
        this->min = this->min < p.y ? this->min : p.y;
        this->max = this->max > p.y ? this->max : p.y;
            
        if (this->max - this->min > 2 * this->error) {
            this->yield();
            this->prev_end = this->curr_end;

            this->min = p.y;
            this->max = p.y;
            this->value = p.y;
        }
        else {
            this->value = (this->max + this->min) / 2;
        }

        this->curr_end = p.x;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
    }

    void Compression::finalize() {
        this->yield();
    }

    void Compression::compress(Univariate* data) {
        if (this->mode == "midrange") this->__midrange(data);
        else if (this->mode == "mean") this->__mean(data);
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        obj->put((short)(this->curr_end - this->prev_end));
        obj->put((float) this->value);

        return obj;
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

        unsigned short length = compress_data->getShort();
        float value = compress_data->getFloat();

        for (int i=this->index + 1; i<=this->index + length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(value));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(value));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        this->index += length;
        return base_obj;
    }
    // End: decompression
};