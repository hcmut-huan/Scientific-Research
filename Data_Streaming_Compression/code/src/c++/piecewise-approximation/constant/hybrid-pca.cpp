#include "piecewise-approximation/constant.hpp"

namespace HybridPCA {

    Window::~Window() {
        for (Point2D* p : this->data) {
            delete p;
        }
    }

    int Window::size() {
        return this->data.size();
    }

    void Window::append(Point2D& p) {
        this->data.push_back(new Point2D(p.x, p.y));
        this->min = this->min < p.y ? this->min : p.y;
        this->max = this->max > p.y ? this->max : p.y;
    }


    int Buffer::size() {
        return this->windows.size();
    }

    void Buffer::pop() {
        this->windows.pop_back();
    }

    void Buffer::clear() {
        this->min = INFINITY;
        this->max = -INFINITY;
        for (Window* window : this->windows) 
            delete window;

        this->windows.clear();
    }

    double Buffer::value() {
        return (this->max + this->min) / 2;
    }

    void Buffer::append(Window* window) {
        this->windows.push_back(window);   
        this->min = this->min < window->min ? this->min : window->min;
        this->max = this->max > window->max ? this->max : window->max;
    }

    bool Buffer::is_appendable(Window* window, float bound) {
        float n_min = this->min < window->min ? this->min : window->min;
        float n_max = this->max > window->max ? this->max : window->max;

        return (n_max - n_min) <= 2*bound;
    }

    

    // Begin: compression
    void Compression::__pmc() {
        double min = INFINITY;
        double max = -INFINITY;
        double value = 0;
        
        int length = 0;
        for (Point2D* p : window->data) {
            min = min < p->y ? min : p->y;
            max = max > p->y ? max : p->y;
            
            if (max - min > 2 * this->error) {
                this->length = length;
                this->value = value;
                this->yield();
                min = p->y;
                max = p->y;
                length = 0;
            }
            value = (max + min) / 2;
            
            length++;
        }

        if (length > 0) {
            this->length = length;
            this->value = value;
            this->yield();
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->w_size = atoi(params[1]);
        this->n_window = atoi(params[2]);
        this->window = new Window;
    }

    void Compression::finalize() {
        if (this->buffer.is_appendable(this->window, this->error)) {
            this->buffer.append(this->window);
            this->length = (this->buffer.size() - 1) * this->w_size + this->window->size();
            this->value = this->buffer.value();
            this->yield();
            this->buffer.clear();
        }
        else {
            if (this->buffer.size() > 0) {
                this->length = this->buffer.size() * this->w_size;
                this->value = this->buffer.value();
                this->yield();
                this->buffer.clear();
            }
    
            if ((this->window->max - this->window->min) > 2 * this->error) {
                this->__pmc();
            }

            delete this->window;
            this->buffer.clear();
        }
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->index++, data->get_value());
        this->window->append(p);

        if (this->window->size() == this->w_size) {
            if (this->buffer.is_appendable(this->window, this->error)) {
                this->buffer.append(this->window);
                if (this->buffer.size() == this->n_window) {
                    this->length = this->buffer.size() * this->w_size;
                    this->value = this->buffer.value();
                    this->yield();
                    this->buffer.clear();
                }
            }
            else {
                if (this->buffer.size() > 0) {
                    this->length = this->buffer.size() * this->w_size;
                    this->value = this->buffer.value();
                    this->yield();
                    this->buffer.clear();
                }

                if ((this->window->max - this->window->min) > 2 * this->error) {
                    this->__pmc();
                    delete this->window;
                }
                else {
                    this->buffer.append(this->window);
                }
            }
            this->window = new Window;
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        obj->put((short) this->length);
        obj->put((float) this->value);

        return obj;
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

        unsigned short length = this->compress_data->getShort();
        float value = this->compress_data->getFloat();

        for (int i=0; i<length; i++) {
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

        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};