#include "piecewise-approximation/linear.hpp"

namespace OptimalPLA {

    bool Compression::__yield_10_bytes(BinObj* obj) {
        Point2D p = Line::intersection(*this->u_line, *this->l_line);
        float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
        float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;

        short embedded = this->length | (0 << 13);
        obj->put((short) embedded);
        obj->put((float) slope);
        obj->put((float) intercept);

        return true;
    }

    bool Compression::__yield_8_bytes(BinObj* obj) {
        Point2D p = Line::intersection(*this->u_line, *this->l_line);
        double u_left = this->u_line->get_intercept();
        double l_left = this->l_line->get_intercept();
        double u_right = this->u_line->subs(this->length);
        double l_right = this->l_line->subs(this->length);
        double u_root = - this->u_line->get_intercept() / this->u_line->get_slope();
        double l_root = - this->l_line->get_intercept() / this->l_line->get_slope();

        if (std::abs((int) u_right - (int) l_right) >= 1) {
            double value = (u_right > l_right) ? std::floor(u_right) : std::floor(l_right);
            Line line = Line::line(p, Point2D(this->length, value));

            if (std::abs(value) < 32000) {
                short embedded = this->length | (1 << 13);
                obj->put((short) embedded);
                obj->put((float) line.get_slope());
                obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                return true;
            }
        }
        
        else if (std::abs((int) u_left - (int) l_left) >= 1) {
            double value = (u_left > l_left) ? std::floor(u_left) : std::floor(l_left);
            Line line = Line::line(p, Point2D(0, value));

            if (std::abs(value) < 32000) {
                short embedded = this->length | (2 << 13);
                obj->put((short) embedded);
                obj->put((float) line.get_slope());
                obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                return true;
            }
        }

        else if (std::abs((int) u_root - (int) l_root) >= 1) {
            double root = 0;

            if (this->u_line->get_slope() * this->l_line->get_slope() < 0) {
                if (p.y >= 0) root = std::floor(u_root);
                else root = std::ceil(u_root);
            }
            else {
                root = u_root > l_root ? std::floor(u_root) : std::floor(l_root);
            }

            if (std::abs(root) < 32000 && root != this->length) {
                Point2D p2(root, 0);
                Line line = Line::line(p2, p);
                short embedded = this->length | (3 << 13);

                obj->put((short) embedded);
                obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) root)));
                obj->put((float) line.subs(this->length));

                return true;
            }
        }

        return false;
    }

    bool Compression::__yield_6_bytes(BinObj* obj) {
        if (this->u_cvx.size() < 2 || this->l_cvx.size() < 2) return false;

        Line u_line = Line::line(this->l_cvx.at(0), this->l_cvx.at(this->l_cvx.size()-1)); 
        for (int i=0; i<this->l_cvx.size(); i++) {
            Point2D p = this->l_cvx.at(i);
            if (u_line.subs(p.x) > p.y) u_line = Line::line(u_line.get_slope(), p);
        }
        
        Line l_line = Line::line(this->u_cvx.at(0), this->u_cvx.at(this->u_cvx.size()-1)); 
        for (int i=0; i<this->u_cvx.size(); i++) {
            Point2D p = this->u_cvx.at(i);
            if (l_line.subs(p.x) < p.y) l_line = Line::line(l_line.get_slope(), p);
        }

        Point2D p = Line::intersection(u_line, l_line);
        if (p.x >= 0 && p.x <= this->length) {
            if (p.x >= Line::intersection(*this->u_line, *this->l_line).x) {
                Point2D p1 = Line::intersection(u_line, *this->u_line);
                Point2D p2 = Line::intersection(l_line, *this->l_line);
                u_line = Line::line(p1, Point2D(this->length, l_line.subs(this->length)));
                l_line = Line::line(p2, Point2D(this->length, u_line.subs(this->length)));  
                p = Line::intersection(u_line, l_line);   
            }
            else {
                Point2D p1 = Line::intersection(u_line, *this->l_line);
                Point2D p2 = Line::intersection(l_line, *this->u_line);
                u_line = Line::line(p1, Point2D(0, l_line.get_intercept()));
                l_line = Line::line(p2, Point2D(0, u_line.get_intercept()));    
                p = Line::intersection(u_line, l_line);    
            }
        }

        if (p.x < 0 || p.x > this->length) {
            double u_right = u_line.subs(this->length);
            double l_right = l_line.subs(this->length);
            double u_left = u_line.get_intercept();
            double l_left = l_line.get_intercept();
            double u_root = - u_line.get_intercept() / u_line.get_slope();
            double l_root = - l_line.get_intercept() / l_line.get_slope();

            if (u_line.get_slope() * l_line.get_slope() > 0) {
                if (std::abs((int) u_left - (int) l_left) >= 1 
                    && std::abs((int) u_right - (int) l_right) >= 1) {
                        double value_1 = u_left > l_left ? std::floor(u_left) : std::floor(l_left);
                        double value_2 = u_right > l_right ? std::floor(u_right) : std::floor(l_right);;
                        short embedded = this->length | (4 << 13);

                        obj->put((short) embedded);
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value_1)));
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value_2)));

                        return true;
                    }

                else if (std::abs((int) u_right - (int) l_right) >= 1
                    && std::abs((int) u_root - (int) l_root) >= 1) {
                        double root = 70000;
                        double f_root = u_root > l_root ? std::floor(u_root) : std::floor(l_root);
                        double c_root = u_root < l_root ? std::ceil(u_root) : std::ceil(l_root);

                        if (p.x > this->length) root = c_root;
                        else if (f_root >= p.x) root = f_root;
                        else if (c_root >= p.x) root = c_root;

                        if (std::abs(root) < 32000 && root <= 0) {
                            double value = u_right > l_right ? std::floor(u_right) : std::floor(l_right);;
                            short embedded = this->length | (5 << 13);

                            obj->put((short) embedded);
                            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) root)));
                            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                            return true;
                        }
                    }

                else if (std::abs((int) u_left - (int) l_left) >= 1
                    && std::abs((int) u_root - (int) l_root) >= 1) {
                        double root = 70000;
                        double f_root = u_root > l_root ? std::floor(u_root) : std::floor(l_root);
                        double c_root = u_root < l_root ? std::ceil(u_root) : std::ceil(l_root);

                        if (p.x < 0) root = f_root;
                        else if (f_root <= p.x) root = f_root;
                        else if (c_root <= p.x) root = c_root;

                        if (std::abs(root) < 32000 && root >= this->length) {
                            double value = u_left > l_left ? std::floor(u_left) : std::floor(l_left);;
                            short embedded = this->length | (6 << 13);

                            obj->put((short) embedded);
                            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) root)));
                            obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                            return true;
                        }
                    }
            }
            else {             
                if (std::abs((int) u_left - (int) l_left) >= 1 
                    && std::abs((int) u_right - (int) l_right) >= 1) {
                        double value_1 = u_left > l_left ? std::floor(u_left) : std::floor(l_left);
                        double value_2 = u_right > l_right ? std::floor(u_right) : std::floor(l_right);;
                        short embedded = this->length | (4 << 13);

                        obj->put((short) embedded);
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value_1)));
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value_2)));

                        return true;
                    }

                else if (std::abs((int) u_right - (int) l_right) >= 1) {
                    double root = 70000;

                    if (p.x > this->length) {
                        if (l_root > p.x) root = std::floor(u_root);
                        else root = std::floor(l_root);
                    }
                    else {
                        if (l_root < p.x) root = std::ceil(u_root);
                        else root = std::ceil(l_root);

                        if (root < p.x) root = 70000;
                    }

                    if (std::abs(root) < 32000 && root <= 0) {
                        double value = u_right > l_right ? std::floor(u_right) : std::floor(l_right);;
                        short embedded = this->length | (5 << 13);

                        obj->put((short) embedded);
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) root)));
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                        return true;
                    }
                }
                    
                else if (std::abs((int) u_left - (int) l_left) >= 1) {
                    double root = 70000;

                    if (p.x > this->length) {
                        if (l_root > p.x) root = std::floor(u_root);
                        else root = std::floor(l_root);

                        if (root > p.x) root = 70000;
                    }
                    else {
                        if (l_root < p.x) root = std::ceil(u_root);
                        else root = std::ceil(l_root);
                    }
                    
                    if (std::abs(root) < 32000 && root >= this->length) {
                        double value = u_left > l_left ? std::floor(u_left) : std::floor(l_left);;
                        short embedded = this->length | (6 << 13);

                        obj->put((short) embedded);
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) root)));
                        obj->put(VariableByteEncoding::encode(ZigZagEncoding::encode((short) value)));

                        return true;
                    }   
                }
            }
        }

        return false;
    }

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

        // if (this->__yield_6_bytes(obj)) {return obj;}
        // else if (this->__yield_8_bytes(obj)) {return obj;}
        // else {
        //     this->__yield_10_bytes(obj); 
        //     return obj; 
        // }

        Point2D p = Line::intersection(*this->u_line, *this->l_line);
        float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
        float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;

        obj->put(VariableByteEncoding::encode(this->length));
        obj->put((float) slope);
        obj->put((float) intercept);

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
            if (this->l_line->subs(p.x) > p.y + this->error || p.y - this->error > this->u_line->subs(p.x) || this->length > 8000) {
                this->length--;
                this->yield();
                
                delete this->u_line; this->u_line = nullptr;
                delete this->l_line; this->l_line = nullptr;
                this->u_cvx.clear(); this->l_cvx.clear();

                this->pivot->x = 0; this->pivot->y = p.y;
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

        // unsigned short embedded = compress_data->getShort();
        // int flag = embedded >> 13;
        // int length = embedded & (0xffff >> 3);
        // float slp = 0;
        // float intercept = 0;

        // if (flag == 0) {
        //     slp = compress_data->getFloat();
        //     intercept = compress_data->getFloat();
        // }
        // else if (flag == 1) {
        //     slp = compress_data->getFloat();
        //     short value = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     Line line = Line::line(slp, Point2D(length, value));
        //     intercept = line.get_intercept();
        // }
        // else if (flag == 2) {
        //     slp = compress_data->getFloat();
        //     short value = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     Line line = Line::line(slp, Point2D(0, value));
        //     intercept = line.get_intercept();
        // }
        // else if (flag == 3) {
        //     short root = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     double value = compress_data->getFloat();
        //     Line line = Line::line(Point2D(root, 0), Point2D(length, value));
            
        //     slp = line.get_slope();
        //     intercept = line.get_intercept();
        // }
        // else if (flag == 4) {
        //     short value_1 = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     short value_2 = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     Line line = Line::line(Point2D(0, value_1), Point2D(length, value_2));

        //     slp = line.get_slope();
        //     intercept = line.get_intercept();
        // }
        // else if (flag == 5) {
        //     short root = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     short value = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     Line line = Line::line(Point2D(root, 0), Point2D(length, value));

        //     slp = line.get_slope();
        //     intercept = line.get_intercept();
        // }
        // else if (flag == 6) {
        //     short root = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     short value = ZigZagEncoding::decode((unsigned short) VariableByteEncoding::decode(compress_data));
        //     Line line = Line::line(Point2D(root, 0), Point2D(0, value));

        //     slp = line.get_slope();
        //     intercept = line.get_intercept();
        // }

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