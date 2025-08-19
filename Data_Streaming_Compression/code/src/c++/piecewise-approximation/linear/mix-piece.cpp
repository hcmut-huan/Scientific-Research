#include "piecewise-approximation/linear.hpp"

namespace MixPiece {
    // Begin: Material 
    Segment::Segment(long t, double a, double b) {
        this->t = t;
        this->a = a;
        this->b = b;
    }

    Segment::Segment(long t, double aMin, double aMax, double b) {
        this->t = t;
        this->aMin = aMin;
        this->aMax = aMax;
        this->b = b;
        this->a = (aMin + aMax) / 2;
    }

    BBlock::BBlock(double b) {
        this->b = b;
        this->blocks.push_back(Block());
    }

    BBlock::BBlock(double b, double a_u, double a_l, long t) {
        this->b = b;
        this->blocks.push_back(Block(a_u, a_l, t));
    }

    void BBlock::push_back(double a_u, double a_l, long t) {
        this->blocks.push_back(Block(a_u, a_l, t));
    }

    Segment BBlock::pop_back() {
        Segment result(this->blocks.back().t[0], this->blocks.back().a_l, 
            this->blocks.back().a_u, b);

        this->blocks.pop_back();
        return result;
    }

    void BBlock::sort() {
        for (Block& block : this->blocks) {
            std::sort(block.t.begin(), block.t.end());
        }
    }

    bool BBlock::is_intersect(double a_u, double a_l) {
        return a_l <= this->blocks.back().a_u && a_u >= this->blocks.back().a_l;
    }

    void BBlock::intersect(double a_u, double a_l, long t) {
        if (this->blocks.back().a_u > a_u) this->blocks.back().a_u = a_u;
        if (this->blocks.back().a_l < a_l) this->blocks.back().a_l = a_l;

        this->blocks.back().t.push_back(t);
    }

    BBlock::Block::Block() {
        this->a_u = INFINITY;
        this->a_l = -INFINITY;
    }

    BBlock::Block::Block(double a_u, double a_l, long t) {
        this->a_u = a_u;
        this->a_l = a_l;
        this->t.push_back(t);
    }

    ABlock::ABlock() {
        this->a_u = INFINITY;
        this->a_l = -INFINITY;
    }

    ABlock::ABlock(double b, double a_u, double a_l, long t) {
        this->a_u = a_u;
        this->a_l = a_l;

        this->blocks.push_back(Block(b, t));
    }

    void ABlock::sort() {
        std::sort(this->blocks.begin(), this->blocks.end(), 
            [](const Block& a, const Block& b) { return a.t < b.t; });
    }

    bool ABlock::is_intersect(double a_u, double a_l) {
        return a_l <= this->a_u && a_u >= this->a_l;
    } 

    void ABlock::intersect(double b, double a_u, double a_l, long t) {
        if (this->a_u > a_u) this->a_u = a_u;
        if (this->a_l < a_l) this->a_l = a_l;

        this->blocks.push_back(Block(b, t));
    }

    ABlock::Block::Block(double b, long t) {
        this->b = b;
        this->t = t;
    }
    // End: Material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->n_segment = atoi(params[1]);
    }

    void Compression::finalize() {
        if (this->flag > 0) {
            this->intervals[this->b_1].push_back(Segment(this->index - 1, this->slp_l_1, this->slp_u_1, this->b_1));
        }
        else {
            this->intervals[this->b_2].push_back(Segment(this->index - 1, this->slp_l_2, this->slp_u_2, this->b_2));
        }

        this->__group();
        this->intervals.clear();
    }

    void Compression::__serializeBblock(BinObj* obj) {
        obj->put((int) this->BBlocks.size());
        for (BBlock& BBlock : this->BBlocks) {
            BBlock.sort();
            float b = BBlock.b;
            obj->put((float) b);
            obj->put((int) BBlock.blocks.size());
            for (BBlock::Block& block : BBlock.blocks) {
                float a = (block.a_l + block.a_u) / 2;
                obj->put((float) a);
                obj->put((int) block.t.size());
                long pivot = 0;
                for (int i=0; i<block.t.size(); i++) {
                    long time = block.t[i];
                    obj->put(VariableByteEncoding::encode(time - pivot));
                    pivot = time;
                }
            }
        }
    }

    void Compression::__serializeAblock(BinObj* obj) {
        obj->put((int) this->ABlocks.size());
        for (ABlock& ABlock : ABlocks) {
            ABlock.sort();
            float a = (ABlock.a_u + ABlock.a_l) / 2;
            obj->put((float) a);
            obj->put((int) ABlock.blocks.size());
            long pivot = 0;
            for (ABlock::Block& block : ABlock.blocks) {
                obj->put((float) block.b);
                obj->put(VariableByteEncoding::encode(block.t - pivot));
                pivot = block.t;
            }
        }
    }

    void Compression::__serializeRblock(BinObj* obj) {
        std::sort(this->r_blocks.begin(), this->r_blocks.end(), 
            [](const Segment& a, const Segment& b)
            { return a.t < b.t; });
        
        obj->put((int) this->r_blocks.size());
        long pivot = 0;
        for (Segment& block : this->r_blocks) {
            obj->put((float) ((block.aMax + block.aMin) / 2));
            obj->put((float) block.b);
            obj->put(VariableByteEncoding::encode(block.t - pivot));
            pivot = block.t;
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;        
        this->__serializeBblock(obj);
        this->__serializeAblock(obj);
        this->__serializeRblock(obj);

        return obj;
    }

    void Compression::__group() {
        std::vector<Segment> ungrouped;
        
        for (auto it = this->intervals.begin(); it != this->intervals.end(); it++) {
            double b = it->first;
            
            std::vector<Segment> segments = it->second;
            std::sort(segments.begin(), segments.end(), 
                [](const Segment& a, const Segment& b){ return a.aMin < b.aMin; });

            BBlock group(b);
            for (Segment& segment : segments) {
                if (group.is_intersect(segment.aMax, segment.aMin)) {
                    group.intersect(segment.aMax, segment.aMin, segment.t);
                }
                else if (group.blocks.back().t.size() > 1) {
                    group.push_back(segment.aMax, segment.aMin, segment.t);
                }
                else {
                    ungrouped.push_back(group.pop_back());
                    group.push_back(segment.aMax, segment.aMin, segment.t);
                }
            }

            if (group.blocks.back().t.size() <= 1) {
                ungrouped.push_back(group.pop_back());
            }
            
            if (group.blocks.size() > 0) this->BBlocks.push_back(group);
        }

        std::sort(ungrouped.begin(), ungrouped.end(), 
            [](const Segment& a, const Segment& b)
            { return a.aMin < b.aMin; });

        ABlock group;
        for (Segment& segment : ungrouped) {
            if (group.is_intersect(segment.aMax, segment.aMin)) {
                group.intersect(segment.b, segment.aMax, segment.aMin, segment.t);
            }
            else if (group.blocks.size() > 1) {
                this->ABlocks.push_back(group);
                group = ABlock(segment.b, segment.aMax, segment.aMin, segment.t); 
            }
            else {
                this->r_blocks.push_back(Segment(group.blocks[0].t, group.a_l, group.a_u, group.blocks[0].b));
                group = ABlock(segment.b, segment.aMax, segment.aMin, segment.t); 
            }
        }

        if (group.blocks.size() > 1) {
            this->ABlocks.push_back(group);
        }
        else {
            this->r_blocks.push_back(Segment(
                group.blocks[0].t, group.a_l, group.a_u, group.blocks[0].b));
        }

        this->yield();
        this->BBlocks.clear();
        this->ABlocks.clear();
        this->r_blocks.clear();
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->index++, data->get_value());

        if (this->index == 1) {
            this->t_s = p.x;
            this->b_1 = std::floor(p.y / this->error) * this->error;
            this->b_2 = std::ceil(p.y / this->error) * this->error;
        }
        else {
            if (p.y > this->slp_u_1 * (p.x - this->t_s) + this->b_1 + this->error || 
                p.y < this->slp_l_1 * (p.x - this->t_s) + this->b_1 - this->error) {
                this->floor_flag = false;
            }
            if (p.y > this->slp_u_2 * (p.x - this->t_s) + this->b_2 + this->error || 
                p.y < this->slp_l_2 * (p.x - this->t_s) + this->b_2 - this->error) {
                this->ceil_flag = false;
            }

            if (this->floor_flag) this->flag++;
            if (this->ceil_flag) this->flag--;

            if (!this->floor_flag && !this->ceil_flag) {
                if (this->flag > 0) {
                    this->intervals[this->b_1].push_back(Segment(p.x - 1, this->slp_l_1, this->slp_u_1, this->b_1));
                }
                else {
                    this->intervals[this->b_2].push_back(Segment(p.x - 1, this->slp_l_2, this->slp_u_2, this->b_2));
                }

                this->flag = 0;
                this->floor_flag = true; 
                this->ceil_flag = true;
                this->curr_segment++;

                this->t_s = p.x;
                this->b_1 = std::floor(p.y / this->error) * this->error;
                this->b_2 = std::ceil(p.y / this->error) * this->error;
                this->slp_u_1 = INFINITY; this->slp_u_2 = INFINITY;
                this->slp_l_1 = -INFINITY; this->slp_l_2 = -INFINITY;
            }

            if (p.y < this->slp_u_1 * (p.x - this->t_s) + this->b_1 - this->error) {
                this->slp_u_1 = (p.y + this->error - this->b_1) / (p.x - this->t_s);
            }
            if (p.y > this->slp_l_1 * (p.x - this->t_s) + this->b_1 + this->error) {
                this->slp_l_1 = (p.y - this->error - this->b_1) / (p.x - this->t_s);
            }
            if (p.y < this->slp_u_2 * (p.x - this->t_s) + this->b_2 - this->error) {
                this->slp_u_2 = (p.y + this->error - this->b_2) / (p.x - this->t_s);
            }
            if (p.y > this->slp_l_2 * (p.x - this->t_s) + this->b_2 + this->error) {
                this->slp_l_2 = (p.y - this->error - this->b_2) / (p.x - this->t_s);
            }

            if (this->curr_segment == this->n_segment) {
                this->curr_segment = 0;
                this->__group();
                this->intervals.clear();
            }
        }                
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        // Do nothing
    }

    void Decompression::finalize() {

    }

    std::pair<CSVObj*, CSVObj*> Decompression::__decompress_segment(Segment& segment) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        Point2D p(this->start, segment.b);
        Line line = Line::line(segment.a, p);
        
        while (this->start <= segment.t) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + this->interval * start));
                base_obj->pushData(std::to_string(start * line.get_slope() + line.get_intercept()));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + this->interval * start));
                obj->pushData(std::to_string(start * line.get_slope() + line.get_intercept()));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }

            this->start++;
        }

        return std::make_pair(base_obj, prev_obj);
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;
        std::vector<Segment> segments;

        // Decompress part 1
        int BBlocks = compress_data->getInt();
        while (BBlocks-- > 0) {
            float b = compress_data->getFloat();
            int ABlocks = compress_data->getInt();
            while (ABlocks-- > 0) {
                float a = compress_data->getFloat();
                int blocks = compress_data->getInt();
                long pivot = 0;
                while (blocks-- > 0) {
                    long time = pivot + VariableByteEncoding::decode(compress_data);
                    segments.push_back(Segment(time, a, b));
                    pivot = time;
                }
            }
        }

        // Decompress part 2
        int ABlocks = compress_data->getInt();
        while (ABlocks-- > 0) {
            float a = compress_data->getFloat();
            int blocks = compress_data->getInt();
            long pivot = 0;
            while (blocks-- > 0) {
                float b = compress_data->getFloat();
                long time = pivot + VariableByteEncoding::decode(compress_data);
                segments.push_back(Segment(time, a, b));
                pivot = time;
            }
        }

        // Decompress part 3
        int r_blocks = compress_data->getInt();
        long pivot = 0;
        while (r_blocks-- > 0) {
            float a = compress_data->getFloat();
            float b = compress_data->getFloat();
            long time = pivot + VariableByteEncoding::decode(compress_data);
            segments.push_back(Segment(time, a, b));
            pivot = time;
        }

        // Decompress segments
        std::sort(segments.begin(), segments.end(), 
            [](const Segment& a, const Segment& b) { return a.t < b.t; });


        for (Segment& segment : segments) {
            std::pair<CSVObj*, CSVObj*> head_tail = __decompress_segment(segment);
            if (base_obj == nullptr) {
                base_obj = head_tail.first;
            }
            else {
                prev_obj->setNext(head_tail.first);
            }
            
            prev_obj = head_tail.second;
        }

        return base_obj;
    }
    // End: decompression
};