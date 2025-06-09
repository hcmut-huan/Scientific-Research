#ifndef TIMESERIES_HPP
#define TIMESERIES_HPP

#include <vector>

// Source for manipulating time series data
class Data {};

class Univariate : public Data {
    private:
        std::time_t time;
        double value;
    
    public:
        Univariate(std::time_t time, double val) {
            this->time = time;
            this->value = val;
        }

        std::time_t get_time() const {
            return this->time;
        }

        std::tm* timestamp() const {
            return std::localtime(&this->time);
        }

        double get_value() const {
            return this->value;
        }
};

class TimeSeries {
    private:
        int _index = 0;
        std::vector<Data*> series;

    public:
        void push(Data* data) {
            this->series.push_back(data);
        }

        Data* get(int i) {
            return this->series.at(i);
        }

        Data* next() {
            return this->series.at(this->_index++);
        }

        bool hasNext() {
            return this->_index < this->series.size();
        }

        void reset() {
            this->_index = 0;
        }

        int size() {
            return this->series.size();
        }

        void finalize() {
            for (Data* data : this->series) {
                delete data;
            }
        }
};

#endif