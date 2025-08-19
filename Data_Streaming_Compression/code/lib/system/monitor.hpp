#ifndef SYSTEM_MONITOR_HPP
#define SYSTEM_MONITOR_HPP

#include <string>
#include <fstream>
#include <chrono>
#include <unistd.h>
#include <thread>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std::chrono;

class Monitor {
    private:
        bool flag;
        long page_size;
        std::thread task;

        Monitor() {
            this->flag = false;
            this->page_size = sysconf(_SC_PAGE_SIZE);
        }

        void __monitor(std::string output) {
            std::fstream file(output, std::ios::out);
            file << "user_cpu_time,system_cpu_time,vsz,rss\n";

            while (this->flag) {
                // dont know why but monitoring cpu time with ruuage seem more accurate
                rusage usage;
                getrusage(RUSAGE_SELF, &usage);  

                std::string data;
                std::ifstream ifs("/proc/self/stat", std::ios_base::in);

                // ignore
                ifs >> data >> data >> data >> data >> data >> data >> data >> data >> data >> data >> data 
                >> data >> data >> data >> data >> data >> data >> data >> data >> data >> data >> data; 
                
                // get necessary data
                ifs >> data; std::string vsz = data;
                ifs >> data; std::string rss = std::to_string(std::stol(data)*this->page_size);
                std::string utime = std::to_string(usage.ru_utime.tv_sec*1000000 + usage.ru_utime.tv_usec);
                std::string stime = std::to_string(usage.ru_stime.tv_sec*1000000 + usage.ru_stime.tv_usec);

                file << utime << "," << stime << "," << vsz << "," << rss + "\n";
            }

            file.close();
        }
    
    public:
        static Monitor instance;
        
        void start(std::string output) {
            this->flag = true;
            this->task = std::thread(&Monitor::__monitor, this, output);
        }

        void stop() {
            this->flag = false;
            this->task.join();
        }

};

class Clock {
    private:
        bool flag = true;
        int _counter = 0;
        double _avg_duration = 0; 
        long _max_duration = -1;

        high_resolution_clock::time_point _start_pivot;
        high_resolution_clock::time_point _tick_pivot;
        
    public:
        void start() {
            this->_counter = 0;
            this->_avg_duration = 0;
            this->_max_duration = -1;
            this->_start_pivot = high_resolution_clock::now();
            this->_tick_pivot = this->_start_pivot;
        }

        void tick() {
            if (this->flag) {
                this->_tick_pivot = high_resolution_clock::now();
            }
            else {
                high_resolution_clock::time_point curr = high_resolution_clock::now();
                long duration = duration_cast<nanoseconds>(curr - this->_tick_pivot).count();

                this->_avg_duration = ((this->_counter) * this->_avg_duration + duration) / ((double) (this->_counter + 1));
                this->_max_duration = duration > this->_max_duration ? duration : this->_max_duration;

                this->_counter++;
                this->_tick_pivot = curr;
            }

            this->flag = !this->flag;
        }

        double getAvgDuration() {
            return this->_avg_duration;
        }

        long getMaxDuration() {
            return this->_max_duration;
        }

        long stop() {
            return duration_cast<nanoseconds>(high_resolution_clock::now() - this->_start_pivot).count();
        }
};

#endif