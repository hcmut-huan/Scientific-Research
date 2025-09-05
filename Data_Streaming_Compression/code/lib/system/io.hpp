#ifndef SYSTEM_IO_HPP
#define SYSTEM_IO_HPP

#include <string>
#include <fstream>
#include <vector>
#include <sstream>

using Byte = unsigned char;
const char HEX[16] = {'0', '1', '2', '3', '4', '5', '6', '7',
                      '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

class IOObj {
    public:
        virtual IOObj* getNext() = 0;
        virtual void setNext(IOObj* obj) = 0;
        virtual std::string toStr() = 0;

        static void clear(IOObj* head) {
            IOObj* curr = head;
            IOObj* next = nullptr;

            while (curr != nullptr) {
                next = curr->getNext();
                delete curr;
                curr = next;
            }
        }
};

class BinObj : public IOObj {
    private:
        std::vector<Byte> byte_vector;
        BinObj* next = nullptr;

    public:
        void setNext(IOObj* obj) override {
            this->next = (BinObj*) obj;
        }

        IOObj* getNext() override {
            return this->next;
        }

        std::string toStr() override {
            BinObj* curr = this;
            std::string str = "";
            
            while (curr != nullptr) {
                for (Byte byte : this->byte_vector) {
                    str += HEX[(byte & 0xf0) >> 4];
                    str += HEX[(byte & 0x0f) >> 0];
                }

                curr = (BinObj*) curr->getNext();
            }

            return str;
        }

        void put(Byte data) {
            this->byte_vector.push_back(data);
        }

        void put(std::vector<Byte> data) {
            for (Byte& byte : data) {
                this->put(byte);
            }
        }

        Byte getByte() {
            Byte data = this->byte_vector[0];
            this->byte_vector.erase(this->byte_vector.begin());
            return data;
        }

        std::vector<Byte> getBytes(int n, bool last) {
            if (last) {
                std::vector<Byte> bytes(this->byte_vector.end() - n, this->byte_vector.end());
                this->byte_vector.resize(this->byte_vector.size() - n);
                
                return bytes;
            }
            else {
                std::vector<Byte> bytes(this->byte_vector.begin(), this->byte_vector.begin() + n);
                this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + n);

                return bytes;
            }
        }

        void put(short data) {
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        short getShort() {
            short data = 0;
            data += (short) this->byte_vector[0] << 8;
            data += (short) this->byte_vector[1] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 2);
            return data;
        }

        void put(int data) {
            this->byte_vector.push_back((data >> 24) & 0xff);
            this->byte_vector.push_back((data >> 16) & 0xff);
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        int getInt() {
            int data = 0;
            data += (int) this->byte_vector[0] << 24;
            data += (int) this->byte_vector[1] << 16;
            data += (int) this->byte_vector[2] << 8;
            data += (int) this->byte_vector[3] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 4);
            return data;
        }

        void put(long data) {
            this->byte_vector.push_back((data >> 56) & 0xff);
            this->byte_vector.push_back((data >> 48) & 0xff);
            this->byte_vector.push_back((data >> 40) & 0xff);
            this->byte_vector.push_back((data >> 32) & 0xff);
            this->byte_vector.push_back((data >> 24) & 0xff);
            this->byte_vector.push_back((data >> 16) & 0xff);
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        long getLong() {
            long data = 0;
            data += (long) this->byte_vector[0] << 56;
            data += (long) this->byte_vector[1] << 48;
            data += (long) this->byte_vector[2] << 40;
            data += (long) this->byte_vector[3] << 32;
            data += (long) this->byte_vector[4] << 24;
            data += (long) this->byte_vector[5] << 16;
            data += (long) this->byte_vector[6] << 8;
            data += (long) this->byte_vector[7] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 8);
            return data;
        }

        void put(float data) {
            Byte* ptr = (Byte*) &data;
            for (int i=0; i< sizeof(float); i++) {
                this->byte_vector.push_back(ptr[i]);
            }
        }

        float getFloat() {
            float data = 0;
            *((Byte*)(&data) + 3) = this->byte_vector[3];
            *((Byte*)(&data) + 2) = this->byte_vector[2];
            *((Byte*)(&data) + 1) = this->byte_vector[1];
            *((Byte*)(&data) + 0) = this->byte_vector[0];

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 4);
            return data;
        }

        void put(double data) {
            Byte* ptr = (Byte*) &data;
            for (int i = 0; i < sizeof(double); i++) {
                this->byte_vector.push_back(ptr[i]);
            }
        }

        double getDouble() {
            double data = 0;
            *((Byte*)(&data) + 7) = this->byte_vector[7];
            *((Byte*)(&data) + 6) = this->byte_vector[6];
            *((Byte*)(&data) + 5) = this->byte_vector[5];
            *((Byte*)(&data) + 4) = this->byte_vector[4];
            *((Byte*)(&data) + 3) = this->byte_vector[3];
            *((Byte*)(&data) + 2) = this->byte_vector[2];
            *((Byte*)(&data) + 1) = this->byte_vector[1];
            *((Byte*)(&data) + 0) = this->byte_vector[0];

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 8);
            return data;
        }

        void put(std::string data) {
            for (unsigned char c : data) {
                this->byte_vector.push_back(c);
            }
        }

        std::string getString(int length) {
            std::string data = "";
            for (int i = 0; i < length; i++) {
                data += (char) this->byte_vector[i];
            }

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + length);
            return data;
        }

        int getSize() {
            return this->byte_vector.size();
        }

        Byte* serialize() {
            return &this->byte_vector[0];
        }
};

class CSVObj : public IOObj {
    private:
        std::vector<std::string> data;
        CSVObj* next = nullptr;

    public:
        void setNext(IOObj* obj) override {
            this->next = (CSVObj*) obj;
        }

        IOObj* getNext() override {
            return this->next;
        }

        std::string toStr() override {
            std::string str = "";
            for (std::string ele : this->data) {
                str += "," + ele;
            }

            return str.substr(1);
        }

        std::vector<std::string> getData() {
            return this->data;
        }

        std::string getData(int i) {
            return this->data.at(i);
        }

        void pushData(std::string data) {
            this->data.push_back(data);
        }
};

class VariableByteEncoding {
    public:
        static std::vector<Byte> encode(long number) {
            std::vector<Byte> bytes;
            do {
                Byte byte = number % 128; // get last 7 bits
                bytes.insert(bytes.begin(), byte);
                number /= 128;
            } while (number > 0);

            // Set the continuation bit (MSB = 1) on the last byte
            bytes.back() |= 0x80;

            return bytes;
        }

        static long decode(BinObj* obj) {
            long number = 0;

            while (true) {
                Byte byte = obj->getByte();
                if (byte & 0x80) { // this is the last byte (MSB = 1)
                    number = number * 128 + (byte & 0x7F);
                    break;
                }
                number = number * 128 + byte;
            }

            return number;
        }
};

class ZigZagEncoding {
    public:
        static unsigned short encode(short number) {
            return (number << 1) ^ (number >> 15);
        }

        static unsigned int encode(int number) {
            return (number << 1) ^ (number >> 31);
        }

        static short decode(unsigned short number) {
            return (number >> 1) ^ -(number & 1);
        }

        static int decode(unsigned int number) {
            return (number >> 1) ^ -(number & 1);
        }
};

class BatchIO {
    public:
        static long write(std::string filename, IOObj* obj, bool append = false) {
            long count = 0;
            std::ofstream file;

            if (!append) file.open(filename);
            else file.open(filename, std::ios::app);

            while (obj != nullptr) {
                file << obj->toStr();
                obj = obj->getNext();
                if (obj != nullptr) {
                    file << "\n";
                }
                count++;
            }

            file.close();

            return count;
        }

        static CSVObj* readCSV(std::string filename) {
            std::ifstream file(filename);

            std::string line;
            CSVObj* head = nullptr;
            CSVObj* current = nullptr;

            while(std::getline(file, line)) {
                CSVObj* data = new CSVObj;
                std::stringstream ss(line);
                std::string ele;

                while(std::getline(ss, ele, ',')) {
                    data->pushData(ele);
                }

                if (head == nullptr) {
                    head = data;
                    current = head;
                }
                else {
                    current->setNext(data);
                    current = data;
                }
            }

            file.close();
            return head;
        }
};

class IterIO {
    private:
        std::fstream file;
        bool mode;

    public:
        IterIO(std::string filename, bool read, bool binary = false, bool append = false) {
            auto flag = read ? std::ios::in : std::ios::out;
            if (binary) flag = flag | std::ios::binary;
            if (append) flag = flag | std::ios::app;
            
            this->file.open(filename, flag);
        }

        void close() {
            if (this->file.is_open()) {
                this->file.close();
            }
        }

        long write(std::string data, bool endline = true) {
            this->file << data;
            if (endline) {
                this->file << "\n";
            }

            return 1;
        }

        long write(IOObj* obj, bool endline = true) {
            long count = 0;

            while (obj != nullptr) {
                this->file << obj->toStr();
                if (endline) {
                    this->file << "\n";
                }

                obj = obj->getNext();
                count++;
            }

            return count;
        }

        long writeBin(BinObj* obj) {
            long count = 0;
            
            while (obj != nullptr) {
                this->file.write((char*) obj->serialize(), obj->getSize());
                obj = (BinObj*) obj->getNext();
                count++;
            }

            return count;
        }

        CSVObj* readCSV() {
            std::string ele;
            std::string line;
            CSVObj* data = nullptr;
            std::stringstream ss(line);

            if (std::getline(this->file, line)) {
                data = new CSVObj;
                while(std::getline(ss, ele, ',')) {
                    data->pushData(ele);
                }
            }

            return data;
        }

        BinObj* readBin() {
            Byte* byte = new Byte;
            BinObj* obj = new BinObj;

            while (true) {
                this->file.read((char*) byte, sizeof(Byte));
                if (this->file.eof()) break;
                else obj->put(*byte);
            }

            return obj;
        }
};

#endif