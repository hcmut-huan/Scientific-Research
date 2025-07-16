#ifndef ALGEBRAIC_CONVEX_HPP
#define ALGEBRAIC_CONVEX_HPP

#include <vector>
#include "algebraic/function.hpp"

class UpperHull {
    private:
        std::vector<Point2D> points;

    public:
        void append(Point2D point) {
            this->points.push_back(point);

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D& p1 = point;
                Point2D p2 = this->points.at(size-2);
                Point2D p3 = this->points.at(size-3);

                Line line = Line::line(p1, p3);
                if (std::abs(line.subs(p2.x) - p2.y) < 0.00001 || line.subs(p2.x) > p2.y) {
                    this->points.erase(this->points.begin()+size-2);
                }
                else {
                    break;
                }
            }
        }

        void append_backward(Point2D point) {
            this->points.insert(this->points.begin(), point);

            while (this->points.size() > 2) {
                Point2D& p1 = point;
                Point2D p2 = this->points.at(1);
                Point2D p3 = this->points.at(2);

                Line line = Line::line(p1, p3);
                if (std::abs(line.subs(p2.x) - p2.y) < 0.00001 || line.subs(p2.x) > p2.y) {
                    this->points.erase(this->points.begin()+1);
                }
                else {
                    break;
                }
            }
        }

        Point2D at(int i) {
            return this->points.at(i);
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        Point2D pop() {
            Point2D p = this->points.back();
            this->points.pop_back();

            return p;
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }

        void erase_from_end(int index) {
            this->points.erase(this->points.begin() + index, this->points.end());
        }
};

class LowerHull {
    private:
        std::vector<Point2D> points;

    public:
        void append(Point2D point) {
            this->points.push_back(point);

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D p1 = this->points.at(size-1);
                Point2D p2 = this->points.at(size-2);
                Point2D p3 = this->points.at(size-3);

                Line line = Line::line(p1, p3);
                if (std::abs(p2.y - line.subs(p2.x)) < 0.00001 || p2.y > line.subs(p2.x)) {
                    this->points.erase(this->points.begin()+size-2);
                }
                else {
                    break;
                }
            }
        }

        void append_backward(Point2D point) {
            this->points.insert(this->points.begin(), point);

            while (this->points.size() > 2) {
                Point2D p1 = this->points.at(0);
                Point2D p2 = this->points.at(1);
                Point2D p3 = this->points.at(2);

                Line line = Line::line(p1, p3);
                if (std::abs(p2.y - line.subs(p2.x)) < 0.00001 || p2.y > line.subs(p2.x)) {
                    this->points.erase(this->points.begin()+1);
                }
                else {
                    break;
                }
            }
        }

        Point2D at(int i) {
            return this->points.at(i);
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        Point2D pop() {
            Point2D p = this->points.back();
            this->points.pop_back();

            return p;
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }

        void erase_from_end(int index) {
            this->points.erase(this->points.begin() + index, this->points.end());
        }
};

class ConvexHull {
    private:
        LowerHull lower;
        UpperHull upper;

    public:
        void append(Point2D point) {
            this->lower.append(point);
            this->upper.append(point);
        }

        Point2D at(int i) {
            if (i < this->lower.size()) {
                return this->lower.at(i);
            }
            else {
                return this->upper.at(i-this->lower.size()+1);
            }
        }

        Point2D at_cw(int i) {
            if (i < this->upper.size()) {
                return this->upper.at(i);
            }
            else {
                return this->lower.at(this->size()-i);
            }
        }

        Point2D at_ccw(int i) {
            if (i < this->lower.size()) {
                return this->lower.at(i);
            }
            else {
                return this->upper.at(this->size()-i);
            }
        }

        int rightmost_index_cw() {
            return this->upper.size() - 1;
        }

        int rightmost_index_ccw() {
            return this->lower.size() - 1;
        }

        int size() {
            if (this->lower.size() + this->upper.size() == 0) return 0;
            else if (this->lower.size() + this->upper.size() == 2) return 1;
            else return this->lower.size() + this->upper.size() - 2;
        }

        void clear() {
            this->upper.clear();
            this->lower.clear();
        }
};

#endif