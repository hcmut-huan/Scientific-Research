#ifndef ALGEBRAIC_CONVEX_H
#define ALGEBRAIC_CONVEX_H

#include <vector>
#include "algebraic/function.h"

class UpperConvexHull {
    private:
        std::vector<Point2D> points;

    public:
        void append(const Point2D& point) {
            this->points.push_back(point);

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D p1 = this->points.at(size-1);
                Point2D p2 = this->points.at(size-2);
                Point2D p3 = this->points.at(size-3);

                Line* line = Line::line(&p1, &p3);
                if (line->substitute(p2.x) >= p2.y) {
                    this->points.erase(this->points.begin()+size-2);
                    delete line;
                }
                else {
                    delete line;
                    break;
                }
            }
        
        }

        Point2D* at(int i) {
            return &this->points.at(i);
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }
};

class LowerConvexHull {
    private:
        std::vector<Point2D> points;

    public:
        void append(const Point2D& point) {
            this->points.push_back(point);

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D p1 = this->points.at(size-1);
                Point2D p2 = this->points.at(size-2);
                Point2D p3 = this->points.at(size-3);

                Line* line = Line::line(&p1, &p3);
                if (line->substitute(p2.x) <= p2.y) {
                    this->points.erase(this->points.begin()+size-2);
                    delete line;
                }
                else {
                    delete line;
                    break;
                }
            }
        
        }

        Point2D* at(int i) {
            return &this->points.at(i);
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }
};

#endif