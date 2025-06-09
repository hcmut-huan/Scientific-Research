// #include "piecewise-approximation/linear.hpp"

// namespace SlideFilter {
    
//     Clock clock;

//     void __yield(BinObj* obj, double index, float value) {
//         obj->put(index);
//         obj->put(value);
//     }

//     Line __fit(std::vector<Point2D>& segment, Line& u, Line& l) {
//         if (u.get_slope() - l.get_slope() < 0.0000001) {
//             return Line(u.get_slope(), u.get_intercept());
//         }
//         else {
//             Point2D p = Line::intersection(u, l);
//             double A_num = 0, A_den = 0;

//             for (Point2D& data : segment) {
//                 A_num += (data.y-p.y)*(data.x-p.x);
//                 A_den += (data.x-p.x)*(data.x-p.x);
//             }

//             double A_ig = A_num / A_den;
//             double temp = A_ig > u.get_slope() ? u.get_slope() : A_ig;
//             double a_ig = temp > l.get_slope() ? temp : l.get_slope();
//             double b_ig = p.y - a_ig * p.x;

//             return Line(a_ig, b_ig);
//         }   
//     }

//     // check if [a, b] interval exists and return [a, b]
//     // ** important: 
//     //      Something went wrong need to be revisited in both implementation and original paper
//     std::pair<double, double> interval_exist(Line& prev_g, Line& prev_u, Line& prev_l, Line& u_line, Line& l_line, double t_j_k) {
//         std::pair<double, double> interval = std::make_pair(-1, -1);

//         if (u_line.get_slope() - l_line.get_slope() > 0.0000001) {
//             Point2D p = Line::intersection(u_line, l_line);

//             if (prev_g.subs(p.x)-p.y>0) {
//                 // point below g^{k-1}
//                 double f_i_k = Line::intersection(l_line, prev_g).x;
//                 if (l_line.subs(t_j_k)-prev_l.subs(t_j_k)>0) {
//                     // find hyperplane s
//                     Point2D s1(t_j_k, prev_l.subs(t_j_k));
//                     Point2D s2 = Line::intersection(u_line, l_line);
//                     Line s = Line::line(s1, s2);

//                     double c_i = Line::intersection(prev_g, u_line).x;     // c_i^k
//                     double d_i = Line::intersection(prev_g, s).x;     // d_i^k
//                     double max_c_vs_d = (c_i > d_i) ? c_i : d_i;

//                     if (max_c_vs_d < t_j_k && f_i_k < t_j_k) {
//                         interval = max_c_vs_d > f_i_k ? std::make_pair(f_i_k, max_c_vs_d) : std::make_pair(max_c_vs_d, f_i_k);
//                     }
//                 }
                
//             }
//             else if (p.y-prev_g.subs(p.x)>0) {
//                 // point above g^{k-1}
//                 double f_i_k = Line::intersection(u_line, prev_g).x;                
//                 if (prev_u.subs(t_j_k)-u_line.subs(t_j_k)>0) {
//                     // find hyperplane q
//                     Point2D q1(t_j_k, prev_u.subs(t_j_k));
//                     Point2D q2 = Line::intersection(u_line, l_line);
//                     Line q = Line::line(q1, q2);

//                     double c_i = Line::intersection(prev_g, l_line).x;   // c_i^k'
//                     double d_i = Line::intersection(prev_g, q).x;   // d_i^k'
//                     double max_c_vs_d = (c_i > d_i) ? c_i : d_i;

//                     if (max_c_vs_d < t_j_k && f_i_k < t_j_k) {
//                         interval = max_c_vs_d > f_i_k ? std::make_pair(f_i_k, max_c_vs_d) : std::make_pair(max_c_vs_d, f_i_k);
//                     }
//                 }
                
//             }
//         }
        
//         return interval;
//     }

//     void compress(TimeSeries& timeseries, float bound, std::string output) {
//         clock.start();
//         IterIO outputFile(output, false);
//         BinObj* compress_data = new BinObj;

//         Univariate* d1 = (Univariate*) timeseries.next();
//         compress_data->put(d1->get_time());
//         Point2D p1(0, d1->get_value());
        
//         Univariate* d2 = (Univariate*) timeseries.next();
//         Point2D p2(1, d2->get_value());
        
//         Line u_line = Line::line(Point2D(p1.x, p1.y-bound), Point2D(p2.x, p2.y+bound));
//         Line l_line = Line::line(Point2D(p1.x, p1.y+bound), Point2D(p2.x, p2.y-bound));
//         ConvexHull cvx; cvx.append(p1); cvx.append(p2);
        
//         Line prev_u(0, 0); Line prev_l(0, 0); Line prev_g(0, 0); long prev_i = 0;
//         std::vector<Point2D> segment = {p1, p2};
//         long index = 2;
//         bool first = true;
//         while (timeseries.hasNext()) {
//             Point2D p(index, ((Univariate*) timeseries.next())->get_value());
            
//             if (p.y-u_line.subs(p.x)>bound || l_line.subs(p.x)-p.y>bound) {
//                 if (!first) {
//                     std::pair<double, double> interval = interval_exist(prev_g, prev_u, prev_l, u_line, l_line, prev_i);
                    
//                     if (interval.first != -1) {
//                         Point2D z = Line::intersection(u_line, l_line);
//                         if (prev_g.subs(z.x)-z.y>0) {
//                             u_line = Line::line(z, Point2D(interval.first, prev_g.subs(interval.first)));
//                             l_line = Line::line(z, Point2D(interval.second, prev_g.subs(interval.second)));
//                         }
//                         else if(z.y-prev_g.subs(z.x)>0) {
//                             u_line = Line::line(z, Point2D(interval.second, prev_g.subs(interval.second)));
//                             l_line = Line::line(z, Point2D(interval.first, prev_g.subs(interval.first)));
//                         }
//                     }

//                     Line g_k = __fit(segment, u_line, l_line);
//                     if (interval.first != -1) {
//                         Point2D knot = Line::intersection(g_k, prev_g);
//                         __yield(compress_data, (double) knot.x, knot.y);
//                     }
//                     else {
//                         Point2D p1 = Point2D(prev_i, prev_g.subs(prev_i));
//                         Point2D p2 = Point2D(prev_i+1, g_k.subs(prev_i+1));

//                         __yield(compress_data, (double) -p1.x, p1.y);
//                         __yield(compress_data, (double) p2.x, p2.y);
//                     }
//                     prev_g = g_k;
//                 }
//                 else {
//                     first = false;
//                     Line g_k = __fit(segment, u_line, l_line);
//                     __yield(compress_data, (double) 0.0, g_k.subs(0));
//                     prev_g = g_k;
//                 }

//                 if (timeseries.hasNext()) {
//                     prev_u = u_line, prev_l = l_line; prev_i = index-1;

//                     Point2D next_p(++index, ((Univariate*) timeseries.next())->get_value());
//                     u_line = Line::line(Point2D(p.x, p.y-bound), Point2D(next_p.x, next_p.y+bound));
//                     l_line = Line::line(Point2D(p.x, p.y+bound), Point2D(next_p.x, next_p.y-bound));

//                     index++;
//                     segment = {p, next_p};
//                     cvx.clear(); cvx.append(p); cvx.append(next_p);
//                 }
//             }
//             else {
//                 cvx.append(p);
//                 if (p.y-l_line.subs(p.x)>bound) {
//                     for (int i=0; i<cvx.size(); i++) {
//                         Point2D cvx_p = cvx.at(i);
//                         if (cvx_p.x != p.x) {
//                             Line l = Line::line(Point2D(cvx_p.x, cvx_p.y+bound), Point2D(p.x, p.y-bound));
//                             l_line = l.get_slope() > l_line.get_slope() ? l: l_line;
//                         }    
//                     }
//                 }
//                 if (u_line.subs(p.x)-p.y>bound) {
//                     for (int i=0; i<cvx.size(); i++) {
//                         Point2D cvx_p = cvx.at(i);
//                         if (cvx_p.x != p.x) {
//                             Line l = Line::line(Point2D(cvx_p.x, cvx_p.y-bound), Point2D(p.x, p.y+bound));
//                             u_line = l.get_slope() < u_line.get_slope() ? l: u_line;
//                         }
//                     }
//                 }
                
//                 index++;
//                 segment.push_back(p);
//             }
//         }

//         outputFile.writeBin(compress_data);
//         outputFile.close();
//         delete compress_data;
        
//         clock.tick();
//         double avg_time = clock.getAvgDuration() / timeseries.size();

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each data point (ns): " << avg_time << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each data point (ns): " + std::to_string(avg_time));
//         timeFile.close();
//     }

    
//     void __decompress_segment(IterIO& file, int interval, time_t basetime, Point2D p1, Point2D p2) {
//         if (p2.x < 0) {
//             Line l = Line::line(p1, Point2D(-p2.x, p2.y));
//             long start = std::ceil(p1.x);
//             float end = -p2.x;
//             while (start <= end) {
//                 CSVObj obj;
//                 obj.pushData(std::to_string(basetime + interval * start));
//                 obj.pushData(std::to_string(start * l.get_slope() + l.get_intercept()));
//                 file.write(&obj);
//                 start++;
//             }
//         }
//         else {
//             Line l = Line::line(p1, p2);
//             long start = std::ceil(p1.x);
//             double end = p2.x;
//             while (start < end) {
//                 CSVObj obj;
//                 obj.pushData(std::to_string(basetime + interval * start));
//                 obj.pushData(std::to_string(start * l.get_slope() + l.get_intercept()));
//                 file.write(&obj);
//                 start++;
//             }
//         }
//     }

//     void decompress(std::string input, std::string output, int interval) {
//         IterIO inputFile(input, true, true);
//         IterIO outputFile(output, false);
//         BinObj* compress_data = inputFile.readBin();

//         time_t basetime = compress_data->getLong();
//         float s_index = compress_data->getDouble();
//         float s_value = compress_data->getFloat();
//         Point2D start(s_index, s_value);

//         clock.start();
//         while (compress_data->getSize() != 0) {
//             float e_index = compress_data->getDouble();
//             float e_value = compress_data->getFloat();
//             Point2D end = Point2D(e_index, e_value);

//             if (start.x >= 0) {
//                 __decompress_segment(outputFile, interval, basetime, start, end);
//             }
//             start = end;

        
//             clock.tick();
//         }

//         delete compress_data;
//         inputFile.close();
//         outputFile.close();

//         // Profile average latency
//         std::cout << std::fixed << "Time taken for each segment (ns): " << clock.getAvgDuration() << "\n";
//         IterIO timeFile(output+".time", false);
//         timeFile.write("Time taken for each segment (ns): " + std::to_string(clock.getAvgDuration()));
//         timeFile.close();
//     }
    
// };

