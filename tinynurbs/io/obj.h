#ifndef IO_H
#define IO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "glm/glm.hpp"
#include "../util/util.h"
#include "../util/array2.h"

namespace tinynurbs {

template <typename T>
bool readOBJ(const std::string &filename, unsigned int &deg_u, unsigned int &deg_v,
             std::vector<T> &knots_u, std::vector<T> &knots_v,
             array2<glm::vec<3, T>> &ctrlPts, array2<T> &weights, bool &rational) {
    T uknot_min = 0, uknot_max = 1;
    T vknot_min = 0, vknot_max = 1;

    std::vector<glm::vec<3, T>> ctrl_pts_buf;
    std::vector<T> weights_buf;
    std::vector<int> indices;
    std::vector<T> temp_uknots;
    std::vector<T> temp_vknots;

    std::string start, token, sline;
    std::istringstream ssline;

    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("File not found: " + filename);
        return false;
    }

    struct ToParse {
        bool deg, cstype, surf, parm;
    };

    ToParse parsed;

    while (std::getline(file, sline)) {
        if (sline.size() == 0) {
            break;
        }
        ssline.str(sline);
        ssline >> start;
        if (start == "v") {
            std::vector<double> four_coords;
            four_coords.resize(4);
            four_coords[3] = 1.0;
            int index = 0;
            while (ssline && index <= 3) {
                ssline >> four_coords[index++];
            }
            ctrl_pts_buf.emplace_back(four_coords[0], four_coords[1], four_coords[2]);
            weights_buf.push_back(four_coords[3]);
        }
        else if (start == "cstype") {
            std::string token1;
            ssline >> token1;
            if (token1 == "bspline") {
                rational = false;
                parsed.cstype = true;
            }
            else if (token1 == "rat") {
                std::string token2;
                ssline >> token2;
                if (token2 == "bspline") {
                    rational = true;
                    parsed.cstype = true;
                }
            }
        }
        else if (start == "deg") {
            ssline >> deg_u >> deg_v;
            parsed.deg = true;
        }
        else if (start == "surf") {
            ssline >> uknot_min >> uknot_max >> vknot_min >> vknot_max;
            while (ssline >> token) {
                if (token == "\\") {
                    ssline.clear();
                    getline(file, sline);
                    ssline.str(sline);
                }
                else {
                    indices.push_back(std::stof(token));
                }
            }
            parsed.surf = true;
        }
        else if (start == "parm") {
            ssline >> start;
            if (start == "u") {
                while (ssline >> token) {
                    if (token == "\\") {
                        ssline.clear();
                        std::getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp_uknots.push_back(std::stof(token));
                    }
                }
            }
            else if (start == "v") {
                while (ssline >> token) {
                    if (token == "\\") {
                        ssline.clear();
                        std::getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp_vknots.push_back(std::stof(token));
                    }
                }
            }
            parsed.parm = true;
        }
        else if (start == "end") {
            break;
        }
        ssline.clear();
    }
    file.close();

    // Check if necessary data was available in file
    if (!parsed.cstype) {
        throw std::runtime_error("'cstype bspline / cstype rat bspline' line missing in file");
    }
    if (!parsed.deg) {
        throw std::runtime_error("'deg' line missing/incomplete in file");
    }
    if (!parsed.parm) {
        throw std::runtime_error("'parm' line missing/incomplete in file");
    }

    int num_knots_u = temp_uknots.size();
    int num_knots_v = temp_vknots.size();
    int num_cp_u = num_knots_u - deg_u - 1;
    int num_cp_v = num_knots_v - deg_v - 1;

    ctrlPts.resize(num_cp_u, num_cp_v);
    weights.resize(num_cp_u, num_cp_v);
    size_t num = 0;
    for (int j = 0; j < num_cp_v; ++j) {
        for (int i = 0; i < num_cp_u; ++i) {
            assert(i < ctrlPts.rows() && j < ctrlPts.cols());
            ctrlPts(i, j) = ctrl_pts_buf[indices[num] - 1];
            weights(i, j) = weights_buf[indices[num] - 1];
            ++num;
        }
    }

    knots_u = temp_uknots;
    knots_v = temp_vknots;

    return true;
}

template <typename T>
void saveOBJ(const std::string &filename, unsigned int deg_u, unsigned int deg_v, const std::vector<T>& knots_u, const std::vector<T>& knots_v,
             const array2<glm::vec<3, T>> &ctrlPts, const array2<T> &weights, bool rational) {

    using std::endl;
    std::ofstream fout(filename);

    if (ctrlPts.rows() == 0 || ctrlPts.cols() == 0) {
        return;
    }

    for (int j = 0; j < ctrlPts.cols(); j++) {
        for (int i = 0; i < ctrlPts.rows(); i++) {
            fout << "v " << ctrlPts(i, j).x << " " << ctrlPts(i, j).y << " " << ctrlPts(i, j).z << " " << weights(i, j) << endl;
        }
    }

    int nknots_u = knots_u.size();
    int nknots_v = knots_v.size();

    int nCpU = ctrlPts.rows();
    int nCpV = ctrlPts.cols();

    if (!rational) {
        fout << "cstype bspline" << endl;
    }
    else {
        fout << "cstype rat bspline" << endl;
    }
    fout << "deg " << deg_u << " " << deg_v << endl << "surf ";
    fout << knots_u[deg_u] << " " << knots_u[nknots_u - deg_u - 1] << " "
         << knots_v[deg_v] << " " << knots_v[nknots_v - deg_v - 1];
    for (int i = 0; i < nCpU*nCpV; i++) {
        fout << " " << i + 1;
    }
    fout << endl << "parm u";
    for (auto knot : knots_u) {
        fout << " " << knot;
    }
    fout << endl << "parm v";
    for (auto knot : knots_v) {
        fout << " " << knot;
    }
    fout << endl << "end";
    fout.close();
}

} // namespace tinynurbs

#endif // IO_H
