/*
@file tinynurbs/io/obj.h
@author Pradeep Kumar Jayaraman <pradeep.pyro@gmail.com>

Read/write NURBS surface data from OBJ files

Use of this source code is governed by a BSD-style license that can be found in
the LICENSE.txt file.
*/

#ifndef IO_H
#define IO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "glm/glm.hpp"
#include "util.h"
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
    for (auto idx : indices) {
        size_t idxm1 = idx - 1;
        size_t i = idxm1 / num_cp_v, j = idxm1 % num_cp_v;
        ctrlPts(i, j) = ctrl_pts_buf[idxm1];
        weights(i, j) = weights_buf[idxm1];
    }

    knots_u = temp_uknots;
    knots_v = temp_vknots;

    return true;
}

template <typename T>
bool readOBJ(const std::string &filename, unsigned int &deg_u, unsigned int &deg_v,
             std::vector<T> &knots_u, std::vector<T> &knots_v,
             std::vector<std::vector<glm::vec<3, T>>> &ctrlPts, std::vector<std::vector<T>> &weights, bool &rational) {
    int num = 1;
    double uknot_min = 0, uknot_max = 1;
    double vknot_min = 0, vknot_max = 1;

    std::vector<glm::vec<3, T>> Posi;
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
    }

    deg_u = 0, deg_v = 0;

    struct ToParse {
        bool deg, cstype, surf, parm;
    };

    ToParse to_parse;

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
                to_parse.cstype = true;
            }
            else if (token1 == "rat") {
                std::string token2;
                ssline >> token2;
                if (token2 == "bspline") {
                    rational = true;
                    to_parse.cstype = true;
                }
            }
        }
        else if (start == "deg") {
            ssline >> deg_u >> deg_v;
            to_parse.deg = true;
        }
        else if (start == "surf") {
            ssline >> uknot_min >> uknot_max >> vknot_min >> vknot_max;
            while (ssline >> token) {
                if (token == "\\") {
                    ssline.clear();
                    std::getline(file, sline);
                    ssline.str(sline);
                }
                else {
                    indices.push_back(std::stof(token));
                }
            }
            to_parse.surf = true;
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
            to_parse.parm = true;
        }
        else if (start == "end") {
            break;
        }
        ssline.clear();
    }
    file.close();

    // Check if necessary data was available in file
    if (!to_parse.cstype) {
        throw std::runtime_error("'cstype bspline / cstype rat bspline' line missing in file");
    }
    if (!to_parse.deg) {
        throw std::runtime_error("'deg' line missing/incomplete in file");
    }
    if (!to_parse.parm) {
        throw std::runtime_error("'parm' line missing/incomplete in file");
    }

    for (int i = 0; i < (int)indices.size(); ++i) {
        Posi.push_back(ctrl_pts_buf[(int)indices[i] - 1]);
    }
    int nknots_u = temp_uknots.size();
    int nknots_v = temp_vknots.size();
    int NIPts = (int)Posi.size();
    int nCtrlPtsU = nknots_u - deg_u - 1;
    int nCtrlPtsV = nknots_v - deg_v - 1;

    /*if (NIPts != nCtrlPtsU*nCtrlPtsV) {
        cerr << "Invalid relation between knots, degree and control points: " <<
             NIPts << " != " << nCtrlPtsU*nCtrlPtsV << endl;
        return false;
    }*/

    ctrlPts.resize(nCtrlPtsU);
    for (int i = 0; i < nCtrlPtsU; i++) {
        ctrlPts[i].resize(nCtrlPtsV);
    }
    num = 0;
    for (int j = 0; j < nCtrlPtsV; j++) {
        for (int i = 0; i < nCtrlPtsU; i++) {
            ctrlPts[i][j] = Posi[num];
            num++;
        }
    }

    weights.resize(nCtrlPtsU);
    for (int i = 0; i < nCtrlPtsU; i++) {
        weights[i].resize(nCtrlPtsV);
    }
    num = 0;
    for (int j = 0; j < nCtrlPtsV; j++) {
        for (int i = 0; i < nCtrlPtsU; i++) {
            weights[i][j] = weights_buf[num];
            num++;
        }
    }

    knots_u.reserve(nknots_u);
    knots_v.reserve(nknots_v);
    auto mnmxu = std::minmax_element(temp_uknots.begin(), temp_uknots.end());
    for (int i = 0; i < nknots_u; ++i) {
        knots_u.push_back(util::mapToRange(temp_uknots[i], *mnmxu.first, *mnmxu.second, 0.0, 1.0));
    }
    auto mnmxv = std::minmax_element(temp_vknots.begin(), temp_vknots.end());
    for (int j = 0; j < nknots_v; ++j) {
        knots_v.push_back(util::mapToRange(temp_vknots[j], *mnmxv.first, *mnmxv.second, 0.0, 1.0));
    }

    return true;
}


/*
static Surface<nd, T> fromOBJ(const std::string &filename) {
    unsigned int deg_u, deg_v;
    std::vector<double> knots_u, knots_v;
    array2<glm::dvec3> ctrl_pts;
    array2<double> weights;
    bool rational;
    readOBJ(filename, deg_u, deg_v, knots_u, knots_v, ctrl_pts, weights, rational);
    return Surface<nd, T>(deg_u, deg_v, knots_u, knots_v, ctrl_pts);
}
*/

template <typename T>
void saveOBJ(const std::string &filename, unsigned int deg_u, unsigned int deg_v, const std::vector<T>& knots_u, const std::vector<T>& knots_v,
             const std::vector<std::vector<glm::vec<3, T>>> &ctrlPts, bool rational) {

    using std::endl;
    std::ofstream fout(filename);

    if (ctrlPts.empty()) {
        return;
    }

    for (int j = 0; j < ctrlPts[0].size(); j++) {
        for (int i = 0; i < ctrlPts.size(); i++) {
            fout << "v " << ctrlPts[i][j].x << " " << ctrlPts[i][j].y << " " << ctrlPts[i][j].z << endl;
        }
    }

    int nknots_u = knots_u.size();
    int nknots_v = knots_v.size();

    int nCpU = ctrlPts.size();
    int nCpV = ctrlPts[0].size();

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
