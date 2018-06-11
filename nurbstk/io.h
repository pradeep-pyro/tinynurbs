#ifndef IO_H
#define IO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "glm/glm.hpp"
#include "util.h"

using namespace std;

template <typename T>
T mapToRange(T val, T old_min, T old_max, T new_min, T new_max) {
    T old_range = old_max - old_min;
    T new_range = new_max - new_min;
    return (((val - old_min) * new_range) / old_range) + new_min;
}

template <typename T>
bool readOBJ(const string &filename, unsigned int &degU, unsigned int &degV,
             vector<T> &knotsU, vector<T> &knotsV,
             nurbstk::util::array2<glm::vec<3, T>> &ctrlPts, nurbstk::util::array2<T> &weights, bool &rational) {
    T uknot_min = 0, uknot_max = 1;
    T vknot_min = 0, vknot_max = 1;

    vector<glm::vec<3, T>> ctrl_pts_buf;
    vector<T> weights_buf;
    vector<int> indices;
    vector<T> temp_uknots;
    vector<T> temp_vknots;

    string start, token;
    string sline;
    istringstream ssline;

    ifstream file(filename);
    if (!file) {
        cerr << "File not found!" << endl;
        return false;
    }

    struct ToParse {
        bool deg, cstype, surf, parm;
    };

    ToParse parsed;

    while (getline(file, sline)) {
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
            ssline >> degU >> degV;
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
                        getline(file, sline);
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
                        getline(file, sline);
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
    int num_cp_u = num_knots_u - degU - 1;
    int num_cp_v = num_knots_v - degV - 1;

    ctrlPts.resize(num_cp_u, num_cp_v);
    weights.resize(num_cp_u, num_cp_v);
    for (auto idx : indices) {
        size_t idxm1 = idx - 1;
        size_t i = idxm1 / num_cp_v, j = idxm1 % num_cp_v;
        ctrlPts(i, j) = ctrl_pts_buf[idxm1];
        weights(i, j) = weights_buf[idxm1];
    }

    knotsU = temp_uknots;
    knotsV = temp_vknots;

    return true;
}

template <typename T>
bool readOBJ(string filename, unsigned int &degU, unsigned int &degV,
             vector<T> &knotsU, vector<T> &knotsV,
             vector<vector<glm::vec<3, T>>> &ctrlPts, vector<vector<T>> &weights, bool &rational) {
    int num = 1;
    double uknot_min = 0, uknot_max = 1;
    double vknot_min = 0, vknot_max = 1;

    vector<glm::vec<3, T>> Posi;
    vector<glm::vec<3, T>> ctrl_pts_buf;
    vector<T> weights_buf;
    vector<int> indices;
    vector<T> temp_uknots;
    vector<T> temp_vknots;

    string start, token;
    string sline;
    istringstream ssline;

    ifstream file(filename);
    if (!file) {
        cerr << "File not found!" << endl;
        return false;
    }

    degU = 0, degV = 0;

    struct ToParse {
        bool deg, cstype, surf, parm;
    };

    ToParse to_parse;

    while (getline(file, sline)) {
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
            ssline >> degU >> degV;
            to_parse.deg = true;
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
            to_parse.surf = true;
        }
        else if (start == "parm") {
            ssline >> start;
            if (start == "u") {
                while (ssline >> token) {
                    if (token == "\\") {
                        ssline.clear();
                        getline(file, sline);
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
                        getline(file, sline);
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
    int nKnotsU = temp_uknots.size();
    int nKnotsV = temp_vknots.size();
    int NIPts = (int)Posi.size();
    int nCtrlPtsU = nKnotsU - degU - 1;
    int nCtrlPtsV = nKnotsV - degV - 1;

    if (NIPts != nCtrlPtsU*nCtrlPtsV) {
        cerr << "Invalid relation between knots, degree and control points: " <<
             NIPts << " != " << nCtrlPtsU*nCtrlPtsV << endl;
        return false;
    }

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

    knotsU.reserve(nKnotsU);
    knotsV.reserve(nKnotsV);
    auto mnmxu = std::minmax_element(temp_uknots.begin(), temp_uknots.end());
    for (int i = 0; i < nKnotsU; ++i) {
        knotsU.push_back(mapToRange(temp_uknots[i], *mnmxu.first, *mnmxu.second, 0.0, 1.0));
    }
    auto mnmxv = std::minmax_element(temp_vknots.begin(), temp_vknots.end());
    for (int j = 0; j < nKnotsV; ++j) {
        knotsV.push_back(mapToRange(temp_vknots[j], *mnmxv.first, *mnmxv.second, 0.0, 1.0));
    }

    return true;
}

template <typename T>
void saveOBJ(const std::string &filename, unsigned int degU, unsigned int degV, const std::vector<T>& knotsU, const std::vector<T>& knotsV,
             const std::vector<std::vector<glm::vec<3, T>>> &ctrlPts, bool rational) {

    ofstream fout(filename);

    if (ctrlPts.empty()) {
        return;
    }

    for (int j = 0; j < ctrlPts[0].size(); j++) {
        for (int i = 0; i < ctrlPts.size(); i++) {
            fout << "v " << ctrlPts[i][j].x << " " << ctrlPts[i][j].y << " " << ctrlPts[i][j].z << endl;
        }
    }

    int nKnotsU = knotsU.size();
    int nKnotsV = knotsV.size();

    int nCpU = ctrlPts.size();
    int nCpV = ctrlPts[0].size();

    if (!rational) {
        fout << "cstype bspline" << endl;
    }
    else {
        fout << "cstype rat bspline" << endl;
    }
    fout << "deg " << degU << " " << degV << endl << "surf ";
    fout << knotsU[degU] << " " << knotsU[nKnotsU - degU - 1] << " "
         << knotsV[degV] << " " << knotsV[nKnotsV - degV - 1];
    for (int i = 0; i < nCpU*nCpV; i++) {
        fout << " " << i + 1;
    }
    fout << endl << "parm u";
    for (auto knot : knotsU) {
        fout << " " << knot;
    }
    fout << endl << "parm v";
    for (auto knot : knotsV) {
        fout << " " << knot;
    }
    fout << endl << "end";
    fout.close();
}

#endif // IO_H
