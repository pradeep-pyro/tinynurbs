#ifndef IO_H
#define IO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;
/*
bool loadOBJ(string filename, unsigned int &degU, unsigned int &degV,
    vector<double> &knotsU, vector<double> &knotsV,
    vector<vector<glm::dvec3>> &ctrlPts) {
    int k = 0, num = 1;
    double uknot_min = 0, uknot_max = 1;
    double vknot_min = 0, vknot_max = 1;
    double temp;
    glm::vec3 pt;

    vector<glm::vec3> Posi;
    vector<glm::vec3> temp_Posi;
    vector<int> Posi_num;
    vector<double> temp_uknots;
    vector<double> temp_vknots;

    string start, token;
    string sline;
    istringstream ssline;

    ifstream file(filename);
    if (!file) {
        cerr << "File not found!" << endl;
        return false;
    }

    degU = 0, degV = 0;
    int nKnotsU = 0, nKnotsV = 0;

    while (!file.eof()) {
        getline(file, sline);
        if (sline.size() == 0) {
            break;
        }
        ssline.str(sline);
        ssline >> start;
        if (start == "v") {
            ssline >> pt[0] >> pt[1] >> pt[2];
            temp_Posi.push_back(pt);
        }
        else if (start == "deg")
            ssline >> degU >> degV;
        else if (start == "surf") {
            ssline >> uknot_min >> uknot_max >> vknot_min >> vknot_max;
            while (!ssline.eof()) {
                ssline >> token;
                if (token == "\\") {
                    ssline.clear();
                    getline(file, sline);
                    ssline.str(sline);
                }
                else {
                    num = atof(token.c_str());
                    Posi_num.push_back(num);
                }
            }
        }
        else if (start == "parm") {
            ssline >> start;
            if (start == "u") {
                while (!ssline.eof()) {
                    ssline >> token;
                    if (token == "\\") {
                        ssline.clear();
                        getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp = atof(token.c_str());
                        temp_uknots.push_back(temp);
                        ++nKnotsU;
                    }
                }
            }
            else if (start == "v") {
                while (!ssline.eof()) {
                    ssline >> token;
                    if (token == "\\") {
                        ssline.clear();
                        getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp = atof(token.c_str());
                        temp_vknots.push_back(temp);
                        ++nKnotsV;
                    }
                }
            }
        }
        else if (start == "end") {
            break;
        }
        ssline.clear();
    }
    file.close();

    cout << temp_uknots[0] << endl;

    for (int i = 0; i < (int)Posi_num.size(); i++) {
        Posi.push_back(temp_Posi[(int)Posi_num[i] - 1]);
    }

    int NIPts = (int)Posi.size();
    int nCtrlPtsU = nKnotsU - (degU + 1);
    int nCtrlPtsV = nKnotsV - (degV + 1);

    if (NIPts != nCtrlPtsU*nCtrlPtsV) {
        cerr << "Invalid OBJ: " << endl;
        cerr << "Degree: " << degU << ", " << degV << endl;
        cerr << "Knots: " << nKnotsU << ", " << nKnotsV << endl;
        cerr << "Control points: " << nCtrlPtsU << ", " << nCtrlPtsV << endl;
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
    for (int j = 0; j < nCtrlPtsV; j++) {
        for (int i = 0; i < nCtrlPtsU; i++) {
            ctrlPts[i][j] = Posi[num];
            num++;
        }
    }
    knotsU.reserve(nKnotsU);
    knotsV.reserve(nKnotsV);
    for (int i = 0; i < nKnotsU; ++i) {
        knotsU.push_back((temp_uknots[i] - uknot_min) / (uknot_max - uknot_min));
    }
    for (int j = 0; j < nKnotsV; ++j) {
        knotsV.push_back((temp_vknots[j] - vknot_min) / (vknot_max - vknot_min));
    }

    return true;
    //MaxMinMidPoint();
}

bool loadOBJ(string filename, unsigned int &degU, unsigned int &degV,
    vector<double> &knotsU, vector<double> &knotsV,
    vector<vector<glm::dvec3>> &ctrlPts) {
    int k = 0, num = 1;
    double uknot_min = 0, uknot_max = 1;
    double vknot_min = 0, vknot_max = 1;
    double temp;
    glm::vec3 pt;

    vector<glm::vec3> Posi;
    vector<glm::vec3> temp_Posi;
    vector<int> Posi_num;
    string start, token;
    string sline;
    istringstream ssline;

    ifstream file(filename);
    if (!file) {
        cerr << "File not found!" << endl;
        return false;
    }

    degU = 0, degV = 0;

    while (!file.eof()) {
        getline(file, sline);
        if (sline.size() == 0) {
            break;
        }
        ssline.str(sline);
        ssline >> start;
        if (start == "v") {
            ssline >> pt[0] >> pt[1] >> pt[2];
            temp_Posi.push_back(pt);
        }
        else if (start == "deg")
            ssline >> degU >> degV;
        else if (start == "surf") {
            ssline >> uknot_min >> uknot_max >> vknot_min >> vknot_max;
            while (!ssline.eof()) {
                ssline >> token;
                if (token == "\\") {
                    ssline.clear();
                    getline(file, sline);
                    ssline.str(sline);
                }
                else {
                    num = atof(token.c_str());
                    Posi_num.push_back(num);
                }
            }
        }
        else if (start == "parm") {
            ssline >> start;
            if (start == "u") {
                while (!ssline.eof()) {
                    ssline >> token;
                    if (token == "\\") {
                        ssline.clear();
                        getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp = atof(token.c_str());
                        knotsU.push_back(temp);
                    }
                }
            }
            else if (start == "v") {
                while (!ssline.eof()) {
                    ssline >> token;
                    if (token == "\\") {
                        ssline.clear();
                        getline(file, sline);
                        ssline.str(sline);
                    }
                    else {
                        temp = atof(token.c_str());
                        knotsV.push_back(temp);
                    }
                }
            }
        }
        else if (start == "end") {
            break;
        }
        ssline.clear();
    }
    file.close();

    for (int index : Posi_num) {
        Posi.push_back(temp_Posi[index - 1]);
    }

    cout << Posi_num.size() << "  " << Posi.size() << endl;

    int NIPts = (int)Posi.size();
    int nCtrlPtsU = knotsU.size() - (degU + 1);
    int nCtrlPtsV = knotsV.size() - (degV + 1);

    if (NIPts != nCtrlPtsU*nCtrlPtsV) {
    cerr << "Invalid OBJ: " << endl;
    cerr << "NIPts: " << NIPts << " != " << nCtrlPtsU << " * " << nCtrlPtsV << endl;
    cerr << "Degree: " << degU << ", " << degV << endl;
    cerr << "Knots: " << knotsU.size() << ", " << knotsV.size() << endl;
    for (auto val : knotsU) {
        cout << val << " ";
    }
    cout << endl;
    cerr << "Control points: " << nCtrlPtsU << " * " << nCtrlPtsV << " = " << NIPts << endl;
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
    for (int j = 0; j < nCtrlPtsV; j++) {
        for (int i = 0; i < nCtrlPtsU; i++) {
            ctrlPts[i][j] = Posi[num];
            num++;
        }
    }

    for (double &val : knotsU) {
        val = (val - uknot_min) / (uknot_max - uknot_min);
    }
    for (double &val : knotsV) {
        val = (val - vknot_min) / (vknot_max - vknot_min);
    }

    return true;
}
*/

template <typename T>
T mapToRange(T val, T old_min, T old_max, T new_min, T new_max) {
    T old_range = old_max - old_min;
    T new_range = new_max - new_min;
    return (((val - old_min) * new_range) / old_range) + new_min;
}

bool readOBJ(string filename, unsigned int &degU, unsigned int &degV,
             vector<double> &knotsU, vector<double> &knotsV,
             vector<vector<glm::dvec3>> &ctrlPts, bool &rational) {
    int num = 1;
    double uknot_min = 0, uknot_max = 1;
    double vknot_min = 0, vknot_max = 1;

    vector<glm::vec3> Posi;
    vector<glm::vec3> ctrl_pts_buf;
    vector<int> indices;
    vector<double> temp_uknots;
    vector<double> temp_vknots;

    string start, token;
    string sline;
    istringstream ssline;

    ifstream file(filename);
    if (!file) {
        cerr << "File not found!" << endl;
        return false;
    }

    degU = 0, degV = 0;

    bool cstype_parsed = false;

    while (getline(file, sline)) {
        if (sline.size() == 0) {
            break;
        }
        ssline.str(sline);
        ssline >> start;
        if (start == "v") {
            glm::vec4 pt;
            pt[3] = 1.0;
            for (int i = 0; i < 4 && !ssline; ++i) {
                ssline >> pt[i];
            }
            /*
             * glm::vec3 pt;
            ssline >> pt[0] >> pt[1] >> pt[2];
            */
            ctrl_pts_buf.push_back(pt);
        }
        else if (start == "cstype") {
            std::string token1;
            ssline >> token1;
            if (token1 == "bspline") {
                rational = false;
                cstype_parsed = true;
            }
            else if (token == "rat") {
                std::string token2;
                ssline >> token2;
                if (token2 == "bspline") {
                    rational = true;
                    cstype_parsed = true;
                }
            }
        }
        else if (start == "deg") {
            ssline >> degU >> degV;
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
                        temp_uknots.push_back(atof(token.c_str()));
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
                        temp_vknots.push_back(atof(token.c_str()));
                    }
                }
            }
        }
        else if (start == "end") {
            break;
        }
        ssline.clear();
    }
    file.close();

    if (!cstype_parsed) {
        std::cerr << "'cstype bspline / cstype rat bspline' line missing in file" << std::endl;
        return false;
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
    //MaxMinMidPoint();
}

void saveOBJ(std::string filename, unsigned int degU, unsigned int degV, const std::vector<double>& knotsU, const std::vector<double>& knotsV,
             const std::vector<std::vector<glm::dvec3>>& ctrlPts) {

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

    fout << "cstype bspline" << endl << "deg " << degU << " " << degV << endl << "surf ";
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
