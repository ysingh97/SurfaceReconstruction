#ifndef Included_NameModel_H

#include<vector>
#include<string>

class XYZ {
    public:

    XYZ() {
       x = 0;
       y = 0;
       z = 0;
       normal = {0, 0, 0};
    }

    XYZ(double X, double Y, double Z, int I, int J, int K) {
        x = X;
        y = Y;
        z = Z;
        ijk.push_back(I);
        ijk.push_back(J);
        ijk.push_back(K);
        normal = {0, 0, 0};
    }

    std::string toString() {
        std::string stringRep;
        stringRep += std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
        return stringRep;
    }

    double x;
    double y;
    double z;
    int id;
    std::vector<int> ijk;
    std::vector<double> normal;


};

#define Included_NameModel_H
#endif