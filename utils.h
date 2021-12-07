#ifndef UTILS_H

#include <vector>
#include <array>
#include "XYZ.h"

std::vector<double> normalize(std::vector<double> vec);

std::vector<double> calculateNormal(std::array<XYZ, 3> points);

double distance(XYZ point, XYZ otherPoint);

std::string normalToString(std::vector<double> normal);

extern double modifier;

#define UTILS_H
#endif