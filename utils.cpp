/*
* utils.cpp: Contains utility functions such as normal vector calculations
*/

#include <vector>
#include <cmath>
#include <array>
#include <float.h>
#include "XYZ.h"

std::vector<double> normalize(std::vector<double> vec) {
	double length = std::sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
	std::vector<double> normalized(3);
	if (length > 0) {
		normalized[0] = -vec[0]/length;
		normalized[1] = -vec[1]/length;
		normalized[2] = -vec[2]/length;
	}
	return normalized;
}

std::vector<double> calculateNormal(std::array<XYZ, 3> points) {
	XYZ p1 = points[0];
	XYZ p2 = points[1];
	XYZ p3 = points[2];

	// printf("P1: %f, %f, %f\n", p1.x, p1.y, p1.z);
	// printf("P2: %f, %f, %f\n", p2.x, p2.y, p2.z);
	// printf("P3: %f, %f, %f\n", p3.x, p3.y, p3.z);

	std::vector<double> u{p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};

	// printf("u: %f, %f, %f\n", u[0], u[1], u[2]);

	std::vector<double> v{p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
	// printf("v: %f, %f, %f\n", v[0], v[1], v[2]);

	std::vector<double> n(3);
	n[0] = u[1]*v[2] - u[2]*v[1];
	n[1] = u[2]*v[0] - u[0]*v[2];
	n[2] = u[0]*v[1] - u[1]*v[0];

	// printf("n: %f, %f, %f\n", n[0], n[1], n[2]);

	std::vector<double> normalN =  normalize(n);

	return normalN;
}

double distance(XYZ point, XYZ otherPoint) {
	return sqrt((point.x - otherPoint.x)*(point.x - otherPoint.x)
	+ (point.y - otherPoint.y) * (point.y - otherPoint.y)
	+ (point.z - otherPoint.z)*(point.z - otherPoint.z));
}