/*
* implicitFunctions.cpp: contains various kernel functions used to generate scalar field
*/

#include "XYZ.h"
#include "utils.h"
#include <cmath>
#include <iostream>

//scaling factor for radial basis function R
double modifier = 2;

double sphere(XYZ point, XYZ poi) {
	XYZ center;
	// center.x = 1;
	// center.y = 1;
	center.z = 0;
	double value = sqrt(pow(point.x - center.x, 2) + pow(point.y - center.y, 2) + pow(point.z - center.z, 2));
	// printf("Value: %f\n", value);
	return value;  
}

double gaussian(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double a = .5;
	double value = r > 0 ? value = exp(-pow(a, 2) * pow(r, 2)) : 0;
	// printf("Value: %f\n", value);
	return value;
}

double cauchy(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double s = .5;
	double value = r > 0 ? 1/pow(1 + pow(s, 2) * pow(r, 2), 2) : 0;
	// printf("Value: %f\n", value);
	return value;
}

double inverse(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double value = r > 0 ? 1.0/r : 0;
	// printf("Value: %f\n", value);
	return value;
}

double inverseSquared(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double value = r > 0 ? 1.0/pow(r, 2) : 0;
	// printf("Value: %f\n", value);
	return value;
}

double metaballs(XYZ point, XYZ poi) {
	double r1 = distance(point, poi)*2;
	
	double value = 0;

	if (r1 >= 0 && r1 <= 1/3.0) {
		value += 1 - 3 * r1 * r1;
	} else if (r1 > 1/3.0 && r1 <= 1) {
		value += 1.5 * (1 - r1)*(1 - r1);
	}

	// printf("distance in function: %f\n", r1);
	return value;
}

double modifiedMetaballs(XYZ point, XYZ poi) {
	double r1 = distance(point, poi)*2;
	r1 *= modifier;
	
	double value = 0;

	if (r1 >= 0 && r1 <= 1) {
		value += 1 - .5 * r1 * r1;
	} else if (r1 > 1 && r1 <= 2) {
		value += .5 * (2 - r1)*(2 - r1);
	}

	// printf("distance in function: %f\n", r1);
	return value;
}

std::vector<double> metaballGradient(XYZ point, XYZ center) {
    // printf("in metaball gradients\n");
    double r = distance(point, center);

	std::vector<double> gradient = {0, 0, 0};
	if (r >= 0 && r <= 1/3.0) {
		gradient.clear();
		double gX = -6*(point.x - center.x);
		double gY = -6*(point.y - center.y);
		double gZ = -6*(point.z - center.z);
		gradient = {gX, gY, gZ};
	} else if (r > 1/3.0 && r <= 1) {
		gradient.clear();
		double gX = (point.x - center.x) * (3 - (3/r));
		double gY = (point.y - center.y) * (3 - (3/r));
		double gZ = (point.z - center.z) * (3 - (3/r));
		gradient = {gX, gY, gZ};
	}
    
    return gradient;
}

std::vector<double> modifiedMetaballGradient(XYZ point, XYZ center) {
	double r = distance(point, center);

	std::vector<double> gradient = {0, 0, 0};
	if (r >= 0 && r <= 1) {
		gradient.clear();
		// double gX = (-3*(point.x - center.x)*(1 - r))/r;
		// double gY = (-3*(point.y - center.y)*(1 - r))/r;
		// double gZ = (-3*(point.z - center.z)*(1 - r))/r;
		double gX = -1 * modifier * modifier*(point.x - center.x);
		double gY = -1 * modifier * modifier*(point.y - center.y);
		double gZ = -1 * modifier * modifier*(point.z - center.z);
		gradient = {gX, gY, gZ};
	} else if (r > 1 && r <= 2) {
		gradient.clear();
		double gX = (point.x - center.x) * (modifier * modifier - ((2 * modifier)/r));
		double gY = (point.y - center.y) * (modifier * modifier - ((2 * modifier)/r));
		double gZ = (point.z - center.z) * (modifier * modifier - ((2 * modifier)/r));
		gradient = {gX, gY, gZ};
	}
    
    return gradient;
}

double softObjects(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double value = r <= 1 ? 1 - (4/9) * pow(r, 6) + (17/9) * pow(r, 4) - (22/9) * pow(r, 2) : 0;
	// printf("Value: %f\n", value);
	return value;
}

double quartic(XYZ point, XYZ poi) {
	double r = distance(point, poi);
	double value = r <= 1 ? pow(1 - pow(r, 2), 2) : 0;
	// printf("Value: %f\n", value);
	return value;
}