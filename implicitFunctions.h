#ifndef IMPLICITFUNCTIONS_H

#include "XYZ.h"

double sphere(XYZ point, XYZ poi);

double gaussian(XYZ point, XYZ poi);

double cauchy(XYZ point, XYZ poi);

double inverse(XYZ point, XYZ poi);

double inverseSquared(XYZ point, XYZ poi);

double metaballs(XYZ point, XYZ poi);

double modifiedMetaballs(XYZ point, XYZ poi);

double softObjects(XYZ point, XYZ poi);

double quartic(XYZ point, XYZ poi);

std::vector<double> metaballGradient(XYZ point, XYZ center);

std::vector<double> modifiedMetaballGradient(XYZ point, XYZ center);

#define IMPLICITFUNCTIONS_H
#endif
