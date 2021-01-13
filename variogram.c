// Original File: Variogram.h
// Author: Rodolfo Mora Z.
// Date: Jun 8, 2016
// Purpose: Multiple semivariogram models for Kriging Implementation
// https://github.com/CNCA-CeNAT/kriging/blob/master/src/Variogram.h

#ifndef VARIOGRAM_H
#define VARIOGRAM_H

#define __USE_MISC

//Standard and Math
#include <stdlib.h>
#include "math.h"
#include "constants.h"
#include "variogram.h"

/**
 * @brief Spherical Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double sphericalVariogram(double H, double C0, double CX, double A)
{
	if (H > A)
	{
		return CX;
	}
	else
	{
		return C0 + (CX - C0) * (((3 * H) / (2 * A)) - (H * H * H) / (2 * A * A * A));
	}
}

/**
 * @brief Exponential Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double exponentialVariogram(double H, double C0, double CX, double A)
{
	return C0 + (CX - C0) * (1 - exp(-H / A));
}

/**
 * @brief Gaussian Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double gaussianVariogram(double H, double C0, double CX, double A)
{
	double e = (H / A);
	e = e * e;
	return C0 + (CX - C0) * (1 - exp(-e));
}

/**
 * @brief Wave Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double waveVariogram(double H, double C0, double CX, double A)
{
	double e = H / A;
	return C0 + (CX - C0) * (1 - (sin(e) / e));
}

/**
 * @brief Rational Quadratic Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double rationalQuadraticVariogram(double H, double C0, double CX, double A)
{
	double e = (H * H) / (A * A);
	return C0 + (CX - C0) * (e / (1 + e));
}

/**
 * @brief Circular Semivariogram Model based on http://spatial-analyst.net/ILWIS/htm/ilwisapp/sec/semivar_models_sec.htm
 * @param H distance between points to covariate
 * @param C0 nugget
 * @param CX sill: C0+C
 * @param A range
 * @return The output of the semivariogram determines how much impact will have the sample in the stimation of the predictants value
 */
double circularVariogram(double H, double C0, double CX, double A)
{
	if (H > A)
	{
		return CX;
	}
	else
	{
		double e = H / A;
		double p = 2 / M_PI;
		double r = sqrt(1 - e * e);
		return C0 + (CX - C0) * (1 - p * acos(e) + p * e * r);
	}
}
#endif