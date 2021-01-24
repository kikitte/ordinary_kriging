#ifndef VARIOGRAM_H
#define VARIOGRAM_H

typedef double (*VariogramFunction)(double H, double C0, double CX, double A);

struct VariogramModel
{
  double NUGGET; //NUGGET
  double SILL;   //SILL (C0+C)
  double RANGE;  //RANGE (Max distance to consider v(a)=SILL)
  char *TYPE;    //Type of variogram to use
};

#define VARIOGRAM_MODEL_SPHERICAL "Spherical"
#define VARIOGRAM_MODEL_EXPONENTIAL "Exponential"
#define VARIOGRAM_MODEL_GAUSSIAN "Gaussian"
#define VARIOGRAM_MODEL_WAVE "Wave"
#define VARIOGRAM_MODEL_RATIONAL_QUADRATIC "RationalQuadratic"
#define VARIOGRAM_MODEL_CIRCLULAR "Circular"

double sphericalVariogram(double H, double C0, double CX, double A);
double exponentialVariogram(double H, double C0, double CX, double A);
double gaussianVariogram(double H, double C0, double CX, double A);
double waveVariogram(double H, double C0, double CX, double A);
double rationalQuadraticVariogram(double H, double C0, double CX, double A);
double circularVariogram(double H, double C0, double CX, double A);

VariogramFunction getFunctionForVariogramType(const char *TYPE);

#endif
