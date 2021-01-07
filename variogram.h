#ifndef VARIOGRAM_H
#define VARIOGRAM_H

struct VariogramModel
{
  double C0; //NUGGET
  double CX; //SILL (C0+C)
  double A;  //RANGE (Max distance to consider v(a)=SILL)
  int VAR;   //Type of variogram to use
};

#define VARIOGRAM_MODEL_LINEAR 1
#define VARIOGRAM_MODEL_SPHERICAL 2
#define VARIOGRAM_MODEL_EXPONENTIAL 3
#define VARIOGRAM_MODEL_GAUSSIAN 4

double sphericalVariogram(double H, double C0, double CX, double A);
double exponentialVariogram(double H, double C0, double CX, double A);
double gaussianVariogram(double H, double C0, double CX, double A);
double waveVariogram(double H, double C0, double CX, double A);
double rationalQuadraticVariogram(double H, double C0, double CX, double A);
double circularVariogram(double H, double C0, double CX, double A);

#endif
