#ifndef CALCULATE_FUNCTION_H
#define CALCULATE_FUNCTION_H
//------------------------------------------------
#include "TypeDef.h"
//------------------------------------------------
//  マクロ定義(Macro definition)
//------------------------------------------------

//------------------------------------------------
//  型定義(Type definition)
//------------------------------------------------

//------------------------------------------------
//  プロトタイプ宣言(Prototype declaration)
//------------------------------------------------

int check_parallel(double x0, double y0, double c, double d, int aspect_ratio, double ParticleDistance);
int check_range(double x_ol,int aspect_ratio, double x_s, double x_e, double x1, double x2);
int check_cross(double Lx, int aspect_ratio, double xa, double ya,double xb, double yb,double xc, double yc,double xd, double yd);
int GetRandom(int min,int max);
linklist_constant *calcLinklistParameter(double ParticleDistance, double Re_forNumberDensity, double MIN_X, double MAX_X, double MIN_Y, double MAX_Y);
void makeBucket(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist);
void calNumberDensity2D(int NumberOfParticles, double Re_forNumberDensity, int *ParticleType, double *NumberDensity, double *Position,
                              int *BucketFirst, int *Nextof, linklist_constant *linklist);
double calPotential(int NumberOfParticles, double Re_forNumberDensity, int *ParticleType, double *NumberDensity, double *Position,
                              int *BucketFirst, int *Nextof, linklist_constant *linklist);
double calPotential2(int NumberOfParticles, double *NumberDensity, double n0);


//------------------------------------------------
#endif