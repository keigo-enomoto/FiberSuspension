#ifndef CALCULATE_PRESSURE_H
#define CALCULATE_PRESSURE_H
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
void calPressure(petsc_object *petsc,int NumberOfParticles, double *Pressure,
                 int *ParticleType, double *NumberDensity, double *Position, int *BoundaryCondition, int *FlagForCheckingBoundaryCondition,
                 simulation_parameters *parameters, Gradient_constant *constant,int *BucketFirst, int *Nextof, linklist_constant *linklist);

//------------------------------------------------
#endif