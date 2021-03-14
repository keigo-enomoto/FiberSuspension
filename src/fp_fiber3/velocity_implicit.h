#ifndef VELOCITY_IMPLICIT_H
#define VELOCITY_IMPLICIT_H
//------------------------------------------------
#include <petscksp.h>

//------------------------------------------------
//  マクロ定義(Macro definition)
//------------------------------------------------

//------------------------------------------------
//  型定義(Type definition)
//------------------------------------------------

//------------------------------------------------
//  プロトタイプ宣言(Prototype declaration)
//------------------------------------------------
void calGravity_and_Viscosity_i(petsc_object *petsc,int NumberOfParticles, double *Velocity,int *ParticleType, double *Position,
                 simulation_parameters *parameters, Gradient_constant *constant,int *BucketFirst, int *Nextof, linklist_constant *linklist);

//------------------------------------------------
#endif