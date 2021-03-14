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
Gradient_constant *calConstantParameter(int DIM, double ParticleDistance,int NumberOfParticles, double *Position, int *ParticleType, int ParticleNumber_lambda);
linklist_constant *calcLinklistParameter(double ParticleDistance, double Re_forLaplacian, range_sim *range_sim);
void makeBucket(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist);
void calGravity_x( double NumberOfParticles, int *ParticleType, double *Acceleration ,double Gravity);
void calViscosity( int NumberOfParticles,
                    int *ParticleType,double *Position,double *Velocity,double *Acceleration,int *BucketFirst, int *Nextof,
                    simulation_parameters *parameters, Gradient_constant *constant, linklist_constant *linklist);
void moveParticle( int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration );
void collision(int NumberOfParticles, int *ParticleType, double *Position, double *Velocity,int *BucketFirst, int *Nextof,
                 simulation_parameters *parameters, linklist_constant *linklist);
void setMinimumPressure(int NumberOfParticles, double Re_forGradient,int *ParticleType, double *MinimumPressure,
                         double *Pressure, double *Position, int *BucketFirst, int *Nextof, linklist_constant *linklist);
void calPressureGradient(int NumberOfParticles,int *ParticleType, double *Position,double *Acceleration,double *Pressure, double *MinimumPressure,
                         simulation_parameters *parameters, Gradient_constant *constant, int *BucketFirst, int *Nextof, linklist_constant *linklist);
void moveParticleUsingPressureGradient_andPBC(int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist);
double calDT(int NumberOfParticles, double *Velocity, int *ParticleType, simulation_parameters *parameters);
void makeBucket_withIndex(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist,int *BucketIndex);

// for couette
void moveParticle_couette( int NumberOfParticles, simulation_parameters *parameters, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist);
void moveParticle_couette_reverse( int NumberOfParticles, simulation_parameters *parameters, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist);
void moveParticle2_PBC_couette(int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist, simulation_parameters *parameters);
void setWallSpeed(int np, double wall_speed, double ParticleDistance,double *Position, int *ParticleType, double *Velocity);
void setWallSpeed_reverse(int np, double wall_speed, double ParticleDistance,double *Position, int *ParticleType, double *Velocity);
void defWallSpeed(simulation_parameters *parameters, int ny, double ParticleDistance);
void setCouetteSpeed(int NumberOfParticle, int ny, double *Position, int *ParticleType, double *Velocity, simulation_parameters *parameters);

// for calculating stress
void output_WallStress2D(int NumberOfParticles, double Time, int *ParticleType, double *Position, double *Acceleration, simulation_parameters *parameters);
void make_output_stress(simulation_parameters *parameters);
void calViscosity_stress( int NumberOfParticles,
                    int *ParticleType,double *Position,double *Velocity,double *Acceleration,int *BucketFirst, int *Nextof,
                    simulation_parameters *parameters, Gradient_constant *constant, linklist_constant *linklist);
void calPressureGradient_stress(int NumberOfParticles,int *ParticleType, double *Position,double *Acceleration,double *Pressure, double *MinimumPressure,
                         simulation_parameters *parameters, Gradient_constant *constant, int *BucketFirst, int *Nextof, linklist_constant *linklist);

//------------------------------------------------
#endif