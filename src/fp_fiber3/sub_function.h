#ifndef SUB_FUNCTION_H
#define SUB_FUNCTION_H
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

simulation_parameters *read_parameters_from_file(const char *input);
void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, double *NumberDensity, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir);
void output_parameter(simulation_parameters *parameters, Gradient_constant *constant,linklist_constant *linklist, int NumberOfParticles);
void output_couette(int n, double *Position, double *Velocity, int *ParticleType, simulation_parameters *parameters, double Time, int iTimeStep, clock_t start);
void output_poiseuille(int n, double *Position, double *Velocity, int *ParticleType, simulation_parameters *parameters, double Time, int iTimeStep, clock_t start);
void writeData_inVtuFormat_tmp(int DIM,int *ParticleType, double *Position,int NumberOfParticles);
void output_analysis(int n, double ParticleDistance,int iTimeStep, simulation_parameters *parameters, int n_fluid, int n_solid, double Time, int ny, clock_t start);
void check_flow_type(simulation_parameters *parameters);
void writeData_inVtuFormat_fiber_wall(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, double *NumberDensity, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir);
void output_coordinate(int n_fiber, int aspect_ratio,int NumberOfParticles, int n_fluid, linklist_constant *linklist,
double ParticleDistance, simulation_parameters *parameters, int *ParticleType, double *Position);
void output_continue(int n_fiber, int aspect_ratio,int NumberOfParticles, int n_fluid, linklist_constant *linklist,
double ParticleDistance, simulation_parameters *parameters, int *ParticleType, double *Position, double *Velocity, double Time);
#ifdef _OPENMP
void print_omp();
#endif

//------------------------------------------------
#endif