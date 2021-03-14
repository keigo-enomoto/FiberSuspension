#ifndef CALCULATE_SOLID_H
#define CALCULATE_SOLID_H
//------------------------------------------------
// #include <petscksp.h>
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
void calFiberRestoreForce(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, double *Acceleration, simulation_parameters *parameters, linklist_constant *linklist
, double *fiber_l_ave, double *fiber_theta_ave);
void output_tanphi_thin(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist);
// void output_FiberLength(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,linklist_constant *linklist, double Time);
void output_monitor_fiber(int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,linklist_constant *linklist, double Time,double *fiber_l_ave, double *fiber_theta_ave);
void output_fiber_coordinate(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters, double Time);
void make_output_file_fiber(simulation_parameters *parameters, int n_fiber);
void output_unit_thin(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist);


void calFiberRestoreForce_thick(int l_f, int a_f, int n_fiber, int *fiber_consist, double *Position, double *Acceleration, simulation_parameters *parameters, linklist_constant *linklist
, double *fiber_l_ave, double *fiber_theta_ave);
void output_tanphi_thick(int l_f, int a_f, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist);

//------------------------------------------------
#endif