/*************************************************************************************
------------MPS (Moving Particle Simulation) Program----------------
-----------------------------------  Keigo Enomoto (2020)-----------

using in cluster
implicit calculation in viscous term

[difference from main.c]
- calDT();
- calGravity_and_Viscosity_i();

2D or 3D simulation and apply PBC in three axis.(if z-axis cell is below 3, this program will be killed by calcLinklistParameter() )

[ The differences from fp_fiber ] labnote vol.2 p175,176
coupling scheme is different from fp_fiber

In this simulation, the fluid particle's stress(Du/Dt) is caluculated with all neighboring particle
But fiber particle is exerted fiber-fiber spring force and the revers force from fluid

Du/Dt_f = F^s + F^b + F^h (s : stretching, b : Bending, h : hydrodynamics)
F^h = - [ -1/^rho \nabla P_m + \nu \nabla^2 u_m ]_f (m:fluid, f:fiber)
    = - (the force exerting fluid particle by solid particle )



we can simulate poiseuille and couette
[rule in input file]
- If you compute poiseuille flow, let "Gravity" > 0.0 and "wall_speed" = 0.0
- If you compute couette flow, let "Gravity" = 0.0 and "wall_speed" > 0.0

|            | poiseuille | couette |
|:----------:|:----------:|:-------:|
|  Gravity   |    > 0.0   |  = 0.0  |
| wall_speed |    = 0.0   |  > 0.0  |

In this file, the fiber must be thin( n_solid = n_fiber * aspect_ratio )

output :
- input_parameters
- output_*.dat (couette of poiseuille)
- particle*.vtu
- output_analysis.txt
- output_tanphi.dat
- output_FiberMonitor.dat
- output_center.dat

how to execute

```command line
  path/to/mps input_cluster.txt
```
input_cluster.txt : input file
**************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <petscksp.h>

#include "TypeDef.h"
#include "CalculateFunction.h"
#include "parameter.h"
#include "CalculatePressure.h"
#include "velocity_implicit.h"
#include "sub_function.h"
#include "CalculateSolid.h"


/*-----------------------------------------------------------------------*/
// -- for PETSc intitialize --
#define INPUT_FILE "input_cluster.txt"

// for PETSc intitialize
static char help[] = "Use PETSc function when solving poisson equation.\n\n";

/*declare function------------------------------------------------------------------------------------------------------------------*/
double calDT_courant(int NumberOfParticles, double *Velocity, int *ParticleType,simulation_parameters *parameters);

/*  beginning of main ********************************************************************************************************************/

int main( int argc, char** argv )
{
// definition array ----------------------------------------------------------------------------------------
int    *ParticleType;      // DUMMY_WALL(3) or WALL(2) or FLUID(0) or GHOST(-1)
double *Position, *Velocity, *Acceleration;
double *Pressure;
double *NumberDensity;     // n
int    *BoundaryCondition; // the flag of whether we apply boundary condition for the particle
double *MinimumPressure;   // the minimum pressure of each particle in re(effect radius)
int    *FlagForCheckingBoundaryCondition; // Flag of whether there is a particle set that satisfies the boundary condition
int    *BucketFirst;       // the label of first particle in the bucket
int    *BucketLast;        // the lable of last particle in the bucket
int    *Nextof;            // the lable of next particle in the bucket
int    *BucketIndex;

// for PETSc program
petsc_object   petsc;
PC             pc;  // preconditoner
PetscErrorCode ierr;

//set parameters ---------------------------------------------------------------------------------------------
const char *input;
char input_temp[256];
simulation_parameters *parameters; // hold input paramters
Gradient_constant *constant;       // hold constant value(Re,lambda etc)
linklist_constant *linklist;       // hold parameter for cell list
range_sim range_sim={0.0};         // temporaly hold the information about simulatinon area
clock_t start,end;

FILE *fp_input;
int ParticleNumber_lambda = 0;     // for calculate lambda
int NumberOfParticles = 0 ;
int n_fluid = 0;                   // number of fluid particles and solid particles
int n_solid = 0;

int pt = 0;                     // for substitute ParticleType from input dat file
double t_x=0.0,t_y=0.0,t_z=0.0; // for substitute Position     from input dat file

int iTimeStep = 0;       //the times of roop
int TimeStepfordata = 1; //for writeData_inVtuFormat()
double Time = 0.0 ;      //current time
// double two_norm = 0.0; in calPressure
int ny = 0; // for calc wall speed from strain rate

int n_fiber = 0;
int aspect_ratio = 0;
int *fiber_consist;
int tf = 0; // for substitute fiber_consist
double *fiber_l_ave;     // for output average distance between fiber particle
double *fiber_theta_ave; // for output average angle of fiber particle

//read parameters (referring to "parameter.h")----------------------------------------------------------------

strcpy(input_temp,argv[1]);
// strcpy(input_temp,INPUT_FILE);
input = input_temp;

parameters = read_parameters_from_file(input); //parameters is pointa variable

//Get Number of particle -----------------------------------------------------------------------------------

fp_input = fopen(parameters->Init_place, "r");
fscanf(fp_input,"%d\n",&NumberOfParticles);
fscanf(fp_input,"%d\n",&n_fluid);
fscanf(fp_input,"%d\n",&n_fiber);
fscanf(fp_input,"%d\n",&aspect_ratio);

Acceleration = (double *)malloc(3*NumberOfParticles *sizeof(double));
ParticleType = (int    *)malloc(  NumberOfParticles * sizeof(int));
Position     = (double *)malloc(3*NumberOfParticles * sizeof(double));
Velocity     = (double *)malloc(3*NumberOfParticles * sizeof(double));
Pressure     = (double *)malloc(  NumberOfParticles * sizeof(double));
NumberDensity = (double *)malloc(NumberOfParticles * sizeof(double));
BoundaryCondition =(int *)malloc(NumberOfParticles * sizeof(int));
MinimumPressure   = (double *)malloc(NumberOfParticles * sizeof(double));
FlagForCheckingBoundaryCondition = (int *)malloc(NumberOfParticles * sizeof(int));
fiber_consist = (int *)malloc(n_fiber * aspect_ratio * sizeof(int));
fiber_l_ave = (double *)malloc(n_fiber * sizeof(double));
fiber_theta_ave = (double *)malloc(n_fiber * sizeof(double));

// initial arrangement---------------------------------------------------------------------------------------------------------------------
fscanf(fp_input,"MIN_X : %lf, MAX_X : %lf\n",&range_sim.MIN_X,&range_sim.MAX_X);
fscanf(fp_input,"MIN_Y : %lf, MAX_Y : %lf\n",&range_sim.MIN_Y,&range_sim.MAX_Y);
fscanf(fp_input,"MIN_Z : %lf, MAX_Z : %lf\n",&range_sim.MIN_Z,&range_sim.MAX_Z);

for(int i=0;i<NumberOfParticles;i++) {

  fscanf(fp_input,"%d %lf %lf %lf\n",&pt,&t_x,&t_y,&t_z);
  ParticleType[i] = pt;
  Position[i*3] = t_x;
  Position[i*3+1] = t_y;
  Position[i*3+2] = t_z;
}

fscanf(fp_input,"lambda_particle : %d\n",&ParticleNumber_lambda);

for (int i = 0; i < aspect_ratio * n_fiber ; i++){
  fscanf(fp_input,"%d\n",&tf);
  fiber_consist[i] = tf ;
}
fclose(fp_input);

n_solid = n_fiber * aspect_ratio;
// calculate constant ---------------------------------------------------------------------------------------------------------------------

constant = calConstantParameter(parameters->DIM, parameters->ParticleDistance,NumberOfParticles,Position,ParticleType,ParticleNumber_lambda);

// prepare the linklist ------------------------------------------------------------------------------------------------------------------

linklist = calcLinklistParameter(parameters->ParticleDistance,constant->Re_forLaplacian,&range_sim);
BucketFirst = (int *)malloc(sizeof(int) * (linklist->nBxyz));	//バケットに格納された先頭の粒子の番号
BucketLast  = (int *)malloc(sizeof(int) * (linklist->nBxyz));	//バケットに格納された最後尾の粒子の番号
Nextof      = (int *)malloc(sizeof(int) * NumberOfParticles);	//同じバケット内の次の粒子の番号
BucketIndex = (int *)malloc(sizeof(int) * NumberOfParticles);

// define wall speed
ny = (int)( (linklist->Ly)/(parameters->ParticleDistance) ); //the number of particles in y axis
defWallSpeed(parameters,ny,parameters->ParticleDistance) ;   // finish setting *paramters

// output
output_parameter(parameters,constant,linklist,NumberOfParticles);
// prepare output for fibers
if( aspect_ratio>0 ) make_output_file_fiber(parameters,n_fiber);

// set PETSc objects ---------------------------------------------------------------------------------------------------------------------

ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
// create vector
ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,NumberOfParticles,NULL,&petsc.b);CHKERRQ(ierr);
ierr = VecDuplicate(petsc.b,&petsc.x);CHKERRQ(ierr);
// set solver type
ierr = KSPCreate(PETSC_COMM_SELF,&petsc.ksp);CHKERRQ(ierr);
ierr = KSPGetPC(petsc.ksp,&pc);CHKERRQ(ierr);
ierr = PCSetType(pc,PCICC);CHKERRQ(ierr);
ierr = KSPSetType(petsc.ksp,KSPCG);CHKERRQ(ierr);

/*========================================================================================================================================
========================= Main roop =====================================================================================================*/

//initialize array
for(int i=0; i<NumberOfParticles; i++){
  Velocity[3*i] = 0.0; Velocity[3*i+1] = 0.0; Velocity[3*i+2] = 0.0;
  Acceleration[3*i] = 0.0; Acceleration[3*i+1] = 0.0; Acceleration[3*i+2] = 0.0;
  Pressure[i] = 0.0;
  BoundaryCondition[i] = 0; // INNER_PARTICLE
  NumberDensity[i] = 0.0;
  MinimumPressure[i] = 0.0;
  FlagForCheckingBoundaryCondition[i] = 0; // substitute invaild integral on purpose
  BucketIndex[i] = -1;
}

// writeData_inVtuFormat_tmp(parameters->DIM,ParticleType,Position,NumberOfParticles);

// If both wall_speed and Gravity are zero, we have to kill this program above
check_flow_type(parameters);
// fix wall speed
if (parameters->wall_speed != 0.0){
  setWallSpeed(NumberOfParticles,parameters->wall_speed,parameters->ParticleDistance,Position,ParticleType,Velocity);
}
// set initial speed
//setCouetteSpeed(NumberOfParticles,ny,Position,ParticleType,Velocity,parameters);
// output vtu
if (parameters->VTU_flag == 1 || n_fiber > 1){
  makeBucket_withIndex(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist,BucketIndex);
  writeData_inVtuFormat(parameters->DIM,ParticleType,Position,Velocity,Pressure,NumberDensity,BoundaryCondition,NumberOfParticles,0,BucketIndex,parameters->output_dir);
}else if(parameters-> VTU_flag == 2){
  makeBucket_withIndex(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist,BucketIndex);
  writeData_inVtuFormat_fiber_wall(parameters->DIM,ParticleType,Position,Velocity,Pressure,NumberDensity,BoundaryCondition,NumberOfParticles,0,BucketIndex,parameters->output_dir);
}

start = clock();
while(1){
  makeBucket(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist);
  calGravity_and_Viscosity_i(&petsc,NumberOfParticles,Velocity,ParticleType,Position,parameters,constant,BucketFirst,Nextof,linklist);
  //fix wall speed
  if(parameters->wall_speed != 0.0){
    setWallSpeed(NumberOfParticles,parameters->wall_speed,parameters->ParticleDistance,Position,ParticleType,Velocity);
  }
  if(aspect_ratio>0) calFiberRestoreForce(aspect_ratio,n_fiber,fiber_consist,Position,Acceleration,parameters,linklist,fiber_l_ave,fiber_theta_ave);
  //move particle
  // if(parameters->Gravity != 0.0){
  //   moveParticle(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration);
  // }else{
    moveParticle_couette(NumberOfParticles,parameters,ParticleType,Position,Velocity,Acceleration,linklist);
  // }
  collision(NumberOfParticles,ParticleType,Position,Velocity,BucketFirst,Nextof,parameters,linklist);
  calPressure(&petsc,NumberOfParticles,Pressure,ParticleType,NumberDensity,Position,BoundaryCondition,FlagForCheckingBoundaryCondition,parameters,constant,BucketFirst,Nextof,linklist);
  setMinimumPressure(NumberOfParticles,constant->Re_forGradient,ParticleType,MinimumPressure,Pressure,Position,BucketFirst,Nextof,linklist);
  calPressureGradient(NumberOfParticles,ParticleType,Position,Acceleration,Pressure,MinimumPressure,parameters,constant,BucketFirst,Nextof,linklist);
  moveParticleUsingPressureGradient_andPBC(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration,linklist);
  iTimeStep++;
  Time += parameters->DT;
  parameters->DT = calDT_courant(NumberOfParticles,Velocity,ParticleType,parameters);

  // output information about fiber
  if(aspect_ratio>0){
    if ( iTimeStep % ((int)(parameters->OUTPUT_INTERVAL/10)) == 0 ){
      output_tanphi_thin(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time,linklist);
    }else if (parameters->OUTPUT_INTERVAL < 11){ // output every time
      output_tanphi_thin(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time,linklist);
    }
  }

  // output data
  if( (iTimeStep % parameters->OUTPUT_INTERVAL) == 0 ){
    end = clock();
    output_analysis(NumberOfParticles,parameters->ParticleDistance,iTimeStep,parameters,n_fluid,n_solid,Time,ny,start);
    if (parameters->VTU_flag == 1){
      writeData_inVtuFormat(parameters->DIM,ParticleType,Position,Velocity,Pressure,NumberDensity,BoundaryCondition,NumberOfParticles,TimeStepfordata,BucketIndex,parameters->output_dir);
    }else if(parameters->VTU_flag == 2){
      writeData_inVtuFormat_fiber_wall(parameters->DIM,ParticleType,Position,Velocity,Pressure,NumberDensity,BoundaryCondition,NumberOfParticles,TimeStepfordata,BucketIndex,parameters->output_dir);
    }
    // if(parameters->Gravity!=0.0){
    //   output_poiseuille(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
    // }else{
      output_couette(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
    // }
    if( aspect_ratio >0 ) {
      output_monitor_fiber(n_fiber,fiber_consist,Position,parameters,linklist,Time,fiber_l_ave,fiber_theta_ave);
      output_fiber_coordinate(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time);
    }
    TimeStepfordata ++ ;
  }
  if( Time >= parameters->FINISH_TIME ){break;}
}
/*==================== End of main roop ================================================================================================
========================================================================================================================================*/
// if(parameters->Gravity!=0.0){
//   output_poiseuille(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
// }else{
  output_couette(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
// }
output_analysis(NumberOfParticles,parameters->ParticleDistance,iTimeStep,parameters,n_fluid,n_solid,Time,ny,start);


free(parameters); free(parameters->output_dir); free(parameters->Init_place);
free(constant); free(linklist);
free(ParticleType);
free(Position); free(Velocity); free(Acceleration);
free(Pressure); free(NumberDensity);
free(BoundaryCondition); free(MinimumPressure);
free(FlagForCheckingBoundaryCondition);
free(BucketFirst); free(BucketLast); free(Nextof);free(BucketIndex);
free(fiber_consist);
free(fiber_l_ave); free(fiber_theta_ave);

ierr = KSPDestroy(&petsc.ksp);CHKERRQ(ierr);
ierr = VecDestroy(&petsc.x);CHKERRQ(ierr);
ierr = VecDestroy(&petsc.b);CHKERRQ(ierr);
ierr = PetscFinalize();

printf("*** End ***\n");
return 0;
}
/*  end of main  ***********************************************************************************************************************/

double calDT_courant(int NumberOfParticles, double *Velocity, int *ParticleType,simulation_parameters *parameters){
  const double Cmax = 0.2; //max value of Courant number
  double dt;
  double MaxVelocity = parameters->wall_speed;
  double tempVelocity= 0.0;

  for(int i=0; i<NumberOfParticles; i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      tempVelocity = Velocity[3*i]*Velocity[3*i] + Velocity[3*i+1]*Velocity[3*i+1] + Velocity[3*i+2]*Velocity[3*i+2];
      if(tempVelocity > MaxVelocity){
        MaxVelocity = tempVelocity;
      }
    }
  }
  MaxVelocity = sqrt(MaxVelocity);
  // if (MaxVelocity <= (parameters->FINISH_TIME)/1000000 ){ return parameters->DT ;}
  // let Cmax = 0.2
  dt = Cmax * (parameters->ParticleDistance) / MaxVelocity;

  // if dt is too large, the simulation will be finished at one step
  // To avoid this situation, dt must be limited.
  if ( parameters->FINISH_TIME/10.0 < dt ){ return (parameters->FINISH_TIME)/(100); }
  else{ return dt; }
}
