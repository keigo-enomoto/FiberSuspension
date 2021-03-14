/*************************************************************************************
------------MPS (Moving Particle Simulation) Program----------------
-----------------------------------  Keigo Enomoto (2021)-----------


calculating strain_rate = 0.0
and if pressure is 0, this program is killed
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
#include "sub_function.h"
#include "CalculateSolid.h"

/*-----------------------------------------------------------------------*/
// -- for PETSc intitialize --
#define INPUT_FILE "input_cluster.txt"

// for PETSc intitialize
static char help[] = "Use PETSc function when solving poisson equation.\n\n";

/*  beginning of main ********************************************************************************************************************/

int main( int argc, char** argv )
{
// definition array ----------------------------------------------------------------------------------------
int    *ParticleType;      // DUMMY_WALL(3) or WALL(2) or FLUID(0) or GHOST(-1)
double *Position, *Velocity, *Acceleration; // for wall particle, Acceleration is \rho * Du/Dt
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
// int l_f = 0;     // long length
// int a_f = 0;     // short length
int aspect_ratio = 0 ;
int *fiber_consist;
int tf = 0; // for substitute fiber_consist
double *fiber_l_ave;     // for output average distance between fiber particle
double *fiber_theta_ave; // for output average angle of fiber particle

double sum_vel; // for finish determination
int end_flg = 0;

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
// fscanf(fp_input,"%d\n",&l_f);
// fscanf(fp_input,"%d\n",&a_f);
// aspect_ratio = l_f / a_f ;

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
// defWallSpeed(parameters,ny,parameters->ParticleDistance) ;   // finish setting *paramters
parameters->wall_speed = 0.0;

// output
output_parameter(parameters,constant,linklist,NumberOfParticles);
// prepare output for fibers and stress
if( aspect_ratio>0 ) make_output_file_fiber(parameters,n_fiber);
make_output_stress(parameters);

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
========================= Main roop ==================================================================================================*/

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
// check_flow_type(parameters);

// fix wall speed
setWallSpeed(NumberOfParticles,parameters->wall_speed,parameters->ParticleDistance,Position,ParticleType,Velocity);

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
  calViscosity(NumberOfParticles,ParticleType,Position,Velocity,Acceleration,BucketFirst,Nextof,parameters,constant,linklist);
  // calViscosity_stress(NumberOfParticles,ParticleType,Position,Velocity,Acceleration,BucketFirst,Nextof,parameters,constant,linklist);

  if(aspect_ratio>0) calFiberRestoreForce(aspect_ratio,n_fiber,fiber_consist,Position,Acceleration,parameters,linklist,fiber_l_ave,fiber_theta_ave);
  //move particle
  // if(parameters->Gravity != 0.0){
    moveParticle(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration);
  // }else{
    // moveParticle_couette(NumberOfParticles,parameters,ParticleType,Position,Velocity,Acceleration,linklist);
  // }
  collision(NumberOfParticles,ParticleType,Position,Velocity,BucketFirst,Nextof,parameters,linklist);
  calPressure(&petsc,NumberOfParticles,Pressure,ParticleType,NumberDensity,Position,BoundaryCondition,FlagForCheckingBoundaryCondition,parameters,constant,BucketFirst,Nextof,linklist);
  setMinimumPressure(NumberOfParticles,constant->Re_forGradient,ParticleType,MinimumPressure,Pressure,Position,BucketFirst,Nextof,linklist);
  calPressureGradient(NumberOfParticles,ParticleType,Position,Acceleration,Pressure,MinimumPressure,parameters,constant,BucketFirst,Nextof,linklist);
  // calPressureGradient_stress(NumberOfParticles,ParticleType,Position,Acceleration,Pressure,MinimumPressure,parameters,constant,BucketFirst,Nextof,linklist);
  moveParticleUsingPressureGradient_andPBC(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration,linklist);
  iTimeStep++;
  Time += parameters->DT;
  parameters->DT = calDT(NumberOfParticles,Velocity,ParticleType,parameters);

  // output information about fiber
  // if(aspect_ratio>0){
  //   if (parameters->OUTPUT_INTERVAL < 11){ // output every time
  //     output_tanphi_thin(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time,linklist);
  //   }else if ( iTimeStep % ((int)(parameters->OUTPUT_INTERVAL/10)) == 0 ){
  //     output_tanphi_thin(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time,linklist);
  //   }
  // }

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
      // output_couette(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
    // }
    if( aspect_ratio >0 ) {
      output_monitor_fiber(n_fiber,fiber_consist,Position,parameters,linklist,Time,fiber_l_ave,fiber_theta_ave);
      // output_fiber_coordinate(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time);
    }
    output_coordinate(n_fiber,aspect_ratio,NumberOfParticles,n_fluid,linklist,parameters->ParticleDistance,parameters,ParticleType,Position);
    TimeStepfordata ++ ;
  }

  sum_vel = 0.0;
  for(int i=0;i<NumberOfParticles;i++){
    sum_vel += Velocity[3*i] * Velocity[3*i] + Velocity[3*i+1] * Velocity[3*i+1];
  }
  if( sum_vel < n_fluid ){ //average velocity is below 1.0
    end_flg += 1;
  }else{
    end_flg = 0;
  }
  if( Time >= parameters->FINISH_TIME || end_flg == 5 ){break;}
}
/*==================== End of main roop ================================================================================================
========================================================================================================================================*/
// if(parameters->Gravity!=0.0){
//   output_poiseuille(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
// }else{
  // output_couette(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
// }
output_analysis(NumberOfParticles,parameters->ParticleDistance,iTimeStep,parameters,n_fluid,n_solid,Time,ny,start);
writeData_inVtuFormat(parameters->DIM,ParticleType,Position,Velocity,Pressure,NumberDensity,BoundaryCondition,NumberOfParticles,TimeStepfordata,BucketIndex,parameters->output_dir);
output_fiber_coordinate(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time);
output_coordinate(n_fiber,aspect_ratio,NumberOfParticles,n_fluid,linklist,parameters->ParticleDistance,parameters,ParticleType,Position);
output_unit_thin(aspect_ratio,n_fiber,fiber_consist,Position,parameters,Time,linklist);


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
/*  end of main  ******************************************************************************************************************/
