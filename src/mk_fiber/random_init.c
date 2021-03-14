/*************************************************************************************
------------MPS (Moving Particle Simulation) Program----------------
-----------------------------------  Keigo Enomoto (2020)-----------

give initial velocity

output : input
**************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <petscksp.h>

#include "TypeDef.h"
#include "CalculateFunction.h"
// #include "parameter.h"
// #include "CalculatePressure.h"

#define M_PI 3.141593
/*for parameter.c--------------------------------------------------------*/

#define BUFFER_SIZE 0x100
#define INPUT_DIRECTORY "./../mk_particle/poiseuille/"
/*-----------------------------------------------------------------------*/

enum{
    ON = 1,
    OFF = 0
};

// for PETSc intitialize
static char help[] = "Use PETSc function when solving poisson equation.\n\n";
/*declare function------------------------------------------------------------------------------------------------------------------*/
simulation_parameters *read_parameters_from_file(const char *input);
void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir);
void output_parameter(simulation_parameters *parameters, Gradient_constant *constant,linklist_constant *linklist, int NumberOfParticles);
static void init_velosity(int pn,double *Velocity,int *ParticleType);
void output_place(int NumberOfParticles, int ParticleNumber_lambda, double ParticleDistance,double *Position,int *ParticleType, linklist_constant *linklist);
static void output_velocity(int pn, double *vl, int *ParticleType, char *output_dir);
/*  beginning of main ********************************************************************************************************************/

int main( int argc, char** argv )
{
// definition array ----------------------------------------------------------------------------------------
int    *ParticleType;      //DUMMY_WALL(3) or WALL(2) or FLUID(0) or GHOST(-1)
double *Position, *Velocity, *Acceleration;
double *Pressure;
double *NumberDensity;     //n
int    *BoundaryCondition; //the flag of whether we apply boundary condition for the particle
double *MinimumPressure;   //the minimum pressure of each particle in re(effect radius)
int    *FlagForCheckingBoundaryCondition; //Flag of whether there is a particle set that satisfies the boundary condition
int    *BucketFirst;       //the label of first particle in the bucket
int    *BucketLast;        //the lable of last particle in the bucket
int    *Nextof;            //the lable of next particle in the bucket
int    *BucketIndex;

double sum_pre;
double sum_vl;
int end_flg = 0;

// for PETSc program
petsc_object   petsc;
PC             pc;  //preconditoner
PetscErrorCode ierr;

//set parameters --------------------------------------------------------------------------------------------
const char *input;
char input_temp[256];
simulation_parameters *parameters; //hold input paramters
Gradient_constant *constant;
linklist_constant *linklist;
range_sim range_sim={0.0};
clock_t start,end;

FILE *fp_input;
int ParticleNumber_lambda = 0;
int NumberOfParticles = 0 ;

int iTimeStep = 0;       //the times of roop
int TimeStepfordata = 1; //for writeData_inVtuFormat()
double Time = 0.0 ;      //current time

//read parameters (referring to "parameter.h")----------------------------------------------------------------

//strcpy(input_temp,INPUT_DIRECTORY);
//strcat(input_temp,argv[1]);
strcpy(input_temp,argv[1]);
input = input_temp;
//input = "input_random.txt";

printf("input file %s\n",input);
parameters = read_parameters_from_file(input); //parameters is pointa variable

//Get Number of particle -----------------------------------------------------------------------------------

fp_input = fopen(parameters->Init_place, "r");
fscanf(fp_input,"%d\n",&NumberOfParticles);
printf("NumberOfparticle : %d \n",NumberOfParticles);

Acceleration = (double *)malloc(3*NumberOfParticles *sizeof(double));
ParticleType = (int    *)malloc(  NumberOfParticles * sizeof(int));
Position     = (double *)malloc(3*NumberOfParticles * sizeof(double));
Velocity     = (double *)malloc(3*NumberOfParticles * sizeof(double));
Pressure     = (double *)malloc(  NumberOfParticles * sizeof(double));
NumberDensity = (double *)malloc(NumberOfParticles * sizeof(double));
BoundaryCondition =(int *)malloc(NumberOfParticles * sizeof(int));
MinimumPressure   = (double *)malloc(NumberOfParticles * sizeof(double));
FlagForCheckingBoundaryCondition = (int *)malloc(NumberOfParticles * sizeof(int));

// initial arrangement---------------------------------------------------------------------------------------------------------------------
fscanf(fp_input,"MIN_X : %lf, MAX_X : %lf\n",&range_sim.MIN_X,&range_sim.MAX_X);
fscanf(fp_input,"MIN_Y : %lf, MAX_Y : %lf\n",&range_sim.MIN_Y,&range_sim.MAX_Y);
fscanf(fp_input,"MIN_Z : %lf, MAX_Z : %lf\n",&range_sim.MIN_Z,&range_sim.MAX_Z);

for(int i=0;i<NumberOfParticles;i++) {
		fscanf(fp_input,"%d %lf %lf %lf\n",&ParticleType[i],&Position[i*3],&Position[i*3+1],&Position[i*3+2]);
}
fscanf(fp_input,"lambda_particle : %d",&ParticleNumber_lambda);
fclose(fp_input);

// calculate constant ---------------------------------------------------------------------------------------------------------------------

constant = calConstantParameter(parameters->DIM, parameters->ParticleDistance,NumberOfParticles,Position,ParticleType,ParticleNumber_lambda);

// prepare the linklist ------------------------------------------------------------------------------------------------------------------
linklist = calcLinklistParameter(parameters->ParticleDistance,constant->Re_forLaplacian,&range_sim);
BucketFirst = (int *)malloc(sizeof(int) * (linklist->nBxyz));	//バケットに格納された先頭の粒子の番号
BucketLast  = (int *)malloc(sizeof(int) * (linklist->nBxyz));	//バケットに格納された最後尾の粒子の番号
Nextof      = (int *)malloc(sizeof(int) * NumberOfParticles);	//同じバケット内の次の粒子の番号
BucketIndex = (int *)malloc(sizeof(int) * NumberOfParticles);

output_parameter(parameters,constant,linklist,NumberOfParticles);

// set PETSc objects ---------------------------------------------------------------------------------------------------------------------

ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
//create vector
ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,NumberOfParticles,NULL,&petsc.b);CHKERRQ(ierr);
ierr = VecDuplicate(petsc.b,&petsc.x);CHKERRQ(ierr);
//set solver type
ierr = KSPCreate(PETSC_COMM_SELF,&petsc.ksp);CHKERRQ(ierr);
ierr = KSPGetPC(petsc.ksp,&pc);CHKERRQ(ierr);
ierr = PCSetType(pc,PCICC);CHKERRQ(ierr);
ierr = KSPSetType(petsc.ksp,KSPCG);CHKERRQ(ierr);

/*========================================================================================================================================
========================= Main roop ======================================================================================================*/
printf("*** Start ***\n");

//initialize array for writeData_inVtuFormat
for(int i=0; i<NumberOfParticles; i++){
  //Velocity[3*i] = 0.0; Velocity[3*i+1] = 0.0; Velocity[3*i+2] = 0.0;
  Pressure[i] = 0.0;
  BoundaryCondition[i] = 0; //INNER_PARTICLE
  BucketIndex[i] = -1;
}
init_velosity(NumberOfParticles,Velocity,ParticleType);
output_velocity(NumberOfParticles,Velocity,ParticleType,parameters->output_dir);
//makeBucket_withIndex(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist,BucketIndex);
writeData_inVtuFormat(parameters->DIM,ParticleType,Position,Velocity,Pressure,BoundaryCondition,NumberOfParticles,0,BucketIndex,parameters->output_dir);

start = clock();
while(1){
  makeBucket(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist);
  //if ( parameters->GRAVITY_SWITCH == ON){
  //  calGravity_x(NumberOfParticles,ParticleType,Acceleration);
  //}
  calViscosity(NumberOfParticles,ParticleType,Position,Velocity,Acceleration,BucketFirst,Nextof,parameters,constant,linklist);
  moveParticle(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration);
  collision(NumberOfParticles,ParticleType,Position,Velocity,BucketFirst,Nextof,parameters,linklist);
  calPressure(&petsc,NumberOfParticles,Pressure,ParticleType,NumberDensity,Position,BoundaryCondition,FlagForCheckingBoundaryCondition,parameters,constant,BucketFirst,Nextof,linklist);
  setMinimumPressure(NumberOfParticles,constant->Re_forGradient,ParticleType,MinimumPressure,Pressure,Position,BucketFirst,Nextof,linklist);
  calPressureGradient(NumberOfParticles,ParticleType,Position,Acceleration,Pressure,MinimumPressure,parameters,constant,BucketFirst,Nextof,linklist);
  moveParticleUsingPressureGradient_andPBC(NumberOfParticles,parameters->DT,ParticleType,Position,Velocity,Acceleration,linklist);
  iTimeStep++;
  Time += parameters->DT;
  parameters->DT = calDT(NumberOfParticles,Velocity,parameters);
  if( (iTimeStep % parameters->OUTPUT_INTERVAL) == 0 ){
    end = clock();
    //printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, NumberOfParticles);
    //printf("  past time: %lf sec   DT : %lf\n",(double)(end - start)/ CLOCKS_PER_SEC,parameters->DT);
    if((iTimeStep % parameters->OUTPUT_INTERVAL*3) == 0){
      writeData_inVtuFormat(parameters->DIM,ParticleType,Position,Velocity,Pressure,BoundaryCondition,NumberOfParticles,TimeStepfordata,BucketIndex,parameters->output_dir);
      TimeStepfordata ++ ;
    }
    //output_poiseuille(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);

    sum_pre = 0.0; sum_vl=0.0;
    for(int i=0;i<NumberOfParticles;i++){
      sum_pre += Pressure[i];
      sum_vl += Velocity[3*i] + Velocity[3*i+1];
    }
    if(sum_vl == 0.0 && sum_pre == 0.0 ){
      end_flg += 1;
    }else{
      end_flg = 0;
    }
  }

 if( Time >= parameters->FINISH_TIME || end_flg > 5 ){break;}
}
/*==================== End of main roop ================================================================================================
========================================================================================================================================*/
//output_poiseuille(NumberOfParticles,Position,Velocity,ParticleType,parameters,Time,iTimeStep,start);
output_place(NumberOfParticles,ParticleNumber_lambda,parameters->ParticleDistance,Position,ParticleType,linklist);
printf("*** End ***");

free(parameters); free(parameters->output_dir); free(parameters->Init_place);
free(constant); free(linklist);
free(ParticleType);
free(Position); free(Velocity); free(Acceleration);
free(Pressure); free(NumberDensity);
free(BoundaryCondition); free(MinimumPressure);
free(FlagForCheckingBoundaryCondition);
free(BucketFirst); free(BucketLast); free(Nextof);free(BucketIndex);

ierr = KSPDestroy(&petsc.ksp);CHKERRQ(ierr);
ierr = VecDestroy(&petsc.x);CHKERRQ(ierr);
ierr = VecDestroy(&petsc.b);CHKERRQ(ierr);
ierr = PetscFinalize();

return 0;
}
/*  end of main  ***********************************************************************************************************************/


/**************************************************************
  *These two function is for reading input paramters          *
***************************************************************/

simulation_parameters *read_parameters_from_file(const char *input)
{
    FILE *fp;
    simulation_parameters *parameters;
    char output_dir[BUFFER_SIZE];
    char Init_place[BUFFER_SIZE];
    parameters = malloc(sizeof(simulation_parameters));

    input_parameter input_parameters[] =
        {{"DIM","%d",&(parameters->DIM)},
         {"ParticleDistance","%lf",&(parameters->ParticleDistance)},
         {"FluidDensity","%lf",&(parameters->FluidDensity)},
         {"KinematicViscosity","%lf",&(parameters->KinematicViscosity)},
         {"GRAVITY_SWITCH","%d",&(parameters->GRAVITY_SWITCH)},
         {"DT","%lf",&(parameters->DT)},
         {"FINISH_TIME","%lf",&(parameters->FINISH_TIME)},
         {"OUTPUT_INTERVAL","%d",&(parameters->OUTPUT_INTERVAL)},
         {"output_dir","%s",&output_dir},
         {"Init_place","%s",&Init_place}
        };

    fp = fopen(input,"rt");

    read_input_parameters(fp,input_parameters);

    parameters->output_dir = malloc((strlen(output_dir) + 1) * sizeof(char));
    strcpy(parameters->output_dir,output_dir);

    parameters->Init_place = malloc((strlen(Init_place) + 1) * sizeof(char));
    strcpy(parameters->Init_place,Init_place);

    fclose(fp);

    return parameters;
}

//-----------------------------------------------------------------------------------------------------
void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024]; //This is array to have the word of filename (fileName[0] = p, fileName[1] = a)
  if(DIM == 2){
    sprintf(fileName, "%s/particle_%04d.vtu",output_dir ,FileNumber); //sprintf is the function to write the word to array
  }else if(DIM == 3){
    sprintf(fileName, "%s/particle_%04d.vtu", output_dir,FileNumber); //sprintf is the function to write the word to array
  }

  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",Position[i*3],Position[i*3+1],Position[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",ParticleType[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( Velocity[i*3]*Velocity[i*3] + Velocity[i*3+1]*Velocity[i*3+1] + Velocity[i*3+2]*Velocity[i*3+2] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='boundarycondition' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",BoundaryCondition[i]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='BucketIndex' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",BucketIndex[i]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  //点(粒子)のデータなら'1'が入る
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


/*---------------------------------------------------------------------------------------------------*/
void output_parameter(simulation_parameters *parameters, Gradient_constant *constant,linklist_constant *linklist, int NumberOfParticles){
  FILE *fp;
  char fileName[1024];

  if(parameters->DIM == 2){
    sprintf(fileName, "%s/input_parameters.txt",parameters->output_dir); //sprintf is the function to write the word to array
  }else if(parameters->DIM == 3){
    sprintf(fileName, "%s/input_parameters.txt",parameters->output_dir); //sprintf is the function to write the word to array
  }
  fp = fopen(fileName,"w");

  fprintf(fp,"*******input parameter *********\n");
  fprintf(fp,"DIM : %d \n",parameters -> DIM);
  fprintf(fp,"ParticleDistance : %lf \n",parameters -> ParticleDistance);
  fprintf(fp,"FluidDensity : %lf \n",parameters -> FluidDensity);
  fprintf(fp,"KinematicViscosity : %lf \n",parameters -> KinematicViscosity);
  fprintf(fp,"GRAVITY_SWITCH : %d (ON -> 1, OFF -> 0)\n",parameters -> GRAVITY_SWITCH);
  fprintf(fp,"DT : %lf \n",parameters -> DT);
  fprintf(fp,"FINISH_TIME : %lf \n",parameters -> FINISH_TIME);
  fprintf(fp,"OUTPUT_INTERVAL : %d \n",parameters -> OUTPUT_INTERVAL);
  fprintf(fp,"output directory : %s \n",parameters -> output_dir);
  fprintf(fp,"Initial placement : %s \n",parameters -> Init_place);
  fprintf(fp,"********************************\n");

  fprintf(fp,"*******constant parameter *********\n");
  fprintf(fp,"Re_forNumberDensity %lf \n",constant -> Re_forNumberDensity);
  fprintf(fp,"Re_forGradient : %lf \n",constant -> Re_forGradient);
  fprintf(fp,"Re_forLaplacian : %lf \n",constant -> Re_forLaplacian);
  fprintf(fp,"N0_forNumebrDensity: %lf \n",constant -> N0_forNumberDensity);
  fprintf(fp,"N0_forGradiant: %lf\n",constant -> N0_forGradient);
  fprintf(fp,"N0_forLaplacian : %lf \n",constant -> N0_forLaplacian);
  fprintf(fp,"Lambda : %lf \n",constant -> Lambda);
  fprintf(fp,"***********************************\n");

  fprintf(fp,"******** linklist parameters ***********\n");
  fprintf(fp,"DB  : %lf\n",linklist->DB) ;
  fprintf(fp,"nBx : %d\n",linklist->nBx) ;
  fprintf(fp,"nBy : %d\n",linklist->nBy) ;
  fprintf(fp,"nBz : %d\n",linklist->nBz) ;
  fprintf(fp,"nBxy : %d\n",linklist->nBxy) ;
  fprintf(fp,"nBxyz : %d\n",linklist->nBxyz) ;
  fprintf(fp,"\n");
  fprintf(fp,"MIN_X : %lf \n",linklist->MIN_X) ;
  fprintf(fp,"MIN_Y : %lf \n",linklist->MIN_Y) ;
  fprintf(fp,"MIN_Z : %lf \n",linklist->MIN_Z) ;
  fprintf(fp,"MAX_X : %lf \n",linklist->MAX_X) ;
  fprintf(fp,"MAX_Y : %lf \n",linklist->MAX_Y) ;
  fprintf(fp,"MAX_Z : %lf \n",linklist->MAX_Z) ;
  fprintf(fp,"****************************************\n");



  fclose(fp);
}
/*-----------------------------------------------------------------------------------------------*/

void output_place(int NumberOfParticles, int ParticleNumber_lambda, double ParticleDistance,double *Position,int *ParticleType, linklist_constant *linklist){

    FILE *fp_output;
    char outputfile[256];
    sprintf(outputfile,"poiseuille_%d_ramdom.dat",NumberOfParticles);
    fp_output = fopen(outputfile,"w");
    printf("output file : %s \n",outputfile);

	fprintf(fp_output,"%d\n",NumberOfParticles);
    fprintf(fp_output,"MIN_X : %lf, MAX_X : %lf\n",linklist->MIN_X,linklist->MAX_X);
    fprintf(fp_output,"MIN_Y : %lf, MAX_Y : %lf\n",linklist->MIN_Y,linklist->MAX_Y);
    fprintf(fp_output,"MIN_Z : %lf, MAX_Z : %lf\n",-6.0*ParticleDistance,6.0*ParticleDistance);
    for(int i=0; i<NumberOfParticles; i++){
        fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[i],Position[i*3],Position[i*3+1],Position[i*3+2]);
    }
	fprintf(fp_output,"lambda_particle : %d",ParticleNumber_lambda);
    fclose(fp_output);
}




/*---------------------------------------------------------------------------*/
/*
static void apply_rondom_move_PBC(int NumberOfParticles, double *Position, int *ParticleType, double ParticleDistance, linklist_constant *linklist){

  double rc1,rc2,rc3; //rondom coefficient
  int i;
  srand(1);
  // 速度の設定
  for(i = 0; i < NumberOfParticles; i++){
    if(ParticleType[i] == FLUID){
      rc1 = (double)rand()/RAND_MAX -0.5 ;   //-0.5 < rc < 0.5
      rc2 = (double)rand()/RAND_MAX -0.5 ;
      rc3 = (double)rand()/RAND_MAX -0.5 ;
      Position[3*i  ] += rc1 * ParticleDistance;
      Position[3*i+1] += rc2 * ParticleDistance;
      Position[3*i+2] += rc3 * ParticleDistance;

      //x
      if(Position[i*3]>=linklist->MAX_X)        {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3  ]<= linklist->MIN_X){Position[i*3]+=linklist->Lx;}

      //y
      if(Position[i*3+1]>=linklist->MAX_Y)      {Position[i*3+1]-=linklist->Ly;}
      else if(Position[i*3+1]<= linklist->MIN_Y){Position[i*3+1]+=linklist->Ly;}

      //z
      if(Position[i*3+2]>=linklist->MAX_Z)      {Position[i*3+2]-=linklist->Lz;}
      else if(Position[i*3+2]<= linklist->MIN_Z){Position[i*3+2]+=linklist->Lz;}
    }
  }

}*/

static void init_velosity(int pn,double *Velocity,int *ParticleType)
{
  int i,j;
  double u0,u1;
  double sum[3];
  int n_fluid=0;

  srand(1);
  /* 速度の設定 */
  for(i = 0; i < pn; i++){
    if(ParticleType[i] == FLUID){
      for(j=0;j<2;j++){
        u0 = (double)rand()/RAND_MAX;
        u1 = (double)rand()/RAND_MAX;
        Velocity[3*i+j] = 0.8 * sqrt(-2*log(u0))*cos(2*M_PI*u1);
      }
      Velocity[3*i+2] = 0.0;
      n_fluid++;
    }else{
      Velocity[3*i] = 0.0; Velocity[3*i+1] = 0.0; Velocity[3*i+2] = 0.0;
    }
  }

  //initialize sum
  for(i = 0; i < 3; i++){
      sum[i] = 0;
  }
  // 速度の和の計算
  for(i = 0; i < pn; i++){
    sum[0] += Velocity[3*i];  //vx->sum[0], vy->sum[1], vz->sum[2]
    sum[1] += Velocity[3*i+1];
    sum[2] += Velocity[3*i+2];
  }

  // 速度の和をゼロになるようにする
  for(i = 0; i < 3; i++){
      sum[i] /= n_fluid;
  }
  for(i = 0; i < pn; i++){
    if(ParticleType[i] == FLUID){
      Velocity[3*i  ] -= sum[0];
      Velocity[3*i+1] -= sum[1];
      Velocity[3*i+2] -= sum[2];
    }
  }

  /* 速度の和がゼロになっているかどうかの確認 */
  for(i = 0; i < 3; i++){
      sum[i] = 0;
  }

  for(i = 0; i < pn; i++){
    sum[0] += Velocity[3*i];  //vx->sum[0], vy->sum[1], vz->sum[2]
    sum[1] += Velocity[3*i+1];
    sum[2] += Velocity[3*i+2];
  }

  printf("sum of initial speed x:%f y:%f z:%f\n",sum[0],sum[1],sum[2]);

}

static void output_velocity(int pn, double *vl, int *ParticleType, char *output_dir){
  FILE *fp_output;
  char outputfile[256];
  sprintf(outputfile, "%s/velocity.txt",output_dir);
  fp_output = fopen(outputfile,"w");
  printf("output file : %s \n",outputfile);

  for(int i=0; i<pn; i++){
      fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[i],vl[i*3],vl[i*3+1],vl[i*3+2]);
  }
  fclose(fp_output);
}
