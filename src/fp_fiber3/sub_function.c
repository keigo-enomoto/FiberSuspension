#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "TypeDef.h"
#include "parameter.h"
#define BUFFER_SIZE 0x100


/**************************************************************
  *These two function is for reading input paramters          *
***************************************************************/

simulation_parameters *read_parameters_from_file(const char *input)
{
    FILE *fp;
    simulation_parameters *parameters;
    double Rho_s ;
    double Nu_f, Nu_s;
    char output_dir[BUFFER_SIZE];
    char Init_place[BUFFER_SIZE];
    parameters = (simulation_parameters *)malloc(sizeof(simulation_parameters));

    input_parameter input_parameters[] =
        {{"DIM","%d",&(parameters->DIM)},
         {"Rho_s","%lf",&Rho_s},
         {"Nu_f","%lf",&Nu_f},
         {"Nu_s","%lf",&Nu_s},
         {"strain_rate","%lf",&(parameters->strain_rate)},
         {"Gravity","%lf",&(parameters->Gravity)},
         {"DT","%lf",&(parameters->DT)},
         {"FINISH_TIME","%lf",&(parameters->FINISH_TIME)},
         {"OUTPUT_INTERVAL","%d",&(parameters->OUTPUT_INTERVAL)},
         {"output_dir","%s",&output_dir},
         {"Init_place","%s",&Init_place},
         {"ks","%lf",&(parameters->ks)},
         {"kb","%lf",&(parameters->kb)},
         {"VTU_flag","%d",&(parameters->VTU_flag)},
        };


    fp = fopen(input,"rt");

    read_input_parameters(fp,input_parameters);

    // parameters->wall_speed was defined in defWallSpeed() (CalculateFunction.c)

    parameters->ParticleDistance = 1.0; // if change, check calFiberRestoreForce()
    parameters->Rho[FLUID] = 1.0 ;
    parameters->Rho[SOLID] = Rho_s;
    parameters->Rho[WALL] = 1.0;
    parameters->Rho[DUMMY_WALL] = 1.0;
    if(parameters->strain_rate != 0.0){
      parameters->Re_inv[FLUID] = Nu_f / parameters->strain_rate; // l_0(ParticleDistance) is omitted because l_0 = 1.0
      parameters->Re_inv[SOLID] = Nu_s / parameters->strain_rate; // l_0(ParticleDistance) is omitted because l_0 = 1.0
      parameters->Re_inv[WALL]  = Nu_f / parameters->strain_rate; // l_0(ParticleDistance) is omitted because l_0 = 1.0
      parameters->Re_inv[DUMMY_WALL]  = Nu_f / parameters->strain_rate; // l_0(ParticleDistance) is omitted because l_0 = 1.0
    }else{
      parameters->Re_inv[FLUID] = Nu_f ; // if strain_rate == 0.0, strain_rate = 1.0 (wall speed is 0.0)
      parameters->Re_inv[SOLID] = Nu_s ;
      parameters->Re_inv[WALL]  = Nu_f ;
      parameters->Re_inv[DUMMY_WALL]  = Nu_f ;
    }

    parameters->output_dir = malloc((strlen(output_dir) + 1) * sizeof(char));
    strcpy(parameters->output_dir,output_dir);

    parameters->Init_place = malloc((strlen(Init_place) + 1) * sizeof(char));
    strcpy(parameters->Init_place,Init_place);

    fclose(fp);

    return parameters;
}


// -----------------------------------------------------------------------------

void output_parameter(simulation_parameters *parameters, Gradient_constant *constant,linklist_constant *linklist, int NumberOfParticles){
  FILE *fp;
  char fileName[1024];

  sprintf(fileName, "%s/input_parameters.txt",parameters->output_dir); //sprintf is the function to write the word to array
  fp = fopen(fileName,"w");

    fprintf(fp,"*******input parameter *********\n");
    fprintf(fp,"DIM : %d\n",(parameters->DIM)) ;
    fprintf(fp,"ParticleDistance : %lf\n",(parameters->ParticleDistance)) ;
    fprintf(fp,"Rho_fluid  : %lf\n",(parameters->Rho[FLUID])) ;
    fprintf(fp,"Rho_solid  : %lf\n",(parameters->Rho[SOLID])) ;
    fprintf(fp,"Rho_wall  : %lf\n",(parameters->Rho[WALL])) ;
    fprintf(fp,"Re_inv fluid  : %lf\n",(parameters->Re_inv[FLUID])) ;
    fprintf(fp,"Re_inv solid  : %lf\n",(parameters->Re_inv[SOLID])) ;
    fprintf(fp,"Re_inv wall : %lf\n",(parameters->Re_inv[WALL])) ;
    fprintf(fp,"strain_rate : %lf\n",(parameters->strain_rate)) ;
    fprintf(fp,"wall_speed : %lf\n",(parameters->wall_speed)) ;
    fprintf(fp,"Gravity : %lf\n",(parameters->Gravity)) ;
    fprintf(fp,"DT : %lf\n",(parameters->DT)) ;
    fprintf(fp,"FINISH_TIME : %lf\n",(parameters->FINISH_TIME)) ;
    fprintf(fp,"OUTPUT_INTERVAL : %d\n",(parameters->OUTPUT_INTERVAL)) ;
    fprintf(fp,"output_dir : %s\n",parameters->output_dir) ;
    fprintf(fp,"Init_place : %s\n",parameters->Init_place) ;
    fprintf(fp,"ks : %lf\n",(parameters->ks)) ;
    fprintf(fp,"kb : %lf\n",(parameters->kb)) ;
    fprintf(fp,"VTU_flag : %d\n",(parameters->VTU_flag)) ;
    fprintf(fp,"\n");

  fprintf(fp,"*******constant parameter *********\n");
  fprintf(fp,"Re_forNumberDensity %lf \n",constant -> Re_forNumberDensity);
  fprintf(fp,"Re_forGradient : %lf \n",constant -> Re_forGradient);
  fprintf(fp,"Re_forLaplacian : %lf \n",constant -> Re_forLaplacian);
  fprintf(fp,"N0_forNumebrDensity: %lf \n",constant -> N0_forNumberDensity);
  fprintf(fp,"N0_forGradiant: %lf\n",constant -> N0_forGradient);
  fprintf(fp,"N0_forLaplacian : %lf \n",constant -> N0_forLaplacian);
  fprintf(fp,"Lambda : %lf \n",constant -> Lambda);
  fprintf(fp,"\n");

  fprintf(fp,"******** linklist parameters ***********\n");
  fprintf(fp,"DB  : %lf\n",linklist->DB) ;
  fprintf(fp,"nBx : %d\n",linklist->nBx) ;
  fprintf(fp,"nBy : %d\n",linklist->nBy) ;
  fprintf(fp,"nBz : %d\n",linklist->nBz) ;
  fprintf(fp,"nBxy : %d\n",linklist->nBxy) ;
  fprintf(fp,"nBxyz : %d\n",linklist->nBxyz) ;
  fprintf(fp,"LX : %lf \n",linklist->Lx) ;
  fprintf(fp,"LY : %lf \n",linklist->Ly) ;
  fprintf(fp,"LZ : %lf \n",linklist->Lz) ;
  fprintf(fp,"MIN_X : %lf \n",linklist->MIN_X) ;
  fprintf(fp,"MIN_Y : %lf \n",linklist->MIN_Y) ;
  fprintf(fp,"MIN_Z : %lf \n",linklist->MIN_Z) ;
  fprintf(fp,"MAX_X : %lf \n",linklist->MAX_X) ;
  fprintf(fp,"MAX_Y : %lf \n",linklist->MAX_Y) ;
  fprintf(fp,"MAX_Z : %lf \n",linklist->MAX_Z) ;

  fclose(fp);
}

//---------------------------------------------------------------------------------------------------
/**************************
 * For analysis efficiency
 **************************/

void output_analysis(int n, double ParticleDistance,int iTimeStep, simulation_parameters *parameters, int n_fluid, int n_solid, double Time, int ny, clock_t start){
  FILE *fp;
  char outputfile[1024];
  double width = (ny - 5) * ParticleDistance ;
  clock_t end = clock();

  sprintf(outputfile, "%s/output_analysis.txt",parameters->output_dir); //sprintf is the function to write the word to array
  fp = fopen(outputfile,"w");
  fprintf(fp,"NumberOfParticle : %d\n",n);
  fprintf(fp,"IterationTime : %d\n",iTimeStep);
  fprintf(fp,"n_fluid : %d\n",n_fluid);
  fprintf(fp,"n_solid : %d\n",n_solid);
  fprintf(fp,"width_y : %lf\n",width);
  fprintf(fp,"Past_Time : %lf \n",Time);
  fprintf(fp,"Past_Time(real) : %lf",(double)(end - start)/ CLOCKS_PER_SEC);
  fclose(fp);
}

//-----------------------------------------------------------------------------------------------
/******************************
 * Check poiseuille or Couette
 ******************************/

void check_flow_type(simulation_parameters *parameters){
  if (parameters->wall_speed == 0.0 && parameters->Gravity==0.0){
    fprintf(stderr,"either Gravity or wall_speed must be nonzero value \n");
    exit(3);
  }else if (parameters->wall_speed != 0.0 && parameters->Gravity!=0.0){
    fprintf(stderr,"either Gravity or wall_speed must be zero \n");
    exit(3);
  }
}

//-----------------------------------------------------------------------------------------------

/**************************
 * Output velocity profile
 **************************/

// for couette flow
void output_couette(int n, double *Position, double *Velocity, int *ParticleType, simulation_parameters *parameters, double Time, int iTimeStep, clock_t start){
  int i;
  FILE *fp;
  char outputfile[1024];
  int count = 0;
  clock_t end = clock();

  sprintf(outputfile, "%s/output_couette.dat",parameters->output_dir); //sprintf is the function to write the word to array
  fp = fopen(outputfile,"a");
  fprintf(fp,"\n");
  fprintf(fp,"TimeStepNumber: %4d   Time: %lf   NumberOfParticles: %d\n", iTimeStep, Time, n);
  fprintf(fp,"  past time: %lf sec   DT : %lf\n",(double)(end - start)/ CLOCKS_PER_SEC,parameters->DT);
  for(i=0; i<n; i++){
    if(ParticleType[i] == FLUID){
      fprintf(fp,"%d y: %lf , vx : %lf \n",count,Position[3*i+1],Velocity[3*i]);
      count++;
    }
  }
  fclose(fp);
}


void output_poiseuille(int n, double *Position, double *Velocity, int *ParticleType, simulation_parameters *parameters, double Time, int iTimeStep, clock_t start){
  int i;
  FILE *fp;
  char outputfile[1024];
  int count = 0;
  clock_t end = clock();

  sprintf(outputfile, "%s/output_poiseuille.dat",parameters->output_dir); //sprintf is the function to write the word to array
  fp = fopen(outputfile,"a");
  fprintf(fp,"\n");
  fprintf(fp,"TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticles: %d\n", iTimeStep, Time, n);
  fprintf(fp,"  past time: %lf sec   DT : %lf\n",(double)(end - start)/ CLOCKS_PER_SEC,parameters->DT);
  for(i=0; i<n; i++){
    if(ParticleType[i] == FLUID){
      fprintf(fp,"%d y: %lf , vx : %lf \n",count,Position[3*i+1],Velocity[3*i]);
      count++;
    }
  }
  fclose(fp);
}

//------------------------------------------------------------------------------------------------

// for debug
/*
void output_accel(int n, int *ParticleType, double *Acceleration, simulation_parameters *parameters){
  FILE *fp;
  char outputfile[1024];

  sprintf(outputfile, "%s/output_test.txt",parameters->output_dir); //sprintf is the function to write the word to array
  fp = fopen(outputfile,"w");
  for(int i=0; i<n; i++){
    fprintf(fp,"%D %lf %lf %lf \n",ParticleType[i],Acceleration[3*i],Acceleration[3*i+1],Acceleration[3*i+2]);
  }
  fprintf(fp,"---------------------------------------------------\n\n");
  fclose(fp);
}
*/

void output_coordinate(int n_fiber, int aspect_ratio,int NumberOfParticles, int n_fluid, linklist_constant *linklist,
double ParticleDistance, simulation_parameters *parameters, int *ParticleType, double *Position){

	FILE *fp_output;
	char outputfile[256];
	int count = 0;
	// int calc_lambda_particle;
  int nx = (int)(linklist->Lx / ParticleDistance);
  int ny = (int)(linklist->Ly / ParticleDistance);

	// sprintf(outputfile,"%s/f%d-%d_%d_%d-%d_2d.dat",parameters->output_dir,n_fiber,aspect_ratio,NumberOfParticles,nx,ny);
	sprintf(outputfile,"%s/%s",parameters->output_dir,parameters->Init_place);
	fp_output = fopen(outputfile, "w");
	// printf("output file : %s \n",outputfile);
	fprintf(fp_output,"%d\n",NumberOfParticles);
	fprintf(fp_output,"%d\n",n_fluid);
	fprintf(fp_output,"%d\n",n_fiber);
	fprintf(fp_output,"%d\n",aspect_ratio);
  fprintf(fp_output,"MIN_X : %lf, MAX_X : %lf\n",linklist->MIN_X,linklist->MAX_X);
  fprintf(fp_output,"MIN_Y : %lf, MAX_Y : %lf\n",linklist->MIN_Y,linklist->MAX_Y);
  fprintf(fp_output,"MIN_Z : %lf, MAX_Z : %lf\n",-8.0*ParticleDistance,8.0*ParticleDistance);

	count=0;
	// Firstly, fiber partiles are printed
	for(int i=0;i<n_fiber*aspect_ratio;i++){
		// if(i == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[i],Position[i*3],Position[i*3+1],Position[i*3+2]);
		count++;
	}

	// other particle are printed
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iy*nx + ix;
		if(ParticleType[ip]==GHOST || ParticleType[ip]==SOLID) continue;
		// if(ip == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		count++;
	}}

	if(count != NumberOfParticles){
		// printf("don't coresspond count with Number of particles \n");
	}

	// fprintf(fp_output,"lambda_particle : %d\n",calc_lambda_particle);
	fprintf(fp_output,"lambda_particle : -1\n");

	// output fiber consist
	for (int i = 0; i < aspect_ratio * n_fiber ; i++){
		fprintf(fp_output,"%d\n",i);
	}

	fclose(fp_output);



}

// to continue calculate, output coordinate, velocity, Time
void output_continue(int n_fiber, int aspect_ratio,int NumberOfParticles, int n_fluid, linklist_constant *linklist,
double ParticleDistance, simulation_parameters *parameters, int *ParticleType, double *Position, double *Velocity, double Time){

	FILE *fp_output;
	char outputfile[256];
	int count = 0;
	// int calc_lambda_particle;
  int nx = (int)(linklist->Lx / ParticleDistance);
  int ny = (int)(linklist->Ly / ParticleDistance);

	// sprintf(outputfile,"%s/f%d-%d_%d_%d-%d_2d.dat",parameters->output_dir,n_fiber,aspect_ratio,NumberOfParticles,nx,ny);
	sprintf(outputfile,"%s/continue.dat",parameters->output_dir);
	fp_output = fopen(outputfile, "w");
	// printf("output file : %s \n",outputfile);
	fprintf(fp_output,"%d\n",NumberOfParticles);
	fprintf(fp_output,"%d\n",n_fluid);
	fprintf(fp_output,"%d\n",n_fiber);
	fprintf(fp_output,"%d\n",aspect_ratio);
  fprintf(fp_output,"MIN_X : %lf, MAX_X : %lf\n",linklist->MIN_X,linklist->MAX_X);
  fprintf(fp_output,"MIN_Y : %lf, MAX_Y : %lf\n",linklist->MIN_Y,linklist->MAX_Y);
  fprintf(fp_output,"MIN_Z : %lf, MAX_Z : %lf\n",linklist->MIN_Z,linklist->MAX_Z);
	fprintf(fp_output,"%lf\n",Time);


	count=0;

	for(int i=0;i<NumberOfParticles;i++){
		// if(i == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[i],Position[i*3],Position[i*3+1],Position[i*3+2]);
		fprintf(fp_output,"v %lf %lf %lf\n",Velocity[i*3],Velocity[i*3+1],Velocity[i*3+2]);
		count++;
	}

	if(count != NumberOfParticles){
		// printf("don't coresspond count with Number of particles \n");
	}

	// fprintf(fp_output,"lambda_particle : %d\n",calc_lambda_particle);
	fprintf(fp_output,"lambda_particle : -1\n");

	// output fiber consist
	for (int i = 0; i < aspect_ratio * n_fiber ; i++){
		fprintf(fp_output,"%d\n",i);
	}

	fclose(fp_output);



}



/**************************
 * Output data in vtu format
 **************************/


void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, double *NumberDensity, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024]; //This is array to have the word of filename (fileName[0] = p, fileName[1] = a)
  sprintf(fileName, "%s/particle_%05d.vtu",output_dir ,FileNumber); //sprintf is the function to write the word to array


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
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Velocity_arr' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",Velocity[i*3],Velocity[i*3+1],Velocity[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( Velocity[i*3]*Velocity[i*3] + Velocity[i*3+1]*Velocity[i*3+1] + Velocity[i*3+2]*Velocity[i*3+2]);
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='number_density' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)NumberDensity[i]);
  }

  // fprintf(fp,"</DataArray>\n");
  // fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='boundarycondition' format='ascii'>\n");
  // for(i=0;i<NumberOfParticles;i++){
  //   fprintf(fp,"%d\n",BoundaryCondition[i]);
  // }

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
/*****************************
 * output only wall and solid
 ****************************/
void writeData_inVtuFormat_fiber_wall(int DIM,int *ParticleType, double *Position, double *Velocity,double *Pressure, double *NumberDensity, int *BoundaryCondition, int NumberOfParticles, int FileNumber ,int *BucketIndex, char *output_dir){
  int i, j;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024]; //This is array to have the word of filename (fileName[0] = p, fileName[1] = a)
  sprintf(fileName, "%s/particle_%05d.vtu",output_dir ,FileNumber); //sprintf is the function to write the word to array
  int *exclude_fluid ;  // hold index of particle
  int N = 0;            // the number of particles excludeint fluid

  exclude_fluid = (int *)malloc(NumberOfParticles * sizeof(int));

  // check the number of particles
  for ( i = 0; i < NumberOfParticles; i++){
    if(ParticleType[i] != FLUID){
      exclude_fluid[N] = i;
      N ++;
    }
  }

  for( i=N; i<NumberOfParticles; i++){
    exclude_fluid[i] = -1 ;
  }

  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",N,N);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%lf %lf %lf\n",Position[j*3],Position[j*3+1],Position[j*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%d\n",ParticleType[j]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Velocity_arr' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%lf %lf %lf\n",Velocity[j*3],Velocity[j*3+1],Velocity[j*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    absoluteValueOfVelocity=
      sqrt( Velocity[j*3]*Velocity[j*3] + Velocity[j*3+1]*Velocity[j*3+1] + Velocity[j*3+2]*Velocity[j*3+2]);
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%f\n",(float)Pressure[j]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='number_density' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%f\n",(float)NumberDensity[j]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='boundarycondition' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%d\n",BoundaryCondition[j]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='BucketIndex' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%d\n",BucketIndex[j]);
  }

  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"%d\n",i+1);
  }
  //点(粒子)のデータなら'1'が入る
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<N;i++){
    j = exclude_fluid[i] ;
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  free(exclude_fluid);
  fclose(fp);
}

//-------------------------------------------------------------------------------------------------
void writeData_inVtuFormat_tmp(int DIM,int *ParticleType, double *Position,int NumberOfParticles){
  int i;
  FILE *fp;
  char fileName[1024]; //This is array to have the word of filename (fileName[0] = p, fileName[1] = a)
  sprintf(fileName, "particle_tmp.vtu"); //sprintf is the function to write the word to array

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


/**************************
 * for parallel computing
 **************************/
#ifdef _OPENMP

void print_omp(){
  FILE *fp;

  fp = fopen("omp-info.txt","w");
  fprintf(fp,"max number of threads : %d \n",omp_get_max_threads());
  fclose(fp);

}

#endif