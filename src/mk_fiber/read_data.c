/********************************************
 * This is test program.
 * read *.dat file and translate *.vtu file
 * [ how to execute ]
 * ./read fp_256_2d.dat
 *
 * [ output ]
 * particle_*D.vtu
 ********************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <string.h>
//#include <time.h>


void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, int NumberOfParticles, int FileNumber );
//* beginning of main *************************************************************
int main( int argc, char** argv )
{

  FILE *fp;
  int DIM=2;
  int NumberOfParticles,n_fluid;
  int *ParticleType;
  double *Position;
  double MIN_X,MAX_X;
  double MIN_Y,MAX_Y;
  double MIN_Z,MAX_Z;
  int n_fiber,aspect_ratio;

  fp = fopen(argv[1], "r");
	fscanf(fp,"%d\n",&NumberOfParticles);
	printf("nP: %d\n",NumberOfParticles);
	fscanf(fp,"%d\n",&n_fluid);
	printf("n_fluid: %d\n",n_fluid);
  fscanf(fp,"%d\n",&n_fiber);
  printf("n_fiber : %d\n",n_fiber);
  fscanf(fp,"%d\n",&aspect_ratio);
  printf("aspect_ratio : %d\n",aspect_ratio);


  fscanf(fp,"MIN_X : %lf, MAX_X : %lf\n",&MIN_X,&MAX_X);
  fscanf(fp,"MIN_Y : %lf, MAX_Y : %lf\n",&MIN_Y,&MAX_Y);
  fscanf(fp,"MIN_Z : %lf, MAX_Z : %lf\n",&MIN_Z,&MAX_Z);
  printf("MIN_X : %lf, MAX_X : %lf\n",MIN_X,MAX_X);
  printf("MIN_Y : %lf, MAX_Y : %lf\n",MIN_Y,MAX_Y);
  printf("MIN_Z : %lf, MAX_Z : %lf\n",MIN_Z,MAX_Z);

  Position     = (double *)malloc(3 * NumberOfParticles * sizeof(double));
  ParticleType = (int *)malloc(NumberOfParticles * sizeof(int));
	for(int i=0;i<NumberOfParticles;i++) {
		fscanf(fp,"%d %lf %lf %lf\n",&ParticleType[i],&Position[i*3],&Position[i*3+1],&Position[i*3+2]);
	}
	fclose(fp);

  writeData_inVtuFormat(DIM,ParticleType,Position,NumberOfParticles,0);

  return 0;
}
//* end of main *************************************************************

void writeData_inVtuFormat(int DIM,int *ParticleType, double *Position, int NumberOfParticles, int FileNumber ){
  int i;
  FILE *fp;
  const char *fileName; //This is array to have the word of filename (fileName[0] = p, fileName[1] = a)
  if(DIM == 2){
    fileName = "particle_2D.vtu";
  }else if(DIM == 3){
    fileName = "particle_3D.vtu";
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
