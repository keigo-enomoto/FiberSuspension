
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TypeDef.h"

#define PI 3.14159265358979323846264338327950288

/********************************
 *calculate the weight function *
 ********************************/
static double weight(double distance, double re){
    double weight_ij = 0.0;

    if( distance >= re ){
        weight_ij = 0.0;
    }else{
        weight_ij = re / distance - 1 ;
    }

    return weight_ij ;
}

/*************************************************************
 * calculate distance between i th particle and j th particle
 *************************************************************/
/*
static double calDistance_ij(int i, int j, double *Position){
  double xij, yij, zij ;

  xij = Position[j*3  ] - Position[i*3  ];
  yij = Position[j*3+1] - Position[i*3+1];
  zij = Position[j*3+2] - Position[i*3+2];

  return  sqrt((xij*xij) + (yij*yij) + (zij*zij)); //the distance

}
*/


static double calDistance_ij_pbc(int i, int j, double *Position,linklist_constant *linklist){
  double xij, yij, zij ;

  xij = Position[j*3  ] - Position[i*3  ];
  yij = Position[j*3+1] - Position[i*3+1];
  zij = Position[j*3+2] - Position[i*3+2];

  // apply pbc
  if(xij < -linklist->Lx*0.50)  { xij += linklist->Lx; }
  else if(xij>linklist->Lx*0.50){ xij -= linklist->Lx; }

  if(yij < -linklist->Ly*0.50)  { yij += linklist->Ly; }
  else if(yij>linklist->Ly*0.50){ yij -= linklist->Ly; }

  if(zij < -linklist->Lz*0.50)  { zij += linklist->Lz; }
  else if(zij>linklist->Lz*0.50){ zij -= linklist->Lz; }

  return  sqrt((xij*xij) + (yij*yij) + (zij*zij)); //the distance

}

/********************************
 * PBC operation
 *******************************/

// apply pbc for xij , yij, and zij
static void apply_pbc_distance_each_axis(double *xij, double *yij, double *zij, linklist_constant *linklist){

  // apply pbc
  if(*xij < -linklist->Lx*0.50)  { *xij += linklist->Lx; }
  else if(*xij>linklist->Lx*0.50){ *xij -= linklist->Lx; }

  if(*yij < -linklist->Ly*0.50)  { *yij += linklist->Ly; }
  else if(*yij>linklist->Ly*0.50){ *yij -= linklist->Ly; }

  if(*zij < -linklist->Lz*0.50)  { *zij += linklist->Lz; }
  else if(*zij>linklist->Lz*0.50){ *zij -= linklist->Lz; }

}

// if using peliodic boundary, use this function
static int calBucketIndex_pbc(int b_jx, int b_jy, int b_jz, linklist_constant *linklist){
  int tb_jx, tb_jy, tb_jz; //the temporary neighbor bucket's coordinate

  tb_jx = (b_jx + (linklist->nBx)) % (linklist->nBx);
  tb_jy = (b_jy + (linklist->nBy)) % (linklist->nBy);
  tb_jz = (b_jz + (linklist->nBz)) % (linklist->nBz);

  return  tb_jz*(linklist->nBxy) + tb_jy*(linklist->nBx) + tb_jx;
}

/**************************************************************
  *These two function is for calculating constant paramters   *
***************************************************************/

//only use in calcConstantParameter()
static void calNZeroAndLambda( int DIM, double ParticleDistance, int NumberOfParticles, double *Position, int *ParticleType, Gradient_constant *constant, int ParticleNumber_lambda){

  int i=0,j=0; // i : the index number of INNER PARTICLE
  double xij,yij,zij;
  double distance, distance2;

  i = ParticleNumber_lambda;
  constant -> N0_forNumberDensity = 0.0;
  constant -> N0_forGradient      = 0.0;
  constant -> N0_forLaplacian     = 0.0;
  constant -> Lambda              = 0.0;

  if(i == -1){ // use lattice version constant parameter
    if(DIM==2){
      constant -> N0_forNumberDensity = 6.539697;
      constant -> N0_forGradient      = 6.539697;
      constant -> N0_forLaplacian     = 18.976417;
      constant -> Lambda              = 2.448472;
    }else if(DIM == 3){ // re_laplacian = 3.1 * l_0 and lattice version
      constant -> N0_forNumberDensity = 14.418575;
      constant -> N0_forGradient      = 14.418575;
      constant -> N0_forLaplacian     = 54.321209;
      constant -> Lambda              = 3.243422;
    }
  }else{
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      constant -> N0_forNumberDensity += weight(distance, (constant -> Re_forNumberDensity));
      constant -> N0_forGradient      += weight(distance, (constant -> Re_forGradient));
      constant -> N0_forLaplacian     += weight(distance, (constant -> Re_forLaplacian));
      constant -> Lambda              += distance2 * weight(distance, (constant -> Re_forLaplacian));
    }
    (constant -> Lambda) = (constant -> Lambda)/(constant -> N0_forLaplacian);
  }

}

/*--------------------------------------------------------------------------------------------------------*/

Gradient_constant *calConstantParameter(int DIM, double ParticleDistance,int NumberOfParticles, double *Position, int *ParticleType, int ParticleNumber_lambda){

    Gradient_constant *constant;
    constant = malloc(sizeof(Gradient_constant));

    constant -> Re_forNumberDensity  = 2.1 * ParticleDistance;
    constant -> Re_forGradient       = 2.1 * ParticleDistance;
    if(DIM == 2){
      constant -> Re_forLaplacian      = 3.1 * ParticleDistance;
    }else if(DIM == 3){
      constant -> Re_forLaplacian      = 2.1 * ParticleDistance;
    }else{
      fprintf(stderr,"*** invaild Dimension : %d ***\n",DIM);
      exit(3);
    }
    calNZeroAndLambda(DIM, ParticleDistance,NumberOfParticles,Position,ParticleType,constant,ParticleNumber_lambda);
    return constant;
}


/********************************************
 * These 2 function are for linklist        *
 ********************************************/
linklist_constant *calcLinklistParameter(double ParticleDistance, double Re_forLaplacian, range_sim *range_sim){
  const double Courant_Number = 0.2 ; //クーラン数
  linklist_constant *linklist;
  double temp_DBx = 0.0;
  int temp_divide_number = 0; // x方向について
  linklist = malloc(sizeof(linklist_constant));

  linklist -> MIN_X = range_sim->MIN_X ;//+ (-4.0) * ParticleDistance
  linklist -> MIN_Y = range_sim->MIN_Y ;//+ (-4.0) * ParticleDistance
  linklist -> MIN_Z = range_sim->MIN_Z ;//+ (-4.0) * ParticleDistance
  linklist -> MAX_X = range_sim->MAX_X ;//+ ParticleDistance * 3.0;
  linklist -> MAX_Y = range_sim->MAX_Y ;//+ ParticleDistance * 30;
  linklist -> MAX_Z = range_sim->MAX_Z ;//+ ParticleDistance * 3.0;
  linklist -> Lx    = linklist -> MAX_X - linklist -> MIN_X ;
  linklist -> Ly    = linklist -> MAX_Y - linklist -> MIN_Y ;
  linklist -> Lz    = linklist -> MAX_Z - linklist -> MIN_Z ;

  temp_DBx = Re_forLaplacian * (1.0 + Courant_Number);
  temp_divide_number = (int)(linklist->Lx/temp_DBx);      // use Lx (Since x direction is peliodic boundary, Lx must be `int * DB` )
  linklist -> DB    = linklist->Lx/(temp_divide_number);  // DB > temp_DBx --> safety (labnote vol.2 p108) if divide_number decrease, DB will increase
  linklist -> nBx   =(int)((linklist->Lx)/(linklist->DB));
  linklist -> nBy   =(int)((linklist->Ly)/(linklist->DB));
  linklist -> nBz   =(int)((linklist->Lz)/(linklist->DB));
  linklist -> nBxy  =(linklist->nBx) * (linklist->nBy);
  linklist -> nBxyz = (linklist->nBxy) * (linklist->nBz);

  //printf("tempnumber : %d\n",temp_divide_number);
  if(linklist->nBx < 3|| linklist->nBy < 3|| linklist->nBz < 3){
    fprintf(stderr,"******** WARNING ***********\n");
    fprintf(stderr,"check the simulation size \n");
    fprintf(stderr,"DB size : %lf \n",linklist->DB);
    exit(1);
  }

  return linklist;
}

/*-------------------------------------------------------------------------------------------*/
void makeBucket(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist){
  int i = 0;
  int ix, iy, iz; //the coordinate of bucket
  int ib ;        //index of bucket
  int temp_BucketLast;
  int nBxy = linklist->nBxy; int nBxyz = linklist->nBxyz;
  int nBx = linklist->nBx; double DBinv = 1/(linklist->DB);
  double MIN_X = linklist->MIN_X;
  double MIN_Y = linklist->MIN_Y;
  double MIN_Z = linklist->MIN_Z;

  //initialize these three array
  for(i=0; i<nBxyz; i++){BucketFirst[i] = -1; BucketLast[i] = -1;}
  for(i=0; i<NumberOfParticles ; i++){Nextof[i] = -1;}


  //set Bucket
  for(int i=0;i<NumberOfParticles;i++){
		if(ParticleType[i] == GHOST)continue;
    ix = (int)((Position[i*3  ] - MIN_X)*DBinv);
    iy = (int)((Position[i*3+1] - MIN_Y)*DBinv);
    iz = (int)((Position[i*3+2] - MIN_Z)*DBinv);
    ib = iz*nBxy + iy*nBx + ix;  // bucket index (0 <= ib <= nBxy)

    /*
    //remove ghost particle (water collapse)
    if((ix < 1) || (iy < 1) || (iz < 1) || (ix >= linklist->nBx) || (iy >= linklist->nBy) || (iz >= linklist->nBz)){
      ParticleType[i] = GHOST;
      continue;
    }*/

    temp_BucketLast = BucketLast[ib]; //temporarily save the label of last particle in the bucket
    BucketLast[ib] = i;
		if(temp_BucketLast == -1){
      BucketFirst[ib] = i;
    }else{
      Nextof[temp_BucketLast] = i;
    }
	}

}

void makeBucket_withIndex(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist,int *BucketIndex){
  int i = 0;
  int ix, iy, iz; //the coordinate of bucket
  int ib ;        //index of bucket
  int temp_BucketLast;
  int nBxy = linklist->nBxy; int nBxyz = linklist->nBxyz;
  int nBx = linklist->nBx; double DBinv = 1/(linklist->DB);
  double MIN_X = linklist->MIN_X;
  double MIN_Y = linklist->MIN_Y;
  double MIN_Z = linklist->MIN_Z;

  //initialize these three array
  for(i=0; i<nBxyz; i++){BucketFirst[i] = -1; BucketLast[i] = -1;}
  for(i=0; i<NumberOfParticles ; i++){Nextof[i] = -1;}


  //set Bucket
  for(int i=0;i<NumberOfParticles;i++){
		if(ParticleType[i] == GHOST)continue;
    ix = (int)((Position[i*3  ] - MIN_X)*DBinv);
    iy = (int)((Position[i*3+1] - MIN_Y)*DBinv);
    iz = (int)((Position[i*3+2] - MIN_Z)*DBinv);
    ib = iz*nBxy + iy*nBx + ix;  // bucket index (0 <= ib <= nBxy)
    BucketIndex[i] = ib;

    /*
    //remove ghost particle (water collapse)
    if((ix < 1) || (iy < 1) || (iz < 1) || (ix >= linklist->nBx) || (iy >= linklist->nBy) || (iz >= linklist->nBz)){
      ParticleType[i] = GHOST;
      continue;
    }*/

    temp_BucketLast = BucketLast[ib]; //temporarily save the label of last particle in the bucket
    BucketLast[ib] = i;
		if(temp_BucketLast == -1){
      BucketFirst[ib] = i;
    }else{
      Nextof[temp_BucketLast] = i;
    }
	}

}

/***************************************************
 * These 3 function are for velosity update to u^* *
 ***************************************************/

void calGravity_x( double NumberOfParticles, int *ParticleType, double *Acceleration ,double Gravity){
  int i=0;
  const double  G_X = Gravity;
  const double  G_Y = 0.0 ;
  const double  G_Z = 0.0 ;

// #pragma omp parallel for
  for(i=0; i<NumberOfParticles; i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Acceleration[i*3  ]=G_X;
      Acceleration[i*3+1]=G_Y;
      Acceleration[i*3+2]=G_Z;
    }else{
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;
    }
  }
}
//----------------------------------------------------------------------------
/**************************************************
 * Calculate Viscosity
 ******************************************************/
// ==========================================================================================================================

// cell list in previous version
void calViscosity( int NumberOfParticles,
                    int *ParticleType,double *Position,double *Velocity,double *Acceleration,int *BucketFirst, int *Nextof,
                    simulation_parameters *parameters, Gradient_constant *constant, linklist_constant *linklist){
  int i,j;
  double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z; //the inner contents of summention
  double distance ;
  double w;//weight
  double a;//coefficient
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //index for jth roop
  int ib,jb ;        //index of bucket
  double Re_forLaplacian = constant->Re_forLaplacian;
  int type_i = GHOST, type_j = GHOST;

  a = (2.0*(parameters->DIM))/((constant->N0_forLaplacian)*(constant->Lambda));

// don't need reduction for viscosity_term_x etc Accerelation
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,type_i,viscosityTerm_x,viscosityTerm_y,viscosityTerm_z,b_jz,b_jy,b_jx,jb,j,type_j,distance,w)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divided to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      type_i = ParticleType[i];
      if(type_i == FLUID || type_i == SOLID){
      viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;  viscosityTerm_z = 0.0;
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            type_j = ParticleType[j];
            if((j!=i) && (type_j!=GHOST)){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<Re_forLaplacian){
                w =  weight(distance, Re_forLaplacian);
                viscosityTerm_x +=(Velocity[j*3  ]-Velocity[i*3  ])*w;
                viscosityTerm_y +=(Velocity[j*3+1]-Velocity[i*3+1])*w;
                viscosityTerm_z +=(Velocity[j*3+2]-Velocity[i*3+2])*w;
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        Acceleration[i*3  ] += viscosityTerm_x * (parameters->Re_inv[type_i]) * a;
        Acceleration[i*3+1] += viscosityTerm_y * (parameters->Re_inv[type_i]) * a;
        Acceleration[i*3+2] += viscosityTerm_z * (parameters->Re_inv[type_i]) * a; //use '+=' because we also have to condisder acceleration by gravity
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}
}


// calculate stress for wall particle -------------------------------------------------------------------

void calViscosity_stress( int NumberOfParticles,
                    int *ParticleType,double *Position,double *Velocity,double *Acceleration,int *BucketFirst, int *Nextof,
                    simulation_parameters *parameters, Gradient_constant *constant, linklist_constant *linklist){
  int i,j;
  double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z; //the inner contents of summention
  double distance ;
  double w;//weight
  double a, a2;//coefficient
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //index for jth roop
  int ib,jb ;        //index of bucket
  double Re_forLaplacian = constant->Re_forLaplacian;
  int type_i = GHOST, type_j = GHOST;

  a = (2.0*(parameters->DIM))/((constant->N0_forLaplacian)*(constant->Lambda));

// don't need reduction for viscosity_term_x etc Accerelation
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,type_i,viscosityTerm_x,viscosityTerm_y,viscosityTerm_z,b_jz,b_jy,b_jx,jb,j,type_j,distance,w)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divide to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      type_i = ParticleType[i];
      if(type_i == FLUID || type_i == SOLID){
      viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;  viscosityTerm_z = 0.0;
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            type_j = ParticleType[j];
            if((j!=i) && (type_j!=GHOST)){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<Re_forLaplacian){
                w =  weight(distance, Re_forLaplacian);
                viscosityTerm_x +=(Velocity[j*3  ]-Velocity[i*3  ])*w;
                viscosityTerm_y +=(Velocity[j*3+1]-Velocity[i*3+1])*w;
                viscosityTerm_z +=(Velocity[j*3+2]-Velocity[i*3+2])*w;

                if(type_j == WALL){ // stress for wall
                  a2 = (parameters->Re_inv[type_i]) * a * (parameters->Rho[type_i]);
                  Acceleration[j*3  ] -= (Velocity[j*3  ]-Velocity[i*3  ])* w * a2;
                  Acceleration[j*3+1] -= (Velocity[j*3+1]-Velocity[i*3+1])* w * a2;
                  Acceleration[j*3+2] -= (Velocity[j*3+2]-Velocity[i*3+2])* w * a2;
                }
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        Acceleration[i*3  ] += viscosityTerm_x * (parameters->Re_inv[type_i]) * a;
        Acceleration[i*3+1] += viscosityTerm_y * (parameters->Re_inv[type_i]) * a;
        Acceleration[i*3+2] += viscosityTerm_z * (parameters->Re_inv[type_i]) * a; //use '+=' because we also have to condisder acceleration by gravity
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}
}

// ==========================================================================================================================

/*********************************
 * Calculate temporary velocity  *
 *********************************/

// don't apply PBC here (the result of calculation does'nt change [labnote vol2. p4])
void moveParticle( int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration){

  for(int i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID ){
      Velocity[i*3  ] += Acceleration[i*3  ]*DT;
      Velocity[i*3+1] += Acceleration[i*3+1]*DT;
      Velocity[i*3+2] += Acceleration[i*3+2]*DT;

      Position[i*3  ] += Velocity[i*3  ]*DT;
      Position[i*3+1] += Velocity[i*3+1]*DT;
      Position[i*3+2] += Velocity[i*3+2]*DT;

     // initialize acceleration
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;
    }
  }

}

//move both upper wall and lower wall
void moveParticle_couette( int NumberOfParticles, simulation_parameters *parameters, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist){

// wall speed was setted in "setWALLSPEED()"
// reset wall velocity because of error of ICCG ( in case of implicit calculation )
// the wall velocity direction is x-axis

  for(int i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Velocity[i*3  ] += Acceleration[i*3  ]*parameters->DT;
      Velocity[i*3+1] += Acceleration[i*3+1]*parameters->DT;
      Velocity[i*3+2] += Acceleration[i*3+2]*parameters->DT;

      Position[i*3  ] += Velocity[i*3  ]*parameters->DT;
      Position[i*3+1] += Velocity[i*3+1]*parameters->DT;
      Position[i*3+2] += Velocity[i*3+2]*parameters->DT;

      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;

    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] > parameters->ParticleDistance ){ //upper wall
      Velocity[i*3] = parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] < parameters->ParticleDistance ){ //lower wall
      Velocity[i*3] = -parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }
  }
}


//move both upper wall and lower wall in reverse direction
void moveParticle_couette_reverse( int NumberOfParticles, simulation_parameters *parameters, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist){

// wall speed was setted in "setWALLSPEED()"
// reset wall velocity because of error of ICCG ( in case of implicit calculation )
// the wall velocity direction is x-axis

  for(int i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Velocity[i*3  ] += Acceleration[i*3  ]*parameters->DT;
      Velocity[i*3+1] += Acceleration[i*3+1]*parameters->DT;
      Velocity[i*3+2] += Acceleration[i*3+2]*parameters->DT;

      Position[i*3  ] += Velocity[i*3  ]*parameters->DT;
      Position[i*3+1] += Velocity[i*3+1]*parameters->DT;
      Position[i*3+2] += Velocity[i*3+2]*parameters->DT;

      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;

    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] > parameters->ParticleDistance ){ //upper wall
      Velocity[i*3] = -parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] < parameters->ParticleDistance ){ //lower wall
      Velocity[i*3] = parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }
  }
}


/*********************************************
 * This function is for collision processing *
 *********************************************/


void collision(int NumberOfParticles, int *ParticleType, double *Position, double *Velocity,int *BucketFirst, int *Nextof,
                 simulation_parameters *parameters, linklist_constant *linklist){
  int    i,j;
  double xij=0.0, yij=0.0, zij=0.0; //the distance between two particles in each axis
  double distance,distance2;
  double forceDT; /* forceDT is the impulse of collision between particles */
  double mi, mj;
  double velocity_ix, velocity_iy, velocity_iz;
  const double e = 0.2; //#define COEFFICIENT_OF_RESTITUTION
  const double dt = parameters -> DT;
  const double collisionDistance2 = 0.5*(parameters -> ParticleDistance) * 0.5*(parameters -> ParticleDistance);
  double *VelocityAfterCollision;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //index for jth roop
  int ib,jb ;        //index of bucket

  VelocityAfterCollision = (double *)malloc(3 * NumberOfParticles * sizeof(double));

  for(i=0;i<3*NumberOfParticles;i++){
    VelocityAfterCollision[i] = Velocity[i];
  }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,mi,velocity_ix,velocity_iy,velocity_iz,b_jz,b_jy,b_jx,jb,j,mj,xij,yij,zij,distance2,distance,forceDT)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divide to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
        mi = parameters -> Rho[ParticleType[i]];
        velocity_ix = Velocity[i*3  ]; velocity_iy = Velocity[i*3+1]; velocity_iz = Velocity[i*3+2];
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if((j!=i) && (ParticleType[j]!=GHOST)){
              mj = parameters -> Rho[ParticleType[j]];
              xij = Position[j*3  ] - Position[i*3  ];
              yij = Position[j*3+1] - Position[i*3+1];
              zij = Position[j*3+2] - Position[i*3+2];
              apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist); //PBC
              distance2 = (xij*xij) + (yij*yij) + (zij*zij); //use xij etc after this
              if(distance2<collisionDistance2){
                distance = sqrt(distance2);
                forceDT = (velocity_ix-Velocity[j*3  ])*(xij/distance)
                        +(velocity_iy-Velocity[j*3+1])*(yij/distance)
                        +(velocity_iz-Velocity[j*3+2])*(zij/distance);
                if(forceDT > 0.0){
                  forceDT *= (1.0+e)*mi*mj/(mi+mj);
                  velocity_ix -= (forceDT/mi)*(xij/distance);
                  velocity_iy -= (forceDT/mi)*(yij/distance);
                  velocity_iz -= (forceDT/mi)*(zij/distance);
                  //if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
                }
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        VelocityAfterCollision[i*3  ] = velocity_ix;
        VelocityAfterCollision[i*3+1] = velocity_iy;
        VelocityAfterCollision[i*3+2] = velocity_iz;
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

  for(i=0;i<NumberOfParticles;i++){
  if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Position[i*3  ] += (VelocityAfterCollision[i*3  ]-Velocity[i*3  ])*dt;
      Position[i*3+1] += (VelocityAfterCollision[i*3+1]-Velocity[i*3+1])*dt;
      Position[i*3+2] += (VelocityAfterCollision[i*3+2]-Velocity[i*3+2])*dt;
      Velocity[i*3  ] = VelocityAfterCollision[i*3  ];
      Velocity[i*3+1] = VelocityAfterCollision[i*3+1];
      Velocity[i*3+2] = VelocityAfterCollision[i*3+2];
  }}

  free(VelocityAfterCollision);
}

/***************************************************
 * calc Minimum pressure to calc pressure gradient *
 ***************************************************/

void setMinimumPressure(int NumberOfParticles, double Re_forGradient,int *ParticleType, double *MinimumPressure,
                         double *Pressure, double *Position, int *BucketFirst, int *Nextof, linklist_constant *linklist){
  double distance;
  int i,j;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,b_jz,b_jy,b_jx,jb,j,distance)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divided to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(ParticleType[i]!=GHOST && ParticleType[i]!=DUMMY_WALL){
        MinimumPressure[i]=Pressure[i];
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if((j!=i) && (ParticleType[j]!=GHOST) && (ParticleType[j]!=DUMMY_WALL)){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              if( (distance < Re_forGradient) && (MinimumPressure[i] > Pressure[j]) ){
                MinimumPressure[i] = Pressure[j];
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

}
/**********************************************
 * calculate pressure gradient                *
 **********************************************/

void calPressureGradient(int NumberOfParticles,int *ParticleType, double *Position,double *Acceleration,double *Pressure, double *MinimumPressure,
                         simulation_parameters *parameters, Gradient_constant *constant, int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int    i,j;
  double gradient_x, gradient_y, gradient_z;
  double xij, yij, zij,distance, distance2;
  double a,w,pij;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket
  double Re_forGradient = constant->Re_forGradient;
  int type_i;

  a =(parameters->DIM) / (constant->N0_forGradient);

// dynamic
// if change .c --> .cpp you can define some variable in roop
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,type_i,gradient_x,gradient_y,gradient_z,b_jz,b_jy,b_jx,jb,j,xij,yij,zij,distance2,distance,w,pij)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divide to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      type_i = ParticleType[i];
      if(type_i == FLUID || type_i == SOLID){
        gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if((j!=i) && (ParticleType[j]!=GHOST) && (ParticleType[j]!=DUMMY_WALL)){
              xij = Position[j*3  ] - Position[i*3  ];
              yij = Position[j*3+1] - Position[i*3+1];
              zij = Position[j*3+2] - Position[i*3+2];
              apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);
              distance2 = (xij*xij) + (yij*yij) + (zij*zij);
              distance = sqrt(distance2);
              if(distance<Re_forGradient){
                w =  weight(distance, Re_forGradient);
                pij = (Pressure[j] - MinimumPressure[i])/distance2;
                gradient_x += xij*pij*w;
                gradient_y += yij*pij*w;
                gradient_z += zij*pij*w;
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        Acceleration[i*3  ]= (-1.0)*a*gradient_x/(parameters->Rho[type_i]);
        Acceleration[i*3+1]= (-1.0)*a*gradient_y/(parameters->Rho[type_i]);
        Acceleration[i*3+2]= (-1.0)*a*gradient_z/(parameters->Rho[type_i]);
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

}


// calculate stress version -----------------------------------------------------------------------
void calPressureGradient_stress(int NumberOfParticles,int *ParticleType, double *Position,double *Acceleration,double *Pressure, double *MinimumPressure,
                         simulation_parameters *parameters, Gradient_constant *constant, int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int    i,j;
  double gradient_x, gradient_y, gradient_z;
  double xij, yij, zij,distance, distance2;
  double a,w,pij,a2;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket
  double Re_forGradient = constant->Re_forGradient;
  int type_i;

  a =(parameters->DIM) / (constant->N0_forGradient);

// dynamic
// if change .c --> .cpp you can define some variable in roop
#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,type_i,gradient_x,gradient_y,gradient_z,b_jz,b_jy,b_jx,jb,j,xij,yij,zij,distance2,distance,w,pij,a2)
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){ // ix is divide to multi threads
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
#else
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
#endif
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      type_i = ParticleType[i];
      if(type_i == FLUID || type_i == SOLID){
        gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if((j!=i) && (ParticleType[j]!=GHOST) && (ParticleType[j]!=DUMMY_WALL)){
              xij = Position[j*3  ] - Position[i*3  ];
              yij = Position[j*3+1] - Position[i*3+1];
              zij = Position[j*3+2] - Position[i*3+2];
              apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);
              distance2 = (xij*xij) + (yij*yij) + (zij*zij);
              distance = sqrt(distance2);
              if(distance<Re_forGradient){
                w =  weight(distance, Re_forGradient);
                pij = (Pressure[j] - MinimumPressure[i])/distance2;
                gradient_x += xij*pij*w;
                gradient_y += yij*pij*w;
                gradient_z += zij*pij*w;

                if (ParticleType[j] == WALL){ // stress for wall particle
                  a2 = w * pij * a;
                  Acceleration[j*3  ] += a2*xij; // - \rho * [ - 1 / rho * nabla P]
                  Acceleration[j*3+1] += a2*yij;
                  Acceleration[j*3+2] += a2*zij;

                }
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        Acceleration[i*3  ]= (-1.0)*a*gradient_x/(parameters->Rho[type_i]);
        Acceleration[i*3+1]= (-1.0)*a*gradient_y/(parameters->Rho[type_i]);
        Acceleration[i*3+2]= (-1.0)*a*gradient_z/(parameters->Rho[type_i]);
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

}

//-------------------------------------------------------------------

/************************************************
 *  update the position using pressure gradient
 ************************************************/

void moveParticleUsingPressureGradient_andPBC(int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist){
  int i;

// #pragma omp parallel for
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Velocity[i*3  ] +=Acceleration[i*3  ]*DT;
      Velocity[i*3+1] +=Acceleration[i*3+1]*DT;
      Velocity[i*3+2] +=Acceleration[i*3+2]*DT;

      Position[i*3  ] +=Acceleration[i*3  ]*DT*DT;
      Position[i*3+1] +=Acceleration[i*3+1]*DT*DT;
      Position[i*3+2] +=Acceleration[i*3+2]*DT*DT;

      //x
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}

      //y
      if(Position[i*3+1]>=linklist->MAX_Y)      {Position[i*3+1]-=linklist->Ly;}
      else if(Position[i*3+1]< linklist->MIN_Y){Position[i*3+1]+=linklist->Ly;}

      //z
      if(Position[i*3+2]>=linklist->MAX_Z)      {Position[i*3+2]-=linklist->Lz;}
      else if(Position[i*3+2]< linklist->MIN_Z){Position[i*3+2]+=linklist->Lz;}
    }

    // Acceleration  for all particles are initialized here
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  }
}
// ----------------------------------------------------------------------------------------------------

// update wall position here
// Don't use it with moveParticle_couette
void moveParticle2_PBC_couette(int NumberOfParticles, double DT, int *ParticleType,double *Position, double *Velocity, double *Acceleration, linklist_constant *linklist, simulation_parameters *parameters){
  int i;

// #pragma omp parallel for
  for(i=0;i<NumberOfParticles;i++){

    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Velocity[i*3  ] +=Acceleration[i*3  ]*DT;
      Velocity[i*3+1] +=Acceleration[i*3+1]*DT;
      Velocity[i*3+2] +=Acceleration[i*3+2]*DT;

      Position[i*3  ] +=Acceleration[i*3  ]*DT*DT;
      Position[i*3+1] +=Acceleration[i*3+1]*DT*DT;
      Position[i*3+2] +=Acceleration[i*3+2]*DT*DT;

      //x
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}

      //y
      if(Position[i*3+1]>=linklist->MAX_Y)      {Position[i*3+1]-=linklist->Ly;}
      else if(Position[i*3+1]< linklist->MIN_Y){Position[i*3+1]+=linklist->Ly;}

      //z
      if(Position[i*3+2]>=linklist->MAX_Z)      {Position[i*3+2]-=linklist->Lz;}
      else if(Position[i*3+2]< linklist->MIN_Z){Position[i*3+2]+=linklist->Lz;}

    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] > parameters->ParticleDistance ){ //upper wall
      Velocity[i*3] = parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] < parameters->ParticleDistance ){ //lower wall
      Velocity[i*3] = -parameters->wall_speed;
      Position[i*3] += Velocity[i*3] *parameters->DT;
      if(Position[i*3]>=linklist->MAX_X)      {Position[i*3]-=linklist->Lx;}
      else if(Position[i*3]< linklist->MIN_X){Position[i*3]+=linklist->Lx;}
    }

      // Acceleration for all particle are initialized here
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;


  }
}

//--------------------------------------------------------------------

/*****************************************************
 * calculate Courant Number (update DT every time )
 *****************************************************/

double calDT(int NumberOfParticles, double *Velocity, int *ParticleType, simulation_parameters *parameters){
  const double Cmax = 0.2; //max value of Courant number
  const double dmax = 0.2; //max value of diffusion number
  double dt;
  double MaxVelocity = parameters->wall_speed;
  int MaxIndex = -1 ;
  double tempVelocity= 0.0;
  double dt_diff;                        //consider diffusion number limit
  double dt_s = PI/sqrt(parameters->ks); // 0.5T = pi * sqrt(m/k)

  for(int i=0; i<NumberOfParticles; i++){
    if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      tempVelocity = Velocity[3*i]*Velocity[3*i] + Velocity[3*i+1]*Velocity[3*i+1] + Velocity[3*i+2]*Velocity[3*i+2];
      if(tempVelocity > MaxVelocity){
        MaxVelocity = tempVelocity;
        MaxIndex = i;
      }
    }
  }
  MaxVelocity = sqrt(MaxVelocity);
  // Error handling
  // Dt is defined by Dt * Maxspeed = 0.2 * l_0
  // the  simulation size is l_0 * 10 ~ 1000, then Maxspeed > wallspeed * 100 ~ 1000 --> calcualtion will be destroyed
  if( parameters->wall_speed && MaxVelocity > parameters->wall_speed * 10000 ){
    fprintf(stderr,"max speed is 10000 times larger than wall speed \n");
    fprintf(stderr,"MaxIndex : %d , MaxVelocity : %lf, (vx: %lf vy: %lf, vz: %lf), PartileType : %d ( 0 -> FLUID , 1-> SOLID )\n",MaxIndex,MaxVelocity,Velocity[3*MaxIndex],Velocity[3*MaxIndex+1],Velocity[3*MaxIndex+2],ParticleType[MaxIndex]);
    exit(2);
  }

  // let Cmax = 0.2
  if(MaxVelocity != 0.0){
    dt = Cmax * (parameters->ParticleDistance) / MaxVelocity;
  }else{
    dt = parameters->DT;
  }

  // let d_i = 0.2
  dt_diff = dmax * (parameters->ParticleDistance) * (parameters->ParticleDistance) / (parameters->Re_inv[FLUID]);
  if (dt > dt_diff){dt = dt_diff;}
  if (dt > dt_s   ){dt = dt_s;}
  return dt;
}


/***********************
 * for couette wall
 ********************/

void defWallSpeed(simulation_parameters *parameters, int ny, double ParticleDistance){

  double width = (ny - 5) * ParticleDistance ;
  parameters->wall_speed = 0.5 * width * (parameters->strain_rate) ;
  // parameters->wall_speed = 0.0;

}


void setWallSpeed(int np, double wall_speed, double ParticleDistance,double *Position, int *ParticleType, double *Velocity){
  int i;
  for( i=0; i<np; i++){
    //all wall
    if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] > ParticleDistance ){ //upper wall
      Velocity[i*3] = wall_speed;
    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] < ParticleDistance ){ //lower wall
      Velocity[i*3] = -(1.0) * wall_speed;
    }
  }
}

void setWallSpeed_reverse(int np, double wall_speed, double ParticleDistance,double *Position, int *ParticleType, double *Velocity){
  int i;
  for( i=0; i<np; i++){
    //all wall
    if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] > ParticleDistance ){ //upper wall
      Velocity[i*3] = -(1.0) * wall_speed;
    }else if ( (ParticleType[i] == WALL || ParticleType[i] == DUMMY_WALL) && Position[i*3+1] < ParticleDistance ){ //lower wall
      Velocity[i*3] = wall_speed;
    }
  }
}

void setCouetteSpeed(int NumberOfParticle, int ny, double *Position, int *ParticleType, double *Velocity, simulation_parameters *parameters){

  int i;
  double width_half = 0.5 * (ny - 5) * parameters->ParticleDistance ;

  for( i=0; i<NumberOfParticle; i++){
    if (ParticleType[i] == FLUID || ParticleType[i] == SOLID){
      Velocity[i*3+1] = (Position[3*i+1] - width_half) * parameters->strain_rate ;
    }
  }
}


/***********************************
 * calculate wall stress
 ***********************************/

// stress is calculated in
void output_WallStress2D(int NumberOfParticles, double Time, int *ParticleType, double *Position, double *Acceleration, simulation_parameters *parameters){

  int i ;
  double stress_u[2] = {0.0,0.0};
  double stress_d[2] = {0.0,0.0};
  FILE *fp;
  char outputfile[256];
  int Nu=0, Nd=0; // number of wall particle

  for(i=0; i<NumberOfParticles; i++){
    if(ParticleType[i] == WALL){
      if (Position[i*3+1] > parameters->ParticleDistance){ // upper wall
          stress_u[0] += Acceleration[3*i];
          stress_u[1] += Acceleration[3*i+1];
          Nu++;
      }else if (Position[i*3+1] < parameters->ParticleDistance){ // lower wall
          stress_d[0] += Acceleration[3*i];
          stress_d[1] += Acceleration[3*i+1];
          Nd++;
      }
    }
  }

  sprintf(outputfile, "%s/output_wallstress.dat",parameters->output_dir);
  fp = fopen(outputfile,"a");

  fprintf(fp,"%lf  ",Time);
  fprintf(fp,"%lf  %lf  ",stress_u[0]/Nu, stress_u[1]/Nu);
  fprintf(fp,"%lf  %lf\n",stress_d[0]/Nd, stress_d[1]/Nd);

  fclose(fp);

}

void make_output_stress(simulation_parameters *parameters){

  FILE *fp;
  char outputfile[1024];

  sprintf(outputfile, "%s/output_wallstress.dat",parameters->output_dir);
  fp = fopen(outputfile,"w");

  fprintf(fp,"Time  ");

  if (parameters->DIM == 2){
    fprintf(fp,"s_u_x  s_u_y  ");
    fprintf(fp,"s_d_x  s_d_y\n");
  }else{
    fprintf(fp,"s_u_x  s_u_y  s_u_z  ");
    fprintf(fp,"s_d_x  s_d_y  s_d_z\n");
  }

  fclose(fp);

}



/**********************************************
 * check the number of neighbor particles      *
 ***********************************************/
//PBC Done 9/24
/*
void Check_Neighbor_Particles(int NumberOfParticles, double Re_forLaplacian, double *Position, int *ParticleType,
                              int *BucketFirst, int *Nextof,linklist_constant *linklist,int Timestep){

  int i,j;
  int *NumOfNeighbor;
  int b_ix,b_iy,b_iz,ib;
  int b_jx,b_jy,b_jz,jb;
  double distance;
  int Max_count = 0;
  FILE *fp;

  fp = fopen("NumberOfNeighborList.txt","a");
  if(Timestep == 0){
    fprintf(fp,"-----------------------------------------------------------");
    fprintf(fp,"Number of neighbor particles \n");
    fprintf(fp,"Number of particles : %d \n",NumberOfParticles);
  }

  NumOfNeighbor = (int *)malloc(NumberOfParticles * sizeof(int));

  //cell list
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      NumOfNeighbor[i] = 0;
      if(ParticleType[i] == FLUID || ParticleType[i] == SOLID || ParticleType[i] == WALL){
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          //jb = b_jz*(linklist->nBxy) + b_jy*(linklist->nBx) + b_jx; //the index of j th bucket
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if(i != j && (ParticleType[j] == FLUID || ParticleType[i] == SOLID || ParticleType[j] == WALL)){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance < Re_forLaplacian) NumOfNeighbor[i] += 1;
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

  for(i=0; i<NumberOfParticles; i++){
    if(Max_count < NumOfNeighbor[i]){
      Max_count = NumOfNeighbor[i];
    }
  }
  fprintf(fp,"%d : %d\n",Timestep,Max_count);

  fclose(fp);

}*/
