#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TypeDef.h"
#include "petsc_solver.h"

/*******************************
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
/*
// apply pbc for xij , yij, and zij
static void apply_pbc_distance_each_axis(double *xij, double *yij, double *zij, linklist_constant *linklist){

  // apply pbc
  if(*xij < -linklist->Lx*0.50)  { *xij += linklist->Lx; }
  else if(*xij>linklist->Lx*0.50){ *xij -= linklist->Lx; }

  if(*yij < -linklist->Ly*0.50)  { *yij += linklist->Ly; }
  else if(*yij>linklist->Ly*0.50){ *yij -= linklist->Ly; }

  if(*zij < -linklist->Lz*0.50)  { *zij += linklist->Lz; }
  else if(*zij>linklist->Lz*0.50){ *zij -= linklist->Lz; }

}*/

// if using peliodic boundary, use this function
static int calBucketIndex_pbc(int b_jx, int b_jy, int b_jz, linklist_constant *linklist){
  int tb_jx, tb_jy, tb_jz; //the temporary neighbor bucket's coordinate

  tb_jx = (b_jx + (linklist->nBx)) % (linklist->nBx);
  tb_jy = (b_jy + (linklist->nBy)) % (linklist->nBy);
  tb_jz = (b_jz + (linklist->nBz)) % (linklist->nBz);

  return  tb_jz*(linklist->nBxy) + tb_jy*(linklist->nBx) + tb_jx;
}

//----------------------------------------------------------------------------
/**************************
 * initialize COO matrix
***************************/

static void initialize_matrix(int n, int N_neighbor, double *diag, double *val, int *col, int *row, double *u_x, double *u_y, double *u_z){
  int i;
  for(i=0;i<n;i++){
    u_x[i] = 0.0;
    u_y[i] = 0.0;
    u_z[i] = 0.0;
    diag[i]= 0.0;
    val[i] = 0.0;
    col[i] = -1;
    row[i] = -1;
  }
  for(i=n;i<n*N_neighbor;i++){
    val[i] = 0.0;
    col[i] = -1;
    row[i] = -1;
  }
}

//----------------------------------------------------------------------------
/***************
 * ser matrix
****************/

static void setSourceTerm(double *Velocity, double *b_x, double *b_y, double *b_z, int n){
    int i;
    for(i=0; i<n; i++){
        b_x[i] = Velocity[3*i];
        b_y[i] = Velocity[3*i+1];
        b_z[i] = Velocity[3*i+2];
    }
}

static void addGravity(int n, double DT, double *b_x, double *b_y, double *b_z, int *ParticleType,double gravity){
    int i;

    for(i=0; i<n; i++){
      if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
        b_x[i] += DT * gravity;
        //b_y[i] += DT * gravity;
      }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------
static int setMatrix_returnNonzeroNumber(int NumberOfParticles,int *ParticleType, double *Position,simulation_parameters *parameters, Gradient_constant *constant,
               int *BucketFirst, int *Nextof, linklist_constant *linklist,double *val, int *col, int *row, double *diag, double *b_x, double *b_y, double *b_z, double *Velocity){

  double distance;
  double coefficientIJ,a;
  int    i,j;
  int index_coo; //index for COO matrix(val,col,row)
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;           //index of bucket

  index_coo = 0;
  a = parameters->DT*2.0*(parameters->DIM)/(constant->N0_forLaplacian*constant->Lambda);
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          //jb = b_jz*(linklist->nBxy) + b_jy*(linklist->nBx) + b_jx; //the index of j th bucket
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( (j!=i) && (ParticleType[j]==FLUID || ParticleType[j] == SOLID) ){
              distance  = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<constant->Re_forLaplacian){
                coefficientIJ = a * weight(distance, constant->Re_forLaplacian) * parameters->Re_inv[ParticleType[i]];
                val[index_coo] = -(1.0) * coefficientIJ;
                col[index_coo] = j;
                row[index_coo] = i; //insert value and position in COO matrix
                index_coo++;  //update index for COO matrix
                diag[i] += coefficientIJ;
              }
            }else if( (j!=i) && (ParticleType[j]==WALL) ){ // add wall coefficient to right side b vector (if poiseuille unnecessary)
              distance  = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<constant->Re_forLaplacian){
                coefficientIJ = a * weight(distance, constant->Re_forLaplacian) * parameters->Re_inv[ParticleType[i]];
                b_x[i] += coefficientIJ * Velocity[3*j];
                b_y[i] += coefficientIJ * Velocity[3*j+1];
                b_z[i] += coefficientIJ * Velocity[3*j+2];
                diag[i] += coefficientIJ;
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
      }
      diag[i] += 1.0;
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

  return index_coo;

}

/*************************************
* solution vector --> Velocity array
**************************************/
static void putSolution(int NumberOfParticles, double *Velocity, double *u_x, double *u_y, double *u_z){
  int i ;

  for(i=0; i<NumberOfParticles; i++){
    Velocity[3*i  ] = u_x[i];
    Velocity[3*i+1] = u_y[i];
    Velocity[3*i+2] = u_z[i];
  }

}

//-----------------------------------------------------------------------

/******************
 * main parts
 ******************/

void calGravity_and_Viscosity_i(petsc_object *petsc,int NumberOfParticles, double *Velocity,int *ParticleType, double *Position,
                 simulation_parameters *parameters, Gradient_constant *constant,int *BucketFirst, int *Nextof, linklist_constant *linklist){

  int n_neighbor, n_nonzero;
  double *b_x,*b_y,*b_z ;
  double *u_x, *u_y, *u_z; //solution of matrix
  double *diag,*val;
  int *col,*row;

  if(parameters->DIM == 2){
      n_neighbor = 48;
  }else{
      n_neighbor = 144;
  }
  b_x = (double *)malloc(NumberOfParticles * sizeof(double));
  b_y = (double *)malloc(NumberOfParticles * sizeof(double));
  b_z = (double *)malloc(NumberOfParticles * sizeof(double));
  u_x = (double *)malloc(NumberOfParticles * sizeof(double));
  u_y = (double *)malloc(NumberOfParticles * sizeof(double));
  u_z = (double *)malloc(NumberOfParticles * sizeof(double));
  diag= (double *)malloc(NumberOfParticles * sizeof(double));
  val = (double *)malloc(NumberOfParticles * n_neighbor * sizeof(double));
  col = (int    *)malloc(NumberOfParticles * n_neighbor * sizeof(int));
  row = (int    *)malloc(NumberOfParticles * n_neighbor * sizeof(int));
  initialize_matrix(NumberOfParticles,n_neighbor,diag,val,col,row,u_x,u_y,u_z);

  setSourceTerm(Velocity,b_x,b_y,b_z,NumberOfParticles);
  if ( parameters->Gravity != 0.0 ){
    addGravity(NumberOfParticles,parameters->DT,b_x,b_y,b_z,ParticleType,parameters->Gravity);
  }
  n_nonzero = setMatrix_returnNonzeroNumber(NumberOfParticles,ParticleType,Position,parameters,constant,BucketFirst,Nextof,linklist,val,col,row,diag,b_x,b_y,b_z,Velocity);

  // PETSc_Solver_coo(petsc,u_x,b_x,NumberOfParticles,n_nonzero,diag,val,col,row);
  // PETSc_Solver_coo(petsc,u_y,b_y,NumberOfParticles,n_nonzero,diag,val,col,row);
  // if(parameters->DIM==3) {PETSc_Solver_coo(petsc,u_z,b_z,NumberOfParticles,n_nonzero,diag,val,col,row);}

  PETSc_Solver_coo_xyz(petsc,u_x,u_y,u_z,b_x,b_y,b_z,NumberOfParticles,n_nonzero,diag,val,col,row,parameters->DIM);
  putSolution(NumberOfParticles,Velocity,u_x,u_y,u_z);

  free(b_x); free(b_y); free(b_z); free(u_x); free(u_y); free(u_z);
  free(val); free(col); free(row); free(diag);
}
