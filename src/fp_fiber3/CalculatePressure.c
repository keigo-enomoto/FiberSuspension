
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <petscksp.h>


#include "TypeDef.h"
#include "petsc_solver.h"

/*****for boundary condition**********/
enum{
     GHOST_OR_DUMMY = -1,
     SURFACE_PARTICLE = 1,
     INNER_PARTICLE = 0
};

/*
// for FlagForCheckingBoundaryCondition (exceptional process) --> for free surface flow
enum{
    DIRICHLET_BOUNDARY_IS_NOT_CONNECTED = 2,
    DIRICHLET_BOUNDARY_IS_CONNECTED = 3,
    DIRICHLET_BOUNDARY_IS_CHECKED = 4
};
*/

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
/*
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

//-------------------------------------------------------------------------------------------------------------

/***************************
 * calculate number density
 ***************************/
static void calNumberDensity(int NumberOfParticles, double Re_forNumberDensity, int *ParticleType, double *NumberDensity, double *Position,
                              int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int    i,j;
  double distance;
  double w;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket

#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(ib,i,b_jz,b_jy,b_jx,jb,j,distance,w)
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
      NumberDensity[i] = 0.0;
      if(ParticleType[i] != GHOST){
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( (j!=i) && (ParticleType[j]!=GHOST) ){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              w =  weight(distance, Re_forNumberDensity);
              NumberDensity[i] += w;
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

//--------------------------------------------------------------------------------------------------
/**************************
 * set Boundary condition *
 **************************/

/*
static void setBoundaryCondition( int NumberOfParticles, double N0_forNumberDensity, int *ParticleType, double *NumberDensity, int *BoundaryCondition){
  int i;
  const double beta = 0.97; // = THRESHOLD_RATIO_OF_NUMBER_DENSITY

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ){ //Don't calculate Pressure
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if( NumberDensity[i] < beta * N0_forNumberDensity ){ //fix the pressure to 0
      BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      BoundaryCondition[i]=INNER_PARTICLE; //calculate pressure (FLUID and WALL)
    }
  }
}
*/

static int setBoundaryCondition_returnNumberOfInnerPartiles( int NumberOfParticles, double N0_forNumberDensity,
              int *ParticleType, double *NumberDensity, int *BoundaryCondition, int *order_inner, int *index_inner){
  int i;
  // const double beta = 0.97; // = THRESHOLD_RATIO_OF_NUMBER_DENSITY
  int n_inner = 0;

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ){ //Don't calculate Pressure
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    // }else if( NumberDensity[i] < beta * N0_forNumberDensity ){ //fix the pressure to 0
    //   BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      BoundaryCondition[i]=INNER_PARTICLE; //calculate pressure (FLUID or SOLID or WALL)
      order_inner[i] = n_inner;
      index_inner[n_inner] = i;
      n_inner ++;
    }
  }

  return n_inner;
}

//------------------------------------------------------------------------------------------------------

/**************************************
 * initialize 5 array for inner matrix
 **************************************/
static void initialize_inner_matrix(int n_inner, int N_neighbor, double *diag,double *val, int *col, int *row){
  int i;
  for(i=0;i<n_inner;i++){
    diag[i]= 0.0;
    val[i] = 0.0;
    col[i] = -1;
    row[i] = -1;
    // nnz[i] = -1;
  }
  for(i=n_inner;i<n_inner*N_neighbor;i++){
    val[i] = 0.0;
    col[i] = -1;
    row[i] = -1;
  }

}
//--------------------------------------------------------------------------------------------------------

/**********************
 * set source term (b)
 **********************/
/*
static void setSourceTerm(double N0_forNumberDensity, int NumberOfParticles, double DT,
                   double *SourceTerm,int *ParticleType, int *BoundaryCondition, double *NumberDensity){
  int i;
  const double gamma = 0.2; //RELAXATION_COEFFICIENT_FOR_PRESSURE generally 0.2 (Labnote p123)
                      // for stable calculation. gamma = 1.0 -> completely match the theory equation

  for(i=0;i<NumberOfParticles;i++){
    SourceTerm[i]=0.0;
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ) continue;
    if(BoundaryCondition[i]==INNER_PARTICLE){
      SourceTerm[i] = gamma * (1.0/(DT*DT))*((NumberDensity[i]-N0_forNumberDensity)/N0_forNumberDensity);
    }else if(BoundaryCondition[i]==SURFACE_PARTICLE){
      SourceTerm[i]=0.0;
    }
  }
}
*/

// return array of Source Term
static double *makeSourceTerm(double N0_forNumberDensity, int n_inner, int *index_inner,double DT,
                              int *ParticleType, int *BoundaryCondition, double *NumberDensity,int NumberOfParticles){
  int i_inner;  //index of inner matrix
  int i_origin; //index of original particle label
  double *SourceTerm;
  const double gamma = 0.2; //RELAXATION_COEFFICIENT_FOR_PRESSURE generally 0.2 (Labnote vol.1 p123)
                            // for stable calculation. gamma = 1.0 -> completely match the theory equation

  SourceTerm = (double *)malloc(NumberOfParticles * sizeof(double));
  for(i_inner=0;i_inner<n_inner;i_inner++){
    i_origin = index_inner[i_inner];
    SourceTerm[i_inner]=0.0;
    if(ParticleType[i_origin]==GHOST || ParticleType[i_origin]== DUMMY_WALL ) continue;
    if(BoundaryCondition[i_origin]==INNER_PARTICLE){
      SourceTerm[i_inner] = gamma * (1.0/(DT*DT))*((NumberDensity[i_origin]-N0_forNumberDensity)/N0_forNumberDensity);
    }else if(BoundaryCondition[i_origin]==SURFACE_PARTICLE){
      SourceTerm[i_inner]=0.0;
    }
  }
  // ignore range
  for(i_inner=n_inner; i_inner<NumberOfParticles; i_inner++){
    SourceTerm[i_inner] = 0.0;
  }
  return SourceTerm;
}

//--------------------------------------------------------------------------------------------------------
/*****************************************************
  *These four function is to set Matrix              *
 *****************************************************/
/*
//PBC Done 9/24
static void checkBoundaryCondition(int NumberOfParticles,int *ParticleType,double *Position, double Re_forLaplacian,
                            int *BoundaryCondition,int *FlagForCheckingBoundaryCondition,int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int i,j,count;
  double distance;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket

  for(i=0;i<NumberOfParticles;i++){
    if (BoundaryCondition[i]==GHOST_OR_DUMMY){
      FlagForCheckingBoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if (BoundaryCondition[i]==SURFACE_PARTICLE){
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else if (BoundaryCondition[i]==INNER_PARTICLE){
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED; //Inner Particle
    }else{
      printf("The boundary condition of %d th particle is invalid.\n",i);
    }
  }

  do{
    count = 0;
    for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
    for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
    for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
      ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
      i = BucketFirst[ib];
      if(i == -1) continue;
      for(;;){
        if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_CONNECTED){
          for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
          for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
          for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
            //jb = b_jz*(linklist->nBxy) + b_jy*(linklist->nBx) + b_jx; //the index of j th bucket
            jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
            j = BucketFirst[jb];
            if(j == -1) continue;
            for(;;){ //beginning of neighbor list of ith particle
              if(FlagForCheckingBoundaryCondition[j]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
              if((j!=i) && (ParticleType[j]!=GHOST) && (ParticleType[j]!= DUMMY_WALL)){
                distance = calDistance_ij_pbc(i,j,Position,linklist);
                if(distance<Re_forLaplacian){
                FlagForCheckingBoundaryCondition[j]=DIRICHLET_BOUNDARY_IS_CONNECTED;
                }
              }}
              j = Nextof[j];
              if(j==-1) break;
            }//end of neighbor list of ith particle
          }}}
          FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CHECKED;
          count++;
        }
        i = Nextof[i];
        if(i == -1) break;
      }
    }}}
  }while (count!=0);// This rootin will be continued until all particles become "DIRICHLET_BOUNDARY_IS_CONNECTED".

*/
  /*
  for(i=0;i<NumberOfParticles;i++){
    if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }*/
/*
}

//-------------------------------------------------------------------------------------------------------------------

//for petsc
static void increaseDiagonalTerm(int n_inner, int *FlagForCheckingBoundaryCondition, double *diag, int *index_inner){
  int i;
  int i_origin;//original label of i th inner particle

  for(i=0;i<n_inner;i++) {
    i_origin = index_inner[i];
    if(FlagForCheckingBoundaryCondition[i_origin] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      diag[i] = 2.0 * diag[i];
    }
  }
}
//---------------------------------------------------------------------------------------------------------------------


static void exceptionalProcessingForBoundaryCondition(int NumberOfParticles,int *ParticleType,double *Position, double Re_forLaplacian,
                                               int *BoundaryCondition,int *FlagForCheckingBoundaryCondition, int *BucketFirst, int *Nextof,
                                               linklist_constant *linklist, int n_inner, double *diag, int *index_inner){
  // If tere is no Dirichlet boundary condition on the fluid,
  //increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions.
  checkBoundaryCondition(NumberOfParticles,ParticleType,Position,Re_forLaplacian,BoundaryCondition,FlagForCheckingBoundaryCondition,BucketFirst,Nextof,linklist);
  increaseDiagonalTerm(n_inner,FlagForCheckingBoundaryCondition,diag,index_inner);
}
*/
//--------------------------------------------------------------------------------------------------------------------------
/*********************
 * set matrix        *
 *********************/

// labnote vol.2 p148 --------------------------------------------------------------------
static int setMatrix_returnNonzeroNumber(double N0_forLaplacian,int NumberOfParticles,double Lambda,double Re_forLaplacian,
               int *ParticleType, int *FlagForCheckingBoundaryCondition,int *BoundaryCondition,double *Position,simulation_parameters *parameters,
               int *BucketFirst, int *Nextof, linklist_constant *linklist,double *val, int *col, int *row, double *diag, int *order_inner,int *index_inner, int n_inner){
  const double COMPRESSIBILITY = 0.45E-9 ; // [1/Pa]

  double distance;
  double coefficientIJ,a;
  int    i,j;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;           //index of bucket
  int i_in, j_in;//index for matrix of inner particles (k * k)
  int index_coo; //index for COO matrix(val,col,row) <- don't match to index_inner[]
  // int pre_index_coo = 0;//for nnz array

  index_coo = 0;
  a = 2.0*(parameters->DIM)/(N0_forLaplacian*Lambda);

// use index_coo --> be careful !!
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(BoundaryCondition[i] == INNER_PARTICLE){ // ------------------------------------------------
        i_in = order_inner[i];      //set index for inner matrix
        // pre_index_coo = index_coo;
        for(b_jz=b_iz;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          if( (b_jx==b_ix-1) && (b_jz==b_iz) ) continue;
          if( (b_jx==b_ix) && (b_jy==b_iy-1) && (b_jz==b_iz) ) continue;
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( (j!=i) && (BoundaryCondition[j]!=GHOST_OR_DUMMY) && ( ib!=jb || (ib==jb && i < j))){
              distance  = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<Re_forLaplacian){
                coefficientIJ = a * weight(distance, Re_forLaplacian) ;
                if(BoundaryCondition[j] == INNER_PARTICLE){
                  j_in = order_inner[j];    //set index for inner matrix
                  // a_ij
                  val[index_coo] = -(1.0) * coefficientIJ / (parameters->Rho[ParticleType[i]]);
                  col[index_coo] = j_in;
                  row[index_coo] = i_in; //insert value and position in COO matrix
                  index_coo++;  //update index for COO matrix
                  // a_ji
                  val[index_coo] = -(1.0) * coefficientIJ / (parameters->Rho[ParticleType[j]]);
                  col[index_coo] = i_in;
                  row[index_coo] = j_in;
                  index_coo++;
                  diag[j_in] += coefficientIJ / (parameters->Rho[ParticleType[j]]);
                }
                diag[i_in] += coefficientIJ / (parameters->Rho[ParticleType[i]]);
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        diag[i_in] += (COMPRESSIBILITY)/((parameters->DT)*(parameters->DT)); //allow compressibility for stability
        // nnz[i_in] = index_coo - pre_index_coo;
      }else if (BoundaryCondition[i] == SURFACE_PARTICLE){ // --------------------------------------------------------
        i_in = order_inner[i];              // set index for inner matrix
        for(b_jz=b_iz;b_jz<=b_iz+1;b_jz++){ // roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          if( (b_jx==b_ix-1) && (b_jz==b_iz) ) continue;
          if( (b_jx==b_ix) && (b_jy==b_iy-1) && (b_jz==b_iz) ) continue;
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( (j!=i) && (BoundaryCondition[j]!=GHOST_OR_DUMMY) && ( ib!=jb || (ib==jb && i < j))){
              distance  = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<Re_forLaplacian){
                coefficientIJ = a * weight(distance, Re_forLaplacian) ;
                if(BoundaryCondition[j] == INNER_PARTICLE){
                  j_in = order_inner[j];    //set index for inner matrix
                  diag[j_in] += coefficientIJ / (parameters->Rho[ParticleType[j]]);
                }
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

  //exceptionalProcessingForBoundaryCondition(NumberOfParticles,ParticleType,Position,Re_forLaplacian,BoundaryCondition,FlagForCheckingBoundaryCondition,BucketFirst,Nextof,linklist,n_inner,diag,index_inner);

  return index_coo;
}

/*
// for petsc --------------------------------------------------------------------
// previous cell list ----------------
static int setMatrix_returnNonzeroNumber(double N0_forLaplacian,int NumberOfParticles,double Lambda,double Re_forLaplacian,
               int *ParticleType, int *FlagForCheckingBoundaryCondition,int *BoundaryCondition,double *Position,simulation_parameters *parameters,
               int *BucketFirst, int *Nextof, linklist_constant *linklist,double *val, int *col, int *row, double *diag, int *order_inner,int *index_inner, int n_inner){
  const double COMPRESSIBILITY = 0.45E-9 ; // [1/Pa]

  double distance;
  double coefficientIJ,a;
  int    i,j;
  int b_ix, b_iy, b_iz; //the coordinate of bucket
  int b_jx, b_jy, b_jz; //the neighbor bucket's coordinate
  int ib,jb ;           //index of bucket
  int i_in, j_in;//index for matrix of inner particles (k * k)
  int index_coo; //index for COO matrix(val,col,row) <- don't match to index_inner[]
  // int pre_index_coo = 0;//for nnz array

  index_coo = 0;
  a = 2.0*(parameters->DIM)/(N0_forLaplacian*Lambda);
  for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(BoundaryCondition[i] == INNER_PARTICLE){
        i_in = order_inner[i];      //set index for inner matrix
        // pre_index_coo = index_coo;
        for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          //jb = b_jz*(linklist->nBxy) + b_jy*(linklist->nBx) + b_jx; //the index of j th bucket
          jb = calBucketIndex_pbc(b_jx,b_jy,b_jz,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( (j!=i) && (BoundaryCondition[j]!=GHOST_OR_DUMMY) ){
              distance  = calDistance_ij_pbc(i,j,Position,linklist);
              if(distance<Re_forLaplacian){
                coefficientIJ = a * weight(distance, Re_forLaplacian) / (parameters->Rho[ParticleType[i]]);
                if(BoundaryCondition[j] == INNER_PARTICLE){
                  j_in = order_inner[j];    //set index for inner matrix
                  val[index_coo] = -(1.0) * coefficientIJ;
                  col[index_coo] = j_in;
                  row[index_coo] = i_in; //insert value and position in COO matrix
                  index_coo++;  //update index for COO matrix
                }
                diag[i_in] += coefficientIJ;
              }
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }}}
        diag[i_in] += (COMPRESSIBILITY)/((parameters->DT)*(parameters->DT)); //allow compressibility for stability
        // nnz[i_in] = index_coo - pre_index_coo;
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }}}

  //exceptionalProcessingForBoundaryCondition(NumberOfParticles,ParticleType,Position,Re_forLaplacian,BoundaryCondition,FlagForCheckingBoundaryCondition,BucketFirst,Nextof,linklist,n_inner,diag,index_inner);

  return index_coo;
}
*/

//------------------------------------------------------------------------------------------------------------
/********************************************
 * solve linear equation  --> petsc_solver.c
 ********************************************/

//-------------------------------------------------------------------------
/**********************************
 * remove negative perssure       *
 **********************************/
/*
static void removeNegativePressure(int NumberOfParticles, double *Pressure){
  int i;

  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0)Pressure[i]=0.0;
  }
}*/

// when using LinearSolver()
static void removeNegativePressure_and_applyBoundaryCondition(int NumberOfParticles, double *Pressure, int *BoundaryCondition){
  int i;

// #pragma omp parallel for
  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0 || BoundaryCondition[i] == GHOST_OR_DUMMY || BoundaryCondition[i] == SURFACE_PARTICLE){
      Pressure[i]=0.0;
    }
  }
}

//-------------------------------------------
/*********************
 * useful function
 *********************/
/*
static void print_array_int(int *array,int n){
    for(int i=0; i<n; i++){
        printf("%d, ",array[i]);
    }
    printf("\n");
}

static void print_array_double(double *array,int n){
    for(int i=0; i<n; i++){
        printf("%lf, ",array[i]);
    }
    printf("\n");
}
*/

/*
static double return_two_norm(double *SourceTerm,int n){
  int i;
  double norm = 0.0;

  for( i=0; i<n; i++ ){
    norm += SourceTerm[i] * SourceTerm[i] ;
  }
  return sqrt(norm);
}
*/





//=========================================================================================================================
/*************************************
 * This is the main funciton
 *************************************/

void calPressure(petsc_object *petsc,int NumberOfParticles, double *Pressure,
                 int *ParticleType, double *NumberDensity, double *Position, int *BoundaryCondition, int *FlagForCheckingBoundaryCondition,
                 simulation_parameters *parameters, Gradient_constant *constant,int *BucketFirst, int *Nextof, linklist_constant *linklist){

  const double Re_forNumberDensity = constant->Re_forNumberDensity;
  const double Re_forLaplacian = constant->Re_forLaplacian;
  const double N0_forNumberDensity = constant -> N0_forNumberDensity;
  const double N0_forLaplacian = constant -> N0_forLaplacian;
  const double Lambda = constant -> Lambda;
  double *SourceTerm; //[n_inner]

  int n_inner = 1;    //number of inner particles
  int n_nonzero;      //number of nonzero value in inner matrix
  int *order_inner;   //An array that holds an index of each particle number in the inner matrix
  int *index_inner;   // array that holds index of inner particles
  double *val,*diag;  // Coefficient Matrix
  int *col, *row;     // Coefficient Matrix
  // int *nnz;           // for PETSc
  int n_neighbor;     // Number of neighbor particles

  if(parameters->DIM == 2){
    n_neighbor = 48;
  }else{
    n_neighbor = 144;
  }

  index_inner = (int *)malloc(NumberOfParticles * sizeof(int));
  order_inner = (int *)malloc(NumberOfParticles * sizeof(int));
  // initialize two array
  for(int i=0; i<NumberOfParticles; i++){
    order_inner[i] = -1;
    index_inner[i] = -1;
  }

//-----------------------------------------------------------------------------------------------------------------------
//set Boundary condition
  calNumberDensity(NumberOfParticles,Re_forNumberDensity,ParticleType,NumberDensity,Position,BucketFirst,Nextof,linklist);
  n_inner = setBoundaryCondition_returnNumberOfInnerPartiles(NumberOfParticles,N0_forNumberDensity,ParticleType,NumberDensity,BoundaryCondition,order_inner,index_inner);

//allocate coeffient matrix and source term
  diag =(double *)malloc(n_inner * sizeof(double));
  val = (double *)malloc(n_inner * n_neighbor * sizeof(double)); //only inner particles
  col = (int    *)malloc(n_inner * n_neighbor * sizeof(int));
  row = (int    *)malloc(n_inner * n_neighbor * sizeof(int));
  // nnz = (int    *)malloc(n_inner * sizeof(int));
  initialize_inner_matrix(n_inner,n_neighbor,diag,val,col,row); //initialize val,col,row,nnz

//calculate
  SourceTerm = makeSourceTerm(N0_forNumberDensity,n_inner,index_inner,parameters->DT,ParticleType,BoundaryCondition,NumberDensity,NumberOfParticles);
  n_nonzero = setMatrix_returnNonzeroNumber(N0_forLaplacian,NumberOfParticles,Lambda,Re_forLaplacian,ParticleType,FlagForCheckingBoundaryCondition,BoundaryCondition,Position,parameters,BucketFirst,Nextof,linklist,val,col,row,diag,order_inner,index_inner,n_inner);
  PETSc_Solver(petsc,Pressure,SourceTerm,NumberOfParticles,n_inner,n_nonzero,index_inner,order_inner,diag,val,col,row);
  removeNegativePressure_and_applyBoundaryCondition(NumberOfParticles,Pressure,BoundaryCondition);

  free(SourceTerm);
  free(order_inner);free(index_inner);
  free(val); free(col); free(row);
  free(diag);
  // free(nnz);
}

