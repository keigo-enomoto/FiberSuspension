
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"


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

// 2D version
static int calBucketIndex_pbc(int b_jx, int b_jy, linklist_constant *linklist){
  int tb_jx, tb_jy; //tb_jz; the temporary neighbor bucket's coordinate

  tb_jx = (b_jx + (linklist->nBx)) % (linklist->nBx);
  tb_jy = (b_jy + (linklist->nBy)) % (linklist->nBy);
//   tb_jz = (b_jz + (linklist->nBz)) % (linklist->nBz);

//   return  tb_jz*(linklist->nBxy) + tb_jy*(linklist->nBx) + tb_jx;
  return  tb_jy*(linklist->nBx) + tb_jx;
}


// 2D version
static double calDistance_ij_pbc(int i, int j, double *Position,linklist_constant *linklist){
  double xij, yij; // zij ;

  xij = Position[j*3  ] - Position[i*3  ];
  yij = Position[j*3+1] - Position[i*3+1];
//   zij = Position[j*3+2] - Position[i*3+2];

  // apply pbc
  if(xij < -linklist->Lx*0.50)  { xij += linklist->Lx; }
  else if(xij>linklist->Lx*0.50){ xij -= linklist->Lx; }

  if(yij < -linklist->Ly*0.50)  { yij += linklist->Ly; }
  else if(yij>linklist->Ly*0.50){ yij -= linklist->Ly; }

//   if(zij < -linklist->Lz*0.50)  { zij += linklist->Lz; }
//   else if(zij>linklist->Lz*0.50){ zij -= linklist->Lz; }

//   return  sqrt((xij*xij) + (yij*yij) + (zij*zij)); //the distance
  return  sqrt((xij*xij) + (yij*yij)); //the distance

}



/***************************
 * check overlap of fiber
 ***************************/

int check_parallel(double x0, double y0, double c, double d, int aspect_ratio, double ParticleDistance){

    // distance to j th fiber
    double r = fabs( c*x0 - y0 + d) / sqrt( c*c + 1);

    // don't overlap
    if( r > 0.1 * ParticleDistance){
        return 0;
    }else{
        return 1;
    }
}


// if overlap, return 1, else return 0
int check_range(double x_ol,int aspect_ratio, double x_s, double x_e, double x1, double x2){

    int flag1 = 0; // flag for making fiber
    int flag2 = 0; // flag for existing fiber

    // check x_ol is in x range of current fiber
    if(x_s < x_e){
        if( x_s <= x_ol && x_ol <= x_e ) flag1 = 1;
    }else{
        if( x_e <= x_ol && x_ol <= x_s ) flag1 = 1;
    }

    // check x_ol is in x range of existing fiber
    if(x1 < x2){
        if( x1 <= x_ol && x_ol <= x2 ) flag2 = 1;
    }else{
        if( x2 <= x_ol && x_ol <= x1 ) flag2 = 1;
    }


    return flag1 * flag2 ;
}

// example usuage
			// c = tan(theta[jf]);
			// d = Position[3*j+1] - Position[3*j] * c;
			/*
			if( a == c ){
				overlap_flag = check_parallel(Position_tmp[0],Position_tmp[1],c,d,aspect_ratio,ParticleDistance);
				if(overlap_flag == 1) break;
			}else{
				x_ol = (d-b)/(a-c);
				x1 = Position[3*j];
				x2 = Position[3*(j + aspect_ratio - 1)];
				overlap_flag = check_range(x_ol,aspect_ratio,x_s,x_e,x1,x2);
				if(overlap_flag == 1) break;
			}*/


// if overlap, return 1, else return 0
int check_cross(double Lx, int aspect_ratio, double xa, double ya,double xb, double yb,double xc, double yc,double xd, double yd){

  // consider the crossing AB and CD
  // labnote vol.3 p5

  int flag1 = 0; // flag for making fiber
  int flag2 = 0; // flag for existing fiber
  double CA_CD, CB_CD, AC_AB, AD_AB;  // cross product

  //check PBC
  if(xa - xb > aspect_ratio){ xb += Lx; }
  else if (xa - xb < -aspect_ratio){ xa += Lx; }

  if(xc - xd > aspect_ratio){ xd += Lx; }
  else if (xc - xd < -aspect_ratio){ xc += Lx; }

  CA_CD = (xa - xc) * (yd - yc) - (ya - yc) * (xd - xc);
  CB_CD = (xb - xc) * (yd - yc) - (yb - yc) * (xd - xc);

  AC_AB = (xc - xa) * (yb - ya) - (yc - ya) * (xb - xa);
  AD_AB = (xd - xa) * (yb - ya) - (yd - ya) * (xb - xa);

  if (CA_CD * CB_CD < 0) flag1 = 1;
  if (AC_AB * AD_AB < 0) flag2 = 1;

  return flag1 * flag2 ;
}


/*********************
 * get random integer
 *********************/

int GetRandom(int min,int max){
	return min + (int)(rand()*(max-min+1.0)/(1.0+RAND_MAX));
}



/*************************
 * function for main roop
 *************************/


linklist_constant *calcLinklistParameter(double ParticleDistance, double Re_forNumberDensity, double MIN_X, double MAX_X, double MIN_Y, double MAX_Y){
  linklist_constant *linklist;
  int temp_divide_number = 0; // x方向について
  linklist = malloc(sizeof(linklist_constant));

  linklist -> MIN_X = MIN_X ;
  linklist -> MIN_Y = MIN_Y ;
  linklist -> MIN_Z = -2.1 * ParticleDistance; // MIN_Z ;
  linklist -> MAX_X = MAX_X ;
  linklist -> MAX_Y = MAX_Y ;
  linklist -> MAX_Z = 2.1 * 2 * ParticleDistance; // MAX_Z ;
  linklist -> Lx    = linklist -> MAX_X - linklist -> MIN_X ;
  linklist -> Ly    = linklist -> MAX_Y - linklist -> MIN_Y ;
  linklist -> Lz    = linklist -> MAX_Z - linklist -> MIN_Z ;

  temp_divide_number = (int)(linklist->Lx/Re_forNumberDensity);      // use Lx (Since x direction is peliodic boundary, Lx must be `int * DB` )
  linklist -> DB    = linklist->Lx/(temp_divide_number);  // DB > temp_DBx --> safety (labnote vol.2 p108) if divide_number decrease, DB will increase
  linklist -> nBx   =(int)((linklist->Lx)/(linklist->DB));
  linklist -> nBy   =(int)((linklist->Ly)/(linklist->DB));
  // linklist -> nBz   =(int)((linklist->Lz)/(linklist->DB));
  linklist -> nBxy  =(linklist->nBx) * (linklist->nBy);
  // linklist -> nBxyz = (linklist->nBxy) * (linklist->nBz);

  //printf("tempnumber : %d\n",temp_divide_number);
  if(linklist->nBx < 3|| linklist->nBy < 3/*|| linklist->nBz < 3*/){
    fprintf(stderr,"******** WARNING ***********\n");
    fprintf(stderr,"check the simulation size \n");
    fprintf(stderr,"DB size : %lf \n",linklist->DB);
    exit(1);
  }

  return linklist;
}


// ----------------------------------------------------------


void makeBucket(int NumberOfParticles, int *BucketFirst, int *BucketLast, int *Nextof, int *ParticleType, double *Position,linklist_constant *linklist){
  int i = 0;
  int ix, iy; // iz; the coordinate of bucket
  int ib ;        //index of bucket
  int temp_BucketLast;
  int nBxy = linklist->nBxy;
  // int nBxyz = linklist->nBxyz;
  int nBx = linklist->nBx; double DBinv = 1/(linklist->DB);
  double MIN_X = linklist->MIN_X;
  double MIN_Y = linklist->MIN_Y;
  // double MIN_Z = linklist->MIN_Z;

  //initialize these three array
  for(i=0; i<nBxy; i++){BucketFirst[i] = -1; BucketLast[i] = -1;} // be careful about array size
  for(i=0; i<NumberOfParticles ; i++){Nextof[i] = -1;}


  //set Bucket
  for(int i=0;i<NumberOfParticles;i++){
	if(ParticleType[i] == GHOST)continue;
    ix = (int)((Position[i*3  ] - MIN_X)*DBinv);
    iy = (int)((Position[i*3+1] - MIN_Y)*DBinv);
    // iz = (int)((Position[i*3+2] - MIN_Z)*DBinv);

    // ib = iz*nBxy + iy*nBx + ix;  // bucket index (0 <= ib < nBxyz)
    ib = iy*nBx + ix;  // bucket index (0 <= ib <= nBxy)

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

/***************************
 * calculate number density
 ***************************/

void calNumberDensity2D(int NumberOfParticles, double Re_forNumberDensity, int *ParticleType, double *NumberDensity, double *Position,
                              int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int    i,j;
  double distance;
  double w;
  int b_ix, b_iy;// , b_iz  the coordinate of bucket
  int b_jx, b_jy;// , b_jz  the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket


// #pragma omp parallel for private(ib,i,b_jz,b_jy,b_jx,jb,j,distance,w)
//   for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    // ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    ib = b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      NumberDensity[i] = 0.0;
      if(ParticleType[i] != GHOST){
        // for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,linklist);
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
        }} //}
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }} //}

}

/**********************
 * calculate potential
 **********************/

// \sum\sum (n_i - n_j)^2 w(r_ij)

double calPotential(int NumberOfParticles, double Re_forNumberDensity, int *ParticleType, double *NumberDensity, double *Position,
                              int *BucketFirst, int *Nextof, linklist_constant *linklist){
  int    i,j;
  double distance;
  double w;
  int b_ix, b_iy;// , b_iz  the coordinate of bucket
  int b_jx, b_jy;// , b_jz  the neighbor bucket's coordinate
  int ib,jb ;        //index of bucket
  double pot = 0.0;       // total potential
  double ni ;


// #pragma omp parallel for private(ib,i,b_jz,b_jy,b_jx,jb,j,distance,w)
//   for(b_iz=0;b_iz<(linklist->nBz);b_iz++){
  for(b_iy=0;b_iy<(linklist->nBy);b_iy++){
  for(b_ix=0;b_ix<(linklist->nBx);b_ix++){
    // ib = b_iz*(linklist->nBxy) + b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    ib = b_iy*(linklist->nBx) + b_ix; //the index of i th bucket
    i = BucketFirst[ib];
    if(i == -1) continue;
    for(;;){
      if(ParticleType[i] == FLUID || ParticleType[i] == SOLID){
		ni = NumberDensity[i];
        // for(b_jz=b_iz-1;b_jz<=b_iz+1;b_jz++){ //roop in neighbor buckets
        for(b_jy=b_iy-1;b_jy<=b_iy+1;b_jy++){
        for(b_jx=b_ix-1;b_jx<=b_ix+1;b_jx++){
          jb = calBucketIndex_pbc(b_jx,b_jy,linklist);
          j = BucketFirst[jb];
          if(j == -1) continue;
          for(;;){ //beginning of neighbor list of ith particle
            if( ( i<j ) && (ParticleType[j]!=GHOST) ){
              distance = calDistance_ij_pbc(i,j,Position,linklist);
              w =  weight(distance, Re_forNumberDensity);
              pot += (ni - NumberDensity[j]) * (ni - NumberDensity[j]) * w;
            }
            j = Nextof[j];
            if(j==-1) break;
          }//end of neighbor list of ith particle
        }} //}
      }
      i = Nextof[i];
      if(i == -1) break;
    }
  }} //}

	return pot;
}


// U = \sum (n_i - n_0)^2
double calPotential2(int NumberOfParticles, double *NumberDensity, double n0){
  int    i;
  double pot = 0.0;       // total potential


  for(i=0; i<NumberOfParticles; i++){
    pot += (NumberDensity[i] - n0) * (NumberDensity[i] - n0);
  }
	return pot;
}
