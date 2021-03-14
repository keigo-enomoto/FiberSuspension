#ifndef TYPE_DEF_H
#define TYPE_DEF_H
//------------------------------------------------
#include <petscksp.h>
//------------------------------------------------
//  マクロ定義(Macro definition)
//------------------------------------------------
enum{

 GHOST = -1,
 FLUID =  0,
 SOLID = 1,
 WALL  =  2,
 DUMMY_WALL = 3

};

//------------------------------------------------
//  型定義(Type definition)
//------------------------------------------------
/**************************************************************
  *This struct is for input paramters                         *
***************************************************************/

typedef struct simulation_parameters
{
  double ParticleDistance ;
  int DIM ;
  double Rho[4] ;       // Density             : 0 -> fluid, 1 -> solid, 2 -> wall
  double Re_inv[4] ;    // 1.0/(Reynolds number) = nu_alpha / (l_0^2 * strain_rate) : 0 -> fluid, 1 -> solid, 2 -> wall
  double strain_rate ;
  double wall_speed ;
  double Gravity ;
  double DT ;
  double FINISH_TIME ;
  int OUTPUT_INTERVAL ;
  char *Init_place ;
  char *output_dir ;
  double ks;           // coefficient between next solid particle
  double kb;           // coefficient between next next solid particle
  int VTU_flag;        // if 1 , output vtu file. else don't output


    // int DIM ;
    // double ParticleDistance ;
    // double FluidDensity ;
    // double KinematicViscosity ;
    // double Gravity;
    // double DT ;
    // double FINISH_TIME ;
    // int OUTPUT_INTERVAL ;
    // char *Init_place;
    // char *output_dir;
    // double strain_rate;
    // double wall_speed;
    // double ks; // coefficient between next solid particle
    // double kd; // coefficient between next next solid particle

} simulation_parameters;

/**************************************************************
  *This struct is constats for Gradient and Laplacian model   *
***************************************************************/

typedef struct Gradient_constant
{
    double Re_forNumberDensity ;
    double Re_forGradient ;
    double Re_forLaplacian ;
    double N0_forNumberDensity ;
    double N0_forGradient ;
    double N0_forLaplacian ;
    double Lambda ;

} Gradient_constant;

/************************************************
 * This structs is for linklist                 *
 *************************************************/
typedef struct linklist_constant
{
  double DB ; //the one side of the bucket
  int nBx ; //the number of the buckets in x axis
  int nBy ; //the number of the buckets in y axis
  int nBz ; //the number of the buckets in z axis
  int nBxy ;
  int nBxyz ;

  double MIN_X ; //the minimum side of x in analysis area
  double MIN_Y ;
  double MIN_Z ;
  double MAX_X ; //the maximum side of x in analysis area
  double MAX_Y ;
  double MAX_Z ;

  double Lx;
  double Ly;
  double Lz;
} linklist_constant;

/****************************************************************
 * This structs is used when initial particle position is read. *
 ****************************************************************/

typedef struct range_sim{
  double MIN_X;
  double MIN_Y;
  double MIN_Z;
  double MAX_X;
  double MAX_Y;
  double MAX_Z;
}range_sim;

/************************************************
 * This structs is for holding petsc objects     *
 *************************************************/

typedef struct petsc_object{
  KSP         ksp;       /* linear solver context */
  Vec         x, b;       //x-> pressure , b->SourceTerm
} petsc_object;


//------------------------------------------------
//  プロトタイプ宣言(Prototype declaration)
//------------------------------------------------

//------------------------------------------------
#endif