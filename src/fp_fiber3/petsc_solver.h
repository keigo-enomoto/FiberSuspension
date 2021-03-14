#ifndef PETSC_SOLVER_H
#define PETSC_SOLVER_H
//------------------------------------------------
#include <petscksp.h>

//------------------------------------------------
//  マクロ定義(Macro definition)
//------------------------------------------------

//------------------------------------------------
//  型定義(Type definition)
//------------------------------------------------

//------------------------------------------------
//  プロトタイプ宣言(Prototype declaration)
//------------------------------------------------
PetscErrorCode PETSc_Solver(petsc_object *petsc,double *Pressure, double *SourceTerm, int n,
                  int n_inner,int n_nonzero,  int *index_inner, int *order_inner, double *diag, double *val, int *col, int *row);
PetscErrorCode PETSc_Solver_coo(petsc_object *petsc,double *solution, double *SourceTerm, int n,
                  int n_nonzero, double *diag, double *val, int *col, int *row);
PetscErrorCode PETSc_Solver_coo_xyz(petsc_object *petsc,double *u_x,double *u_y,double *u_z, double *b_x,double *b_y,double *b_z, int n,
                  int n_nonzero, double *diag, double *val, int *col, int *row, int DIM) ;

//------------------------------------------------
#endif