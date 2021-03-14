
/**T
   Concepts: KSP^solving a system of linear equations
   Processors: 1
T**/

#include <petscksp.h>

#include "TypeDef.h"

/***********************************************
 * calculate only INNER_PARTICLEs
 ***********************************************/
/*
 static double *makeInnerArray(double *Array, int *index_inner, int n_inner){
   double *InnerArray;
   InnerArray = (double *)malloc(n_inner * sizeof(double));
   for(int i=0; i<n_inner; i++){
      InnerArray[i] = Array[index_inner[i]];
   }
   return InnerArray;
}
*/
static void returnInnerArray(double *Array, double *InnerArray, int *index_inner, int n_inner){
   for(int i=0; i<n_inner; i++){
      Array[index_inner[i]] = InnerArray[i];
   }
}

PetscErrorCode PETSc_Solver(petsc_object *petsc,double *Pressure, double *SourceTerm, int n,
                  int n_inner,int n_nonzero,  int *index_inner, int *order_inner, double *diag, double *val, int *col, int *row)
{
   PetscErrorCode ierr;
   Mat            A ;
   // PetscBool      isSymmetric;
   double *Inner_Pressure;
   int i;
   double one = 1.0;

   //initialize
   Inner_Pressure = (double *)malloc(n * sizeof(double));
   for( i=0; i<n; i++){
      Inner_Pressure[i] = 0.0;
   }

   //Assemble Matrix here----------------------------------------------------------------------------------------
   ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
   ierr = MatSetUp(A);CHKERRQ(ierr);

   //non diag part
   for( i=0; i<n_nonzero; i++){
      ierr = MatSetValues(A,1,&row[i],1,&col[i],&val[i],INSERT_VALUES); CHKERRQ(ierr);
   }
   //diag part
   for( i=0; i<n_inner; i++){
      ierr = MatSetValues(A,1,&i,1,&i,&diag[i],INSERT_VALUES); CHKERRQ(ierr);
   }
   //non fluid part
   for( i=n_inner; i<n; i++){
      ierr = MatSetValues(A,1,&i,1,&i,&one,INSERT_VALUES); CHKERRQ(ierr);
   }
   ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
   // A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner
   ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/*
   printf("----- A in petsc -------\n");
   ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
*/
/*
   ierr = MatIsSymmetric(A,0.0,&isSymmetric);CHKERRQ(ierr);
    if (!isSymmetric) {
      printf("*******************\n");
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: A is non-symmetric \n");CHKERRQ(ierr);
      printf("*******************\n");
   }else{
      printf("*******************\n");
      printf("A is symmetric\n");
      printf("*******************\n");
   }
*/
   //printf("----- A in petsc -------\n");
   //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


   //set vector from array------------------------------------------------------------------------------------------
   ierr = VecPlaceArray(petsc->x,Inner_Pressure);CHKERRQ(ierr);
   ierr = VecPlaceArray(petsc->b,SourceTerm);CHKERRQ(ierr);

   //set solver ---------------------------------------------------------------------------------------------------
   ierr = KSPSetOperators(petsc->ksp,A,A);CHKERRQ(ierr);
   ierr = KSPSetFromOptions(petsc->ksp);CHKERRQ(ierr);
   ierr = KSPSolve(petsc->ksp,petsc->b,petsc->x);CHKERRQ(ierr);
   // ierr = KSPView(petsc->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

   //Put back the PETSc array that belongs in the vector petsc->x--------------------------------------------------
   ierr = VecResetArray(petsc->x);CHKERRQ(ierr);
   ierr = VecResetArray(petsc->b);CHKERRQ(ierr);
   returnInnerArray(Pressure,Inner_Pressure,index_inner,n_inner);

   ierr = MatDestroy(&A);CHKERRQ(ierr);
   free(Inner_Pressure);
   return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------
/******************************
 * Don't assemble inner matrix
 *
 * the wall part in matrix is
 *       a_ij  = 1 (i = j)
 *             = 0 (i != j)
 * and b_i = 0 (i = wall)
 ******************************/

PetscErrorCode PETSc_Solver_coo(petsc_object *petsc,double *solution, double *SourceTerm, int n,
                  int n_nonzero, double *diag, double *val, int *col, int *row)
{
   PetscErrorCode ierr;
   Mat            A ;
   //PetscBool      isSymmetric;
   int i;
   //double one = 1.0;

   //Assemble Matrix here----------------------------------------------------------------------------------------
   ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
   ierr = MatSetUp(A);CHKERRQ(ierr);

   //non diag part
   for( i=0; i<n_nonzero; i++){
      ierr = MatSetValues(A,1,&row[i],1,&col[i],&val[i],INSERT_VALUES); CHKERRQ(ierr);
   }
   //diag part
   for( i=0; i<n; i++){
      ierr = MatSetValues(A,1,&i,1,&i,&diag[i],INSERT_VALUES); CHKERRQ(ierr);
   }

   ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
   // A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner
   ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/*
   printf("----- A for velocity -------\n");
   ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
*/
/*
   ierr = MatIsSymmetric(A,0.0,&isSymmetric);CHKERRQ(ierr);
    if (!isSymmetric) {
      printf("*******************\n");
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: A is non-symmetric \n");CHKERRQ(ierr);
      printf("*******************\n");
   }else{
      printf("*******************\n");
      printf("A is symmetric\n");
      printf("*******************\n");
   }
*/
   //set vector from array------------------------------------------------------------------------------------------
   ierr = VecPlaceArray(petsc->x,solution);CHKERRQ(ierr);
   ierr = VecPlaceArray(petsc->b,SourceTerm);CHKERRQ(ierr);

   // ierr = VecView(petsc->b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

   //set solver ---------------------------------------------------------------------------------------------------
   ierr = KSPSetOperators(petsc->ksp,A,A);CHKERRQ(ierr);
   ierr = KSPSetFromOptions(petsc->ksp);CHKERRQ(ierr);
   ierr = KSPSolve(petsc->ksp,petsc->b,petsc->x);CHKERRQ(ierr);
   //ierr = KSPView(petsc->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

   //Put back the PETSc array that belongs in the vector petsc->x--------------------------------------------------
   ierr = VecResetArray(petsc->x);CHKERRQ(ierr);
   ierr = VecResetArray(petsc->b);CHKERRQ(ierr);

   ierr = MatDestroy(&A);CHKERRQ(ierr);
   return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------

/************************************
 * Don't make A every time.
 * xyz are calculated at one function
 *
 * This function also adapts 2D (skip calculation of v_z)
 *************************************/

PetscErrorCode PETSc_Solver_coo_xyz(petsc_object *petsc,double *u_x,double *u_y,double *u_z, double *b_x,double *b_y,double *b_z, int n,
                  int n_nonzero, double *diag, double *val, int *col, int *row, int DIM)
{
   PetscErrorCode ierr;
   Mat            A ;
   //PetscBool      isSymmetric;
   int i;
   //double one = 1.0;

   //Assemble Matrix here----------------------------------------------------------------------------------------
   ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
   ierr = MatSetFromOptions(A);CHKERRQ(ierr);
   ierr = MatSetUp(A);CHKERRQ(ierr);

   //non diag part
   for( i=0; i<n_nonzero; i++){
      ierr = MatSetValues(A,1,&row[i],1,&col[i],&val[i],INSERT_VALUES); CHKERRQ(ierr);
   }
   //diag part
   for( i=0; i<n; i++){
      ierr = MatSetValues(A,1,&i,1,&i,&diag[i],INSERT_VALUES); CHKERRQ(ierr);
   }

   ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
   // A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner
   ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/*
   printf("----- A for velocity -------\n");
   ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
*/
/*
   ierr = MatIsSymmetric(A,0.0,&isSymmetric);CHKERRQ(ierr);
    if (!isSymmetric) {
      printf("*******************\n");
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: A is non-symmetric \n");CHKERRQ(ierr);
      printf("*******************\n");
   }else{
      printf("*******************\n");
      printf("A is symmetric\n");
      printf("*******************\n");
   }
*/

// calculation will be repeated about x,y,z
// ******** x ********************************************************************
   //set vector from array------------------------------------------------
   ierr = VecPlaceArray(petsc->x,u_x);CHKERRQ(ierr);
   ierr = VecPlaceArray(petsc->b,b_x);CHKERRQ(ierr);

   // ierr = VecView(petsc->b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

   //set solver ---------------------------------------------------------
   ierr = KSPSetOperators(petsc->ksp,A,A);CHKERRQ(ierr);
   ierr = KSPSetFromOptions(petsc->ksp);CHKERRQ(ierr);
   ierr = KSPSolve(petsc->ksp,petsc->b,petsc->x);CHKERRQ(ierr);
   //ierr = KSPView(petsc->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

   //Put back the PETSc array that belongs in the vector petsc->x--------
   ierr = VecResetArray(petsc->x);CHKERRQ(ierr);
   ierr = VecResetArray(petsc->b);CHKERRQ(ierr);


// ******** y **********************************************************************
   //set vector from array------------------------------------------------
   ierr = VecPlaceArray(petsc->x,u_y);CHKERRQ(ierr);
   ierr = VecPlaceArray(petsc->b,b_y);CHKERRQ(ierr);

   // ierr = VecView(petsc->b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

   //set solver ---------------------------------------------------------
   ierr = KSPSetOperators(petsc->ksp,A,A);CHKERRQ(ierr);
   ierr = KSPSetFromOptions(petsc->ksp);CHKERRQ(ierr);
   ierr = KSPSolve(petsc->ksp,petsc->b,petsc->x);CHKERRQ(ierr);
   //ierr = KSPView(petsc->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

   //Put back the PETSc array that belongs in the vector petsc->x--------
   ierr = VecResetArray(petsc->x);CHKERRQ(ierr);
   ierr = VecResetArray(petsc->b);CHKERRQ(ierr);


// ******** z ************************************************************************
   if (DIM == 3){
      //set vector from array------------------------------------------------
      ierr = VecPlaceArray(petsc->x,u_z);CHKERRQ(ierr);
      ierr = VecPlaceArray(petsc->b,b_z);CHKERRQ(ierr);

      // ierr = VecView(petsc->b,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

      //set solver ---------------------------------------------------------
      ierr = KSPSetOperators(petsc->ksp,A,A);CHKERRQ(ierr);
      ierr = KSPSetFromOptions(petsc->ksp);CHKERRQ(ierr);
      ierr = KSPSolve(petsc->ksp,petsc->b,petsc->x);CHKERRQ(ierr);
      //ierr = KSPView(petsc->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

      //Put back the PETSc array that belongs in the vector petsc->x--------
      ierr = VecResetArray(petsc->x);CHKERRQ(ierr);
      ierr = VecResetArray(petsc->b);CHKERRQ(ierr);
   }


   ierr = MatDestroy(&A);CHKERRQ(ierr);
   return 0;
}