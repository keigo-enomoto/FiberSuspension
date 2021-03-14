
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h> // stat() in calTangentPhi()


#include "TypeDef.h"

#define PI 3.14159265358979323846264338327950288


// apply pbc for xij , yij, and zij
static void apply_pbc_distance_each_axis(double *xij, double *yij, double *zij, linklist_constant *linklist){

    // apply pbc
    if(*xij < -linklist->Lx*0.50)  { *xij += linklist->Lx; }
    else if(*xij>linklist->Lx*0.50){ *xij -= linklist->Lx; }

    // fiber particles don't go over y axis period because of wall
    // if(*yij < -linklist->Ly*0.50)  { *yij += linklist->Ly; }
    // else if(*yij>linklist->Ly*0.50){ *yij -= linklist->Ly; }

    if(*zij < -linklist->Lz*0.50)  { *zij += linklist->Lz; }
    else if(*zij>linklist->Lz*0.50){ *zij -= linklist->Lz; }

}
/*
static void checkSolidForce(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, double *Acceleration, simulation_parameters *parameters){
    int jf ;
    int i ,index;

    for ( jf = 0; jf < n_fiber; jf++)
    {
        printf("%d \n",jf);
        for (index = 0; index < aspect_ratio; index++)
        {
            i = fiber_consist[aspect_ratio * jf + index] ;              // give the index of the particle
            printf("%d (%d): %lf %lf %lf \n",index,i,Acceleration[i*3],Acceleration[i*3+1],Acceleration[i*3+2]);
        }

    }

}
*/
/****************************
 * calc fiber restore force ( labnote vol.2 p117 )
 * calc fiber bending force ( labnote vol.2 p120 ~ 122 )
 *
 * int aspect_ratio --> the number of particle in a fiber
 * int n_fiber      --> the number of fibers
 * int fiber_consist --> holding the Particle label of particles in fiber
 ****************************/

void calFiberRestoreForce(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, double *Acceleration, simulation_parameters *parameters, linklist_constant *linklist
, double *fiber_l_ave, double *fiber_theta_ave){
    int jf;     // index for fiber
    int i,j,k ; // index for particle in fiber
    double fs ; // fs = ks * (r_ij - l_0) / rij --> F = fs * (x_j - x_i)
    double xij=0.0, yij=0.0, zij=0.0, r_ij = 0.0 ;
    double xjk=0.0, yjk=0.0, zjk=0.0, r_jk = 0.0 ;
    double cos_j = 0.0 ;  // cos(theta_j)
    // double sin_j = 0.0 ;  // sin(theta_j)
    double theta_j = 0.0; // radian
    double ks = parameters->ks;  double kb = parameters->kb;
    double ai = 0.0 ;      // 1.0 / (rij * sin_j )
    double ak = 0.0 ;      // 1.0 / (rjk * sin_j )
    double nd_s ;
    double nd_b ;
    double f_i[3],f_k[3] ; // force exerted to i or k th particle
    int count_l=0, count_theta=0 ; // for fiber_*_ave

    if(parameters->strain_rate != 0.0){
        nd_s = 1.0 / (parameters->strain_rate * parameters->strain_rate) ; // for nondimension parameter (labnote vol.2 p123)
        nd_b = 1.0 / (parameters->strain_rate * parameters->strain_rate) ; // ParticleDistance is omitted because l_0 = 1.0
    }else{
        nd_s = 1.0;
        nd_b = 1.0;
    }

// Accerelation[i or j or k] isn't refferred by other thread (because parent roop is divided for each fiber )
// This calculation may not be heavy computation. plz reffer to gprof
// #pragma omp parallel for private(count_l,count_theta,i,j,xij,yij,zij,r_ij,fs,k,xjk,yjk,zjk,r_jk,cos_j,theta_j,ai,ak)
    for ( jf=0; jf<n_fiber; jf++){

        count_l = 0;
        count_theta = 0;
        fiber_l_ave[jf] = 0.0;
        fiber_theta_ave[jf] = 0.0;

        // calc the contribution of next fiber particle
        // Since exerting force to both side, the for roop is only one.
        // the most right particle does'nt have right side particle --> index<(aspect_ratio - 1 )
        for(int index=0; index<(aspect_ratio - 1); index++){

        // calc fs (sprint potential for next particle)
            i = fiber_consist[aspect_ratio * jf + index] ;              // give the index of the particle
            j = fiber_consist[aspect_ratio * jf + index + 1] ;          // give the index of next particle
            xij = Position[j*3  ] - Position[i*3  ];
            yij = Position[j*3+1] - Position[i*3+1];
            zij = Position[j*3+2] - Position[i*3+2];
            apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);     //PBC
            r_ij = sqrt( (xij*xij) + (yij*yij) + (zij*zij) ) ;
            fs = nd_s * ks * (r_ij - 1.0) / r_ij ;   // coefficient
            Acceleration[i*3  ] += fs * xij;
            Acceleration[i*3+1] += fs * yij;
            Acceleration[i*3+2] += fs * zij;
            Acceleration[j*3  ] -= fs * xij;                           // calc inverse direction
            Acceleration[j*3+1] -= fs * yij;
            Acceleration[j*3+2] -= fs * zij;

            fiber_l_ave[jf] += r_ij;
            count_l++ ;

        // calc fb ( bending potentail for i and k=(i+2) th particle )
            if (index < (aspect_ratio - 2)){
                k = fiber_consist[aspect_ratio * jf + index + 2] ;       // give the index of next next particle
                xjk = Position[k*3  ] - Position[j*3  ];
                yjk = Position[k*3+1] - Position[j*3+1];
                zjk = Position[k*3+2] - Position[j*3+2];
                apply_pbc_distance_each_axis(&xjk,&yjk,&zjk,linklist);   //PBC
                r_jk = sqrt( (xjk*xjk) + (yjk*yjk) + (zjk*zjk) ) ;
            // calc theta and some coefficient (labnote vol.2 p140)
                cos_j = (- 1.0) * ( (xij * xjk) + (yij * yjk) + (zij * zjk) ) / (r_ij * r_jk) ;
                if ( cos_j < -1.0 ) cos_j = -1.0;
                if ( cos_j > 1.0 )  cos_j = 1.0 ;
                theta_j = acos(cos_j) ;
                ai = nd_b / r_ij ;
                ak = nd_b / r_jk ;

                f_i[0] = - kb * ai * ( xjk/r_jk + cos_j * xij / r_ij );
                f_i[1] = - kb * ai * ( yjk/r_jk + cos_j * yij / r_ij );
                f_i[2] = - kb * ai * ( zjk/r_jk + cos_j * zij / r_ij );
                f_k[0] = kb * ak * ( xij/r_ij + cos_j * xjk / r_jk );
                f_k[1] = kb * ak * ( yij/r_ij + cos_j * yjk / r_jk );
                f_k[2] = kb * ak * ( zij/r_ij + cos_j * zjk / r_jk );

                Acceleration[i*3  ] += f_i[0];  // calc F_i
                Acceleration[i*3+1] += f_i[1];
                Acceleration[i*3+2] += f_i[2];
                Acceleration[j*3  ] -= f_i[0] + f_k[0];  // reaction force (This conserves the center of gravity )
                Acceleration[j*3+1] -= f_i[1] + f_k[1];
                Acceleration[j*3+2] -= f_i[2] + f_k[2];
                Acceleration[k*3  ] += f_k[0];  // calc F_k
                Acceleration[k*3+1] += f_k[1];
                Acceleration[k*3+2] += f_k[2];

                fiber_theta_ave[jf] += theta_j;
                count_theta ++;

            }
        }

        fiber_l_ave[jf] /= count_l;
        fiber_theta_ave[jf] /= count_theta;

    }
    // checkSolidForce(aspect_ratio,n_fiber,fiber_consist,Position,Acceleration,parameters);

}

// ----------------------------------------------------------------------------------------------------

/***********************************************************************************
 * calc "tan \phi" for calculating Jeffery orbit
 *
 * $$
 * tan(\phi) = 1/re tan[- \dot{\gamma}\frac{re}{re^2 + 1} + tan^{-1}(re tan \phi_0)]
 * $$
 *
 ************************************************************************************/

void output_tanphi_thin(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist){

    int jf ;
    int i,j;
    double xij, yij, zij=0;
    double tan_phi = 0.0 ;
    FILE *fp_t;
    char outputfile[1024];
    // struct stat st;
    // int output_flag;

    sprintf(outputfile, "%s/output_tanphi.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_t = fopen(outputfile,"a");

    // // if not exist outputfile
    // if (0 != output_flag) {
    //     fp_trintf(fp_t,"n_fiber : %d\n",n_fiber);
    //     fp_trintf(fp_t,"aspect_ratio : %d\n",aspect_ratio);
    // }

    fprintf(fp_t,"%lf  ",Time);

    for ( jf=0; jf<n_fiber; jf++){
        tan_phi = 0.0;  // iniitalize

        i = fiber_consist[aspect_ratio * jf ] ;                     // give the index of side particle
        j = fiber_consist[aspect_ratio * jf + aspect_ratio - 1] ;   // give the index of other side particle
        xij = Position[j*3  ] - Position[i*3  ];
        yij = Position[j*3+1] - Position[i*3+1];
        apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);
        if ( xij * xij < 1.0E-24 ){
            tan_phi = 1.0E+11 ;
        }else{
            tan_phi = yij / xij ;
        }
        fprintf(fp_t,"%lf  ",tan_phi);
    }
        fprintf(fp_t,"\n");

    fclose(fp_t);
}



/******************************
 * output unit vector component
 ******************************/
void output_unit_thin(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist){

    int jf ;
    int i,j;
    double xij, yij, zij=0;
    double norm = 0.0;
    FILE *fp_u;
    char outputfile[1024];
    // struct stat st;
    // int output_flag;

    sprintf(outputfile, "%s/output_unit.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_u = fopen(outputfile,"a");

    // // if not exist outputfile
    // if (0 != output_flag) {
    //     fp_urintf(fp_u,"n_fiber : %d\n",n_fiber);
    //     fp_urintf(fp_u,"aspect_ratio : %d\n",aspect_ratio);
    // }

    fprintf(fp_u,"%lf  ",Time);

    for ( jf=0; jf<n_fiber; jf++){

        i = fiber_consist[aspect_ratio * jf ] ;                     // give the index of side particle
        j = fiber_consist[aspect_ratio * jf + aspect_ratio - 1] ;   // give the index of other side particle
        xij = Position[j*3  ] - Position[i*3  ];
        yij = Position[j*3+1] - Position[i*3+1];
        // zij = Position[j*3+2] - Position[i*3+2];
        apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);
        norm = (xij * xij) + (yij * yij) + (zij * zij);
        norm = sqrt(norm);
        fprintf(fp_u,"%lf  %lf  ",xij/norm,yij/norm);

    }

        fprintf(fp_u,"\n");

    fclose(fp_u);
}





// need fiber_l_ave[], fiber_theta_ave[] calculate in calFiberRestoreForce() funciton
void output_monitor_fiber(int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,linklist_constant *linklist, double Time,
double *fiber_l_ave, double *fiber_theta_ave){

    int jf ;
    FILE *fp_l;
    char outputfile2[1024];
    double l=0.0, theta=0.0; // average
    // struct stat st;
    // int output_flag;

    sprintf(outputfile2, "%s/output_FiberMonitor.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_l = fopen(outputfile2,"a");

    // // if not exist outputfile
    // if (0 != output_flag) {
    //     fprintf(fp_l,"n_fiber : %d\n",n_fiber);
    //     fprintf(fp_l,"aspect_ratio : %d\n",aspect_ratio);
    // }

    fprintf(fp_l,"%lf  ",Time);
    for ( jf=0; jf<n_fiber; jf++){
        l += fiber_l_ave[jf];
        theta += fiber_theta_ave[jf];
    }
    l /= n_fiber;
    theta /= n_fiber;

    fprintf(fp_l,"%lf  ",l);
    fprintf(fp_l,"%lf  ",theta);
    fprintf(fp_l,"\n");
    fclose(fp_l);
}


// if use many fibers, plz use array holding the center coordinate of each fiber and parallel computing

void output_fiber_coordinate(int aspect_ratio, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters, double Time){

    int jf,i;
    int index;
    double cent_x, cent_y, cent_z; // the center of coordinate
    FILE *fp_c;
    char outputfile3[1024];
    double div = 1.0 / aspect_ratio ;

    sprintf(outputfile3, "%s/output_center.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_c = fopen(outputfile3,"a");

    fprintf(fp_c,"%lf  ",Time);
    for(jf = 0; jf < n_fiber; jf++){
        cent_x = 0.0; cent_y = 0.0; cent_z = 0.0;
        for(index = 0; index<aspect_ratio; index++){
            i = fiber_consist[aspect_ratio*jf + index] ;
            cent_x += Position[3*i];
            cent_y += Position[3*i+1];
            cent_z += Position[3*i+2];
        }
        fprintf(fp_c,"%lf  %lf  %lf  ",cent_x*div,cent_y*div,cent_z*div);
    }
    fprintf(fp_c,"\n");
    fclose(fp_c);

}


// initial operation for output file
void make_output_file_fiber(simulation_parameters *parameters, int n_fiber){
    int i;
    // FILE *fp_t;
    // char outputfile[1024];
    FILE *fp_l;
    char outputfile2[1024];
    FILE *fp_c;
    char outputfile3[1024];
    FILE *fp_u;
    char outputfile4[1024];

    // output_tanphi_thin()
    // sprintf(outputfile, "%s/output_tanphi.dat",parameters->output_dir); //sprintf is the function to write the word to array
    // fp_t = fopen(outputfile,"w");
    // fprintf(fp_t,"Time  ");
    // for ( i = 0; i < n_fiber; i++){
    //     fprintf(fp_t,"theta%d  ",i+1);
    // }
    // fprintf(fp_t,"\n");
    // fclose(fp_t);

    // output_tanphi_thin()
    sprintf(outputfile4, "%s/output_unit.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_u = fopen(outputfile4,"w");
    fprintf(fp_u,"Time  ");
    for ( i = 0; i < n_fiber; i++){
        fprintf(fp_u,"ux%d  uy%d  ",i+1,i+1);
    }
    fprintf(fp_u,"\n");
    fclose(fp_u);

    // output_monitor_fiber()
    sprintf(outputfile2, "%s/output_FiberMonitor.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_l = fopen(outputfile2,"w");
    fprintf(fp_l,"Time  ");
    fprintf(fp_l,"ave_len  ave_theta\n");
    fprintf(fp_l,"\n");
    fclose(fp_l);

    // output_fiber_coordinate()
    sprintf(outputfile3, "%s/output_center.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_c = fopen(outputfile3,"w");
    fprintf(fp_c,"Time  ");
    for ( i = 0; i < n_fiber; i++){
        fprintf(fp_c,"x%d  y%d  z%d  ",i+1, i+1, i+1);
    }
    fprintf(fp_c,"\n");
    fclose(fp_c);
}

// -----------------------------------------------------------------------------------------------------------------------
/*******************************************************************
 *      Code for thick fiber
 * *****************************************************************/

// check labnote vol.2 p185,186
void calFiberRestoreForce_thick(int l_f, int a_f, int n_fiber, int *fiber_consist, double *Position, double *Acceleration, simulation_parameters *parameters, linklist_constant *linklist
, double *fiber_l_ave, double *fiber_theta_ave){
    int jf;     // index for fiber
    int i,j,k ; // index for particle in fiber
    double fs ; // fs = ks * (r_ij - l_0) / rij --> F = fs * (x_j - x_i)
    double xij=0.0, yij=0.0, zij=0.0, r_ij = 0.0 ;
    double xjk=0.0, yjk=0.0, zjk=0.0, r_jk = 0.0 ;
    double cos_j = 0.0 ;  // cos(theta_j)
    // double sin_j = 0.0 ;  // sin(theta_j)
    double theta_j = 0.0; // radian
    double ks = parameters->ks;  double kb = parameters->kb;
    double ai = 0.0 ;      // 1.0 / (rij * sin_j )
    double ak = 0.0 ;      // 1.0 / (rjk * sin_j )
    double nd_s ; // for nondimension parameter (labnote vol.2 p123)
    double nd_b ; // ParticleDistance is omitted because l_0 = 1.0
    double f_i[3],f_k[3] ; // force exerted to i or k th particle
    int count_l=0, count_theta=0 ; // for fiber_*_ave
    int n_in_fiber = l_f * a_f; // the number of particle in a fiber

    if(parameters->strain_rate != 0.0){
        nd_s = 1.0 / (parameters->strain_rate * parameters->strain_rate) ; // for nondimension parameter (labnote vol.2 p123)
        nd_b = 1.0 / (parameters->strain_rate * parameters->strain_rate) ; // ParticleDistance is omitted because l_0 = 1.0
    }else{
        nd_s = 1.0;
        nd_b = 1.0;
    }

    for ( jf=0; jf<n_fiber; jf++){

        count_l = 0;
        count_theta = 0;
        fiber_l_ave[jf] = 0.0;
        fiber_theta_ave[jf] = 0.0;

        // calc the contribution of next fiber particle and next row particle
        // Since exerting force to both side, the for roop is only one for each direction (l_f and a_f).
        // the most right particle does'nt have right side particle --> index<(l_f - 1 )

        for(int irow=0; irow < a_f; irow++ ){ // 0 or 1
        for(int index=0; index<(l_f - 1); index++){

        // calc fs (sprint potential for next particle)
            i = fiber_consist[n_in_fiber * jf + irow * l_f + index] ;              // give the index of the particle
            j = fiber_consist[n_in_fiber * jf + irow * l_f + index + 1] ;          // give the index of next particle
            xij = Position[j*3  ] - Position[i*3  ];
            yij = Position[j*3+1] - Position[i*3+1];
            zij = Position[j*3+2] - Position[i*3+2];
            apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);     //PBC
            r_ij = sqrt( (xij*xij) + (yij*yij) + (zij*zij) ) ;
            fs = nd_s * ks * (r_ij - 1.0) / r_ij ;   // coefficient
            Acceleration[i*3  ] += fs * xij;
            Acceleration[i*3+1] += fs * yij;
            Acceleration[i*3+2] += fs * zij;
            Acceleration[j*3  ] -= fs * xij;                           // calc inverse direction
            Acceleration[j*3+1] -= fs * yij;
            Acceleration[j*3+2] -= fs * zij;

            fiber_l_ave[jf] += r_ij;
            ++ count_l;

        // calc fb ( bending potentail for i and k=(i+2) th particle )
            if (index < (l_f - 2)){
                k = fiber_consist[n_in_fiber * jf + irow * l_f + index + 2] ;       // give the index of next next particle
                xjk = Position[k*3  ] - Position[j*3  ];
                yjk = Position[k*3+1] - Position[j*3+1];
                zjk = Position[k*3+2] - Position[j*3+2];
                apply_pbc_distance_each_axis(&xjk,&yjk,&zjk,linklist);   //PBC
                r_jk = sqrt( (xjk*xjk) + (yjk*yjk) + (zjk*zjk) ) ;
            // calc theta and some coefficient (labnote vol.2 p140)
                cos_j = (- 1.0) * ( (xij * xjk) + (yij * yjk) + (zij * zjk) ) / (r_ij * r_jk) ;
                if ( cos_j < -1.0 ) cos_j = -1.0;
                if ( cos_j > 1.0 )  cos_j = 1.0 ;
                theta_j = acos(cos_j) ;
                ai = nd_b / r_ij ;
                ak = nd_b / r_jk ;

                f_i[0] = - kb * ai * ( xjk/r_jk + cos_j * xij / r_ij );
                f_i[1] = - kb * ai * ( yjk/r_jk + cos_j * yij / r_ij );
                f_i[2] = - kb * ai * ( zjk/r_jk + cos_j * zij / r_ij );
                f_k[0] = kb * ak * ( xij/r_ij + cos_j * xjk / r_jk );
                f_k[1] = kb * ak * ( yij/r_ij + cos_j * yjk / r_jk );
                f_k[2] = kb * ak * ( zij/r_ij + cos_j * zjk / r_jk );

                Acceleration[i*3  ] += f_i[0];  // calc F_i
                Acceleration[i*3+1] += f_i[1];
                Acceleration[i*3+2] += f_i[2];
                Acceleration[j*3  ] -= (f_i[0] + f_k[0]);  // reaction force (This conserves the center of gravity )
                Acceleration[j*3+1] -= (f_i[1] + f_k[1]);
                Acceleration[j*3+2] -= (f_i[2] + f_k[2]);
                Acceleration[k*3  ] += f_k[0];  // calc F_k
                Acceleration[k*3+1] += f_k[1];
                Acceleration[k*3+2] += f_k[2];

                fiber_theta_ave[jf] += theta_j;
                ++ count_theta;

            }

            // calc fs for next row particle (sprint potential for next particle)
            if(irow == 0){
                j = fiber_consist[n_in_fiber*jf + l_f + index] ;  // give the index of next row particle
                xij = Position[j*3  ] - Position[i*3  ];
                yij = Position[j*3+1] - Position[i*3+1];
                zij = Position[j*3+2] - Position[i*3+2];
                apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);     //PBC
                r_ij = sqrt( (xij*xij) + (yij*yij) + (zij*zij) ) ;
                fs = nd_s * ks * (r_ij - 1.0) / r_ij ;   // coefficient
                Acceleration[i*3  ] += fs * xij;
                Acceleration[i*3+1] += fs * yij;
                Acceleration[i*3+2] += fs * zij;
                Acceleration[j*3  ] -= fs * xij;         // calc inverse direction
                Acceleration[j*3+1] -= fs * yij;
                Acceleration[j*3+2] -= fs * zij;
            }
        }

            // side end particle and next row particle
            if(irow == 0){
                i = fiber_consist[n_in_fiber*jf    +    l_f - 1] ;    // give the index of the particle
                j = fiber_consist[n_in_fiber*jf + l_f + l_f - 1] ;  // give the index of next row particle
                xij = Position[j*3  ] - Position[i*3  ];
                yij = Position[j*3+1] - Position[i*3+1];
                zij = Position[j*3+2] - Position[i*3+2];
                apply_pbc_distance_each_axis(&xij,&yij,&zij,linklist);     //PBC
                r_ij = sqrt( (xij*xij) + (yij*yij) + (zij*zij) ) ;
                fs = nd_s * ks * (r_ij - 1.0) / r_ij ;   // coefficient
                Acceleration[i*3  ] += fs * xij;
                Acceleration[i*3+1] += fs * yij;
                Acceleration[i*3+2] += fs * zij;
                Acceleration[j*3  ] -= fs * xij;         // calc inverse direction
                Acceleration[j*3+1] -= fs * yij;
                Acceleration[j*3+2] -= fs * zij;
            }
        }

        fiber_l_ave[jf] /= count_l;
        fiber_theta_ave[jf] /= count_theta;

    }
    // checkSolidForce(aspect_ratio,n_fiber,fiber_consist,Position,Acceleration,parameters);

}

/*****************************
 * calc tanphi
 *****************************/

void output_tanphi_thick(int l_f, int a_f, int n_fiber, int *fiber_consist, double *Position, simulation_parameters *parameters,double Time,linklist_constant *linklist){

    int jf ;
    int i1,i2;
    int j1,j2;
    double xij, yij;
    double tan_phi = 0.0 ;
    FILE *fp_t;
    char outputfile[1024];
    int n_in_fiber = l_f * a_f; // the number of particle in a fiber
    int aspect_ratio = l_f / a_f;
    double xs, ys; // start
    double xe, ye; // end
    // struct stat st;
    // int output_flag;

    sprintf(outputfile, "%s/output_tanphi.dat",parameters->output_dir); //sprintf is the function to write the word to array
    fp_t = fopen(outputfile,"a");

    // // if not exist outputfile
    // if (0 != output_flag) {
    //     fp_trintf(fp_t,"n_fiber : %d\n",n_fiber);
    //     fp_trintf(fp_t,"aspect_ratio : %d\n",aspect_ratio);
    // }

    fprintf(fp_t,"%lf  ",Time);

    for ( jf=0; jf<n_fiber; jf++){
        tan_phi = 0.0;  // iniitalize

        i1 = fiber_consist[ n_in_fiber*jf ] ;                  // give the index of start particle
        i2 = fiber_consist[ n_in_fiber*jf + l_f] ;             // give the index of start particle
        j1 = fiber_consist[ n_in_fiber*jf + l_f - 1] ;         // give the index of end  particle
        j2 = fiber_consist[ n_in_fiber*jf + l_f - 1 + l_f] ;   // give the index of end  particle

        xs = (Position[3*i1] + Position[3*i2]) * 0.5;
        ys = (Position[3*i1+1] + Position[3*i2+1]) * 0.5;
        xe = (Position[3*j1] + Position[3*j2]) * 0.5;
        ye = (Position[3*j1+1] + Position[3*j2+1]) * 0.5;
        xij = xs - xe;
        yij = ys - ye;

        if ( fabs( xij ) < 1.0E-12 ){
            tan_phi = aspect_ratio * 1.0E+12 ;
        }else{
            tan_phi = yij / xij ;
        }
        fprintf(fp_t,"%lf  ",tan_phi);
    }
        fprintf(fp_t,"\n");

    fclose(fp_t);
}
