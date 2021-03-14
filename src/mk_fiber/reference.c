/***********************************************
 * algin n fiber on same y cordinate.
 * the distance between each fibers is dx (= nx / n_fiber)
 * Each fiber's angle is setted in theta = pi
 * Particle Distance = 1.0
 * concentareted regime (x direction : dx, y direction : 2*dx - 2 )
 *
 * [how to execute]
 * ./2d input.txt
 *
 * labnote vol.2 p124
 * labnote vol.2 p187
 *
 *
 * This is copy of dev/mk_particle/src/mk_fn-mn-1_fp2D.c
************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"


int align_fiber(int count, int ny, int nx, int aspect_ratio, int lambda_y, int lambda_x, int *ParticleType, int *fiber_consist);

int main(int argc, char** argv){

    FILE *fp_input, *fp_output;
	const char *input;
	char input_temp[256];
    int DIM;
    double ParticleDistance = 1.0;
    int nx,ny,nxy;
    int NumberOfParticle = 0;
    int count,calc_lambda_particle,lambda_flag,lambda_x,lambda_y;
	double eps;
    int *ParticleType;
    double *Position;
	char outputfile[256];
	double MIN_X, MAX_X;
	double MIN_Y, MAX_Y;
	int n_fluid = 0;

	int aspect_ratio = 0;
	int n_fiber = 1;
	int *fiber_consist;
	// int dx, dy;
	int nfx, nfy;
	double bairitsu; // This is adjust the density of fiber
	int dy;


	// read
	printf("start mk_particle\n");
	strcpy(input_temp,argv[1]);
	input=input_temp;
	fp_input = fopen(input,"r");
	fscanf(fp_input,"DIM %d\n",&DIM);
	fscanf(fp_input,"nx %d\n",&nx);
	fscanf(fp_input,"ny %d\n",&ny);

	fscanf(fp_input,"aspect_ratio %d\n",&aspect_ratio);
	fscanf(fp_input,"bairitsu %lf\n",&bairitsu);
	// fscanf(fp_input,"n_fiber %d\n",&n_fiber);



	nfx = nx / aspect_ratio;
	dy = aspect_ratio / bairitsu ;
	nfy = 1 + 2 * (int)(((ny - 6) - aspect_ratio) * 0.5 / dy );
	n_fiber = nfx * nfy;

	// check condition
	if((ny - 6) < aspect_ratio){
		printf("*** invalid aspect ratio *** \n");
		printf("\n");
		exit(0);
	}
	// else if ( nx % n_fiber != 0){
	// 	printf("*** invalid n_fiber *** \n");
	// 	printf(" nx / n_fiber = %lf \n",(double)nx/n_fiber);
	// 	printf("\n");
	// 	exit(0);
	// }

    nxy = nx*ny;
	// dx = nx / n_fiber;
	// dy = (int)(ny - 5) / 4;

	printf("DIM %d\n",DIM);
	printf("nx:%d ny:%d nxy:%d\n", nx, ny,nxy);
	printf("aspect ratio : %d\n",aspect_ratio);
	printf("n_fiber : %d\n",n_fiber);

	MIN_X = 0.0;
	MAX_X = nx * ParticleDistance ;
	MIN_Y = -2 * ParticleDistance ;
	MAX_Y = (ny - 2) * ParticleDistance ;

	printf("MIN_X %lf\n",MIN_X);
	printf("MAX_X %lf\n",MAX_X);
	printf("MIN_Y %lf\n",MIN_Y);
	printf("MAX_Y %lf\n",MAX_Y);

    // arrange particles ----------------------------------------------------------------------------
	ParticleType = (int*)malloc(sizeof(int)*nxy);
	Position = (double*)malloc(sizeof(double)*nxy*3);
	fiber_consist = (int *)malloc(sizeof(int)*aspect_ratio*n_fiber);

	// arrange all particles
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iy*nx + ix;
		ParticleType[ip] = GHOST;
		Position[ip*3  ] = MIN_X + ParticleDistance*ix;
		Position[ip*3+1] = MIN_Y + ParticleDistance*iy;
		Position[ip*3+2] = ParticleDistance;
	}}

	// set particle type
	eps = 0.1 * ParticleDistance;
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iy*nx + ix;
		if( iy < 2 || iy > ny-3 ){
			ParticleType[ip] = DUMMY_WALL;
			NumberOfParticle++;
		}
		else if( iy == 2 /*|| iy == 3 */|| iy == ny-3 /*|| iy == ny-4*/ ){
			ParticleType[ip] = WALL;
			NumberOfParticle++;
		}
		else{
			ParticleType[ip] = FLUID;
			NumberOfParticle++;
			n_fluid++;
		}
	}}

	printf("NumberOfParticle:     %d\n",NumberOfParticle);

	// specify particle for calculating lambda
	// lambda_x and lambda_y equal to nx,ny of lambda particle
	lambda_x = (int)((MAX_X - MIN_X) *0.5/ParticleDistance);
	lambda_y = (int)((MAX_Y - MIN_Y) *0.5/ParticleDistance);
	lambda_flag = lambda_y*nx + lambda_x;
	printf("x : %d, y : %d, flag : %d\n",lambda_x,lambda_y,lambda_flag);
	if(ParticleType[lambda_flag] != FLUID){
		printf("fail to find particle for calculating lambda \n");
		exit(0);
	}

	// align fiber
	count = 0; // the number of fiber particles

	// labnote vol.2 p189
	for(int ifx=0; ifx<nfx; ifx++){
	for(int ify=0; ify<nfy; ify++){
		int ix ;
		int iy = lambda_y + dy * (ify - (int)(nfy * 0.5)) ;
		if( ify % 2 == 0){
			ix = ifx * aspect_ratio + (int)(aspect_ratio * 0.25); // 2 --> if aspect_ratio is 10
		}else {
			ix = ifx * aspect_ratio + (int)(aspect_ratio * 0.75);
		}
		count = align_fiber(count, ny, nx, aspect_ratio, iy, ix, ParticleType, fiber_consist);
	}}

	if(count != aspect_ratio * n_fiber) {
		printf("don't coresspond count with Number of fiber particles \n");
		exit(0);
	}

    // output -------------------------------------------------------------------------------------

	sprintf(outputfile,"f%d-%d_%d_%d-%d_2d.dat",n_fiber,aspect_ratio,NumberOfParticle,nx,ny);
	fp_output = fopen(outputfile, "w");
	printf("output file : %s \n",outputfile);
	fprintf(fp_output,"%d\n",NumberOfParticle);
	fprintf(fp_output,"%d\n",n_fluid);
	fprintf(fp_output,"%d\n",n_fiber);
	fprintf(fp_output,"%d\n",aspect_ratio);
    fprintf(fp_output,"MIN_X : %lf, MAX_X : %lf\n",MIN_X,MAX_X);
    fprintf(fp_output,"MIN_Y : %lf, MAX_Y : %lf\n",MIN_Y,MAX_Y);
    fprintf(fp_output,"MIN_Z : %lf, MAX_Z : %lf\n",-8.0*ParticleDistance,8.0*ParticleDistance);

	count=0;
	// Firstly, fiber partiles are printed
	for(int i=0;i<n_fiber*aspect_ratio;i++){
		int ip = fiber_consist[i];
		if(ip == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		count++;
	}

	// other particle are printed
	for(int iy=0;iy<ny;iy++){
	for(int ix=0;ix<nx;ix++){
		int ip = iy*nx + ix;
		if(ParticleType[ip]==GHOST || ParticleType[ip]==SOLID) continue;
		if(ip == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[ip],Position[ip*3],Position[ip*3+1],Position[ip*3+2]);
		count++;
	}}

	if(count != NumberOfParticle){
		printf("don't coresspond count with Number of particles \n");
	}

	fprintf(fp_output,"lambda_particle : %d\n",calc_lambda_particle);

	// output fiber consist
	for (int i = 0; i < aspect_ratio * n_fiber ; i++){
		fprintf(fp_output,"%d\n",i);
	}

	fclose(fp_input); /*fclose(fp_input2);*/ fclose(fp_output);
	free(ParticleType);	free(Position); free(fiber_consist);
	printf("end mk_particle\n");
	printf("\n");
	return 0;
}


int align_fiber(int count, int ny, int nx, int aspect_ratio, int lambda_y, int lambda_x, int *ParticleType, int *fiber_consist){
	// the center of fiber is (lamba_x, lambda_y)
	// count means the number of fiber particles

	if ((ny % 2 == 0) && (aspect_ratio % 2 == 0)){
		for (int iy = ( lambda_y - aspect_ratio * 0.5 ); iy < (lambda_y + aspect_ratio * 0.5); iy++){
			int ip = iy * nx + lambda_x ;
			ParticleType[ip] = SOLID ;
			fiber_consist[count] = ip ;
			count++ ;
		}
	}else if((ny % 2 == 1) && (aspect_ratio % 2 == 1)){
		for (int iy = ( lambda_y - (aspect_ratio - 1) * 0.5 ); iy < (lambda_y + (aspect_ratio - 1) * 0.5 + 1); iy++){
			int ip = iy * nx + lambda_x ;
			ParticleType[ip] = SOLID ;
			fiber_consist[count] = ip ;
			count++ ;
		}
	}

	return count;
}