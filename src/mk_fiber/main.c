/**********************************************************
 * algin fiber in random direction and pack fluid particle
 *
 * This procedure is (1) ~ (7) 2.6 p20 in my graduate thesis
***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "TypeDef.h"
#include "CalculateFunction.h"

#define PI 3.14159265358979323846264338327950288


void output(int n_fiber, int aspect_ratio,int NumberOfParticles, int nx, int ny, int n_fluid, linklist_constant *linklist,
double ParticleDistance, int lambda_flag, int *ParticleType, double *Position, int seed);

int main(int argc, char** argv){

    FILE *fp_input;
	const char *input;
	char input_temp[256];
    double ParticleDistance = 1.0;
    int DIM;
    int nx,ny,nxy;
    int NumberOfParticles = 0;
    int calc_lambda_particle,lambda_flag = 0,lambda_x,lambda_y;
    int *ParticleType;
    double *Position;
    double *NewPosition;
	double MIN_X,MAX_X;
	double MIN_Y,MAX_Y;
	double border_u, border_d;  // the border of wall in y axis
	double eps; // epsilon
	int n_fluid = 0; int n_wall = 0; int n_solid = 0; int n_sf=0;
	int index = 0; // index for particle

	int aspect_ratio = 0;
	int n_fiber = 1;
	int seed; // for randomness
	// int *fiber_consist;
	int count_fiber, overlap_flag; // for making fiber alignment
	double rc1,rc2,rct;   // random coefficient
	double *Position_tmp; // temporary fiber position
	double theta;

	double *NumberDensity;
	double *NewNumberDensity;
	double U_new = 0.0;
	double U_old = 0.0;
	linklist_constant *linklist;       // hold parameter for cell list
	int    *BucketFirst;       // the label of first particle in the bucket
	int    *BucketLast;        // the label of last particle in the bucket
	int    *Nextof;            // the label of next particle in the bucket
	double Re_forNumberDensity;
	int size_posi; // size of array for memory copy
	int size_n;    // size of array for memory copy
	int mcstep = 0;    // the index for main roop
	int missstep = 0;  // the times of failure try
	int target;    // chosen random particle
	double xa,xb,xc,xd;  // for checking fiber overlapping
	double ya,yb,yc,yd; // for checking target fiber overlapping

	double n0 = 6.539697; // basis of number density

	int i,j,jf;


// ==================================================================================================
	// read input file
	printf("start mk_particle\n");
	strcpy(input_temp,argv[1]);
	input=input_temp;
	fp_input = fopen(input,"r");
	fscanf(fp_input,"DIM %d\n",&DIM);
	fscanf(fp_input,"nx %d\n",&nx);
	fscanf(fp_input,"ny %d\n",&ny);
    nxy = nx*ny;
	fscanf(fp_input,"aspect_ratio %d\n",&aspect_ratio);
	fscanf(fp_input,"n_fiber %d\n",&n_fiber);
	fscanf(fp_input,"seed %d\n",&seed);
	// check condition
	if((ny - 6) < aspect_ratio){
		printf("*** invalid aspect ratio *** \n");
		printf("\n");
		exit(1);
	}else if ( aspect_ratio * n_fiber > nx * (ny - 6) ){
		printf("*** invalid n_fiber *** \n");
		printf(" n_solid = %d > NumberOfParticles = %d \n",aspect_ratio*n_fiber, nxy);
		printf("\n");
		exit(1);
	}
	printf("DIM %d\n",DIM);
	printf("nx:%d ny:%d nxy:%d\n", nx, ny,nxy);
	printf("aspect ratio : %d\n",aspect_ratio);
	printf("n_fiber : %d\n",n_fiber); printf("seed : %d\n",seed);

	MIN_X = 0.0;
	MAX_X = nx * ParticleDistance ;
	MIN_Y = -2 * ParticleDistance ;
	MAX_Y = (ny - 2) * ParticleDistance ;
	// MIN_Z = -2.1 * ParticleDistance ;
	// MAX_Z = 2.1 * 2 * ParticleDistance ;
	Re_forNumberDensity = 2.1 * ParticleDistance;

	// wall border in y axis
	border_d = 0.0 ;
	border_u = MAX_Y - 3 * ParticleDistance;
	eps = 0.5 * ParticleDistance;

	// number of each type of particle
	n_solid = aspect_ratio * n_fiber;
	n_wall = nx * 6; // 3 + 3. WALL + GHOST
	n_fluid = nxy - n_solid - n_wall;

	printf("MIN_X %lf\n",MIN_X); printf("MAX_X %lf\n",MAX_X);
	printf("MIN_Y %lf\n",MIN_Y); printf("MAX_Y %lf\n",MAX_Y);
	printf("upper wall border %lf\n",border_u);
	printf("lower wall border %lf\n",border_d);

// ==================================================================================================
    // arrange particles
	ParticleType = (int*)malloc(sizeof(int)*nxy);
	Position = (double*)malloc(sizeof(double)*nxy*3);
	// fiber_consist = (int *)malloc(sizeof(int)*aspect_ratio*n_fiber);
	Position_tmp = (double *)malloc(sizeof(double) * aspect_ratio * 3);

	linklist = calcLinklistParameter(ParticleDistance,Re_forNumberDensity,MIN_X,MAX_X,MIN_Y,MAX_Y);

	// arrange fiber particle ---------------------------------------
	// labnote vol.3 p1
	srand(seed);
	count_fiber = 0;
	index = 0;

	while (count_fiber < n_fiber){

		overlap_flag = 0;
		// set x, y, theta
		rc1 = (double)rand()/RAND_MAX ;   //0 < rc < 1.0
		rc2 = (double)rand()/RAND_MAX ;
		rct = (double)rand()/RAND_MAX ;
		Position_tmp[0] = rc1 * (MAX_X - MIN_X) ;
		Position_tmp[1] = rc2 * (border_u - border_d - 2*eps) + eps; // keep distance from wall 0.5*ParticleDistance
		Position_tmp[2] = ParticleDistance;
		theta = rct * PI * 2;

		// printf("%lf  %lf  %lf  \n",rc1,rc2,rct);

		// set fiber particle
		for ( i = 1; i < aspect_ratio; i++){
			Position_tmp[3*i  ] = i * cos(theta) + Position_tmp[0] ; //omit multiple particle distance
			Position_tmp[3*i+1] = i * sin(theta) + Position_tmp[1];
			Position_tmp[3*i+2] = Position_tmp[2];
		}

		xa = Position_tmp[0];
		ya = Position_tmp[1];
    	xb = Position_tmp[3*(aspect_ratio-1)];   // other side of particle
    	yb = Position_tmp[3*(aspect_ratio-1)+1];

		// check end particle whether the position is valid or not considering wall arrangement
		if( yb > border_u - eps || yb < border_d + eps ) continue;
		if( xb >= MAX_X  || xb < MIN_X ) continue;

		// check overlapping with other fiber or wall [ labnote vol.2 p200 ]
		// overlap -> overlap_flag = 1

		for( jf=0; jf<count_fiber; jf++){
			j = aspect_ratio * jf;
			xc = Position[3*j];
			xd = Position[3*(j + aspect_ratio - 1)];
			yc = Position[3*j+1];
			yd = Position[3*(j + aspect_ratio - 1)+1];
			overlap_flag = check_cross(linklist->Lx,aspect_ratio,xa,ya,xb,yb,xc,yc,xd,yd);
			if(overlap_flag == 1) break;
		}

		if(overlap_flag == 1) {
			missstep ++;
			continue;
		}else{
			// register fiber particle
			for(i=0; i<aspect_ratio; i++){
				Position[3*index] = Position_tmp[3*i];
				Position[3*index+1] = Position_tmp[3*i+1];
				Position[3*index+2] = Position_tmp[3*i+2];
				ParticleType[index] = SOLID;
				index ++;
			}
			count_fiber ++;
		}
	}

	printf("miss : %d\n",missstep);

	if(count_fiber != n_fiber) {
		fprintf(stderr,"don't coresspond count with Number of fiber \n");
		fprintf(stderr,"count_fiber : %d\n",count_fiber);
		fprintf(stderr,"n_fiber : %d\n",n_fiber);
		exit(2);
	}
	if(index != n_solid){
		fprintf(stderr,"particle index does'nt correspond to ideal value\n");
		fprintf(stderr,"n_solid  = %d\n",n_solid);
		fprintf(stderr,"index = %d",index);
		exit(2);
	}

	// apply PBC for fiber particles
	for( jf=0; jf<n_fiber; jf++ ){
		for( i=0; i<aspect_ratio; i++ ){
			j = aspect_ratio * jf + i; // particle index
			if(Position[3*j] < MIN_X){
				Position[3*j] += ( MAX_X - MIN_X );
			}else if(Position[3*j] >= MAX_X){
				Position[3*j] -= ( MAX_X - MIN_X );
			}
		}
	}

	// arrange fluid particles -------------------------------------------------
  	// labnote vol.2 p16
	for(i = 0; i < n_fluid; i++){
		rc1 = (double)rand()/RAND_MAX ;   //0 < rc < 1.0
		rc2 = (double)rand()/RAND_MAX ;
		Position[3*index  ] = rc1 * (MAX_X - MIN_X) ;
		Position[3*index+1] = rc2 * (border_u - border_d - 2*eps) + eps; // keep distance from wall 0.5*ParticleDistance
		Position[3*index+2] = ParticleDistance;
		ParticleType[index] = FLUID;
		index ++;
	}

	//arrange wall particles ----------------------------------------------------

	//lower wall
	for(int iy=0; iy<3; iy++){
	for(int ix=0; ix<nx; ix++){
		Position[index*3  ] = MIN_X + ParticleDistance*ix;
		Position[index*3+1] = MIN_Y + ParticleDistance*iy;
		Position[index*3+2] = ParticleDistance;
		if(iy < 2){ ParticleType[index] = DUMMY_WALL; }
		else{       ParticleType[index] = WALL;}
		index ++;
	}}

	//upper wall
	for(int iy=0; iy<3; iy++){
	for(int ix=0; ix<nx; ix++){
		Position[index*3  ] = MIN_X + ParticleDistance*ix;
		Position[index*3+1] = border_u + ParticleDistance*iy;
		Position[index*3+2] = ParticleDistance;
		if(iy == 0){ParticleType[index] = WALL;}
		else{ ParticleType[index] = DUMMY_WALL;}
		index ++;
	}}

	if(index != nxy){
		fprintf(stderr,"*** number of particle does'nt correspond to nxy ***\n");
		fprintf(stderr," index : %d",index);
		exit(2);
	}

	NumberOfParticles = index;
	// printf("NumberOfParticles:     %d\n",NumberOfParticles);

// ==================================================================================================
// main part montecalro


linklist = calcLinklistParameter(ParticleDistance,Re_forNumberDensity,MIN_X,MAX_X,MIN_Y,MAX_Y);

output(n_fiber,aspect_ratio,NumberOfParticles,nx,ny,n_fluid,linklist,ParticleDistance,lambda_flag,ParticleType,Position,seed);


// array size is nxy (2D)
BucketFirst = (int *)malloc(sizeof(int) * (linklist->nBxy));	//バケットに格納された先頭の粒子の番号
BucketLast  = (int *)malloc(sizeof(int) * (linklist->nBxy));	//バケットに格納された最後尾の粒子の番号
Nextof      = (int *)malloc(sizeof(int) * NumberOfParticles);	//同じバケット内の次の粒子の番号
NumberDensity = (double *)malloc(NumberOfParticles * sizeof(double));
NewPosition = (double*)malloc(sizeof(double)*NumberOfParticles*3);
NewNumberDensity = (double *)malloc(NumberOfParticles * sizeof(double));

size_posi = sizeof(double)*NumberOfParticles*3;
size_n = NumberOfParticles * sizeof(double);

// make Bucket
makeBucket(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,Position,linklist);
calNumberDensity2D(NumberOfParticles,Re_forNumberDensity,ParticleType,NumberDensity,Position,BucketFirst,Nextof,linklist);
// U_old = calPotential(NumberOfParticles,Re_forNumberDensity,ParticleType,NumberDensity,Position,BucketFirst,Nextof,linklist);
U_old = calPotential2(NumberOfParticles,NumberDensity,n0);
U_new = U_old ;

printf("potential : %lf \n",U_old);

memcpy(NewPosition,Position,size_posi);        // copy Position -> NewPosition
memcpy(NewNumberDensity,NumberDensity,size_n); // copy NumberDensity -> NewNumberDensity

n_sf = n_fluid+n_solid-1; // index of last fluid particle

missstep = 0;
mcstep = 0;
while ( U_old > 0.1 * NumberOfParticles ){

	overlap_flag = 0;
	target = GetRandom(0,n_sf);   // index of target particle

	// give random displacement
	if ( U_old < NumberOfParticles ){
		rc1 = ( (double)rand()/RAND_MAX - 0.5)*0.1 ;   // 0 < rc < 1.0
		rc2 = ( (double)rand()/RAND_MAX - 0.5)*0.1 ;   // omit multiple Particle Distance
	}else{
		rc1 = (double)rand()/RAND_MAX - 0.5 ;   // 0 < rc < 1.0
		rc2 = (double)rand()/RAND_MAX - 0.5 ;   // omit multiple Particle Distance
	}

	// apply random displacement
	if(ParticleType[target] == FLUID){
		NewPosition[3*target] += rc1;
		NewPosition[3*target+1] += rc2;
		//x
		if(NewPosition[3*target]>=linklist->MAX_X)      {NewPosition[3*target]-=linklist->Lx;}
		else if(NewPosition[3*target]< linklist->MIN_X){NewPosition[3*target]+=linklist->Lx;}
		//y
		if(NewPosition[3*target+1]>=linklist->MAX_Y)      {NewPosition[3*target+1]-=linklist->Ly;}
		else if(NewPosition[3*target+1]< linklist->MIN_Y){NewPosition[3*target+1]+=linklist->Ly;}
		// check wall
		if( NewPosition[3*target+1] > border_u - eps || NewPosition[3*target+1] < border_d + eps ){
			overlap_flag = 1;
		}
	}else if (ParticleType[target] == SOLID){
		jf = (int)(target / aspect_ratio);  // index of fiber
		for(i=0; i<aspect_ratio; i++){
			j = jf * aspect_ratio + i;
			NewPosition[3*j] += rc1;
			NewPosition[3*j+1] += rc2;
			//x
			if(NewPosition[3*j]>=linklist->MAX_X)     {NewPosition[3*j]-=linklist->Lx;}
			else if(NewPosition[3*j]< linklist->MIN_X){NewPosition[3*j]+=linklist->Lx;}
			//y
			if(NewPosition[3*j+1]>=linklist->MAX_Y)     {NewPosition[3*j+1]-=linklist->Ly;}
			else if(NewPosition[3*j+1]< linklist->MIN_Y){NewPosition[3*j+1]+=linklist->Ly;}
		}

		// check fiber overlap
		xa = NewPosition[3* jf*aspect_ratio];
		ya = NewPosition[3* jf*aspect_ratio+1];
    	xb = NewPosition[3*(jf*aspect_ratio+aspect_ratio-1)];   // other side of particle
		yb = NewPosition[3*(jf*aspect_ratio+aspect_ratio-1)+1];

		for( int kf=0; kf<count_fiber; kf++){
			int k = aspect_ratio * kf;
			if (k == j) continue;
			xc = NewPosition[3*k];
			yc = NewPosition[3*k+1];
			xd = NewPosition[3*(k + aspect_ratio - 1)];
			yd = NewPosition[3*(k + aspect_ratio - 1)+1];
			overlap_flag = check_cross(linklist->Lx,aspect_ratio,xa,ya,xb,yb,xc,yc,xd,yd);
			if(overlap_flag == 1) break;
		}

		// check whether edge particle is in wall region or not
		if( NewPosition[3*j+1] > border_u - eps || NewPosition[3*j+1] < border_d + eps ){
			overlap_flag = 1;
		}else if( NewPosition[3*( jf*aspect_ratio ) +1] > border_u - eps || NewPosition[3*( jf*aspect_ratio ) +1] < border_d + eps ){
			overlap_flag = 1;
		}
	}

	// calc potential
	if(overlap_flag == 0){
		makeBucket(NumberOfParticles,BucketFirst,BucketLast,Nextof,ParticleType,NewPosition,linklist);
		calNumberDensity2D(NumberOfParticles,Re_forNumberDensity,ParticleType,NewNumberDensity,NewPosition,BucketFirst,Nextof,linklist);
		// U_new = calPotential(NumberOfParticles,Re_forNumberDensity,ParticleType,NewNumberDensity,NewPosition,BucketFirst,Nextof,linklist);
		U_new = calPotential2(NumberOfParticles,NewNumberDensity,n0);
	}

	if(U_new < U_old){ // success
		memcpy(Position,NewPosition,size_posi);        // copy NewPosition -> Position
		memcpy(NumberDensity,NewNumberDensity,size_n); // copy NewNumberDensity -> NumberDensity
		mcstep ++;
		if(mcstep % 1000 == 0){
			printf("mcstep is %d, U_new : %lf \n",mcstep,U_new);
			output(n_fiber,aspect_ratio,NumberOfParticles,nx,ny,n_fluid,linklist,ParticleDistance,lambda_flag,ParticleType,Position,seed);
		}
		U_old = U_new;
	}else { // failure
		memcpy(NewPosition,Position,size_posi);        // copy Position -> NewPosition
		memcpy(NewNumberDensity,NumberDensity,size_n); // copy NumberDensity -> NewNumberDensity
		missstep ++;
	}

	// check fiber overlap

}

printf("mcstep : %d\n",mcstep);
printf("miss step : %d\n",missstep);
printf("U_new : %lf\n",U_new);

output(n_fiber,aspect_ratio,NumberOfParticles,nx,ny,n_fluid,linklist,ParticleDistance,lambda_flag,ParticleType,Position,seed);


	fclose(fp_input); /*fclose(fp_input2);*/
	free(ParticleType);	free(Position); // free(fiber_consist);
	free(Position_tmp);
	free(BucketFirst); free(BucketLast); free(Nextof);
	free(NumberDensity);

	// printf("end mk_particle\n");
	// printf("\n");
	return 0;
}


void output(int n_fiber, int aspect_ratio,int NumberOfParticles, int nx, int ny, int n_fluid, linklist_constant *linklist,
double ParticleDistance, int lambda_flag, int *ParticleType, double *Position,int seed){

	FILE *fp_output;
	char outputfile[256];
	int count = 0;
	int calc_lambda_particle;




	sprintf(outputfile,"%df%d-%d_%d_%d-%d_2d.dat",seed,n_fiber,aspect_ratio,NumberOfParticles,nx,ny);
	fp_output = fopen(outputfile, "w");
	// printf("output file : %s \n",outputfile);
	fprintf(fp_output,"%d\n",NumberOfParticles);
	fprintf(fp_output,"%d\n",n_fluid);
	fprintf(fp_output,"%d\n",n_fiber);
	fprintf(fp_output,"%d\n",aspect_ratio);
    fprintf(fp_output,"MIN_X : %lf, MAX_X : %lf\n",linklist->MIN_X,linklist->MAX_X);
    fprintf(fp_output,"MIN_Y : %lf, MAX_Y : %lf\n",linklist->MIN_Y,linklist->MAX_Y);
    fprintf(fp_output,"MIN_Z : %lf, MAX_Z : %lf\n",-8.0*ParticleDistance,8.0*ParticleDistance);

	count=0;
	// Firstly, fiber partiles are printed
	for(int i=0;i<n_fiber*aspect_ratio;i++){
		if(i == lambda_flag) calc_lambda_particle = count ;
		fprintf(fp_output,"%d %lf %lf %lf\n",ParticleType[i],Position[i*3],Position[i*3+1],Position[i*3+2]);
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

	if(count != NumberOfParticles){
		// printf("don't coresspond count with Number of particles \n");
	}

	// fprintf(fp_output,"lambda_particle : %d\n",calc_lambda_particle);
	fprintf(fp_output,"lambda_particle : -1\n");

	// output fiber consist
	for (int i = 0; i < aspect_ratio * n_fiber ; i++){
		fprintf(fp_output,"%d\n",i);
	}

	fclose(fp_output);

}

