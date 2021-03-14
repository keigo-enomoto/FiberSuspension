#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <string.h>
//#include <time.h>

int GetRandom(int min,int max);

//* beginning of main *************************************************************
int main( int argc, char** argv ){

  int i=0,j;
  double rc1,rc2;

/*
  // You can change random seed using this srand(seed)
  srand(256);
  // labnote vol.2 p16
  for(i = 0; i < 5; i++){
    rc1 = (double)rand()/RAND_MAX ;   //0 < rc < 1.0
    rc2 = (double)rand()/RAND_MAX ;
    printf("%lf %lf \n",rc1,rc2);
  }

  printf("\n");
  // You can change random seed using this srand(seed)
  // srand(256);
  // labnote vol.2 p16
  for(i = 0; i < 5; i++){
    rc1 = (double)rand()/RAND_MAX ;   //0 < rc < 1.0
    rc2 = (double)rand()/RAND_MAX ;
    printf("%lf %lf \n",rc1,rc2);
  }
*/

int array[300] = {0};
srand(256);
for(int i=0; i<100000000; i++){
  j = GetRandom(0,299);
  array[j] += 1;
}

for(int i=0; i<300; i++){
  int n = (int)(array[i]/10000);
  printf("%d ",i);
  for(int j=0; j<n; j++){
  printf("*");
  }
  printf("\n");
}



return 0;
}
//* end of main *************************************************************


int GetRandom(int min,int max)
{
	return min + (int)(rand()*(max-min+1.0)/(1.0+RAND_MAX));
}