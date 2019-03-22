// Jeu de la vie avec sauvegarde de quelques itérations
// compiler avec gcc -O3 -march=native (et -fopenmp si OpenMP souhaité)

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef __linux__
#include <openmpi-x86_64/mpi.h>
#include <memory.h>

#elif  __APPLE__
#include <mpi.h>
#endif

// hauteur et largeur de la matrice
#define HM 1200
#define LM 800

// nombre total d'itérations
#define ITER 10001
// multiple d'itérations à sauvegarder
#define SAUV 1000

#define DIFFTEMPS(a,b) \
(((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)

/* tableau de cellules */
typedef char Tab[HM][LM];

// initialisation du tableau de cellules
void init(Tab);

// calcule une nouveau tableau de cellules à partir de l'ancien
// - paramètres : ancien, nouveau
void calcnouv(Tab, Tab);
void calcnouv_l_i_f(Tab,Tab);

// variables globales : pas de débordement de pile
Tab t1, t2;
Tab tsauvegarde[1+ITER/SAUV];
void exchange_and_calculate(int i);
int size,rank,local_rows_num,next_n,previous_n;
int main(int argc,char**argv)
{
 // struct timeval tv_init, tv_beg, tv_end, tv_save;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  next_n = (rank == size-1)?0:rank+1;
  previous_n = (rank == 0)?size-1:rank-1;
  int* send_count = NULL;
  int *displacement = NULL;
  int reminder = HM%size;
  local_rows_num = HM/size;
  if(rank < reminder)
    local_rows_num += 1; //3,3,,3,2



  //gettimeofday( &tv_init, NULL);

  if(rank == 0){
    init(t2);
    /*
    send_count = malloc(sizeof(int)*size);
    displacement = malloc(sizeof(int)*size);
    int sum = 0;
    for (int i = 0; i < size ; ++i) {
      int tmp = (HM/local_rows_num);
      send_count[i] = (i != size-1) ? tmp :tmp+HM%size;
        if(reminder >0){
          send_count[i]++;
          reminder--;
        }
      send_count[i] *=LM;
      displacement[i] = sum;
      sum += send_count[i];
    }
*/
  }
//  if(rank == size-1){
//    memset(&t2[local_rows_num+1],0,LM);
//  }


 // MPI_Scatterv(t1,send_count,displacement,MPI_CHAR,&t1,local_rows_num,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Scatter(t2,local_rows_num*LM,MPI_CHAR,&t1[1],local_rows_num*LM,MPI_CHAR,0,MPI_COMM_WORLD);

for(int i=0;i<ITER;i++){

  exchange_and_calculate(i);
  if(i % SAUV == 0){
      if(rank == 0){
          printf("--------save(%d)-------\n",i);
         MPI_Gather(&t1[1],local_rows_num*LM,MPI_CHAR,&tsauvegarde[i/SAUV],local_rows_num*LM,MPI_CHAR,0,MPI_COMM_WORLD);

      }
      else
          MPI_Gather(&t1[1],local_rows_num*LM,MPI_CHAR,NULL,local_rows_num*LM,MPI_CHAR,0,MPI_COMM_WORLD);


  }
}

if(rank == 0){
  FILE *f = fopen("test/para.out", "w");
  for(int i=0 ; i<ITER ; i+=SAUV)
  {
    fprintf(f, "------------------ sauvegarde %d ------------------\n", i);
    for(int x=0 ; x<HM ; x++)
    {
      for(int y=0 ; y<LM ; y++)
        fprintf(f, tsauvegarde[i/SAUV][x][y]?"*":" ");
      fprintf(f, "\n");
    }
  }
  fclose(f);
}
/*
    char file[20];
  sprintf(file,"file.%d",rank);
  FILE *f = fopen(file, "w");
  fprintf(f, "------------------ sauvegarde %d ------------------\n", rank);

    for(int x=0 ; x<LM ; x++)
      fprintf(f, t1[0][x]?"*":"-");

    fprintf(f, "\n==========================================\n");
    for(int x=1 ; x<=local_rows_num ; x++)
    {
      for(int y=0 ; y<LM ; y++)
        fprintf(f, t1[x][y]?"*":"-");
      fprintf(f, "\n");
    }
  fprintf(f, "==========================================\n");
  for(int x=0 ; x<LM ; x++)
    fprintf(f, t1[local_rows_num+1][x]?"*":"-");
  fclose(f);

  */

/*
  //gettimeofday( &tv_beg, NULL);
  for(int i=0 ; i<ITER ; i++)
  {
    if( i%2 == 0)
      calcnouv(t1, t2);
    else
      calcnouv(t2, t1);

    if(i%SAUV == 0)
    {
      printf("sauvegarde (%d)\n", i);
      // copie t1 dans tsauvegarde[i/SAUV]
      for(int x=0 ; x<HM ; x++)
        for(int y=0 ; y<LM ; y++)
          tsauvegarde[i/SAUV][x][y] = t1[x][y];
    }
  }
  //gettimeofday( &tv_end, NULL);


  //gettimeofday( &tv_save, NULL);

  //printf("init : %lf s,", DIFFTEMPS(tv_init, tv_beg));
  //printf(" calcul : %lf s,", DIFFTEMPS(tv_beg, tv_end));
  //printf(" sauvegarde : %lf s\n", DIFFTEMPS(tv_end, tv_save));

  */



  free(send_count);
  free(displacement);
  MPI_Finalize();
  return( 0 );
}

void init(Tab t)
{

 // srand(time(0));
#pragma omp prallel for
  for(int i=0 ; i<HM ; i++)
    for(int j=0 ; j<LM ; j++ )
    {
      // t[i][j] = rand()%2;
      // t[i][j] = ((i+j)%3==0)?1:0;
      // t[i][j] = (i==0||j==0||i==h-1||j==l-1)?0:1;
      t[i][j] = 0;
    }

  t[10][10] = 1;
  t[10][11] = 1;
  t[10][12] = 1;
  t[9][12] = 1;
  t[8][11] = 1;

  t[55][50] = 1;
  t[54][51] = 1;
  t[54][52] = 1;
  t[55][53] = 1;
  t[56][50] = 1;
  t[56][51] = 1;
  t[56][52] = 1;
   /*
  t[0][0] = 1;
  t[1][1] = 1;
  t[2][2] = 1;
  t[3][3] = 1;
  t[4][4] = 1;
  t[5][5] = 1;
  t[6][6] = 1;
  t[7][0] = 1;
  t[8][1] = 1;
  t[9][2] = 1;
  t[10][3] = 1;
  t[11][7] = 1;
*/

}

int nbvois(Tab t, int i, int j)
{
  int n=0;
  if( i>0 )
  {  /* i-1 */
    if( j>0 )
      if( t[i-1][j-1] )
        n++;
    if( t[i-1][j] )
        n++;
    if( j<LM-1 )
      if( t[i-1][j+1] )
        n++;
  }
  if( j>0 )
    if( t[i][j-1] )
      n++;
  if( j<LM-1 )
    if( t[i][j+1] )
      n++;
  if( i<HM-1 )
  {  /* i+1 */
    if( j>0 )
      if( t[i+1][j-1] )
        n++;
    if( t[i+1][j] )
        n++;
    if( j<LM-1 )
      if( t[i+1][j+1] )
        n++;
  }
  return( n );
}

void calcnouv(Tab t, Tab n)
{
#pragma omp parallel for
  for(int i=2 ; i < local_rows_num ; i++)
    for(int j=0 ; j<LM ; j++)
    {
      int v = nbvois(t, i, j);
      if(v==3)
        n[i][j] = 1;
      else if(v==2)
        n[i][j] = t[i][j];
      else
        n[i][j] = 0;
    }

}

void calcnouv_l_i_f(Tab t,Tab n){
#pragma omp prallel for
  for(int j=0 ; j<LM ; j++)
  {
    int v = nbvois(t,1, j);
    if(v==3)
      n[1][j] = 1;
    else if(v==2)
      n[1][j] = t[1][j];
    else
      n[1][j] = 0;
  }
#pragma omp parallel for
  for(int j=0 ; j<LM ; j++)
  {
    int v = nbvois(t,local_rows_num, j);
    if(v==3)
      n[local_rows_num][j] = 1;
    else if(v==2)
      n[local_rows_num][j] = t[local_rows_num][j];
    else
      n[local_rows_num][j] = 0;
  }
}

void exchange_and_calculate(int i){
    MPI_Request mpi_send_req1,mpi_send_req2,mpi_recv_req1,mpi_recv_req2;
    if(i%2 == 0){
        if(rank != 0)
            MPI_Isend(&t1[1],LM,MPI_CHAR,previous_n,0,MPI_COMM_WORLD,&mpi_send_req1);
        if(rank != size-1)
            MPI_Isend(&t1[local_rows_num],LM,MPI_CHAR,next_n,0,MPI_COMM_WORLD,&mpi_send_req2);
        if(rank != 0)
            MPI_Irecv(&t1[0],LM,MPI_CHAR,previous_n,0,MPI_COMM_WORLD,&mpi_recv_req1);
        if(rank != size-1)
            MPI_Irecv(&t1[local_rows_num+1],LM,MPI_CHAR,next_n,0,MPI_COMM_WORLD,&mpi_recv_req2);
        calcnouv(t1,t2);

    }else{
        if(rank != 0)
            MPI_Isend(&t2[1],LM,MPI_CHAR,previous_n,0,MPI_COMM_WORLD,&mpi_send_req1);
        if(rank != size-1)
            MPI_Isend(&t2[local_rows_num],LM,MPI_CHAR,next_n,0,MPI_COMM_WORLD,&mpi_send_req2);
        if(rank != 0)
            MPI_Irecv(&t2[0],LM,MPI_CHAR,previous_n,0,MPI_COMM_WORLD,&mpi_recv_req1);
        if(rank != size-1)
            MPI_Irecv(&t2[local_rows_num+1],LM,MPI_CHAR,next_n,0,MPI_COMM_WORLD,&mpi_recv_req2);
        calcnouv(t2,t1);
    }


    if(rank != 0)
        MPI_Wait(&mpi_send_req1,MPI_STATUS_IGNORE);
    if(rank != size-1)
        MPI_Wait(&mpi_send_req2,MPI_STATUS_IGNORE);
    if(rank != 0)
        MPI_Wait(&mpi_recv_req1,MPI_STATUS_IGNORE);
    if(rank != size-1)
        MPI_Wait(&mpi_recv_req2,MPI_STATUS_IGNORE);

    //have previous next...
    //calculate first and last line
    if(i%2 == 0)
    calcnouv_l_i_f(t1,t2);
    else
    calcnouv_l_i_f(t2,t1);



}