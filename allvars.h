/********************************************************************
			       MACROS
********************************************************************/
//Macros generales
#define MAXNPARTICLES 	20000
#define NMAX1 		10
#define NMAX2 		200
#define NMAX3 		1000

//Constantes fisicas
#define GC		1

//Macros de parametros
#define N_1		0
#define N_2		1
#define N_3		2
#define T_STEP		3
#define T_MAX		4

#define TYP_INT		5

#define HOM_DIS		0
#define ARM_DIS		1
#define BAR_DIS		2
#define BULGE		3
#define HALO_ST		4
#define HALO_DK		5

#define COOR_1		12
#define COOR_2		13
#define COOR_3		14

#define VEL_1		15
#define VEL_2		16
#define VEL_3		17

#define MASS		18

//Potential parameters
#define HD_PAR		0
#define AM_PAR		5
#define BM_PAR		10
#define BL_PAR		15
#define SH_PAR		20
#define DH_PAR		25


//Macros de variables
#define X		0
#define Y		1
#define Z		2


/********************************************************************
			     ESTRUCTURAS
********************************************************************/

struct particle{
    double m;
    double r[3], v[3];
    
    double E;
    double Lz;
    
    int label;
    int type;};

struct cluster{
    int N1;
    int N2;
    int N3;
    int N;
    
    double t_snap;
    double t_step;
    double t_end;
    int n_end;
    
    double Ek, Ep;
    
    double pot_conf[NMAX2];
    char filename[NMAX2];
    
    int distros[6];
    
    struct particle prts[MAXNPARTICLES];};
    
struct intg_step{
    double r[3], v[3];};
    
   
/********************************************************************
			      HEADERS
********************************************************************/
#include <stdio.h>
#include <math.h>
#include <proto.h>
#include <stdlib.h>
#include <time.h>