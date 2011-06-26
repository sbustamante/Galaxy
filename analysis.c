#include <allvars.h>

/*****************************************************************
 NAME:       analysis
 FUNCTION:   analysing all particles datas
 INPUTS:     galaxy, parameters and analysis configuration name
 RETURN:     0
 LOG:        31/05/2011    Creation
*****************************************************************/
int analysis( struct cluster galaxy[1],
	      double prmt[],
	      char filename[] )
{
    double orb[NMAX3][6];
    
    orbit( 1500, galaxy[0].filename, orb, galaxy[0].n_end, 1 );
    
    return 0;
}

/*****************************************************************
 NAME:       orbit
 FUNCTION:   calculate a particle orbit
 INPUTS:     particle label, galaxy filename, empty array, number
	     of datas, integer of fileout verification
 RETURN:     0
 LOG:        31/05/2011    Creation
*****************************************************************/
int orbit( int label,
	   char filename[],
	   double orb[NMAX3][6],
	   int N_t, int out )
{
    int i=0;
    struct particle prt;
    
    char filename_fold[NMAX2];
    FILE *file;
    
    for( i=0; i<N_t; i++ )
    {
	//Conversed to plain text
	sprintf( filename_fold, "./datas/%s-%05d", filename, i );
	file = fopen( filename_fold, "r" );
	prt.label = 0;
	//Read datas
	while( label != prt.label && getc( file ) != EOF ){
	    fscanf( file,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d", 
		    &prt.label, 
		    &prt.r[X], &prt.r[Y], &prt.r[Z],
		    &prt.v[X], &prt.v[Y], &prt.v[Z], 
		    &prt.m, &prt.type );}
	fclose( file );
	
	//Write datas
	orb[i][0] = prt.r[X];	orb[i][1] = prt.r[Y];	orb[i][2] = prt.r[Z];
	orb[i][3] = prt.v[X];	orb[i][4] = prt.v[Y];	orb[i][5] = prt.v[Z];	
    }
    //Data out
    if( out == 1 ){
	sprintf( filename_fold, "%05d.orb", label );
	file = fopen( filename_fold, "w" );
	for( i=0; i<N_t; i++ )
	    fprintf( file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i,
		     orb[i][0], orb[i][1], orb[i][2],
		     orb[i][3], orb[i][4], orb[i][5] );
    	fclose( file );}
    
    return 0;
}