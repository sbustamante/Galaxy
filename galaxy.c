#include <allvars.h>

// $./Galaxy.out <Evolution> <Video> <Analysis>
int main( int argc, char *argv[] )
{
    int i=0,j=0,n;
    double t=0;
    double prmt[NMAX2];
    char str[NMAX2];
    
    double zoom;
    char filename[NMAX2];
    FILE *script;

    printf( "\n\n***************** GALAXY SIMULATION *****************\n" );
    
    //Galaxy Initialization --------------------------------------------------------
    struct cluster galaxy[1]; 
    //Read parameters
    if( read_parameters( prmt, "parameters.conf")==1 ) return 0;
    galaxy[0].N1 = (int)prmt[N_1];
    galaxy[0].N2 = (int)prmt[N_2];
    galaxy[0].N3 = (int)prmt[N_3];
    galaxy[0].N = galaxy[0].N1 + galaxy[0].N2 + galaxy[0].N3;
    galaxy[0].distros[HOM_DIS] = (int)prmt[6];
    galaxy[0].distros[ARM_DIS] = (int)prmt[7];
    galaxy[0].distros[BAR_DIS] = (int)prmt[8];
    galaxy[0].distros[BULGE  ] = (int)prmt[9];
    galaxy[0].distros[HALO_ST] = (int)prmt[10];
    galaxy[0].distros[HALO_DK] = (int)prmt[11];
    galaxy[0].t_step = prmt[T_STEP];
    galaxy[0].t_end = prmt[T_MAX];
    galaxy[0].n_end = prmt[T_MAX]/prmt[T_STEP];
    sprintf( galaxy[0].filename, "%s", "Galaxy" );
    //Potentials parameters
    read_parameters( galaxy[0].pot_conf, "potential.conf" );

    //Setting initial conditions----------------------------------------------------
    initial_conditions( galaxy, prmt, "datos.in" );
    
    //Galaxy Evolution--------------------------------------------------------------
    if( atoi(argv[1])==1 ){
    printf("  * Evolución de la galaxia activada.\n");
    system("make cleandatas");
    while( t<prmt[T_MAX] )
    {
	//Current time
	galaxy[0].t_snap = t;

	//Integration step
	integration_step( galaxy, t, prmt[T_STEP], prmt[TYP_INT] );
	
	//Data printing
 	sprintf( filename, "./datas/%s-%05d", galaxy[0].filename, j );
 	data_out( galaxy, filename );
	
	//Video --------------------------
	if( atoi(argv[2])==1 ){
	    if( t==0 ) printf("  * Generación de video activada.\n");
	    //Gnuplot script
	    script = fopen( "script.gpl", "w" );
	    sprintf( str, "set terminal png\nset output '_tmp-%05d.png'\nset nokey\n", j );
	    fprintf( script, "%s", str );
	    sprintf( str, "set title 't = %lf'\n", t );
	    fprintf( script, "%s", str );
	    zoom = 2;
	    sprintf( str, "set xrange [%lf:%lf] \nset yrange [%lf:%lf]\n", -zoom,zoom,-zoom,zoom );
	    fprintf( script, "%s", str );
	    sprintf( str, "plot './datas/%s-%05d' u 2:3 w p 1",galaxy[0].filename, j );
	    fprintf( script, "%s", str );
	    fclose(script);
	    system("gnuplot script.gpl");}
	
	//System step
// 	printf( " t = %lf\n", t );
	t += prmt[T_STEP];
	j++;
    }
    //Deleting temporal files ------------------------------------------------------
    if( atoi(argv[2])==1 ){
	system("ffmpeg -f image2 -i _tmp-%05d.png  video.mpg");
// 	system("convert -delay .1 -loop 0 *.png animacion.gif");
	system("rm -rf *.png");
	system("rm -rf script.gpl");}}
	    
    return 0;
}