#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
/*
*/

/* --------------------- UNCLUDE PARTS --------------------- */
#include"head3mg.c"
#include"load3mg.c"
#include"move3mg.c"
#include"mark3mg.c"
#include"gaus3mg.c"
#include"heat3mg.c"
#include"i3vis_output_paraview.c"
/* --------------------- UNCLUDE PARTS --------------------- */




/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0,f1,fln3;
long int pos0cur0,m0,m1,m2,m3,m4;
/**/
/**/
/**/
/* Load data from input file */
fln3=loadconf()+1;
/**/
/**/
/**/
/* Load from data file */
loader();
/**/
/**/
/**/
/* Reset iterations */
if(maxcyc<0) 
	{
	/* Clear vx,vy,vz,P */
	for(m2=0;m2<nodenum;m2++) 
		{
		vx[m2]=vy[m2]=vz[m2]=pr[m2]=0;
		exx[m2]=eyy[m2]=ezz[m2]=exy[m2]=exz[m2]=eyz[m2]=0;
		sxx[m2]=syy[m2]=szz[m2]=sxy[m2]=sxz[m2]=syz[m2]=0;
		}
	ronurecalc();
	maxcyc=-maxcyc;
	}
/**/
/**/
/**/
/* Output File Cycle */
f1=fln3+filesjob; if (f1>fl0num) f1=fl0num;
if (printmod) printf("\n PROCESSING UP TO %d .prn FILES (from %d to %d) FOR THIS JOB \n",f1-fln3,fln3,f1-1);
for (f0=fln3;f0<f1;f0++)
{
/* Reload Cur Output File Name */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/**/
/* Reload cyc0max_maxxystep_maxtkstep_maxtmstep */
cyc0max=fl0cyc[f0][0];
maxxystep=fl0stp[f0][0];
maxtkstep=fl0stp[f0][1];
maxtmstep=fl0stp[f0][2];
nubeg=fl0stp[f0][3];
nuend=fl0stp[f0][4];
p0koef=fl0stp[f0][5];
p1koef=fl0stp[f0][6];
p2koef=fl0stp[f0][7];
multinum=fl0cyc[f0][1];
/**/
/* General Cycle */
for (n0=0;n0<cyc0max;n0++)
	{
	if (printmod) printf("\n! FILE %s  KRUG ! %d\n",fl1out,n0+1);
	/**/
	/**/
	/**/
	/* Set initial time step */
	timestep=maxtmstep;
	if (printmod) printf("\n !!! MAX VALID TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7);
	/**/
	/**/
	/**/
	/* vX,vY, vZ recalc after Stokes+Contin equation */
	if(movemod)
		{
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		clockbeg=clock();
		if (printmod) printf("\n G, VX, VY, VZ, P ITERATION CYCLE BEGINING (clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
		vpiterate(n0,f0);
		clockend=clock();
		if (printmod) printf("\n G, VX, VY, VZ, P ITERATION CYCLE END (clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		}
	/**/
	/**/
	/**/
	/* Tk recalc after Heat conservation equation */
	if(tempmod && timestep)
		{
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		clockbeg=clock();
		if (printmod) printf("\n T RECALC...(clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
		titerate(n0);
		clockend=clock();
		if (printmod) printf("OK!(clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		}
	/**/
	/**/
	/**/
	/* Move marker */
	if(markmod && timestep)
		{
		clockbeg=clock();
		if (printmod) printf("MARKERS MOVE...(clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
		movemark();
		clockend=clock();
		if (printmod) printf("OK!(clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
		}
	/**/
	/**/
	/**/
	/* ro[],nu[] Recalc */
	if(gridmod)
		{
		clockbeg=clock();
		if (printmod) printf("\n RO, NU, CP etc  RECALC AFTER NEW MARKERS POSITIONS...(clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
		ronurecalc();
		clockend=clock();
		if (printmod) printf("OK!(clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
		}
	/**/
	/**/
	/**/
	/* Increse Timesum */
	timesum+=timestep;
	/* 1year=3.15576*10^7sek */
	if (printmod) printf("\n %e YEARS IN CYCLE     %e YEARS FROM START\n",timestep/3.15576e+7,timesum/3.15576e+7);
	/**/
	/**/
	/**/
	}
/**/
/**/
/**/
/* Print Results */
clockbeg=clock();
if (printmod) printf("\n SAVING RESULTS...(clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
saver(f0+1,n0-1);
if ((f0+1)%50==0) saver_paraview(f0+1,n0-1);
clockend=clock();
if (printmod) printf("OK!(clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
/* Print Results */
/**/
/* End Output file Names Cycle */
}
/* End Program */
return 0;
}
/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */



