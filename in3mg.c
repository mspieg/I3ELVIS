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
#include"mark3mg.c"
#include"move3mg.c"
#include"gaus3mg.c"
#include"heat3mg.c"
#include"i3vis_output_paraview.c"
/* --------------------- UNCLUDE PARTS --------------------- */

FILE *fl1;



/* Formation of Data file for i3mg.c program */
int main()
{
/* Counters */
int n1,n2,n3,mgi,ntype,ynprn=0;
long int m1,m2,m3,m4,m5,m6,m7,mcmax1;
long int mm1,mm2,mm3,wn1[500];
/* Nonstability for markers, m */
int nonstab;
double xnonstab=0,ynonstab=0,znonstab;
/* Distance val */
double x,y,z,dx,dy,dz;
/* Initial Rock type for Grid,Rock bodyes Num, Initial Distribution */
double bx[8],by[8],bz[8];
double k[6][4];
double cnst,koef,koef1,koef2,ival,ival1,eps1[500],wi1[500],mex,mee;
/* T K in Edges of Grid */
double t[8];
long int m10,m11,m20,m21,m30,m31,m1min,m1max,m2min,m2max,m3min,m3max;
long int nshiftx,nshifty,nshiftz;
long int nshiftx1,nshifty1,nshiftz1;
long int nshiftx2,nshifty2,nshiftz2;
long int un1[MAXPOS];
double ui1[MAXPOS];
/**/
/**/
/**/
/* Load configuration from mode.t3c */
loadconf();
/**/
/**/
/**/
/* stop.yn file creation */
fl = fopen("stop.yn","wt");
fprintf(fl,"n \n");
fclose(fl);
/**/
/**/
/**/
/* Open File init.t3c */
fl1 = fopen("init.t3c","rt");
printf("Formation of Initial Condition ...\n");
/**/
/**/
/**/
/* Grid Parameters */
ffscanf1();xnumx=atoi(sa);
ffscanf1();ynumy=atoi(sa);
ffscanf1();znumz=atoi(sa);
ffscanf1();mnumx=atoi(sa);
ffscanf1();mnumy=atoi(sa);
ffscanf1();mnumz=atoi(sa);
ffscanf1();xsize=atof(sa);
ffscanf1();ysize=atof(sa);
ffscanf1();zsize=atof(sa);
ffscanf1();pxinit=atoi(sa); 
ffscanf1();pyinit=atoi(sa); 
ffscanf1();pzinit=atoi(sa); 
ffscanf1();pinit=atof(sa);
ffscanf1();GXKOEF=atof(sa);
ffscanf1();GYKOEF=atof(sa);
ffscanf1();GZKOEF=atof(sa);
ffscanf1();timesum=atof(sa)*3.15576e+7;
ffscanf1();nonstab=atoi(sa);
/* Random Nonstability Read */
ffscanf1();xnonstab=atof(sa)*xsize/((double)(xnumx-1))/((double)(mnumx));
ffscanf1();ynonstab=atof(sa)*ysize/((double)(ynumy-1))/((double)(mnumy));
ffscanf1();znonstab=atof(sa)*zsize/((double)(znumz-1))/((double)(mnumz));
/**/
/* First Calc,Check Grid parameters */
marknum=(xnumx-1)*(ynumy-1)*(znumz-1)*mnumx*mnumy*mnumz;
gridcheck();
/**/
/**/
/**/
/* Set initial gridlines positions */
gx[0]=0;
for (m1=1;m1<xnumx;m1++) gx[m1]=gx[m1-1]+xstpx;
gy[0]=0;
for (m2=1;m2<ynumy;m2++) gy[m2]=gy[m2-1]+ystpy;
gz[0]=0;
for (m3=1;m3<znumz;m3++) gz[m3]=gz[m3-1]+zstpz;
/**/
/*
for (m1=0;m1<xnumx;m1++) {printf("%ld %e",m1,gx[m1]);getchar();}
for (m2=0;m2<ynumy;m2++) {printf("%ld %e",m2,gy[m2]);getchar();}
for (m3=0;m3<znumz;m3++) {printf("%ld %e",m3,gz[m3]);getchar();}
*/
/**/
/**/
/* Load Data from File name */
ffscanf1();
if (sa[0]!='0')
	{
	/* Data from file Y(1)/N(0) */
	ynprn=1; 
	/* Reload File Name */
	for (n1=0;n1<50;n1++) fl1in[n1]=sa[n1];
	ffscanf1();fl1itp=0; if(sa[0] == 'b') fl1itp=1;
	/* Load data from input file */
	loader();
	/**/
	/* Clear BC */
	#pragma omp parallel for shared(bondm) private(m1) schedule(static)
	for(m1=0;m1<MAXPAR;m1++) bondm[m1]=0;
	}
else
	{
	/* Markers Location Set */
	printf("Set markers location \n");
	dx=xstpx/(double)(mnumx);
	dy=ystpy/(double)(mnumy);
	dz=zstpz/(double)(mnumz);
	/* Marker Types add */
	for (m3=0;m3<znumz1;m3++)
	for (m1=0;m1<xnumx1;m1++)
	for (m2=0;m2<ynumy1;m2++)
		{
		/* Cell adres */
		m4=m3*xynumxy1+m1*ynumy1+m2;
		/**/
		/* First marker in cell num */
		mm2=m4*mnumx*mnumy*mnumz;
		/**/
		/* Markers in Cell cycle */
		z=(double)(m3)*zstpz+dz/2.0;
		for (m5=0;m5<mnumz;m5++)
			{
			x=(double)(m1)*xstpx+dx/2.0;
			for (m6=0;m6<mnumx;m6++)
				{
				y=(double)(m2)*ystpy+dy/2.0;
				for (m7=0;m7<mnumy;m7++)
					{
					/* Marker num */
					mm1=mm2+m5*mnumx*mnumy+m6*mnumy+m7;
					/**/
					/* X,Y */
					markx[mm1]=x;
					marky[mm1]=y;
					markz[mm1]=z;
					/**/
					/* Add Y */
					y+=dy;
					}
				/* Add X */
				x+=dx;
				}
			/* Add Z */
			z+=dz;
			}
		}
	/* Nonstability Set on markers */
	if(nonstab)
		{
		printf("Set nonstability on markers\n");
		for (mm1=0;mm1<marknum;mm1++)
			{
			markx[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(xnonstab);
			marky[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(ynonstab);
			markz[mm1]+=(float)(rand() % (nonstab*2+1) - nonstab)/((float)(nonstab))*(float)(znonstab);
			}
		}
	}
/**/
/**/
/**/
/**/
/**/
/**/
/* Output File name */
ffscanf1();
for (n1=0;n1<15;n1++) fl1out[n1]=sa[n1];
ffscanf1(); if(sa[0] == 'b') fl1otp=1;
/**/
/**/
/**/
/* Rock Properties Read */
printf("Read rocks properties \n");
ffscanf1();
while (sa[0]!='~')
	{
	n1=atoi(sa);
	ffscanf1();markn0[n1]=atof(sa);
	ffscanf1();markn1[n1]=atof(sa);
	ffscanf1();marks0[n1]=atof(sa);
	ffscanf1();marks1[n1]=atof(sa);
	ffscanf1();marknu[n1]=atof(sa);
	ffscanf1();markdh[n1]=atof(sa);
	ffscanf1();markdv[n1]=atof(sa);
	ffscanf1();markss[n1]=atof(sa);
	ffscanf1();markmm[n1]=atof(sa);
	ffscanf1();markll[n1]=atof(sa);
	ffscanf1();marka0[n1]=atof(sa);
	ffscanf1();marka1[n1]=atof(sa);
	ffscanf1();markb0[n1]=atof(sa);
	ffscanf1();markb1[n1]=atof(sa);
	ffscanf1();marke0[n1]=atof(sa);
	ffscanf1();marke1[n1]=atof(sa);
	ffscanf1();markro[n1]=atof(sa);
	ffscanf1();markbb[n1]=atof(sa);
	ffscanf1();markaa[n1]=atof(sa);
	ffscanf1();markcp[n1]=atof(sa);
	ffscanf1();markkt[n1]=atof(sa);
	ffscanf1();markkf[n1]=atof(sa);
	ffscanf1();markkp[n1]=atof(sa);
	/* Ht in Wt/kg load Y(k)/N() */
	ffscanf1(); if(sa[0]=='k') markht[n1]=atof(sa+1)*markro[n1]; else markht[n1]=atof(sa);
	ffscanf1();
	rocknum=MAXV(rocknum,n1+1);
	}
/**/
/**/
/**/
/**/
/**/
/**/
/**/
/* Bondary Conditions Read,Set to Multi Grid */
bondnum=0;
printf("Set bondary conditions\n");
ffscanf1();
while (sa[0]!='~')
	{
	/* VAR for bondary conditions: Vx,Vy,T  */
	n1=-1;
	if(sa[0] == 'P') n1=0;
	if(sa[1] == 'x') n1=1;
	if(sa[1] == 'y') n1=2;
	if(sa[1] == 'z') n1=3;
	if(sa[0] == 'T') 
		{
		n1=4;
		n3=0;
		if(sa[1] == 'c' || sa[2] == 'c') n3=1; 
		}
	if(sa[0] == 'X' || sa[1] == 'X') n1=5;
	if(sa[0] == 'Y' || sa[1] == 'Y') n1=6;
	if(sa[0] == 'Z' || sa[1] == 'Y') n1=7;
	if(sa[0] == 'M' || sa[1] == 'M') {n1=8; printf("Set new Cell markers location \n");}
	if(n1 == -1) {printf("Unknown Parameter <%s>",sa); exit(0);}
	/**/
	/* Bond Parameters Read for cur Box */
	/* VAR___m10___m11___m20___m21___m30___m31___Const___Koef___dm1___dm2___dm3 */
	/* Vx    1     x-1   0     0     z     z     0       1.0    0     +1    0   */
	ffscanf1();if (sa[0] == 'x') m10=xnumx-1+atoi(sa+1); else m10=atoi(sa);
	ffscanf1();if (sa[0] == 'x') m11=xnumx-1+atoi(sa+1); else m11=atoi(sa);
	ffscanf1();if (sa[0] == 'y') m20=ynumy-1+atoi(sa+1); else m20=atoi(sa);
	ffscanf1();if (sa[0] == 'y') m21=ynumy-1+atoi(sa+1); else m21=atoi(sa);
	ffscanf1();if (sa[0] == 'z') m30=znumz-1+atoi(sa+1); else m30=atoi(sa);
	ffscanf1();if (sa[0] == 'z') m31=znumz-1+atoi(sa+1); else m31=atoi(sa);
	ffscanf1();cnst=atof(sa);
	ffscanf1();koef=atof(sa);
	nshiftx=nshifty=nshiftz=0;
	if(koef)
		{
		ffscanf1();nshiftx=atoi(sa);
		ffscanf1();nshifty=atoi(sa);
		ffscanf1();nshiftz=atoi(sa);
		}
	ffscanf1();koef1=atof(sa);
	nshiftx1=nshifty1=nshiftz1=0;
	if(koef1)
		{
		ffscanf1();nshiftx1=atoi(sa);
		ffscanf1();nshifty1=atoi(sa);
		ffscanf1();nshiftz1=atoi(sa);
		}
	ffscanf1();koef2=atof(sa);
	nshiftx2=nshifty2=nshiftz2=0;
	if(koef2)
		{
		ffscanf1();nshiftx2=atoi(sa);
		ffscanf1();nshifty2=atoi(sa);
		ffscanf1();nshiftz2=atoi(sa);
		}
	/**/
	/* Next VAR for bondary conditions read */
	ffscanf1();
	/**/
	/**/
	/**/
	/* Conditions set for Grid */
	/* m10,m11 - X min,max Num */
	/* m20,m21 - Y min,max Num */
	/* m30,m31 - Z min,max Num */
	/* n1 - Variable: Vx(0), Vy(1), T(2) */
	/* nshiftx,nshifty,nshiftz - Shifts for Symmetry conditions */
	/**/
	if (n1<5)
	{
	for (mgi=0;mgi<=multinum;mgi++)
		{
		/* Limits for cur grid */
		m1min=(m10-2)/mgs[mgi]+2;
		if(m10<3) m1min=m10;
		if(m10>xnumx-4) m1min=mgx[mgi]-(xnumx-m10);
		m1max=(m11-2)/mgs[mgi]+2;
		if(m11<3) m1max=m11;
		if(m11>xnumx-4) m1max=mgx[mgi]-(xnumx-m11);
		m2min=(m20-2)/mgs[mgi]+2;
		if(m20<3) m2min=m20;
		if(m20>ynumy-4) m2min=mgy[mgi]-(ynumy-m20);
		m2max=(m21-2)/mgs[mgi]+2;
		if(m21<3) m2max=m21;
		if(m21>ynumy-4) m2max=mgy[mgi]-(ynumy-m21);
		m3min=(m30-2)/mgs[mgi]+2;
		if(m30<3) m3min=m30;
		if(m30>znumz-4) m3min=mgz[mgi]-(znumz-m30);
		m3max=(m31-2)/mgs[mgi]+2;
		if(m31<3) m3max=m31;
		if(m31>znumz-4) m3max=mgz[mgi]-(znumz-m31);
/*
printf("%d    %d %d  %ld %ld   %ld %ld     %ld %ld   %ld %ld    %ld %ld   %ld %ld ",n1,multinum,mgi,m10,m11,m1min,m1max,m20,m21,m2min,m2max,m30,m31,m3min,m3max);getchar();     
*/
		/* Grid cycle */
		for (m3=m3min;m3<=m3max;m3++)
		for (m1=m1min;m1<=m1max;m1++)
		for (m2=m2min;m2<=m2max;m2++)
			{
			/* Node num for cur grid */
			m4=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			/* Shifted Node num for cur grid */
			m6=m4+nshiftz*mgxy[mgi]+nshiftx*mgy[mgi]+nshifty;
			/**/
			switch(n1)
				{
				/* P,Vx,Vy,Vz  Bondary equations */
				case 0:
				case 1:
				case 2:
				case 3:
				m5=m4*4+n1;
				bondnum++; bondm[m5]=bondnum;
				bondv1[bondnum][0]=cnst;
				if(nshiftx || nshifty || nshiftz) 
					{
					bondv1[bondnum][1]=koef;
					bondn1[bondnum]=m6*4+n1+1;
					}
				else
					{
					bondv1[bondnum][1]=0;
					bondn1[bondnum]=0;
					}
				break;
				/**/
				/* T bondary equations */
				case 4:
				/* Check distance to the center */
				n2=1;
				/* Cylindrical */
				if(n3==1 && mgi==0)
					{
					n2=0;
					dx=(gx[m1]-(gx[m11]+gx[m10])/2.0)/((gx[m11]-gx[m10])/2.0);
					dz=(gz[m3]-(gz[m31]+gz[m30])/2.0)/((gz[m31]-gz[m30])/2.0);
					ival=pow(dx*dx+dz*dz,0.5);
					if(ival<=1.0) n2=1;
					}
				if(n2==1) 
					{
/*
if(n3==1) {printf("%ld %ldi %ld %e",m1,m2,m3,cnst);getchar();}
*/
					m5=nodenum8+m4;
					bondnum++; bondm[m5]=bondnum;
					bondv1[bondnum][0]=cnst;
					if(nshiftx || nshifty || nshiftz) 
						{
						bondv1[bondnum][1]=koef;
						bondn1[bondnum]=m6+1;
						}
					else
						{
						bondv1[bondnum][1]=0;
						bondn1[bondnum]=0;
						}
					}
				break;
				/**/
				}
			}
		}
	}
	else
	{
	/* Grid cycle */
	for (m3=m30;m3<=m31;m3++)
	for (m1=m10;m1<=m11;m1++)
	for (m2=m20;m2<=m21;m2++)
		{
		switch(n1)
			{
			/**/
			/**/
			/* X coordinates of greedlines definition */
			case 5:
			ival=(double)(m11-m10);
			if(!ival) ival=1.0;
			ival1=(double)(m1-m10);
			gx[m1]=gx[m1-1]+cnst+(koef-cnst)*ival1/ival;
			if(koef1) gx[m1]=gx[m1-1]+exp(log(cnst)+log(koef1/cnst)*ival1/ival);
/*
printf("X  %e %e %e %e %e %e", ival,ival1,cnst,koef,gx[m1-1],gx[m1]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* Y coordinates of greedlines definition */
			case 6:
			ival=(double)(m21-m20);
			if(!ival) ival=1.0;
			ival1=(double)(m2-m20);
			gy[m2]=gy[m2-1]+cnst+(koef-cnst)*ival1/ival;
			if(koef1) gy[m2]=gy[m2-1]+exp(log(cnst)+log(koef1/cnst)*ival1/ival);
/*
printf("Y %ld  %e %e %e %e %e %e %e",m2,ival,ival1,cnst,koef,koef1,gy[m2-1],gy[m2]); getchar();
printf("Y  %e %e %e %e %e %e", ival,ival1,cnst,koef,gy[m2-1],gy[m2]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* Z coordinates of greedlines definition */
			case 7:
			ival=(double)(m31-m30);
			if(!ival) ival=1.0;
			ival1=(double)(m3-m30);
			gz[m3]=gz[m3-1]+cnst+(koef-cnst)*ival1/ival;
			if(koef1) gz[m3]=gz[m3-1]+exp(log(cnst)+log(koef1/cnst)*ival1/ival);
/*
printf("Y  %e %e %e %e %e %e", ival,ival1,cnst,koef,gy[m2-1],gy[m2]); getchar();
*/
			break;
			/**/
			/**/
			/**/
			/* MARKER grid set to cells */
			case 8:
			dx=(gx[m1+1]-gx[m1])/(double)(nshiftx);
			dy=(gy[m2+1]-gy[m2])/(double)(nshifty1);
			dz=(gz[m3+1]-gz[m3])/(double)(nshiftz2);
			/* Marker Types add */
			z=gz[m3]+dz/2.0;
			for (mm3=0;mm3<nshiftz2;mm3++)
				{
				x=gx[m1]+dx/2.0;
				for (mm1=0;mm1<nshiftx;mm1++)
					{
					y=gy[m2]+dy/2.0;
					for (mm2=0;mm2<nshifty1;mm2++)
						{
						/* X,Y */
						markx[marknum]=(float)(x);
						marky[marknum]=(float)(y);
						markz[marknum]=(float)(z);
/*
printf("%ld %ld %ld   %ld %ld %ld   %ld  %e %e %e   %e %e %e ", m1,m2,m3,mm1,mm2,mm3,marknum,koef,koef1,koef2,markx[marknum],marky[marknum],markz[marknum]);getchar();
*/
						/* Random Nonstability Set on markers */
						if(cnst>0 && koef>0)
							{
							/* Random nonstability Set on X */
							markx[marknum]+=(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dx*koef);
/*
printf("X  %ld %ld %ld %ld %e %e %e ", m1,m2,mm1,mm2,dx,dy,(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dx*koef));getchar();
*/
							}
						if(cnst>0 && koef1>0)
							{
							/* Random nonstability Set on Y */
							marky[marknum]+=(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dy*koef1);
/*
printf("Y  %ld %ld %ld %ld %e %e %e ", m1,m2,mm1,mm2,dx,dy,(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dy*koef1));getchar();
*/
							}
						if(cnst>0 && koef2>0)
							{
							/* Random nonstability Set on Y */
							markz[marknum]+=(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dz*koef2);
/*
printf("Y  %ld %ld %ld %ld %e %e %e ", m1,m2,mm1,mm2,dx,dy,(float)(rand() % ((int)(cnst)*2+1) - (int)(cnst))/((float)(cnst))*(float)(dy*koef1));getchar();
*/
							}
/*
printf("%ld %ld %ld   %ld %ld %ld   %ld   %e %e %e ", m1,m2,m3,mm1,mm2,mm3,marknum,markx[marknum],marky[marknum],markz[marknum]);getchar();
*/
						/* Add mark Num */
						marknum++;
						/* Add Y */
						y+=dy;
						}
					/* Add X */
					x+=dx;
					}
				/* Add Z */
				z+=dz;
				}
			break;
			/**/
			/**/
			/**/
			}
		}
	}
	}
/* Restore grid */
gx[0]=gy[0]=gz[0]=0;
gx[xnumx-1]=xsize;
gy[ynumy-1]=ysize;
gz[znumz-1]=zsize;
/* Set bondary Conditions for Grid ------------------ */
/* Count Conditions for Grid ------------------ */
bondnum=0;
for (m1=0;m1<nodenum*10;m1++)
if(bondm[m1])
	{
	bondnum++;
	}
/**/
/**/
/**/
/* Rock bodies set on marker */
printf("Set rocks bodies as markers types \n");
ffscanf1();
while (sa[0]!='~')
	{
	/* Rock type */
	if(sa[0]=='s' || sa[0]=='c') 
		{
		/* Ellipsoid */
		n1=atoi(sa+1);
		ntype=1;
		/* Cylinder */
		if(sa[0]=='c') ntype=2; 
		}
	else
		{
		/* Rectangle */
		n1=atoi(sa);
		ntype=0;
		}
	/* Melt extraction, strain */
	ffscanf1(); mex=atof(sa);
	ffscanf1(); mee=atof(sa);
	/* Read coordinats of Cur boxe */
	/* Z=0     Z=1 */
	/* 0  2    4  6 */
	/* 1  3    5  7 */
	for(n2=0;n2<8;n2++)
		{
		ffscanf1(); if(sa[0]=='m') bx[n2]=atof(sa+1)/xsize; else bx[n2]=atof(sa);
		ffscanf1(); if(sa[0]=='m') by[n2]=atof(sa+1)/ysize; else by[n2]=atof(sa);
		ffscanf1(); if(sa[0]=='m') bz[n2]=atof(sa+1)/zsize; else bz[n2]=atof(sa);
		}
	/**/
	/* Next Rock type read */
	ffscanf1();
	/**/
	/* Making a rectangle */
	if(ntype==0)
		{
		/* Reactangles Bondary Equations Koef  calc */
		/* Left,Right:  A0+A1y+A2z+A3yz=x */
		/* Up,Down:     A0+A1x+A2z+A3xz=y */
		/* Front,Back:  A0+A1y+A2x+A3yx=z */
		/* Z=0     Z=1 */
		/* 0  2    4  6 */
		/* 1  3    5  7 */
		/**/
		/* Left bond - after Edges 0,1,4,5: A0+A1y+A2z+A3yz=x*/
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=by[0]; mat[ 2]=bz[0]; mat[ 3]=by[0]*bz[0]; mat[ 4]=bx[0];
		mat[ 5]=1.0; mat[ 6]=by[1]; mat[ 7]=bz[1]; mat[ 8]=by[1]*bz[1]; mat[ 9]=bx[1];
		mat[10]=1.0; mat[11]=by[4]; mat[12]=bz[4]; mat[13]=by[4]*bz[4]; mat[14]=bx[4];
		mat[15]=1.0; mat[16]=by[5]; mat[17]=bz[5]; mat[18]=by[5]*bz[5]; mat[19]=bx[5];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Left bond Save */
		for(n2=0;n2<4;n2++) k[0][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Right bond - after Edges 3,4,7,8: A0+A1y+A2z+A3yz=x */
		/*         Z=0                         Z=1             */
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=by[2]; mat[ 2]=bz[2]; mat[ 3]=by[2]*bz[2]; mat[ 4]=bx[2];
		mat[ 5]=1.0; mat[ 6]=by[3]; mat[ 7]=bz[3]; mat[ 8]=by[3]*bz[3]; mat[ 9]=bx[3];
		mat[10]=1.0; mat[11]=by[6]; mat[12]=bz[6]; mat[13]=by[6]*bz[6]; mat[14]=bx[6];
		mat[15]=1.0; mat[16]=by[7]; mat[17]=bz[7]; mat[18]=by[7]*bz[7]; mat[19]=bx[7];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Right bond Save */
		for(n2=0;n2<4;n2++) k[1][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Up bond - after Edges 1,3,5,7: A0+A1x+A2z+A3xz=y */
		/*         Z=0                         Z=1             */
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=bx[0]; mat[ 2]=bz[0]; mat[ 3]=bx[0]*bz[0]; mat[ 4]=by[0];
		mat[ 5]=1.0; mat[ 6]=bx[2]; mat[ 7]=bz[2]; mat[ 8]=bx[2]*bz[2]; mat[ 9]=by[2];
		mat[10]=1.0; mat[11]=bx[4]; mat[12]=bz[4]; mat[13]=bx[4]*bz[4]; mat[14]=by[4];
		mat[15]=1.0; mat[16]=bx[6]; mat[17]=bz[6]; mat[18]=bx[6]*bz[6]; mat[19]=by[6];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Up bond Save */
		for(n2=0;n2<4;n2++) k[2][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Down bond - after Edges 2,4,6,8: A0+A1x+A2z+A3xz=y */
		/*         Z=0                         Z=1             */
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=bx[1]; mat[ 2]=bz[1]; mat[ 3]=bx[1]*bz[1]; mat[ 4]=by[1];
		mat[ 5]=1.0; mat[ 6]=bx[3]; mat[ 7]=bz[3]; mat[ 8]=bx[3]*bz[3]; mat[ 9]=by[3];
		mat[10]=1.0; mat[11]=bx[5]; mat[12]=bz[5]; mat[13]=bx[5]*bz[5]; mat[14]=by[5];
		mat[15]=1.0; mat[16]=bx[7]; mat[17]=bz[7]; mat[18]=bx[7]*bz[7]; mat[19]=by[7];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Down bond Save */
		for(n2=0;n2<4;n2++) k[3][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Front bond - after Edges 1,2,3,4: A0+A1y+A2x+A3yx=z */
		/*         Z=0                         Z=1             */
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=by[0]; mat[ 2]=bx[0]; mat[ 3]=by[0]*bx[0]; mat[ 4]=bz[0];
		mat[ 5]=1.0; mat[ 6]=by[1]; mat[ 7]=bx[1]; mat[ 8]=by[1]*bx[1]; mat[ 9]=bz[1];
		mat[10]=1.0; mat[11]=by[2]; mat[12]=bx[2]; mat[13]=by[2]*bx[2]; mat[14]=bz[2];
		mat[15]=1.0; mat[16]=by[3]; mat[17]=bx[3]; mat[18]=by[3]*bx[3]; mat[19]=bz[3];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Front bond Save */
		for(n2=0;n2<4;n2++) k[4][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Back bond - after Edges 5,6,7,8: A0+A1y+A2x+A3yx=z */
		/*         Z=0                         Z=1             */
		/* 1(0,1,2)  3(6, 7, 8)       5(12,13,14)  7(18,19,20) */
		/* 2(3,4,5)  4(9,10,11)       6(15,16,17)  8(21,22,23) */
		/* Clear Matrix */
		matclear(4,mat);
		/* Add Matrix */
		mat[ 0]=1.0; mat[ 1]=by[4]; mat[ 2]=bx[4]; mat[ 3]=by[4]*bx[4]; mat[ 4]=bz[4];
		mat[ 5]=1.0; mat[ 6]=by[5]; mat[ 7]=bx[5]; mat[ 8]=by[5]*bx[5]; mat[ 9]=bz[5];
		mat[10]=1.0; mat[11]=by[6]; mat[12]=bx[6]; mat[13]=by[6]*bx[6]; mat[14]=bz[6];
		mat[15]=1.0; mat[16]=by[7]; mat[17]=bx[7]; mat[18]=by[7]*bx[7]; mat[19]=bz[7];
		/* Solve Matrix - Koef Calc */
		gausmat(4,mat,sol);
		/* Back bond Save */
		for(n2=0;n2<4;n2++) k[5][n2]=sol[n2];
		/**/
		/**/
		/**/
		/* Check Markers for cur reactangle bondary */
		for (mm1=0;mm1<marknum;mm1++)
			{
			/* Relativ coordinats Calc */
			x=markx[mm1]/xsize;
			y=marky[mm1]/ysize;
			z=markz[mm1]/zsize;
			/**/
			/* Left,Right:  A0+A1y+A2z+A3yz=x */
			if(x>=k[0][0]+k[0][1]*y+k[0][2]*z+k[0][3]*y*z && x<=k[1][0]+k[1][1]*y+k[1][2]*z+k[1][3]*y*z)
			/* Up,Down:     A0+A1x+A2z+A3xz=y */
			if(y>=k[2][0]+k[2][1]*x+k[2][2]*z+k[2][3]*x*z && y<=k[3][0]+k[3][1]*x+k[3][2]*z+k[3][3]*x*z)
			/* Front,Back:  A0+A1y+A2x+A3yx=z */
			if(z>=k[4][0]+k[4][1]*y+k[4][2]*x+k[4][3]*y*x && z<=k[5][0]+k[5][1]*y+k[5][2]*x+k[5][3]*y*x)
				{
				/* Rock type set */
				markt[mm1]=n1;
				markex[mm1]=mex;
				marke[mm1]=mee;
				/* Immobile markers */
				if(n1>=100)
					{
					markk[mm1]=markx[mm1];
					markd[mm1]=marky[mm1];
					markw[mm1]=markz[mm1];
					}
				else
					{
					markk[mm1]=0;
					markd[mm1]=0;
					markw[mm1]=0;
					}
				}
			}
		}
	else        
		{
		switch(ntype)
			{
			/* Ellipsoid */
			case 1:
			/* Check Markers for cur reactangle bondary */
			for (mm1=0;mm1<marknum;mm1++)
				{
				/* Relativ coordinats Calc */
				x=markx[mm1]/xsize;
				y=marky[mm1]/ysize;
				z=markz[mm1]/zsize;
				/* Check distance to the first center */
				dx=(x-bx[0])/bx[7];	
				dy=(y-by[0])/by[7];	
				dz=(z-bz[0])/bz[7];	
				ival=pow(dx*dx+dy*dy+dz*dz,0.5);
/*
printf("%e %e %e    %e %e %e   %e %e %e",x,y,z,dx,dy,dz,ival,bx[2],bx[5]);getchar();
*/
				if(ival>=bx[2] && ival<=bx[5])
					{
					/* Rock type set */
					markt[mm1]=n1;
					markex[mm1]=mex;
					marke[mm1]=mee;
					}
				}
			break;
			/**/
			case 2:
			/* Cylindrical */
			/* Check Markers for cur reactangle bondary */
			for (mm1=0;mm1<marknum;mm1++)
				{
				/* Relativ coordinats Calc */
				x=markx[mm1]/xsize;
				y=marky[mm1]/ysize;
				z=markz[mm1]/zsize;
				/* Check distance to the first center */
				dx=(x-bx[0])/bx[7];	
				dy=ABSV((y-by[0])/by[7]);	
				dz=(z-bz[0])/bz[7];	
				ival=pow(dx*dx+dz*dz,0.5);
/*
printf("%e %e %e    %e %e %e   %e %e %e",x,y,z,dx,dy,dz,ival,bx[2],bx[5]);getchar();
*/
				if(ival>=bx[2] && ival<=bx[5] && dy<=1.0)
					{
					/* Rock type set */
					markt[mm1]=n1;
					markex[mm1]=mex;
					marke[mm1]=mee;
					}
				}
			break;
			}
		}
	}
/**/
/**/
/**/
/* Initial distribution of Themperature Read, Set to Grid */
printf("Read distribution of temperature \n");
ffscanf1();
while (sa[0]!='~')
	{
	/* Read Type of the box */
	/* 0 - simple box */
	/* 1 - Box with cooling age */
	ntype=atoi(sa);
	/* Read diagonal coordinats of Cur box */
	/* Z=0     Z=1 */
	/* 0  2    4  6 */
	/* 1  3    5  7 */
	ffscanf1(); if(sa[0]=='m') bx[0]=atof(sa+1); else bx[0]=atof(sa)*xsize;
	ffscanf1(); if(sa[0]=='m') bx[2]=atof(sa+1); else bx[2]=atof(sa)*xsize;
	ffscanf1(); if(sa[0]=='m') by[0]=atof(sa+1); else by[0]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[1]=atof(sa+1); else by[1]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[2]=atof(sa+1); else by[2]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[3]=atof(sa+1); else by[3]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') bz[0]=atof(sa+1); else bz[0]=atof(sa)*zsize;
	ffscanf1(); if(sa[0]=='m') bx[5]=atof(sa+1); else bx[5]=atof(sa)*xsize;
	ffscanf1(); if(sa[0]=='m') bx[7]=atof(sa+1); else bx[7]=atof(sa)*xsize;
	ffscanf1(); if(sa[0]=='m') by[4]=atof(sa+1); else by[4]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[5]=atof(sa+1); else by[5]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[6]=atof(sa+1); else by[6]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') by[7]=atof(sa+1); else by[7]=atof(sa)*ysize;
	ffscanf1(); if(sa[0]=='m') bz[7]=atof(sa+1); else bz[7]=atof(sa)*zsize;
	/**/
	/* Read T K in Edges of Cur boxe */
	for(n1=0;n1<8;n1++)
		{
		ffscanf1();t[n1]=atof(sa);
		}
	/* Read x0 for next T-box */
	ffscanf1();
	/**/
	/* Initial distribution of Temperature set for cur box */
	for (m3=0;m3<znumz;m3++)
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		/* Node Num */
		m4=xynumxy*m3+m1*ynumy+m2;
		/**/
		if(ntype!=2 && ntype!=3)
			{
			/* Relativ dist in Grid */
			z=(gz[m3]-bz[0])/(bz[7]-bz[0]);
			x=(gx[m1]-(bx[0]*(1.0-z)+bx[5]*z))/((bx[2]-bx[0])*(1.0-z)+(bx[7]-bx[5])*z);
			y=(gy[m2]-(by[0]*(1.0-x)*(1.0-z)+by[2]*x*(1.0-z)+by[4]*(1.0-x)*z+by[6]*x*z))/((by[1]-by[0])*(1.0-x)*(1.0-z)+(by[3]-by[2])*x*(1.0-z)+(by[5]-by[4])*(1.0-x)*z+(by[7]-by[6])*x*z);
			/**/
			if(x>=0 && x<=1.0 && y>=0 && y<=1.0 && z>=0 && z<=1.0)
			switch(ntype)
				{
				/* Simple box */
				case 0:
				/* T recalc from eight T values */
				/* Z=0     Z=1 */
				/* 0  2    4  6 */
				/* 1  3    5  7 */
				tk[m4] =(1.0-x)*(1.0-y)*(1.0-z)*t[0];
				tk[m4]+=(1.0-x)*(    y)*(1.0-z)*t[1];
				tk[m4]+=(    x)*(1.0-y)*(1.0-z)*t[2];
				tk[m4]+=(    x)*(    y)*(1.0-z)*t[3];
				tk[m4]+=(1.0-x)*(1.0-y)*(    z)*t[4];
				tk[m4]+=(1.0-x)*(    y)*(    z)*t[5];
				tk[m4]+=(    x)*(1.0-y)*(    z)*t[6];
				tk[m4]+=(    x)*(    y)*(    z)*t[7];
/*
if(tk[m4]<272.0) {printf("AA %ld %ld %ld %e",m1,m2,m3,tk[m4]);getchar();}
*/
				break;
				/**/
				/* 1,4 = Cooling age */
				case 1:
				case 4:
				/* Age interpolation */
				/* Z\X Left  Right */
				/* Back  t4    t5  */
				/* Front t2    t3  */
				/* Z\X Left  Right */
				ival=t[2]*(1.0-x)*(1.0-z)+t[3]*x*(1.0-z)+t[4]*(1.0-x)*z+t[5]*x*z;
				/* Turcotte & Schubert, 2002: Cooling age T=T1-erfc()(T1-T0) */
				ival=(gy[m2]-(by[0]*(1.0-x)*(1.0-z)+by[2]*x*(1.0-z)+by[4]*(1.0-x)*z+by[6]*x*z))/2.0/pow(t[6]*ival*3.15576e+7,0.5);
				tk[m4]=t[1]-(1.0-erf(ival))*(t[1]-t[0])+t[7]*(gy[m2]-by[0]);
				break;
				/**/
				/* 5,6 = Transitional geotherms */
				case 5:
				case 6:
				/* Interpolate T on 4 sides of the box */
				tkcalc1(bx[0],gy[m2],bz[0],eps1,wn1,wi1); t[0]=eps1[2];
				tkcalc1(bx[2],gy[m2],bz[0],eps1,wn1,wi1); t[1]=eps1[2];
				tkcalc1(bx[5],gy[m2],bz[7],eps1,wn1,wi1); t[2]=eps1[2];
				tkcalc1(bx[7],gy[m2],bz[7],eps1,wn1,wi1); t[3]=eps1[2];
				tk[m4]=t[0]*(1.0-x)*(1.0-z)+t[1]*x*(1.0-z)+t[2]*(1.0-x)*z+t[3]*x*z;
/*
{printf("BB %ld %ld %ld %e",m1,m2,m3,tk[m4]);getchar();}
if(tk[m4]<272.0) {printf("BB %ld %ld %ld %e",m1,m2,m3,tk[m4]);getchar();}
if(bx[0]>5e+2){printf("%ld %ld %ld %e %e",m1,m2,m3,gy[m2]-by[0],tk[m4]);getchar();}
*/
				break;
				/**/
				}
			}
		else
			{
			switch(ntype)
				{
				case 2:
				/* Ellipsoid */
				/* Check distance to the center */
				dx=(gx[m1]-bx[0])/bx[7];	
				dy=(gy[m2]-by[0])/by[7];	
				dz=(gz[m3]-bz[0])/bz[7];	
				ival=pow(dx*dx+dy*dy+dz*dz,0.5);
/*
printf("%e %e %e    %e %e %e   %e %e %e",gx[m1],gy[m2],gz[m3],dx,dy,dz,ival,bx[2]/xsize,bx[5]/xsize);getchar();
*/
				if(ival>=bx[2]/xsize && ival<=bx[5]/xsize)
					{
					tk[m4]=t[0]+(ival-bx[2]/xsize)/(bx[5]/xsize-bx[2]/xsize)*(t[1]-t[0]);
/*
printf("%ld %ld %ld \n %e %e %e \n %e %e %e \n %e %e %e \n %e %e %e \n %e %e %e",m1,m2,m3,gx[m1],gy[m2],gz[m3],bx[0],by[0],bz[0],dx,dy,dz,ival,bx[2],by[2],tk[m4],t[0],t[1]);getchar();
*/
					}
				break;
				/**/
				case 3:
				/* Cylindrical */
				/* Check distance to the center */
				dx=(gx[m1]-bx[0])/bx[7];	
				dy=ABSV((gy[m2]-by[0])/by[7]);	
				dz=(gz[m3]-bz[0])/bz[7];	
				ival=pow(dx*dx+dz*dz,0.5);
/*
printf("%e %e %e    %e %e %e   %e %e %e",gx[m1],gy[m2],gz[m3],dx,dy,dz,ival,bx[2]/xsize,bx[5]/xsize);getchar();
*/
				if(ival>=bx[2]/xsize && ival<=bx[5]/xsize && dy<=1.0)
					{
					tk[m4]=t[0]+(ival-bx[2]/xsize)/(bx[5]/xsize-bx[2]/xsize)*(t[1]-t[0]);
/*
printf("%ld %ld %ld \n %e %e %e \n %e %e %e \n %e %e %e \n %e %e %e \n %e %e %e",m1,m2,m3,gx[m1],gy[m2],gz[m3],bx[0],by[0],bz[0],dx,dy,dz,ival,bx[2],by[2],tk[m4],t[0],t[1]);getchar();
*/
					}
				break;
				/**/
				}
			}
		}
	}
/**/
/**/
/**/
/* Close init.t3c file */
fclose(fl1);
/**/
/**/
/**/
/* Reset T for markers */
for (m3=0;m3<marknum;m3++) 
/* Check markers out of grid */
if (markt[m3]<100 && markx[m3]>0 && marky[m3]>0 && markz[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && (double)(markz[m3])<zsize)
	{
	/* Interpolate temperature */
	tkcalc1((double)(markx[m3]),(double)(marky[m3]),(double)(markz[m3]),eps1,wn1,wi1);
	/**/
	/* Reset marker temperature for newly coming markers */
	markk[m3]=(float)(eps1[2]);
/*
if(markk[m3]<272.0) {printf("EEE %ld %d  %e",m3,markt[m3],markk[m3]);getchar();}
{printf("%ld %e",m3,markk[m3]);getchar();}
if(markk[m3]<990.0){printf("%ld %e",m3,markk[m3]);getchar();}
*/
	}
/**/
/**/
/*
for (m3=0;m3<znumz;m3++)
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	m4=xynumxy*m3+m1*ynumy+m2;
	if(tk[m4]<272.0) {printf("CCC %ld %ld %ld %e",m1,m2,m3,tk[m4]);getchar();}
	}
*/
/**/
/**/
/**/
if(ynprn==0)
{
/* ro[],nu[], cp[] kt[] for nodes Recalc after markers */
printf("Recalc nodes properties after markers types \n");
ronurecalc();
/**/
/**/
/**/
/*
for (m3=0;m3<znumz;m3++)
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	m4=xynumxy*m3+m1*ynumy+m2;
	if(tk[m4]<272.0) {printf("DDD %ld %ld %ld %e",m1,m2,m3,tk[m4]);getchar();}
	}
*/
/**/
/**/
/**/
/* Set Initial pressure distribution */
if(GYKOEF>0)
for (m3=1;m3<znumz;m3++)
for (m1=1;m1<xnumx;m1++)
	{
	m4=xynumxy*m3+m1*ynumy+1;
	pr[m4]=pinit;
	for (m2=2;m2<ynumy;m2++)
		{
		m4++;
		pr[m4]=pr[m4-1]+(gy[m2]-gy[m2-2])/2.0*GYKOEF*(ro[m4-1]+ro[m4-ynumy-1]+ro[m4-xynumxy-1]+ro[m4-xynumxy-ynumy-1])/4.0;
		}
	}
/*
*/
/**/
/**/
/**/
/* ro[],nu[], cp[] kt[] for nodes Recalc after markers */
printf("Second recalc nodes properties after markers types \n");
ronurecalc();
/*
*/
}
/**/
/**/
/**/
/* Second Calc,Check Grid parameters */
gridcheck();
/**/
/**/
/**/
/*
for (m3=1;m3<ynumy;m3++) {printf("%ld %e",m3,gy[m3]); getchar();}
*/
/* Save information in output file */
printf("Saving information to output file \n");
printmod=1;
saver_paraview(0,0);
saver(0,0);
/**/
/**/
/**/
printf("OK!\n");
return 0;
}
/* Formation of Data file for i2.c program */
