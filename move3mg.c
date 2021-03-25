/* Left side or Value for Sxx  Equation */ 
/* Sxx=2NU(dVx/dX) */
double sxxcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
/* mgi - cur multigrid level */
{
/* Val buf */
double ival,leftsxx=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
	{
	nueff=nxx[v[7]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSxx7, SIGxx7 */
/**/
/* Return Sxx,Exx val ----------------------------*/
if(ynval==0)
	{
	/* Exx=dVx/dX */
	dvxdx=(vx[v[2]]-vx[v[0]])*xkf;
	leftsxx=dvxdx;
	/**/
	/* Save Exx */
	eps1[0]=leftsxx;
	/**/
	/* Calc Sxx=2Nu*Exx */
	leftsxx=2.0*nueff*leftsxx;
	/**/
	return leftsxx;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSxx7, SIGxx7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxx=2NU(dVx/dX) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=-ynval*nueff*2.0*xkf;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=+ynval*nueff*2.0*xkf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=0;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=0;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=0;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=0;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exx %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxx  Equation */ 



/* Left side or Value for Syy  Equation */ 
/* Syy=2NU(dVy/dY) */
double syycalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsyy=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
	{
	nueff=nxx[v[7]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSyy7, SIGyy7 */
/**/
/* Return Syy,Eyy val ----------------------------*/
if(ynval==0)
	{
	/* Eyy=dVy/dY */
	dvydy=(vy[v[1]]-vy[v[0]])*ykf;
	leftsyy=dvydy;
	/* Save Eyy */
	eps1[0]=leftsyy;
	/**/
	/* Calc Syy=2Nu*Eyy */
	leftsyy=2.0*nueff*leftsyy;
	/**/
	return leftsyy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syy ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Syy=2NU(dVy/dY-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=0;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=0;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=-ynval*nueff*2.0*ykf;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=+ynval*nueff*2.0*ykf;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=0;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=0;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Ezz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Syy  Equation */ 



/* Left side or Value for Szz  Equation */ 
/* Szz=2NU(dVz/dZ) */
double szzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftszz=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
	{
	nueff=nxx[v[7]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/* Return Szz,Ezz val ----------------------------*/
if(ynval==0)
	{
	/* Ezz=dVz/dZ */
	dvzdz=(vz[v[4]]-vz[v[0]])*zkf;
	leftszz=dvzdz;
	/**/
	/* Save Ezz */
	eps1[0]=leftszz;
	/**/
	/* Calc Szz=2Nu*Ezz */
	leftszz=2.0*nueff*leftszz;
	/**/
	return leftszz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Szz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Szz=2NU(dVz/dZ) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=0;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=0;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=0;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=0;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=-ynval*nueff*2.0*zkf;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=+ynval*nueff*2.0*zkf;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Ezz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Szz  Equation */ 





/* Left side or Value for Sxy  Equation */ 
/* Sxy=NU(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsxy=0,nueff,dvxdy,dvydx;
long int v[8];
/* Distances */
double xkf,ykf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+m3*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>3 && m3>2 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
	{
	nueff=nxy[v[3]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[3]]+nu[v[7]])/2.0;
	if(viscmod==1) nueff=exp((log(nu[v[3]])+log(nu[v[7]]))/2.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[3]]+1.0/nu[v[7]])/2.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx3 */
/* Vy1 | Vy3 */
/* EPSxy3, SIGxy3 */
/**/
/* Return Sxy,Exy val ----------------------------*/
if(ynval==0)
	{
	/* Exy=1/2(dVx/dY+dVy/dX) */
	dvxdy=(vx[v[3]]-vx[v[2]])*ykf;
	dvydx=(vy[v[3]]-vy[v[1]])*xkf;
	leftsxy=dvxdy+dvydx;
	/**/
	/* Save Exy */
	eps1[0]=leftsxy;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsxy=2.0*nueff*leftsxy;
	/**/
	return leftsxy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx3 */
/* Vy1 | Vy3 */
/* EPSxy3, SIGxy3 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxy=NU(dVx/dY+dVy/dX) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[2]*4+1;
ui1[un1[0]+1]=-ynval*nueff*2.0*ykf;
un1[un1[0]+2]=v[3]*4+1;
ui1[un1[0]+2]=+ynval*nueff*2.0*ykf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[1]*4+2;
ui1[un1[0]+3]=-ynval*nueff*2.0*xkf;
un1[un1[0]+4]=v[3]*4+2;
ui1[un1[0]+4]=+ynval*nueff*2.0*xkf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exy %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxy  Equation */ 



/* Left side or Value for Sxz  Equation */ 
/* Sxz=NU(dVx/dZ+dVz/dX) */
double sxzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsxz=0,nueff,dvxdz,dvzdx;
long int v[8];
/* Distances */
double xkf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+m2;v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>2 && m3>3 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
	{
	nueff=nxz[v[6]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[6]]+nu[v[7]])/2.0;
	if(viscmod==1) nueff=exp((log(nu[v[6]])+log(nu[v[7]]))/2.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[6]]+1.0/nu[v[7]])/2.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx6 */
/* Vz4 / Vz6 */
/* EPSxz6, SIGxz6 */
/**/
/* Return Sxz,Exz val ----------------------------*/
if(ynval==0)
	{
	/* Exz=1/2(dVx/dZ+dVz/dX) */
	dvxdz=(vx[v[6]]-vx[v[2]])*zkf;
	dvzdx=(vz[v[6]]-vz[v[4]])*xkf;
	leftsxz=dvxdz+dvzdx;
	/**/
	/* Save Exy */
	eps1[0]=leftsxz;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsxz=2.0*nueff*leftsxz;
	/**/
	return leftsxz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx6 */
/* Vz4 / Vz6 */
/* EPSxz6, SIGxz6 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxz=NU(dVx/dZ+dVz/dX) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[2]*4+1;
ui1[un1[0]+1]=-ynval*nueff*2.0*zkf;
un1[un1[0]+2]=v[6]*4+1;
ui1[un1[0]+2]=+ynval*nueff*2.0*zkf;
/* Add Vz with koefficients */
un1[un1[0]+3]=v[4]*4+3;
ui1[un1[0]+3]=-ynval*nueff*2.0*xkf;
un1[un1[0]+4]=v[6]*4+3;
ui1[un1[0]+4]=+ynval*nueff*2.0*xkf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxz  Equation */ 



/* Left side or Value for Syz  Equation */ 
/* Syz=NU(dVy/dZ+dVz/dY) */
double syzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsyz=0,nueff,dvydz,dvzdy;
long int v[8];
/* Distances */
double ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+m1*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0 && m1>2 && m2>3 && m3>3 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
	{
	nueff=nyz[v[5]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[5]]+nu[v[7]])/2.0;
	if(viscmod==1) nueff=exp((log(nu[v[5]])+log(nu[v[7]]))/2.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[5]]+1.0/nu[v[7]])/2.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vy1 | Vy5 */
/* Vz4 / Vz5 */
/* EPSyz5, SIGyz5 */
/**/
/* Return Syz,Eyz val ----------------------------*/
if(ynval==0)
	{
	/* Eyz=1/2(dVy/dZ+dVz/dY) */
	dvydz=(vy[v[5]]-vy[v[1]])*zkf;
	dvzdy=(vz[v[5]]-vz[v[4]])*ykf;
	leftsyz=dvydz+dvzdy;
	/**/
	/* Save Exy */
	eps1[0]=leftsyz;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsyz=2.0*nueff*leftsyz;
	/**/
	return leftsyz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vy1 | Vy5 */
/* Vz4 / Vz5 */
/* EPSyz5, SIGyz5 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Syz=NU(dVy/dZ+dVz/dY) */
/* Add Vy with koefficients */
un1[un1[0]+1]=v[1]*4+2;
ui1[un1[0]+1]=-ynval*nueff*2.0*zkf;
un1[un1[0]+2]=v[5]*4+2;
ui1[un1[0]+2]=+ynval*nueff*2.0*zkf;
/* Add Vz with koefficients */
un1[un1[0]+3]=v[4]*4+3;
ui1[un1[0]+3]=-ynval*nueff*2.0*ykf;
un1[un1[0]+4]=v[5]*4+3;
ui1[un1[0]+4]=+ynval*nueff*2.0*ykf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Eyz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Syz  Equation */ 



/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stoks equation initial form */
/* dSIGxx/dX + dSIGxy/dY + dSIGxz/dZ - dP/dX = - RO*Gx */
double xstokserr(long int m1, long int m2, long int m3, int ynerr, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[8];
/* Err Buf */
double ival,leftx,rightx,nueff;
/* Distances */
double xkf,ykf,zkf,skf,skf1,GXKOEFX=GXKOEF;
/**/
/**/
/**/
/* Distances Calc */
xkf=2.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2]);
zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
rightx  = -GXKOEFX*(ro[v[0]]+ro[v[1]]+ro[v[4]]+ro[v[5]])/4.0;
/**/
/**/
/**/
/* Solution koef calc */
skf=(1.0/xkf+1.0/ykf+1.0/zkf)/3.0;
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P5 - P7 */
/* SIGxx5 - SIGxx7 */
/* SIGxy0 | SIGxy1 */
/* SIGxz0 / SIGxz4 */
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* Effective NU calc */
	if(viscmod==0) nueff=(nu[v[0]-mgy[mgi]]+nu[v[1]-mgy[mgi]]+nu[v[4]-mgy[mgi]]+nu[v[5]-mgy[mgi]]+10.0*(nu[v[0]]+nu[v[1]]+nu[v[4]]+nu[v[5]])+nu[v[2]]+nu[v[3]]+nu[v[6]]+nu[v[7]])/48.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]-mgy[mgi]])+log(nu[v[1]-mgy[mgi]])+log(nu[v[4]-mgy[mgi]])+log(nu[v[5]-mgy[mgi]])+10.0*(log(nu[v[0]])+log(nu[v[1]])+log(nu[v[4]])+log(nu[v[5]]))+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[6]])+log(nu[v[7]]))/48.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]-mgy[mgi]]+1.0/nu[v[1]-mgy[mgi]]+1.0/nu[v[4]-mgy[mgi]]+1.0/nu[v[5]-mgy[mgi]]+10.0*(1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[4]]+1.0/nu[v[5]])+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[6]]+1.0/nu[v[7]])/48.0);
/*
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[4]]+nu[v[5]])/4.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[4]])+log(nu[v[5]]))/4.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[4]]+1.0/nu[v[5]])/4.0);
*/
	if(nueff>nuend) nueff=nuend;
	if(nueff<nubeg) nueff=nubeg;
	skf1=1.0/nueff*(1.0/xkf/xkf+1.0/ykf/ykf+1.0/zkf/zkf)/3.0;
	/**/
	/**/
	/**/
	/* LEFT part of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY + dSIGxz/dZ - dP/dX = - RO*Gx */
	/**/
	/* dSIGxx/dX */
	leftx =(sxx[v[7]]-sxx[v[5]])*xkf;
	/* dSIGxy/dY */
	leftx+=(sxy[v[1]]-sxy[v[0]])*ykf;
	/* dSIGxz/dZ */
	leftx+=(sxz[v[4]]-sxz[v[0]])*zkf;
	/* -dP/dX */
	leftx-=(pr[v[7]]-pr[v[5]])*xkf;
	/**/
/*
printf("Vx1 %ld %ld %ld %e %e %e",m1,m2,m3,leftx,rightx,leftx-rightx);
*/
	/* X-STOKS Error */
	leftx=(leftx-rightx);
	/**/
	eps1[0]=leftx*skf1;
	return leftx;
	}
/**/
/**/
/* Save Right part Save for X-Stokes ---------------------*/
ui1[0]=rightx;
/* Set Initial Num of lines */
un1[0]=17;
/**/
/* Add Coefficients for left parts of X-Stokes ----------------*/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P5 - P7 */
/* SIGxx5 - SIGxx7 */
/* SIGxy0 | SIGxy1 */
/* SIGxz0 / SIGxz4 */
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* -dP/dX */
un1[16]=v[5]*4+0;
ui1[16]=+xkf;
un1[17]=v[7]*4+0;
ui1[17]=-xkf;
/* vx  4  7   vy 8  10  vz  13  15 p 16  17  */
/* 2   1   3     9  11    12  14             */
/*  6  5                                     */
/* dSIGxx/dX */
sxxcalc(m1  ,m2+1,m3+1,-xkf,mgi,eps1,un1,ui1);
un1[ 2]=un1[18]; ui1[ 2] =ui1[18];
un1[ 1]=un1[19]; ui1[ 1] =ui1[19];
un1[ 8]=un1[20]; ui1[ 8] =ui1[20];
un1[ 9]=un1[21]; ui1[ 9] =ui1[21];
un1[12]=un1[22]; ui1[12] =ui1[22];
un1[13]=un1[23]; ui1[13] =ui1[23];
un1[0]=17;
sxxcalc(m1+1,m2+1,m3+1, xkf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[18];
un1[ 3]=un1[19]; ui1[ 3] =ui1[19];
un1[10]=un1[20]; ui1[10] =ui1[20];
un1[11]=un1[21]; ui1[11] =ui1[21];
un1[14]=un1[22]; ui1[14] =ui1[22];
un1[15]=un1[23]; ui1[15] =ui1[23];
/* dSIGxy/dY */
un1[0]=17;
sxycalc(m1  ,m2  ,m3  ,-ykf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[19];
un1[ 4]=un1[18]; ui1[ 4] =ui1[18];
                 ui1[ 8]+=ui1[20];
                 ui1[10]+=ui1[21];
un1[0]=17;
sxycalc(m1  ,m2+1,m3  , ykf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[18];
un1[ 5]=un1[19]; ui1[ 5] =ui1[19];
                 ui1[ 9]+=ui1[20];
                 ui1[11]+=ui1[21];
/* dSIGxz/dZ */
un1[0]=17;
sxzcalc(m1  ,m2  ,m3  ,-zkf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[19];
un1[ 6]=un1[18]; ui1[ 6] =ui1[18];
                 ui1[12]+=ui1[20];
                 ui1[14]+=ui1[21];
un1[0]=17;
sxzcalc(m1  ,m2  ,m3+1, zkf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[18];
un1[ 7]=un1[19]; ui1[ 7] =ui1[19];
                 ui1[13]+=ui1[20];
                 ui1[15]+=ui1[21];
/**/
/* vx  4  7   vy 8  10  vz  13  15 p 16  17  */
/* 2   1   3     9  11    12  14             */
/*  6  5                                     */
/* dDIV(v)/dX */
if(p2koef!=0)
	{
un1[0]=17;
	continadd(m1  ,m2+1,m3+1,-xkf*p2koef,mgi,eps1,un1,ui1);
ui1[ 2]+=ui1[18];
ui1[ 1]+=ui1[19];
ui1[ 8]+=ui1[20];
ui1[ 9]+=ui1[21];
ui1[12]+=ui1[22];
ui1[13]+=ui1[23];
un1[0]=17;
	continadd(m1+1,m2+1,m3+1, xkf*p2koef,mgi,eps1,un1,ui1);
ui1[ 1]+=ui1[18];
ui1[ 3]+=ui1[19];
ui1[10]+=ui1[20];
ui1[11]+=ui1[21];
ui1[14]+=ui1[22];
ui1[15]+=ui1[23];
	}
un1[0]=17;
/**/
/**/
/**/
/*
for(n1=0;n1<=un[0];n1++)printf("Vx %e %ld \n",ui[n1],un[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */


/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stoks equation initial form */
/* dSIGyy/dY + dSIGxy/dX + dSIGyz/dZ - dP/dY = - RO*Gy */
double ystokserr(long int m1, long int m2, long int m3, int ynerr, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[8];
/* Err Buf */
double ival,lefty,righty,nueff;
/* Distances */
double xkf,ykf,zkf,skf,skf1,GYKOEFY=GYKOEF;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1]);
ykf=2.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* RIGHT parts of Y-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
righty  = -GYKOEFY*(ro[v[0]]+ro[v[2]]+ro[v[4]]+ro[v[6]])/4.0;
/**/
/**/
/**/
/* Solution koef calc */
skf=(1.0/xkf+1.0/ykf+1.0/zkf)/3.0;
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P6 | P7 */
/* SIGyy6 | SIGyy7 */
/* SIGxy0 - SIGxy2 */
/* SIGyz0 / SIGyz4 */
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* Effective NU calc */
	if(viscmod==0) nueff=(nu[v[0]-1]+nu[v[2]-1]+nu[v[4]-1]+nu[v[6]-1]+10.0*(nu[v[0]]+nu[v[2]]+nu[v[4]]+nu[v[6]])+nu[v[1]]+nu[v[3]]+nu[v[5]]+nu[v[7]])/48.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]-1])+log(nu[v[2]-1])+log(nu[v[4]-1])+log(nu[v[6]-1])+10.0*(log(nu[v[0]])+log(nu[v[2]])+log(nu[v[4]])+log(nu[v[6]]))+log(nu[v[1]])+log(nu[v[3]])+log(nu[v[5]])+log(nu[v[7]]))/48.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]-1]+1.0/nu[v[2]-1]+1.0/nu[v[4]-1]+1.0/nu[v[6]-1]+10.0*(1.0/nu[v[0]]+1.0/nu[v[2]]+1.0/nu[v[4]]+1.0/nu[v[6]])+1.0/nu[v[1]]+1.0/nu[v[3]]+1.0/nu[v[5]]+1.0/nu[v[7]])/48.0);
/*
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[2]]+nu[v[4]]+nu[v[6]])/4.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[2]])+log(nu[v[4]])+log(nu[v[6]]))/4.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[2]]+1.0/nu[v[4]]+1.0/nu[v[6]])/4.0);
*/
	if(nueff>nuend) nueff=nuend;
	if(nueff<nubeg) nueff=nubeg;
	skf1=1.0/nueff*(1.0/xkf/xkf+1.0/ykf/ykf+1.0/zkf/zkf)/3.0;
	/**/
	/**/
	/**/
	/* LEFT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX + dSIGyz/dZ - dP/dY = - RO*Gy */
	/**/
	/* dSIGyy/dY */
	lefty =(syy[v[7]]-syy[v[6]])*ykf;
	/* dSIGxy/dY */
	lefty+=(sxy[v[2]]-sxy[v[0]])*xkf;
	/* dSIGxz/dZ */
	lefty+=(syz[v[4]]-syz[v[0]])*zkf;
	/* -dP/dX */
	lefty-=(pr[v[7]]-pr[v[6]])*ykf;
	/**/
/*
printf("Vy %ld %ld %ld %e %e %e \n",m1,m2,m3,lefty,righty,lefty-righty);
*/
	/* X-STOKS Error */
	lefty=(lefty-righty);
	/**/
	eps1[0]=lefty*skf1;
	return lefty;
	}
/**/
/**/
/**/
/* Save Right part Save for Y-Stokes ---------------------*/
ui1[0]=righty;
/* Set Initial Num of lines */
un1[0]=17;
/**/
/* Add Coefficients for left parts of Y-Stokes ----------------*/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P6 | P7 */
/* SIGyy6 | SIGyy7 */
/* SIGxy0 - SIGxy2 */
/* SIGyz0 / SIGyz4 */
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* -dP/dX */
un1[16]=v[6]*4+0;
ui1[16]=+ykf;
un1[17]=v[7]*4+0;
ui1[17]=-ykf;
/* vy  4  7   vx 8  10      vz  14   p 16    */
/* 2   1   3     9  11      12  15     17    */
/*  6  5                    13               */
/* dSIGxx/dX */
/* dSIGxx/dX */
syycalc(m1+1,m2  ,m3+1,-ykf,mgi,eps1,un1,ui1);
un1[ 8]=un1[18]; ui1[ 8] =ui1[18];
un1[10]=un1[19]; ui1[10] =ui1[19];
un1[ 4]=un1[20]; ui1[ 4] =ui1[20];
un1[ 1]=un1[21]; ui1[ 1] =ui1[21];
un1[12]=un1[22]; ui1[12] =ui1[22];
un1[14]=un1[23]; ui1[14] =ui1[23];
un1[0]=17;
syycalc(m1+1,m2+1,m3+1, ykf,mgi,eps1,un1,ui1);
un1[ 9]=un1[18]; ui1[ 9] =ui1[18];
un1[11]=un1[19]; ui1[11] =ui1[19];
                 ui1[ 1]+=ui1[20];
un1[ 5]=un1[21]; ui1[ 5] =ui1[21];
un1[13]=un1[22]; ui1[13] =ui1[22];
un1[15]=un1[23]; ui1[15] =ui1[23];
/* dSIGxy/dY */
un1[0]=17;
sxycalc(m1  ,m2  ,m3  ,-xkf,mgi,eps1,un1,ui1);
                 ui1[ 8]+=ui1[18];
                 ui1[ 9]+=ui1[19];
un1[ 2]=un1[20]; ui1[ 2] =ui1[20];
                 ui1[ 1]+=ui1[21];
un1[0]=17;
sxycalc(m1+1,m2  ,m3  , xkf,mgi,eps1,un1,ui1);
                 ui1[10]+=ui1[18];
                 ui1[11]+=ui1[19];
                 ui1[ 1]+=ui1[20];
un1[ 3]=un1[21]; ui1[ 3] =ui1[21];
/* dSIGxz/dZ */
un1[0]=17;
syzcalc(m1  ,m2  ,m3  ,-zkf,mgi,eps1,un1,ui1);
un1[ 6]=un1[18]; ui1[ 6] =ui1[18];
                 ui1[ 1]+=ui1[19];
                 ui1[12]+=ui1[20];
                 ui1[13]+=ui1[21];
un1[0]=17;
syzcalc(m1  ,m2  ,m3+1, zkf,mgi,eps1,un1,ui1);
                 ui1[ 1]+=ui1[18];
un1[ 7]=un1[19]; ui1[ 7] =ui1[19];
                 ui1[14]+=ui1[20];
                 ui1[15]+=ui1[21];
/**/
/* dDIV(v)/dY */
/* vy  4  7   vx 8  10      vz  14   p 16    */
/* 2   1   3     9  11      12  15     17    */
/*  6  5                    13               */
if(p2koef!=0)
	{
un1[0]=17;
	continadd(m1+1,m2  ,m3+1,-ykf*p2koef,mgi,eps1,un1,ui1);
ui1[ 8]+=ui1[18];
ui1[10]+=ui1[19];
ui1[ 4]+=ui1[20];
ui1[ 1]+=ui1[21];
ui1[12]+=ui1[22];
ui1[14]+=ui1[23];
un1[0]=17;
	continadd(m1+1,m2+1,m3+1, ykf*p2koef,mgi,eps1,un1,ui1);
ui1[ 9]+=ui1[18];
ui1[11]+=ui1[19];
ui1[ 1]+=ui1[20];
ui1[ 5]+=ui1[21];
ui1[13]+=ui1[22];
ui1[15]+=ui1[23];
	}
un1[0]=17;
/**/
/* Normalize contents */
/*
for(n1=0;n1<=un[0];n1++) ui[n1]*=skf;
*/
/**/
/**/
/*
for(n1=0;n1<=un[0];n1++)printf("Vy %e %ld \n",ui[n1],un[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */



/* LEFT+Right Side or Err for Z-Stokes Equation */
/* Stoks equation initial form */
/* dSIGzz/dZ + dSIGxz/dX + dSIGyz/dY - dP/dZ = - RO*Gz */
double zstokserr(long int m1, long int m2, long int m3, int ynerr, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[8];
/* Err Buf */
double ival,leftz,rightz,nueff;
/* Distances */
double xkf,ykf,zkf,skf,skf1,GZKOEFZ=GZKOEF;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1]);
ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2]);
zkf=2.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
rightz  = -GZKOEFZ*(ro[v[0]]+ro[v[1]]+ro[v[2]]+ro[v[3]])/4.0;
/**/
/**/
/**/
/* Solution koef calc */
skf=(1.0/xkf+1.0/ykf+1.0/zkf)/3.0;
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P3 | P7 */
/* SIGzz3 / SIGzz7 */
/* SIGxz0 - SIGxz2 */
/* SIGyz0 | SIGyz1 */
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* Effective NU calc */
	if(viscmod==0) nueff=(nu[v[0]-mgxy[mgi]]+nu[v[1]-mgxy[mgi]]+nu[v[2]-mgxy[mgi]]+nu[v[3]-mgxy[mgi]]+10.0*(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]])+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/48.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]-mgxy[mgi]])+log(nu[v[1]-mgxy[mgi]])+log(nu[v[2]-mgxy[mgi]])+log(nu[v[3]-mgxy[mgi]])+10.0*(log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]]))+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/48.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]-mgxy[mgi]]+1.0/nu[v[1]-mgxy[mgi]]+1.0/nu[v[2]-mgxy[mgi]]+1.0/nu[v[3]-mgxy[mgi]]+10.0*(1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]])+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/48.0);
/*
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]])/4.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]]))/4.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]])/4.0);
*/
	if(nueff>nuend) nueff=nuend;
	if(nueff<nubeg) nueff=nubeg;
	skf1=1.0/nueff*(1.0/xkf/xkf+1.0/ykf/ykf+1.0/zkf/zkf)/3.0;
	/**/
	/**/
	/**/
	/* LEFT part of Z-Stokes */
	/* dSIGzz/dZ + dSIGxz/dX + dSIGyz/dY - dP/dZ = - RO*Gz */
	/**/
	/* dSIGyy/dY */
	leftz =(szz[v[7]]-szz[v[3]])*zkf;
	/* dSIGxy/dY */
	leftz+=(sxz[v[2]]-sxz[v[0]])*xkf;
	/* dSIGxz/dZ */
	leftz+=(syz[v[1]]-syz[v[0]])*ykf;
	/* -dP/dX */
	leftz-=(pr[v[7]]-pr[v[3]])*zkf;
	/**/
/*
printf("Vz1 %ld %ld %ld %e %e %e \n",m1,m2,m3,leftz,rightz,leftz-rightz);
*/
	/* Z-STOKS Error */
	leftz=(leftz-rightz);
	/**/
	eps1[0]=leftz*skf1;
	return leftz;
/*
*/
	}
/**/
/**/
/**/
/* Save Right part Save for Z-Stokes ---------------------*/
ui1[0]=rightz;
/* Set Initial Num of lines */
un1[0]=17;
/**/
/* Add Coefficients for left parts of Z-Stokes ----------------*/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* P3 | P7 */
/* SIGzz3 / SIGzz7 */
/* SIGxz0 - SIGxz2 */
/* SIGyz0 | SIGyz1 */
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* -dP/dX */
un1[16]=v[3]*4+0;
ui1[16]=+zkf;
un1[17]=v[7]*4+0;
ui1[17]=-zkf;
/* vz  4  7   vy    10  vx  13  15    p  17  */
/* 2   1   3     8  11    12  14      16     */
/*  6  5         9                           */
/* dSIGxx/dX */
szzcalc(m1+1,m2+1,m3  ,-zkf,mgi,eps1,un1,ui1);
un1[12]=un1[18]; ui1[12] =ui1[18];
un1[14]=un1[19]; ui1[14] =ui1[19];
un1[ 8]=un1[20]; ui1[ 8] =ui1[20];
un1[ 9]=un1[21]; ui1[ 9] =ui1[21];
un1[ 6]=un1[22]; ui1[ 6] =ui1[22];
un1[ 1]=un1[23]; ui1[ 1] =ui1[23];
un1[0]=17;
szzcalc(m1+1,m2+1,m3+1, zkf,mgi,eps1,un1,ui1);
un1[13]=un1[18]; ui1[13] =ui1[18];
un1[15]=un1[19]; ui1[15] =ui1[19];
un1[10]=un1[20]; ui1[10] =ui1[20];
un1[11]=un1[21]; ui1[11] =ui1[21];
                 ui1[ 1]+=ui1[22];
un1[ 7]=un1[23]; ui1[ 7] =ui1[23];
/* dSIGxy/dY */
un1[0]=17;
sxzcalc(m1  ,m2  ,m3  ,-xkf,mgi,eps1,un1,ui1);
                 ui1[12]+=ui1[18];
                 ui1[13]+=ui1[19];
un1[ 2]=un1[20]; ui1[ 2] =ui1[20];
                 ui1[ 1]+=ui1[21];
un1[0]=17;
sxzcalc(m1+1,m2  ,m3  , xkf,mgi,eps1,un1,ui1);
                 ui1[14]+=ui1[18];
                 ui1[15]+=ui1[19];
                 ui1[ 1]+=ui1[20];
un1[ 3]=un1[21]; ui1[ 3] =ui1[21];
/* dSIGxz/dZ */
un1[0]=17;
syzcalc(m1  ,m2  ,m3  ,-ykf,mgi,eps1,un1,ui1);
                 ui1[ 8]+=ui1[18];
                 ui1[10]+=ui1[19];
un1[ 4]=un1[20]; ui1[ 4] =ui1[20];
                 ui1[ 1]+=ui1[21];
un1[0]=17;
syzcalc(m1  ,m2+1,m3  , ykf,mgi,eps1,un1,ui1);
                 ui1[ 9]+=ui1[18];
                 ui1[11]+=ui1[19];
                 ui1[ 1]+=ui1[20];
un1[ 5]=un1[21]; ui1[ 5] =ui1[21];
/**/
/* vz  4  7   vy    10  vx  13  15    p  17  */
/* 2   1   3     8  11    12  14      16     */
/*  6  5         9                           */
/* dDIV(v)/dY */
if(p2koef!=0)
	{
un1[0]=17;
	continadd(m1+1,m2+1,m3  ,-zkf*p2koef,mgi,eps1,un1,ui1);
ui1[12]+=ui1[18];
ui1[14]+=ui1[19];
ui1[ 8]+=ui1[20];
ui1[ 9]+=ui1[21];
ui1[ 6]+=ui1[22];
ui1[ 1]+=ui1[23];
un1[0]=17;
	continadd(m1+1,m2+1,m3+1, zkf*p2koef,mgi,eps1,un1,ui1);
ui1[13]+=ui1[18];
ui1[15]+=ui1[19];
ui1[10]+=ui1[20];
ui1[11]+=ui1[21];
ui1[ 1]+=ui1[22];
ui1[ 7]+=ui1[23];
	}
un1[0]=17;
/**/
/* Normalize contents */
/*
for(n1=0;n1<=un[0];n1++) ui[n1]*=skf;
*/
/**/
/**/
/*
for(n1=0;n1<=un[0];n1++)printf("Vz %e %ld  %e %e %e\n",ui[n1],un[n1],xkf,ykf,zkf);getchar();
*/
return 0;
}
/* Left+Right Side or Err for Z-Stokes Equation */





/* Left side or Err for Continuity Equation  */
/* div(V) = -D(ln(RO))/dt, div(V)=dVx/dX+dVy/dY */
double continerr(long int m1, long int m2, long int m3, int ynerr, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val buf */
double leftc,rightc,nueff;
long int v[8],p[7],m4,m5,m6,m10,m20,m30,m40;
double nubeg1,nuend1,ival,ival1,ival2,ival3,ival4,xkf,ykf,zkf,skf;
int n1,n2;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* -D(ln(RO))/dt add to the right part */
rightc=0;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* Return dVx/dX+dVy/dY+dVz/dZ err ----------------------------*/
if(ynerr==1)
	{
	/* div(V)=dVx/dX+dVy/dY+dVz/dZ */
	leftc=(vx[v[2]]-vx[v[0]])*xkf+(vy[v[1]]-vy[v[0]])*ykf+(vz[v[4]]-vz[v[0]])*zkf;
	skf=(1.0/xkf+1.0/ykf+1.0/zkf)/3.0;
	leftc=(leftc-rightc)*skf;
	return leftc;
/*
printf("Coi1 %ld %ld %ld %e %e %e \n",m1,m2,m3,leftc,rightc,leftc-rightc);getchar();
*/
	}
/**/
/**/
/**/
/* Add continuity equation */
/* Effective NU calc */
/* Basic level */
if(mgi==0 && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
	{
	nueff=nxx[v[7]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Save Right part for Contin ---------------------*/
ui1[0]=rightc;
un1[0]=0;
/**/
/**/
/**/
/* Add Coefficients for left parts of div(V) ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/*  0(P) 1(Vx)  2(Vy) 3(Vz) */
/* Pressure variable */
nd[v[7]]=nueff;
/**/
/* Add Vx with koefficients */
/* Use koef from fre0[], val1[] */
un1[1]=v[0]*4+1;
ui1[1]=-xkf;
un1[2]=v[2]*4+1;
ui1[2]=+xkf;
/* Add Vy with koefficients */
un1[3]=v[0]*4+2;
ui1[3]=-ykf;
un1[4]=v[1]*4+2;
ui1[4]=+ykf;
/* Add Vz with koefficients */
un1[5]=v[0]*4+3;
ui1[5]=-zkf;
un1[6]=v[4]*4+3;
ui1[6]=+zkf;
/* Add total Num of lines */
un1[0]=6;
/*
printf("Cont %ld %ld %ld   %d  %ld  %e",m30,m10,m20,un[0],un[un[0]],ui[un[0]]);getchar();
*/
return leftc;
}
/* Left side or Err for Continuity Equation */



/* Add Continuity Equation  */
/* div(V) = -D(ln(RO))/dt, div(V)=dVx/dX+dVy/dY */
double continadd(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y number */
/* ynval - Val Div(V) Add */
{
/* Val buf */
double rightc;
long int v[8];
double xkf,ykf,zkf,skf,nueff;
int n1,n2;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
/**/
/**/
/**/
/* Effective NU calc */
/* Basic level */
if(mgi==-10 && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
	{
	nueff=nxx[v[7]];
	}
/* Higher levels */
else
	{
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* -D(ln(RO))/dt add to the right part */
rightc=0;
/**/
/**/
/**/
/* Add continuity equation */
/* Save Right part for Contin ---------------------*/
ui1[0]+=rightc*ynval;
/**/
/**/
/**/
/* Add Coefficients for left parts of div(V) ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/*  0(P) 1(Vx)  2(Vy) 3(Vz) */
/* div(V)=dVx/dX+dVy/dY+dVz/dZ */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=-xkf*ynval*nueff;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=+xkf*ynval*nueff;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=-ykf*ynval*nueff;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=+ykf*ynval*nueff;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=-zkf*ynval*nueff;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=+zkf*ynval*nueff;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<3;n1++)printf("Cont %e %d \n",ui[n1],wn[n1]);getchar();
*/
return rightc;
}
/* Add Continuity Equation */




/* Left side or Err for Boundary Condition Equation */ 
/* Ai=CONST+KOEF*An */
double bonderr(long int mcmax, int ynerr, double *eps1, long int *un1, double *ui1)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftb=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr==1)
	{
	/* Add Const */
	leftb=sol0[mcmax]-bondv1[bondm[mcmax]][0];
	/* Add Koef */
	if(bondn1[bondm[mcmax]]) 
		{
		leftb-=bondv1[bondm[mcmax]][1]*sol0[bondn1[bondm[mcmax]]-1];
		}
/*
printf("Bo %ld %e %e \n",mcmax,sol0[mcmax],leftb);
*/
	/**/
	return leftb;
	}
/**/
/**/
/**/
/* Add Y CONST */
un1[0]=1;
ui1[0]=bondv1[bondm[mcmax]][0];
un1[1]=mcmax;
ui1[1]=1.0;
/* Add Y PAR1,PAR2,PAR3 */
if(bondn1[bondm[mcmax]]) 
	{
	un1[0]+=1;
	un1[un1[0]]=bondn1[bondm[mcmax]]-1;
	ui1[un1[0]]=-bondv1[bondm[mcmax]][1];
	}
/*
printf("Bo %ld %e %ld %e \n",mcmax,bondv1[mcmax][0],bondn1[mcmax],bondv1[mcmax][1]);getchar();
printf("Bo %ld %e %e \n",mcmax,sol0[mcmax],leftb);getchar();
for(n1=0;n1<3;n1++)printf("%e %d \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for Boundary Condition Equation */ 



/* Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */
void vpiterate(int m0, int f0)
{
/* Counters */
long int m2prn=0,m1prn=0,m1,m2,m3,m4,m5,m10,m20,m30,m40,mcmax,mcmax0,mcmax1;
double dvxmax,dvymax,dvzmax,maxvxyz,vxmax,vymax,vzmax,minvx,maxvx,minvy,maxvy,minvz,maxvz,minpr,maxpr,minnu,maxnu;
double minvx1[1000],maxvx1[1000],minvy1[1000],maxvy1[1000],minvz1[1000],maxvz1[1000],minpr1[1000],maxpr1[1000],minnu1[1000],maxnu1[1000];
double minvx11,maxvx11,minvy11,maxvy11,minvz11,maxvz11,minpr11,maxpr11,minnu11,maxnu11;
/**/
/* Val buffer */
double ival,ival1,numult1=numult,p00koef=p0koef,p10koef=p1koef,p20koef=p2koef,v00koef=v0koef,nuend0=nuend;
/* Err koef */
double bondsum,bondnum,contsum,contnum,contsumb,contnumb,stoksum,stoksum1,stoksum2,stoknum;
double pbondsum[1000],pbondnum[1000],pcontsum[1000],pcontnum[1000],pcontsumb[1000],pcontnumb[1000],pstoksum[1000],pstoksum1[1000],pstoksum2[1000],pstoknum[1000];
double pbondsum10,pbondnum10,pcontsum10,pcontnum10,pcontsumb10,pcontnumb10,pstoksum10,pstoksum11,pstoksum21,pstoknum10;
int mgi=0,n0,n1,n2,n3,ynstoks=1,multinum0=multinum,multcur=0;
long int v[8];
double nu1=1e+50,nu2=-1e+50;
long int un1[MAXPOS],pos0cur1;
double ui1[MAXPOS],eps1[100];
/**/
/**/
/**/
/* Parallel region ----------------------------------------------------- */
#pragma omp parallel private(n1,m1,m2,m3,mgi,mcmax0,mcmax1,eps1,un1,ui1,minvx11,maxvx11,minvy11,maxvy11,minvz11,maxvz11,minpr11,maxpr11,minnu11,maxnu11) shared(n0,minvx1,maxvx1,minvy1,maxvy1,minvz1,maxvz1,minpr1,maxpr1,minnu1,maxnu1,vx,vy,vz,pr,sxx,syy,szz,sxy,sxz,syz,exx,eyy,ezz,exy,exz,eyz) firstprivate(mgx,mgy,mgz,mgxy,mgp)
{
/*
*/
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Obtain total number of threads */
if (n1==0) n0=omp_get_num_threads();
mgi=0;
/* Define vx,vy,vz,P min,max */
/* Vx,Vy,Vz max-min definition */
minvx11=1e+30;maxvx11=-1e+30;
minvy11=1e+30;maxvy11=-1e+30;
minvz11=1e+30;maxvz11=-1e+30;
minpr11=1e+30;maxpr11=-1e+30;
/* Node  Cycle */
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax0=mcmax1*4;
	/**/	
	/**/	
	/**/	
	/* Reload P */	
	if(m1 && m2 && m3) 
		{
		minpr11=MINV(minpr11,pr[mcmax1]);
		maxpr11=MAXV(maxpr11,pr[mcmax1]);
/*
printf("Pr %ld %ld %ld %e    %e \n",m1,m2,m3,ival,pr[mcmax1]);
printf("Pr %ld %ld %ld %e    %e",m1,m2,m3,ival,pr[mcmax1]); getchar();
*/
		}
	/* Reload Vx */	
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
		minvx11=MINV(minvx11,vx[mcmax1]);
		maxvx11=MAXV(maxvx11,vx[mcmax1]);
/*
printf("Vx %ld %ld %ld    %e",m1,m2,m3,vx[mcmax1]); getchar();
*/
		}
	/* Reload Vy */	
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
		minvy11=MINV(minvy11,vy[mcmax1]);
		maxvy11=MAXV(maxvy11,vy[mcmax1]);
/*
printf("Vy %ld %ld %ld    %e",m1,m2,m3,vy[mcmax1]); getchar();
*/
		}
	/* Reload Vz */	
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
		minvz11=MINV(minvz11,vz[mcmax1]);
		maxvz11=MAXV(maxvz11,vz[mcmax1]);
/*
printf("Vz %ld %ld %ld    %e",m1,m2,m3,vz[mcmax1]); getchar();
*/
		}
	}
minpr1[n1]=minpr11;
maxpr1[n1]=maxpr11;
minvx1[n1]=minvx11;
maxvx1[n1]=maxvx11;
minvy1[n1]=minvy11;
maxvy1[n1]=maxvy11;
minvz1[n1]=minvz11;
maxvz1[n1]=maxvz11;
/* Max Vx,Vy,Vz Diff in Grid Calc */
/**/	
/**/	
/**/	
/* Recalc EPS, SIG Results */
/* Node  Cycle */
minnu11=1e+50;maxnu11=-1e+50;
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/**/	
	if(m1>0 && m2>0 && m3>0)
		{
		/* Sxx,Exx */	
		sxx[mcmax1]=sxxcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exx[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		/**/	
		/* Syy,Eyy */	
		syy[mcmax1]=syycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyy[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		/**/	
		/* Szz,Ezz */	
		szz[mcmax1]=szzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); ezz[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		}
	/**/	
	/* Sxy,Exy */	
	if(m1>0 && m1<mgx[mgi]-1 && m2>0 && m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0];
		}
	/* Sxz,Exz */	
	if(m1>0 && m1<mgx[mgi]-1 && m2<mgy[mgi]-1 && m3>0 && m3<mgz[mgi]-1)
		{
		sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0];
		}
	/* Syz,Eyz */	
	if( m1<mgx[mgi]-1 && m2>0 && m2<mgy[mgi]-1 && m3>0 && m3<mgz[mgi]-1)
		{
		syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0];
		}
	}
minnu1[n1]=minnu11; maxnu1[n1]=maxnu11;
/* End Recalc EPS, SIG Results */
}
/* End parallel */
/* Find maximal values */
minvx=minvx1[0];maxvx=maxvx1[0];
minvy=minvy1[0];maxvy=maxvy1[0];
minvz=minvz1[0];maxvz=maxvz1[0];
minpr=minpr1[0];maxpr=maxpr1[0];
minnu=minnu1[0];maxnu=maxnu1[0];
for(n1=1;n1<n0;n1++)
	{
	minvx=MINV(minvx,minvx1[n1]);
	maxvx=MAXV(maxvx,maxvx1[n1]);
	minvy=MINV(minvy,minvy1[n1]);
	maxvy=MAXV(maxvy,maxvy1[n1]);
	minvz=MINV(minvz,minvz1[n1]);
	maxvz=MAXV(maxvz,maxvz1[n1]);
	minpr=MINV(minpr,minpr1[n1]);
	maxpr=MAXV(maxpr,maxpr1[n1]);
	minnu=MINV(minnu,minnu1[n1]);
	maxnu=MAXV(maxnu,maxnu1[n1]);
	}
dvxmax=maxvx-minvx;
dvymax=maxvy-minvy;
dvzmax=maxvz-minvz;
maxvxyz=pow((dvxmax*dvxmax+dvymax*dvymax+dvzmax*dvzmax)/3.0,0.5);
if (maxvxyz<=0) maxvxyz=1e-8;
if (DIVVMAX>0 && maxvxyz>DIVVMAX) maxvxyz=DIVVMAX;
/* End Define vx,vy,vz,P,NU min,max */
/**/	
/**/	
/**/	
/* Parallel region ----------------------------------------------------- */
#pragma omp parallel private(n1,m1,m2,m3,mgi,mcmax0,eps1,un1,ui1,ival,ival1,pbondsum10,pbondnum10,pcontsum10,pcontnum10,pcontsumb10,pcontnumb10,pstoksum10,pstoksum11,pstoksum21,pstoknum10) shared(bondm,pbondsum,pbondnum,pstoksum,pstoksum1,pstoksum2,pstoknum,pcontsum,pcontnum,pcontsumb,pcontnumb) firstprivate(mgx,mgy,mgz,mgxy,mgp,maxvxyz)
{
/*
*/
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Check Initial Error */
mgi=0;
/* Err koef */
pbondsum10=pbondnum10=0;
pstoksum10=pstoksum11=pstoksum21=pstoknum10=0;
pcontsum10=pcontnum10=0;
pcontsumb10=pcontnumb10=0;
/* Node  Cycle */
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax0=(mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2)*4;
	/**/
	/**/
	/**/
	/* Check Continuity equation for Cells =========================== */
	if(m1 && m2 && m3)
		{
		ival=continerr(m1,m2,m3,1,mgi,eps1,un1,ui1);
		ival/=maxvxyz;
		if(!bondm[mcmax0])
			{
			pcontsum10+=ival*ival;
       	        	pcontnum10+=1.0;
			}
		else
			{
			pcontsumb10+=ival*ival;
       	        	pcontnumb10+=1.0;
			}
		}
	/**/
	/* Check vX-Equations for nodes =========================== */
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
/*
		if(!bondm[mcmax0+1] || (m1>10 && m1<mgx[mgi]-11 && m2>0 && m2<mgy[mgi]-2 && m3>0 && m3<mgz[mgi]-2 && timesum>convend)) 
		if(!bondm[mcmax0+1] || (timesum>convend && bondv1[bondm[mcmax0+1]][0])) 
*/
		if(!bondm[mcmax0+1]) 
			{
			/* Check vX-Stokes */
	                ival=xstokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vx %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vX-Boundary */
			ival=bonderr(mcmax0+1,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	/**/
	/* Check vY-Equations for nodes =========================== */
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
		if(!bondm[mcmax0+2]) 
			{
			/* Check vY-Stokes */
	                ival=ystokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vy %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vY-Boundary */
			ival=bonderr(mcmax0+2,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	/**/
	/* Check vZ-Equations for nodes =========================== */
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
		if(!bondm[mcmax0+3]) 
			{
			/* Check vZ-Stokes */
	                ival=zstokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vz %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vZ-Boundary */
			ival=bonderr(mcmax0+3,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	}
pcontsum[n1]=pcontsum10;
pcontnum[n1]=pcontnum10;
pcontsumb[n1]=pcontsumb10;
pcontnumb[n1]=pcontnumb10;
pstoksum[n1]=pstoksum10;
pstoksum1[n1]=pstoksum11;
pstoksum2[n1]=pstoksum21;
pstoknum[n1]=pstoknum10;
pbondsum[n1]=pbondsum10;
pbondnum[n1]=pbondnum10;
}
/* Find maximum values */
bondsum=pbondsum[0];bondnum=pbondnum[0];
stoksum=pstoksum[0];stoksum1=pstoksum1[0];stoksum2=pstoksum2[0];stoknum=pstoknum[0];
contsum=pcontsum[0];contnum=pcontnum[0];
contsumb=pcontsumb[0];contnumb=pcontnumb[0];
for(n1=1;n1<n0;n1++)
	{
	bondsum+=pbondsum[n1];bondnum+=pbondnum[n1];
	stoksum+=pstoksum[n1];stoksum1+=pstoksum1[n1];stoksum2+=pstoksum2[n1];stoknum+=pstoknum[n1];
	contsum+=pcontsum[n1];contnum+=pcontnum[n1];
	contsumb+=pcontsumb[n1];contnumb+=pcontnumb[n1];
	}
stoksum=pow(stoksum/stoknum,0.5);
stoksum1=pow(stoksum1/stoknum,0.5);
stoksum2=pow(stoksum2/stoknum,0.5);
contsum=pow(contsum/contnum,0.5);
contsumb=pow(contsumb/contnumb,0.5);
bondsum=pow(bondsum/bondnum,0.5);
/**/
/* Print Results */
if (printmod)
	{
	printf("\n FILE %s KRUG %d Cycle %ld TIME %e yr \n",fl0out[f0],m0+1,m2prn,timesum/(365.25*24.0*3600.0));
	printf("v0koef %e v1koef %e p0koef %e p1koef %e p2koef %e \n",v0koef,v1koef,p0koef,p1koef,p2koef);
	printf("X-VELOCITY: min = %e max = %e \n",minvx,maxvx);
	printf("Y-VELOCITY: min = %e max = %e \n",minvy,maxvy);
	printf("Z-VELOCITY: min = %e max = %e \n",minvz,maxvz);
	printf("PRESSURE  : min = %e max = %e \n",minpr,maxpr);
	printf("VISKOSITY : min = %e max = %e \n",minnu,maxnu);
	printf("STOKS LIMIT : num = %e errS= %e errV= %e  errV*errS= %e errVS= %e \n",stoknum,STOKSMIN,STOKSMAX,STOKSMIN*STOKSMAX,STOKSMIN*STOKSMAX);
	printf("STOKS       : num = %e errS= %e errV= %e  errV*errS= %e errVS= %e \n",stoknum,stoksum,stoksum1,stoksum1*stoksum,stoksum2);
	printf("STOKSBC : num = %e err = %e \n",bondnum,bondsum);
	printf("CONT LIMIT  : num = %e err = %e \n",contnum,DIVVMIN);
	printf("CONT        : num = %e err = %e \n",contnum,contsum);
	printf("CONTBC: num = %e err = %e \n",contnumb,contsumb);
	printf("\n");
	}
/**/
/* End Check Initial Error */
/**/
/**/
/**/
/* MAIN ITERATION CYCLE-------------------- */
if(0==0 || stoksum2>STOKSMIN*STOKSMAX || stoksum*stoksum1>STOKSMIN*STOKSMAX || contsum>DIVVMIN)
{
/* Ro, Nu "Restriction" from Finer to Coarcer level */
for(mgi=1;mgi<=multimax;mgi++) ronurestrict(mgi);
/**/
/**/
/**/
/* Clear New Multigrid Solution */
/* P, Vx,Vy */
#pragma omp parallel for shared(sol0,num0) private(m1) firstprivate(multimax,mgp,mgn) schedule(static)
/*
*/
for (m1=0;m1<(mgp[multimax]+mgn[multimax])*4;m1++)
	{
	sol0[m1]=0;
	num0[m1]=0;
	}
/**/
/**/
/**/
if(printmod) printf("\n NU limit %e %e\n",nubeg,nuend);
if(ynstoks)
{
leftnum=0;
rightnum=0;
/* Reload P, Vx, Vy, Vz OLD Results to sol0[] */
/* Node  Cycle */
mgi=0;
#pragma omp parallel for shared(sol0,pr,vx,vy,vz) private(m1,m2,m3,mcmax0,mcmax1) firstprivate(mgi,mgp,mgn,mgy,mgx,mgz,mgxy) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax0=mcmax1*4;
	/**/	
	/**/	
	/**/	
	/* Reload P */	
	if(m1 && m2 && m3) 
		{
		sol0[mcmax0+0]=pr[mcmax1];
/*
printf("Pr %ld %ld %ld    %e",m1,m2,m3,pr[mcmax1]); getchar();
*/
		}
	/* Reload Vx */	
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
		sol0[mcmax0+1]=vx[mcmax1];
/*
printf("Vx %ld %ld %ld    %e",m1,m2,m3,vx[mcmax1]); getchar();
*/
		}
	/* Reload Vy */	
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
		sol0[mcmax0+2]=vy[mcmax1];
/*
printf("Vy %ld %ld %ld    %e",m1,m2,m3,vy[mcmax1]); getchar();
*/
		}
	/* Reload Vz */	
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
		sol0[mcmax0+3]=vz[mcmax1];
/*
printf("Vz %ld %ld %ld    %e",m1,m2,m3,vz[mcmax1]); getchar();
*/
		}
	}
/* End Reload P, Vx, Vy, Vz OLD Results to sol0[] */
/**/
/**/
/**/
/* Add Matrix by vX-vY-Vz-Stokes, Continuity, Boundary, EPS, SIG, Equations */
/* Node  Cycle */
/*
printf("%d %ld %ld %ld %ld %ld %ld\n",mgi,mgx[mgi],mgy[mgi],mgz[mgi],mgxy[mgi],mgp[mgi],mgn[mgi]);getchar();
*/
for(mgi=0;mgi<=multimax;mgi++)
{
#pragma omp parallel for shared(bondm) private(m1,m2,m3,mcmax,mcmax1,pos0cur1,eps1,un1,ui1,ival) firstprivate(mgi,multimax,mgp,mgx,mgy,mgz,mgxy,timesum) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax=mcmax1*4;
	/* Pos P,Vx,Vy,Vz in val1[] */
	pos0cur1=mcmax1*57;
/*
printf("BoP  %ld %ld %ld  %ld %ld %e %ld %e \n",m1,m2,m3,mcmax+0,bondm[mcmax+0],bondv1[mcmax+0][0],bondn1[mcmax+0],bondv1[mcmax+0][1]);
printf("BoX  %ld %ld %ld  %ld %ld %e %ld %e \n",m1,m2,m3,mcmax+1,bondm[mcmax+1],bondv1[mcmax+1][0],bondn1[mcmax+1],bondv1[mcmax+1][1]);
printf("BoY  %ld %ld %ld  %ld %ld %e %ld %e \n",m1,m2,m3,mcmax+2,bondm[mcmax+2],bondv1[mcmax+2][0],bondn1[mcmax+2],bondv1[mcmax+2][1]);
printf("BoZ  %ld %ld %ld  %ld %ld %e %ld %e \n",m1,m2,m3,mcmax+3,bondm[mcmax+3],bondv1[mcmax+3][0],bondn1[mcmax+3],bondv1[mcmax+3][1]);getchar();
printf("P-Vx-Vy-Vz Cycle %ld %ld %ld   %ld",m1,m2,m3,mcmax); getchar();
printf("P-Vx-Vy-Vz Cycle %ld %ld %ld   %ld \n",m1,m2,m3,mcmax);
*/
	/**/
	/**/
	/**/
	/* Add Continuity equation for Cells ------------------------------------------------ */
	if(m1 && m2 && m3)
		{
/*
printf("Pr %ld %ld %ld    %ld",m1,m2,m3,bondm[mcmax+0]); getchar();
printf("Pr %ld %ld %ld    %ld \n",m1,m2,m3,bondm[mcmax+0]);
*/
		/* No Boundary condition for 222-cell for iterative solver */ 
		if(!bondm[mcmax+0]) 
			{
			/* Check Pressure in Cell calculation Y/N */ 
	                continerr(m1,m2,m3,0,mgi,eps1,un1,ui1);
	        	gausmat2(1,mcmax+0,pos0cur1,un1,ui1);
			}
		else
			{
			/* Add P-Boundary */
			bonderr(mcmax+0,0,eps1,un1,ui1);
	        	gausmat2(1,mcmax+0,pos0cur1,un1,ui1);
			}
		}

	/**/
	/**/
	/**/
	/* Add vX-Equations --------------------------------------------------- */
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
/*
printf("Vx %ld %ld %ld    %ld",m1,m2,m3,bondm[mcmax+1]); getchar();
printf("Vx %ld %ld %ld    %ld \n",m1,m2,m3,bondm[mcmax+1]);
printf("XSTOKES %ld",pos0cur1);getchar();
*/
/*
		if(!bondm[mcmax+1] || (mgi==0 && m1>10 && m1<mgx[mgi]-11 && m2>0 && m2<mgy[mgi]-2 && m3>0 && m3<mgz[mgi]-2 && timesum>convend)) 
		if(!bondm[mcmax+1] || (timesum>convend && bondv1[bondm[mcmax+1]][0])) 
*/
		if(!bondm[mcmax+1]) 
			{
			/* Add vX-Stokes */
	                xstokserr(m1,m2,m3,0,mgi,eps1,un1,ui1);
			gausmat2(1,mcmax+1,pos0cur1+6,un1,ui1);
			}
		else
			{
			/* Add vX-Boundary */
			bonderr(mcmax+1,0,eps1,un1,ui1);
			gausmat2(1,mcmax+1,pos0cur1+6,un1,ui1);
			}
		}
	/**/
	/**/
	/**/
	/* Add vY-Equations --------------------------------------------------- */
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
/*
printf("Vy %ld %ld %ld    %ld",m1,m2,m3,bondm[mcmax+2]); getchar();
printf("Vy %ld %ld %ld    %ld \n",m1,m2,m3,bondm[mcmax+2]);
*/
		if(!bondm[mcmax+2]) 
			{
			/* Add vY-Stokes */
        	        ystokserr(m1,m2,m3,0,mgi,eps1,un1,ui1);
			gausmat2(1,mcmax+2,pos0cur1+23,un1,ui1);
			}
		else
			{
			/* Add vY-Boundary */
			bonderr(mcmax+2,0,eps1,un1,ui1);
			gausmat2(1,mcmax+2,pos0cur1+23,un1,ui1);
			}
		}
	/**/
	/**/
	/**/
	/* Add vZ-Equations --------------------------------------------------- */
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
/*
printf("Vz %ld %ld %ld    %ld",m1,m2,m3,bondm[mcmax+3]); getchar();
printf("Vz %ld %ld %ld    %ld \n",m1,m2,m3,bondm[mcmax+3]);
*/
		if(!bondm[mcmax+3]) 
			{
			/* Add vZ-Stokes */
        	        zstokserr(m1,m2,m3,0,mgi,eps1,un1,ui1);
			gausmat2(1,mcmax+3,pos0cur1+40,un1,ui1);
			}
		else
			{
			/* Add vZ-Boundary */
			bonderr(mcmax+3,0,eps1,un1,ui1);
			gausmat2(1,mcmax+3,pos0cur1+40,un1,ui1);
			}
		}
	/**/
	/**/
	/**/
	}
/* End  Add Matrix By vX-vY-Vz-Stokes, Continuity Equations */
/*
if(printmod) printf("mgi = %d \n",mgi);
*/
}
}
mcmax=mgp[multimax]+mgn[multimax];
pos0cur1=mcmax*57;
mcmax1=mcmax*4;
if(printmod) printf("mcmax = %ld mcmax1 = %ld Pos0cur1 = %ld \n",mcmax,mcmax1,pos0cur1);
/**/
/**/
do
{
/* Mulrigrid/Simple iteration cycles */
for(n1=0;n1<multicyc;n1++)
	{
	printf("Iter%d/%d: ",n1+1,multicyc);
	/* Mulrigrid V cycle */
	if(multinum>0)
		{
		/* Up cylcle */
		for(n2=0;n2<multinum;n2++)
			{
			mgi=n2;
/*
			p0koef=multinnp[0][mgi];
			v0koef=multinnv[0][mgi];
*/
			if(multinnn[0][mgi]>0) gausmat4(multinnn[0][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
			vprestrict(mgi+1);
			if(printmod) printf("Up%d/%d ",n2,multinnn[0][mgi]);
			}
		/* Down cycle */
		for(n2=multinum;n2>0;n2--)
			{
			mgi=n2;
/*
			p0koef=multinnp[1][mgi];
			v0koef=multinnv[1][mgi];
*/
			if(multinnn[1][mgi]>0) gausmat4(multinnn[1][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
			vpprolong(mgi,0);
			if(printmod) printf("Down%d/%d ",n2,multinnn[1][mgi]);
			}
		mgi=0;
		if(multinnn[1][mgi]>0) 
			{
			gausmat4(multinnn[1][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
			if(printmod) printf("Down%d/%d ",mgi,multinnn[1][mgi]);
			}
		}
	else
		{
		/* Simple cylcle */
		mgi=0;
		if(multinnn[0][mgi]>0) 
			{
			gausmat4(multinnn[0][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
			if(printmod) printf("SimpleUP%d/%d ",mgi,multinnn[0][mgi]);
			}
		else
			{
			if(multinnn[1][mgi]>0) gausmat4(multinnn[1][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
			if(printmod) printf("SimpleDOWN%d/%d ",mgi,multinnn[1][mgi]);
			}
		}
	if(printmod) printf("\n");
	}
mgi=0;
/*
gausmat4(multinnn[0][mgi],(mgp[mgi]+mgn[mgi])*4-1,mgp[mgi]*4,mgi);
*/
/* End Mulrigrid V cycle */
/**/
/**/
/**/
#pragma omp parallel private(n1,m1,m2,m3,mgi,mcmax0,mcmax1,eps1,un1,ui1,ival,minvx11,maxvx11,minvy11,maxvy11,minvz11,maxvz11,minpr11,maxpr11,minnu11,maxnu11) shared(n0,minvx1,maxvx1,minvy1,maxvy1,minvz1,maxvz1,minpr1,maxpr1,minnu1,maxnu1,vx,vy,vz,pr,sxx,syy,szz,sxy,sxz,syz,exx,eyy,ezz,exy,exz,eyz,sol0) firstprivate(mgp,mgx,mgy,mgz,mgxy,pxinit,pyinit,pzinit,pinit)
{
/*
*/
/* Obtain cur thread number */
n1=omp_get_thread_num();
mgi=0;
/* Define vx,vy,vz,P min,max */
/* Vx,Vy,Vz max-min definition */
minvx11=1e+30;maxvx11=-1e+30;
minvy11=1e+30;maxvy11=-1e+30;
minvz11=1e+30;maxvz11=-1e+30;
minpr11=1e+30;maxpr11=-1e+30;
/* Node  Cycle */
/* Reload P, Vx, Vy, Vz Results */
/* Pressure reset in the first cell */
m1=pxinit;m2=pyinit;m3=pzinit;
mcmax1=(mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2)*4;
ival=pinit-sol0[mcmax1];
/* Node  Cycle */
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax0=mcmax1*4;
	/**/	
	/**/	
	/**/	
	/* Reload P */	
	if(m1 && m2 && m3) 
		{
		pr[mcmax1]=sol0[mcmax0+0]+ival;
		minpr11=MINV(minpr11,pr[mcmax1]);
		maxpr11=MAXV(maxpr11,pr[mcmax1]);
if(pr[mcmax1]>1e+12) {printf("Pressure is too big %e ",pr[mcmax1]); exit(0);}
/*
printf("Pr %ld %ld %ld %e    %e \n",m1,m2,m3,ival,pr[mcmax1]);
printf("Pr %ld %ld %ld %e    %e",m1,m2,m3,ival,pr[mcmax1]); getchar();
*/
		}
	/* Reload Vx */	
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
		vx[mcmax1]=sol0[mcmax0+1];
		minvx11=MINV(minvx11,vx[mcmax1]);
		maxvx11=MAXV(maxvx11,vx[mcmax1]);
/*
printf("Vx %ld %ld %ld    %e",m1,m2,m3,vx[mcmax1]); getchar();
*/
		}
	/* Reload Vy */	
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
		vy[mcmax1]=sol0[mcmax0+2];
		minvy11=MINV(minvy11,vy[mcmax1]);
		maxvy11=MAXV(maxvy11,vy[mcmax1]);
/*
printf("Vy %ld %ld %ld    %e",m1,m2,m3,vy[mcmax1]); getchar();
*/
		}
	/* Reload Vz */	
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
		vz[mcmax1]=sol0[mcmax0+3];
		minvz11=MINV(minvz11,vz[mcmax1]);
		maxvz11=MAXV(maxvz11,vz[mcmax1]);
/*
printf("Vz %ld %ld %ld    %e",m1,m2,m3,vz[mcmax1]); getchar();
*/
		}
	}
minpr1[n1]=minpr11;
maxpr1[n1]=maxpr11;
minvx1[n1]=minvx11;
maxvx1[n1]=maxvx11;
minvy1[n1]=minvy11;
maxvy1[n1]=maxvy11;
minvz1[n1]=minvz11;
maxvz1[n1]=maxvz11;
/* End Reload P, Vx, Vy, Vz Results */
/**/	
/**/	
/**/	
/* Recalc EPS, SIG Results */
/* Node  Cycle */
minnu11=1e+50;maxnu11=-1e+50;
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in vx[], vy[], vz[], pr[], etc. */
	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	/**/	
	/* Sxx,Exx */	
	if(m1>0 && m2>0 && m3>0)
		{
		sxx[mcmax1]=sxxcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exx[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		/**/	
		/* Syy,Eyy */	
		syy[mcmax1]=syycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyy[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		/**/	
		/* Szz,Ezz */	
		szz[mcmax1]=szzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); ezz[mcmax1]=eps1[0];
		/* Min,Max Nu value Calc */
		minnu11=MINV(minnu11,eps1[2]); maxnu11=MAXV(maxnu11,eps1[2]);
		}
	/**/	
	/* Sxy,Exy */	
	if(m1>0 && m1<mgx[mgi]-1 && m2>0 && m2<mgy[mgi]-1 &&  m3<mgz[mgi]-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0];
		}
	/* Sxz,Exz */	
	if(m1>0 && m1<mgx[mgi]-1 && m2<mgy[mgi]-1 && m3>0 && m3<mgz[mgi]-1)
		{
		sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0];
		}
	/* Syz,Eyz */	
	if(m1<mgx[mgi]-1 && m2>0 && m2<mgy[mgi]-1 && m3>0 && m3<mgz[mgi]-1)
		{
		syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0];
		}
	}
minnu1[n1]=minnu11; maxnu1[n1]=maxnu11;
/* End Recalc EPS, SIG Results */
}
/* End Parallel region */
/**/	
/**/	
/**/	
/* Max Vx,Vy,Vz Diff in Grid Calc */
/* Find maximal values */
minvx=minvx1[0];maxvx=maxvx1[0];
minvy=minvy1[0];maxvy=maxvy1[0];
minvz=minvz1[0];maxvz=maxvz1[0];
minpr=minpr1[0];maxpr=maxpr1[0];
minnu=minnu1[0];maxnu=maxnu1[0];
for(n1=1;n1<n0;n1++)
	{
	minvx=MINV(minvx,minvx1[n1]);
	maxvx=MAXV(maxvx,maxvx1[n1]);
	minvy=MINV(minvy,minvy1[n1]);
	maxvy=MAXV(maxvy,maxvy1[n1]);
	minvz=MINV(minvz,minvz1[n1]);
	maxvz=MAXV(maxvz,maxvz1[n1]);
	minpr=MINV(minpr,minpr1[n1]);
	maxpr=MAXV(maxpr,maxpr1[n1]);
	minnu=MINV(minnu,minnu1[n1]);
	maxnu=MAXV(maxnu,maxnu1[n1]);
	}
dvxmax=maxvx-minvx;
dvymax=maxvy-minvy;
dvzmax=maxvz-minvz;
maxvxyz=pow((dvxmax*dvxmax+dvymax*dvymax+dvzmax*dvzmax)/3.0,0.5);
if (maxvxyz<=0) maxvxyz=1e-8;
if (DIVVMAX>0) maxvxyz=DIVVMAX;
/**/
/**/
/**/
/* Parallel region ----------------------------------------------------- */
#pragma omp parallel private(n1,m1,m2,m3,mgi,mcmax0,eps1,un1,ui1,ival,ival1,pbondsum10,pbondnum10,pcontsum10,pcontnum10,pcontsumb10,pcontnumb10,pstoksum10,pstoksum11,pstoksum21,pstoknum10) shared(bondm,pbondsum,pbondnum,pstoksum,pstoksum1,pstoksum2,pstoknum,pcontsum,pcontnum,pcontsumb,pcontnumb) firstprivate(mgx,mgy,mgz,mgxy,maxvxyz)
{
/*
*/
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Check Initial Error */
mgi=0;
/* Err koef */
pbondsum10=pbondnum10=0;
pstoksum10=pstoksum11=pstoksum21=pstoknum10=0;
pcontsum10=pcontnum10=0;
pcontsumb10=pcontnumb10=0;
/* Node  Cycle */
#pragma omp for schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos P,Vx,Vy,Vz in sol0[] */
	mcmax0=(mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2)*4;
	/**/
	/**/
	/**/
	/* Check Continuity equation for Cells =========================== */
	if(m1 && m2 && m3)
		{
		ival=continerr(m1,m2,m3,1,mgi,eps1,un1,ui1);
		ival/=maxvxyz;
		if(!bondm[mcmax0])
			{
			pcontsum10+=ival*ival;
       	        	pcontnum10+=1.0;
			}
		else
			{
			pcontsumb10+=ival*ival;
       	        	pcontnumb10+=1.0;
			}
		}
	/**/
	/* Check vX-Equations for nodes =========================== */
	if(m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
		{
/*
		if(!bondm[mcmax0+1] || (m1>10 && m1<mgx[mgi]-11 && m2>0 && m2<mgy[mgi]-2 && m3>0 && m3<mgz[mgi]-2 && timesum>convend)) 
		if(!bondm[mcmax0+1] || (timesum>convend && bondv1[bondm[mcmax0+1]][0])) 
*/
		if(!bondm[mcmax0+1]) 
			{
			/* Check vX-Stokes */
	                ival=xstokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vx %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vX-Boundary */
			ival=bonderr(mcmax0+1,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	/**/
	/* Check vY-Equations for nodes =========================== */
	if(m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
		{
		if(!bondm[mcmax0+2]) 
			{
			/* Check vY-Stokes */
	                ival=ystokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vy %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vY-Boundary */
			ival=bonderr(mcmax0+2,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	/**/
	/* Check vZ-Equations for nodes =========================== */
	if(m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
		{
		if(!bondm[mcmax0+3]) 
			{
			/* Check vZ-Stokes */
	                ival=zstokserr(m1,m2,m3,1,mgi,eps1,un1,ui1);
			ival1=eps1[0]/maxvxyz;
        	        pstoksum10+=ival*ival;
        	        pstoksum11+=ival1*ival1;
        	        pstoksum21+=ival1*ival1*ival*ival;
                	pstoknum10+=1.0;
/*
printf("Vz %d %ld %ld %ld %e %e ",mgi,m1,m2,m3,maxvxyz,ival);getchar();
*/
			}
		else
			{
			/* Check vZ-Boundary */
			ival=bonderr(mcmax0+3,1,eps1,un1,ui1);
			pbondsum10+=ival*ival;
	                pbondnum10+=1.0;
			}
		}
	}
pcontsum[n1]=pcontsum10;
pcontnum[n1]=pcontnum10;
pcontsumb[n1]=pcontsumb10;
pcontnumb[n1]=pcontnumb10;
pstoksum[n1]=pstoksum10;
pstoksum1[n1]=pstoksum11;
pstoksum2[n1]=pstoksum21;
pstoknum[n1]=pstoknum10;
pbondsum[n1]=pbondsum10;
pbondnum[n1]=pbondnum10;
}
/* Find maximum values */
bondsum=pbondsum[0];bondnum=pbondnum[0];
stoksum=pstoksum[0];stoksum1=pstoksum1[0];stoksum2=pstoksum2[0];stoknum=pstoknum[0];
contsum=pcontsum[0];contnum=pcontnum[0];
contsumb=pcontsumb[0];contnumb=pcontnumb[0];
for(n1=1;n1<n0;n1++)
	{
	bondsum+=pbondsum[n1];bondnum+=pbondnum[n1];
	stoksum+=pstoksum[n1];stoksum1+=pstoksum1[n1];stoksum2+=pstoksum2[n1];stoknum+=pstoknum[n1];
	contsum+=pcontsum[n1];contnum+=pcontnum[n1];
	contsumb+=pcontsumb[n1];contnumb+=pcontnumb[n1];
	}
stoksum=pow(stoksum/stoknum,0.5);
stoksum1=pow(stoksum1/stoknum,0.5);
stoksum2=pow(stoksum2/stoknum,0.5);
contsum=pow(contsum/contnum,0.5);
contsumb=pow(contsumb/contnumb,0.5);
bondsum=pow(bondsum/bondnum,0.5);
/**/
/**/
/**/
/* Print Results */
m1prn=m1prn+1;
m2prn=m2prn+1;
if (printmod && printmod<=m1prn)
	{
	printf("\n FILE %s KRUG %d Cycle %ld CPUs = %d \n",fl0out[f0],m0+1,m2prn,n0);
	printf("v0koef %e v1koef %e p0koef %e p1koef %e p2koef %e \n",v0koef,v1koef,p0koef,p1koef,p2koef);
	printf("X-VELOCITY: min = %e max = %e \n",minvx,maxvx);
	printf("Y-VELOCITY: min = %e max = %e \n",minvy,maxvy);
	printf("Z-VELOCITY: min = %e max = %e \n",minvz,maxvz);
	printf("PRESSURE  : min = %e max = %e \n",minpr,maxpr);
	printf("VISKOSITY : min = %e max = %e \n",minnu,maxnu);
	printf("STOKS LIMIT : num = %e errS= %e errV= %e  errV*errS= %e errVS= %e \n",stoknum,STOKSMIN,STOKSMAX,STOKSMIN*STOKSMAX,STOKSMIN*STOKSMAX);
	printf("STOKS       : num = %e errS= %e errV= %e  errV*errS= %e errVS= %e \n",stoknum,stoksum,stoksum1,stoksum1*stoksum,stoksum2);
	printf("STOKSBC : num = %e err = %e \n",bondnum,bondsum);
	printf("CONT LIMIT  : num = %e err = %e \n",contnum,DIVVMIN);
	printf("CONT        : num = %e err = %e \n",contnum,contsum);
	printf("CONTBC: num = %e err = %e \n",contnumb,contsumb);
	printf("\n");
	m1prn=0;
/*
getchar();
*/
	}
/**/
/* Change multigrid level number and parameters */
if(m2prn>=fl0cyc[f0][2] && multcur==0)
	{
	multinum=fl0cyc[f0][3];
	p0koef=fl0stp[f0][8];
	p1koef=fl0stp[f0][9];
	p2koef=fl0stp[f0][10];
/*
*/
	}
/**/
}
while(stoksum2>STOKSMIN*STOKSMAX || stoksum*stoksum1>STOKSMIN*STOKSMAX || contsum>DIVVMIN);
}
/* End MAIN ITERATION CYCLE-------------------- */
p0koef=p00koef;
p1koef=p10koef;
p2koef=p20koef;
v0koef=v00koef;
multinum=multinum0;
nuend=nuend0;
/**/
/**/
/**/
/* Time step for markers definition */
	{
	vxmax=MAXV(ABSV(minvx),ABSV(maxvx));
	vymax=MAXV(ABSV(minvy),ABSV(maxvy));
	vzmax=MAXV(ABSV(minvz),ABSV(maxvz));
	if (vxmax)
		{
		vxmax=maxxystep/vxmax;
		if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",vxmax/3.15576e+7);
		if(timestep>vxmax) timestep=vxmax;
		}
	if (vymax)
		{
		vymax=maxxystep/vymax;
		if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",vymax/3.15576e+7);
		if(timestep>vymax) timestep=vymax;
		}
	if (vzmax)
		{
		vzmax=maxxystep/vzmax;
		if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vz-MARKER %e YEAR !!!\n",vzmax/3.15576e+7);
		if(timestep>vzmax) timestep=vzmax;
		}
	if(printmod) printf("\n !!!       CURRENT TIME STEP FOR CYCLE %e YEAR !!!\n",timestep/3.15576e+7);
	}
/**/
/**/
/**/
}
/* End Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */



/* Ro, Nu "Restriction" from Finer to Coarser level */
void ronurestrict(int mgi)
/* mgi - Coarser multigrid level number */
{
/* Counters */
long int m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min;
/* Buffers */
double dx,dy,dz,ddx,ddy,ddz,swt,ival2;
/**/
/**/
/**/
/* Clear Old Ro, Nu */
#pragma omp parallel for shared(sol1,nu) private(m1) firstprivate(mgi,mgp,mgn) schedule(static)
/*
*/
for (m1=mgp[mgi];m1<(mgp[mgi]+mgn[mgi]);m1++)
	{
	nu[m1]=0;
	sol1[m1]=0;
	}
/* Add New Ro, Nu */
#pragma omp parallel for shared(sol1,nu) private(m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min,dx,dy,dz,ddx,ddy,ddz,swt,ival2) firstprivate(mgi,mgp,mgx,mgy,mgz,mgxy,mggx,mggy,mggz,viscmod) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi-1];m3++)
for (m1=0;m1<mgx[mgi-1];m1++)
for (m2=0;m2<mgy[mgi-1];m2++)
	{
	/* Pos in ro[] nu[] for the finer level */
	m4=mgp[mgi-1]+m3*mgxy[mgi-1]+m1*mgy[mgi-1]+m2;
	/* Add 8 surround nodes for the coarser level */
	/* Cell indexes */
	m10min=(m1-2)/2+2;
	if(m1<3) m10min=m1;
	if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
	if(m10min>mgx[mgi]-2) m10min=mgx[mgi]-2;
	m20min=(m2-2)/2+2;
	if(m2<3) m20min=m2;
	if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
	if(m20min>mgy[mgi]-2) m20min=mgy[mgi]-2;
	m30min=(m3-2)/2+2;
	if(m3<3) m30min=m3;
	if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
	if(m30min>mgz[mgi]-2) m30min=mgz[mgi]-2;
	/* Cell dimensions */
	dx=mggx[mgi][m10min+1]-mggx[mgi][m10min];
	dy=mggy[mgi][m20min+1]-mggy[mgi][m20min];
	dz=mggz[mgi][m30min+1]-mggz[mgi][m30min];
/*
printf("%ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
	/* Add 8 nodes */
	for (m30=m30min;m30<=m30min+1;m30++)
		{
		/* Normalized Z-distance calc */
		ddz=(mggz[mgi][m30]-mggz[mgi-1][m3])/dz;
		ddz=1.0-ABSV(ddz);
		for (m10=m10min;m10<=m10min+1;m10++)
			{
			/* Normalized X-distance calc */
			ddx=(mggx[mgi][m10]-mggx[mgi-1][m1])/dx;
			ddx=1.0-ABSV(ddx);
			for (m20=m20min;m20<=m20min+1;m20++)
				{
				/* Normalized Y-distance calc */
				ddy=(mggy[mgi][m20]-mggy[mgi-1][m2])/dy;
				ddy=1.0-ABSV(ddy);
				/* Pos in ro[] nu[] for the coarcer level */
				m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
				/* Calc wt */
				swt=ddx*ddy*ddz;
/*
printf("%ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();
*/
				/* Add Ro, Nu for current node */
				if(viscmod==0) ival2=nu[m4]*swt;
				if(viscmod==1) ival2=log(nu[m4])*swt;
				if(viscmod==2) ival2=1.0/nu[m4]*swt;
				nu[m40]+=ival2;
				ival2=ro[m4]*swt;
				/* Add Wt for current node */
				sol1[m40]+=swt;
				}
			}
		}
	}
/* Calc New Ro, Nu */
#pragma omp parallel for shared(sol1,nu) private(m1) firstprivate(mgi,mgp,mgn) schedule(static)
/*
*/
for (m1=mgp[mgi];m1<(mgp[mgi]+mgn[mgi]);m1++)
	{
	nu[m1]/=sol1[m1];
	if(viscmod==1) nu[m1]=exp(nu[m1]);
	if(viscmod==2) nu[m1]=1.0/nu[m1];
/*
printf("%ld %e %e %e",m1,nu[m1],ro[m1],sol1[m1]);getchar();
*/
	sol1[m1]=0;
	}
}
/* End Ro, Nu "Restriction" from Finer to Coarser level */




/* Vx,Vy,Vz,P "Prolongation" from coarser to finer level */
void vpprolong(int mgi, int ynreset)
/* mgi - Coarser multigrid level number */
/* ynreset - Reset solution Y(1)/N(0) */
{
/* Counters */
long int m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min;
/* Buffers */
double dx,dy,dz,ddx,ddy,ddz,ival,swt,ival1,ivalp;
/**/
/**/
/**/
/* Interpolate Add New Vx,Vy,Vz,P */
#pragma omp parallel for shared(sol0,nd) private(m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min,dx,dy,dz,ddx,ddy,ddz,ival,swt,ival1,ivalp) firstprivate(p1koef,v1koef,mgi,ynreset,mgp,mgx,mgy,mgz,mgxy,mggx,mggy,mggz) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi-1];m3++)
for (m1=0;m1<mgx[mgi-1];m1++)
for (m2=0;m2<mgy[mgi-1];m2++)
	{
	/* Pos in vx[] vy[] vz[] pr[] for the finer level */
	m4=mgp[mgi-1]+m3*mgxy[mgi-1]+m1*mgy[mgi-1]+m2;
	/**/
	/**/
	/**/
	/* Vx interpolation */
	if(m2<mgy[mgi-1]-1 && m3<mgz[mgi-1]-1)
		{
		/* Vx Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-2) m10min=mgx[mgi]-2;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-3) m20min=mgy[mgi]-3;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2+1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min+1])/2.0 && m20min>0) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-3) m30min=mgz[mgi]-3;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3+1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min+1])/2.0 && m30min>0) m30min-=1;
		/* Cell dimensions */
		dx= mggx[mgi][m10min+1]-mggx[mgi][m10min];
		dy=(mggy[mgi][m20min+2]-mggy[mgi][m20min])/2.0;
		dz=(mggz[mgi][m30min+2]-mggz[mgi][m30min])/2.0;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset velocity if required */
		if(ynreset) sol0[m4*4+1]=0;
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30+1]-mggz[mgi-1][m3]-mggz[mgi-1][m3+1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]-mggx[mgi-1][m1])/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20+1]-mggy[mgi-1][m2]-mggy[mgi-1][m2+1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					/* Calc wt */
					swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol0[m4*4+1]);getchar();}
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
					/* Add Vx for current node */
/*
					sol0[m4*4+1]+=sol0[m40*4+1]/(nd[m4+1+mgxy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])*(nd[m40+1+mgxy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt*v1koef;
*/
					sol0[m4*4+1]+=sol0[m40*4+1]*swt*v1koef;
					/* Add Wt for current node */
					ival+=swt;
					}
				}
			}
/*
printf("Vx %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vx[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* Vy interpolation */
	if(m1<mgx[mgi-1]-1 && m3<mgz[mgi-1]-1)
		{
		/* Vy Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-3) m10min=mgx[mgi]-3;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1+1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min+1])/2.0 && m10min>0) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-2) m20min=mgy[mgi]-2;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-3) m30min=mgz[mgi]-3;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3+1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min+1])/2.0 && m30min>0) m30min-=1;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+2]-mggx[mgi][m10min])/2.0;
		dy= mggy[mgi][m20min+1]-mggy[mgi][m20min];
		dz=(mggz[mgi][m30min+2]-mggz[mgi][m30min])/2.0;
/*
printf("Vy %ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset velocity if required */
		if(ynreset) sol0[m4*4+2]=0;
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30+1]-mggz[mgi-1][m3]-mggz[mgi-1][m3+1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10+1]-mggx[mgi-1][m1]-mggx[mgi-1][m1+1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]-mggy[mgi-1][m2])/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					/* Calc wt */
					swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol0[m4*4+2]);getchar();}
printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vy[m40]);getchar();
*/
					/* Add Vy for current node */
					sol0[m4*4+2]+=sol0[m40*4+2]*swt*v1koef;
					/* Add Wt for current node */
					ival+=swt;
					}
				}
			}
/*
printf("Vy %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vy[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* Vz interpolation */
	if(m1<mgx[mgi-1]-1 && m2<mgy[mgi-1]-1)
		{
		/* Vz Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-3) m10min=mgx[mgi]-3;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1+1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min+1])/2.0 && m10min>0) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-3) m20min=mgy[mgi]-3;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2+1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min+1])/2.0 && m20min>0) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-2) m30min=mgz[mgi]-2;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+2]-mggx[mgi][m10min])/2.0;
		dy=(mggy[mgi][m20min+2]-mggy[mgi][m20min])/2.0;
		dz= mggz[mgi][m30min+1]-mggz[mgi][m30min];
/*
printf("%ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset velocity if required */
		if(ynreset) sol0[m4*4+3]=0;
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]-mggz[mgi-1][m3])/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10+1]-mggx[mgi-1][m1]-mggx[mgi-1][m1+1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20+1]-mggy[mgi-1][m2]-mggy[mgi-1][m2+1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					/* Calc wt */
					swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol0[m4*4+3]);getchar();}
printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();
*/
					/* Add Vz for current node */
					sol0[m4*4+3]+=sol0[m40*4+3]*swt*v1koef;
					/* Add Wt for current node */
					ival+=swt;
					}
				}
			}
/*
printf("Vz %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vz[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* P  interpolation */
	if(m1 && m2 && m3)
		{
		/* P Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-2) m10min=mgx[mgi]-2;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1-1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min-1])/2.0 && m10min>1) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-2) m20min=mgy[mgi]-2;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2-1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min-1])/2.0 && m20min>1) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-2) m30min=mgz[mgi]-2;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3-1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min-1])/2.0 && m30min>1) m30min-=1;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+1]-mggx[mgi][m10min-1])/2.0;
		dy=(mggy[mgi][m20min+1]-mggy[mgi][m20min-1])/2.0;
		dz=(mggz[mgi][m30min+1]-mggz[mgi][m30min-1])/2.0;
/*
printf("%ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset pressure if required */
		if(ynreset) sol0[m4*4]=0;
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30-1]-mggz[mgi-1][m3]-mggz[mgi-1][m3-1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10-1]-mggx[mgi-1][m1]-mggx[mgi-1][m1-1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20-1]-mggy[mgi-1][m2]-mggy[mgi-1][m2-1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					/* Calc wt */
					swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Pr %ld %ld %ld   %ld %ld %ld %ld  %e %e %e   %e %e %e   %e %e",m1,m2,m3,m10,m20,m30,m40,dx,dy,dz,ddx,ddy,ddz,swt,sol0[m4*4]);getchar();}
printf("Pr %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,pr[m40]);getchar();
*/
					/* Add Pr for current node */
/*
					sol0[m4*4]+=sol0[m40*4]*nd[m4]/nd[m40]*swt*p1koef;
*/
					sol0[m4*4]+=sol0[m40*4]*swt*p1koef;
					/* Add Wt for current node */
					ival+=swt;
					}
				}
			}
/*
printf("Pr %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,pr[m4]);getchar();
*/
		}
	}
/**/
/*
printf("%d   %ld %ld %ld %ld    %e %e",mgi-1,m1,num0[m1],pos0[m1],bondm[m1],sol0[m1],ival);getchar();
*/
/**/
/*
m1=2;m2=2;m3=2;
m1=(mgp[mgi-1]+m3*mgxy[mgi-1]+m1*mgy[mgi-1]+m2)*4;
printf("%d   %ld %ld %ld %ld    %e %e",mgi-1,m1,num0[m1],pos0[m1],bondm[m1],sol0[m1],ival);getchar();
*/
}
/* End Vx,Vy,Vz,P "Prolongation" from coarser to finer level */


/* SOLVE MATRIX BY ITERATIVE METHOD */
int gausmat4(int am, long int mcmax, long int mcmin, int mgi)
/* am - number of iterations */
/* val0[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* mcmax - current line in num0[] */
/* mcmin - Number of iteration cycles */
{
/* Counters */
long int m1,m2,m3,m4,m10,m20,m30,m40;
/* Val Buffer */
double ival,ival1,ival2,ival3;
/**/
/**/
/**/
/* Seudel iteration */
/* Major line cycle */
for (m3=0;m3<am;m3++)
{
#pragma omp parallel for shared(fre1,num0,pos0,sol0,lin0,val1,nd,bondm) private(m30,m40,m1,m2,ival,ival1) firstprivate(m3,mcmin,mcmax,v0koef,p0koef,mgp,mgxy,mgx,mgy,mgz) schedule(static)
for (m30=0;m30<mgz[mgi]*mgx[mgi]*mgy[mgi];m30++)
for (m40=0;m40<4;m40++)
	{
	/* Pos in vx[] vy[] vz[] pr[] for the finer level */
	m1=(mgp[mgi]+m30)*4+m40;
if(num0[m1]) 
	{
	/* Use koef from fre1[], val1[] */
	ival=-fre1[m1];
	ival1=val1[pos0[m1]];
	for (m2=0;m2<num0[m1];m2++) 
		{
		/* Recalc increment */
		ival+=val1[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]];
		}
	if (!m40)
		{
		/* Continuity Eq. */
		if(!bondm[m1]) 
			{
			sol0[m1]-=ival*nd[m1/4]*p0koef;
			}
		else 
			{
			sol0[m1]-=ival;
			}
		}
	else
		{
		/* Stokes Eq. */
		if(!bondm[m1]) 
			{
			sol0[m1]-=ival/ival1*v0koef;
			}
		else 
			{
			sol0[m1]-=ival;
			}
		}
	}
	}
}
/**/
/**/
/**/
/* Residuals calculation cycle */
#pragma omp parallel for shared(fre0,fre1,num0,pos0,sol0,lin0,val1) private(m1,m2,ival) firstprivate(mcmin,mcmax) schedule(static)
for (m1=mcmin;m1<=mcmax;m1++)
if(num0[m1]) 
	{
	/* Use koef from fre1[], val1[] */
	ival=fre1[m1];
/*
printf("A %ld %e %e %e \n",m1,fre1[m1],sol0[m1],ival);
*/
	for (m2=0;m2<num0[m1];m2++) 
		{
		/* Recalc increment */
		ival-=val1[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]];
/*
printf("B %ld %ld %ld   %e %e %e %e \n",m1,m2,lin0[pos0[m1]+m2],val1[pos0[m1]+m2],sol0[lin0[pos0[m1]+m2]],val1[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]],ival);
*/
		}
	fre0[m1]=ival;
/*
if(m1>2600){printf("C %ld %e %e  ",m1,fre1[m1],fre0[m1]); getchar();}
*/
	}
return 0;
}
/* End SOLVE MATRIX BY ITERATIVE METHOD */



/* Vx,Vy,Vz,P residuals "Restriction" from finer to coarser level */
void vprestrict(int mgi)
/* mgi - Coarser multigrid level number */
{
/* Counters */
long int m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min;
/* Buffers */
double dx,dy,dz,ddx,ddy,ddz,ival,swt,ival1,ival2;
/**/
/**/
/**/
/* Clear Old Residuals */
#pragma omp parallel for shared(sol1,fre1) private(m1) firstprivate(mgi,mgp,mgn) schedule(static)
/*
*/
for (m1=mgp[mgi]*4;m1<(mgp[mgi]+mgn[mgi])*4;m1++)
	{
	fre1[m1]=0;
	sol1[m1]=0;
	}
/* Interpolate Add residuals */
#pragma omp parallel for shared(fre0,fre1,sol1,nd,bondm) private(m1,m2,m3,m4,m10,m20,m30,m40,m10min,m20min,m30min,dx,dy,dz,ddx,ddy,ddz,ival,swt,ival1,ival2) firstprivate(p1koef,v1koef,mgi,mgp,mgx,mgy,mgz,mgxy,mggx,mggy,mggz) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi-1];m3++)
for (m1=0;m1<mgx[mgi-1];m1++)
for (m2=0;m2<mgy[mgi-1];m2++)
	{
	/* Pos in vx[] vy[] vz[] pr[] for the finer level */
	m4=mgp[mgi-1]+m3*mgxy[mgi-1]+m1*mgy[mgi-1]+m2;
	/**/
	/**/
	/**/
	/* Vx interpolation */
	if(m2<mgy[mgi-1]-1 && m3<mgz[mgi-1]-1 && !bondm[m4*4+1])
		{
		/* Vx Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-2) m10min=mgx[mgi]-2;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-3) m20min=mgy[mgi]-3;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2+1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min+1])/2.0 && m20min>0) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-3) m30min=mgz[mgi]-3;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3+1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min+1])/2.0 && m30min>0) m30min-=1;
		/* Cell dimensions */
		dx= mggx[mgi][m10min+1]-mggx[mgi][m10min];
		dy=(mggy[mgi][m20min+2]-mggy[mgi][m20min])/2.0;
		dz=(mggz[mgi][m30min+2]-mggz[mgi][m30min])/2.0;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld  %e %e %e %ld",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz,bondm[m4*4+1]);getchar();
*/
		/* Reset velocity if required */
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30+1]-mggz[mgi-1][m3]-mggz[mgi-1][m3+1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]-mggx[mgi-1][m1])/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20+1]-mggy[mgi-1][m2]-mggy[mgi-1][m2+1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					if(!bondm[m40*4+1])
						{
						/* Calc wt */
						swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol1[m40*4+1]);getchar();}
*/
						/* Add Vx for current node */
/*
						fre1[m40*4+1]+=fre0[m4*4+1]/(nd[m4+1+mgxy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])*(nd[m40+1+mgxy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
						fre1[m40*4+1]+=fre0[m4*4+1]*(nd[m4+1+mgxy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])/(nd[m40+1+mgxy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
*/

						ival2=fre0[m4*4+1]*swt;
						fre1[m40*4+1]+=ival2;
						/* Add Wt for current node */
						sol1[m40*4+1]+=swt;
						}
					}
				}
			}
/*
printf("Vx %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vx[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* Vy interpolation */
	if(m1<mgx[mgi-1]-1 && m3<mgz[mgi-1]-1 && !bondm[m4*4+2])
		{
		/* Vy Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-3) m10min=mgx[mgi]-3;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1+1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min+1])/2.0 && m10min>0) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-2) m20min=mgy[mgi]-2;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-3) m30min=mgz[mgi]-3;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3+1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min+1])/2.0 && m30min>0) m30min-=1;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+2]-mggx[mgi][m10min])/2.0;
		dy= mggy[mgi][m20min+1]-mggy[mgi][m20min];
		dz=(mggz[mgi][m30min+2]-mggz[mgi][m30min])/2.0;
/*
printf("Vy %ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset velocity if required */
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30+1]-mggz[mgi-1][m3]-mggz[mgi-1][m3+1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10+1]-mggx[mgi-1][m1]-mggx[mgi-1][m1+1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]-mggy[mgi-1][m2])/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					if(!bondm[m40*4+2])
						{
						/* Calc wt */
						swt=ddx*ddy*ddz;
/*
printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vy[m40]);getchar();
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol1[m40*4+2]);getchar();}
*/
						/* Add Vy for current node */
/*
						fre1[m40*4+2]+=fre0[m4*4+2]/(nd[m4+mgy[mgi-1]+mgxy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])*(nd[m40+mgy[mgi]+mgxy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
						fre1[m40*4+2]+=fre0[m4*4+2]*(nd[m4+mgy[mgi-1]+mgxy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])/(nd[m40+mgy[mgi]+mgxy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
*/
						ival2=fre0[m4*4+2]*swt;
						fre1[m40*4+2]+=ival2;
						/* Add Wt for current node */
						sol1[m40*4+2]+=swt;
						}
					}
				}
			}
/*
printf("Vy %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vy[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* Vz interpolation */
	if(m1<mgx[mgi-1]-1 && m2<mgy[mgi-1]-1 && !bondm[m4*4+3])
		{
		/* Vz Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-3) m10min=mgx[mgi]-3;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1+1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min+1])/2.0 && m10min>0) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-3) m20min=mgy[mgi]-3;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2+1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min+1])/2.0 && m20min>0) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-2) m30min=mgz[mgi]-2;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+2]-mggx[mgi][m10min])/2.0;
		dy=(mggy[mgi][m20min+2]-mggy[mgi][m20min])/2.0;
		dz= mggz[mgi][m30min+1]-mggz[mgi][m30min];
/*
printf("%ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset velocity if required */
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]-mggz[mgi-1][m3])/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10+1]-mggx[mgi-1][m1]-mggx[mgi-1][m1+1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20+1]-mggy[mgi-1][m2]-mggy[mgi-1][m2+1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					if(!bondm[m40*4+3])
						{
						/* Calc wt */
						swt=ddx*ddy*ddz;
/*
printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,sol1[m40*4+3]);getchar();}
*/
						/* Add Vz for current node */
/*
						fre1[m40*4+3]+=fre0[m4*4+3]/(nd[m4+1+mgy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])*(nd[m40+1+mgy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
						fre1[m40*4+3]+=fre0[m4*4+3]*(nd[m4+1+mgy[mgi-1]]+nd[m4+1+mgy[mgi-1]+mgxy[mgi-1]])/(nd[m40+1+mgy[mgi]]+nd[m40+1+mgy[mgi]+mgxy[mgi]])*swt;
*/
						ival2=fre0[m4*4+3]*swt;
						fre1[m40*4+3]+=ival2;
						/* Add Wt for current node */
						sol1[m40*4+3]+=swt;
						}
					}
				}
			}
/*
printf("Vz %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,vz[m4]);getchar();
*/
		}
	/**/
	/**/
	/**/
	/* P  interpolation */
	if(m1 && m2 && m3 && !bondm[m4*4])
		{
		/* P Cell indexes */
		m10min=(m1-2)/2+2;
		if(m1<3) m10min=m1;
		if(m1>mgx[mgi-1]-4) m10min=mgx[mgi]-(mgx[mgi-1]-m1);
		if(m10min>mgx[mgi]-2) m10min=mgx[mgi]-2;
		if((mggx[mgi-1][m1]+mggx[mgi-1][m1-1])/2.0<(mggx[mgi][m10min]+mggx[mgi][m10min-1])/2.0 && m10min>1) m10min-=1;
		m20min=(m2-2)/2+2;
		if(m2<3) m20min=m2;
		if(m2>mgy[mgi-1]-4) m20min=mgy[mgi]-(mgy[mgi-1]-m2);
		if(m20min>mgy[mgi]-2) m20min=mgy[mgi]-2;
		if((mggy[mgi-1][m2]+mggy[mgi-1][m2-1])/2.0<(mggy[mgi][m20min]+mggy[mgi][m20min-1])/2.0 && m20min>1) m20min-=1;
		m30min=(m3-2)/2+2;
		if(m3<3) m30min=m3;
		if(m3>mgz[mgi-1]-4) m30min=mgz[mgi]-(mgz[mgi-1]-m3);
		if(m30min>mgz[mgi]-2) m30min=mgz[mgi]-2;
		if((mggz[mgi-1][m3]+mggz[mgi-1][m3-1])/2.0<(mggz[mgi][m30min]+mggz[mgi][m30min-1])/2.0 && m30min>1) m30min-=1;
		/* Cell dimensions */
		dx=(mggx[mgi][m10min+1]-mggx[mgi][m10min-1])/2.0;
		dy=(mggy[mgi][m20min+1]-mggy[mgi][m20min-1])/2.0;
		dz=(mggz[mgi][m30min+1]-mggz[mgi][m30min-1])/2.0;
/*
printf("%ld %ld %ld   %ld %ld %ld  %e %e %e",m1,m2,m3,m10min,m20min,m30min,dx,dy,dz);getchar();
*/
		/* Reset pressure if required */
		ival=0;
		/* Interpolate from 8 nodes */
		for (m30=m30min;m30<=m30min+1;m30++)
			{
			/* Normalized Z-distance calc */
			ddz=(mggz[mgi][m30]+mggz[mgi][m30-1]-mggz[mgi-1][m3]-mggz[mgi-1][m3-1])/2.0/dz;
			if (m30==m30min) ddz=-ddz;
			ddz=1.0-ddz;
			for (m10=m10min;m10<=m10min+1;m10++)
				{
				/* Normalized X-distance calc */
				ddx=(mggx[mgi][m10]+mggx[mgi][m10-1]-mggx[mgi-1][m1]-mggx[mgi-1][m1-1])/2.0/dx;
				if (m10==m10min) ddx=-ddx;
				ddx=1.0-ddx;
				for (m20=m20min;m20<=m20min+1;m20++)
					{
					/* Normalized Y-distance calc */
					ddy=(mggy[mgi][m20]+mggy[mgi][m20-1]-mggy[mgi-1][m2]-mggy[mgi-1][m2-1])/2.0/dy;
					if (m20==m20min) ddy=-ddy;
					ddy=1.0-ddy;
					/* Pos in ro[] nu[] for the coarcer level */
					m40=mgp[mgi]+m30*mgxy[mgi]+m10*mgy[mgi]+m20;
					if(!bondm[m40*4])
						{
						/* Calc wt */
						swt=ddx*ddy*ddz;
/*
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("Pr %ld %ld %ld   %ld %ld %ld %ld  %e %e %e   %e %e %e   %e %e",m1,m2,m3,m10,m20,m30,m40,dx,dy,dz,ddx,ddy,ddz,swt,sol1[m40*4]);getchar();}
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("%e %e   %e %e",mggy[mgi][m20],mggy[mgi][m20+1],mggy[mgi-1][m2],mggy[mgi-1][m2+1]);getchar();}
if(ddx<0 || ddy<0 || ddz<0 || ddx>1.0 || ddy>1.0 || ddz>1.0) {printf("%e %e",(mggy[mgi][m20]+mggy[mgi][m20+1])/2.0,(mggy[mgi-1][m2]+mggy[mgi-1][m2+1])/2.0);getchar();}
printf("Pr %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,pr[m40]);getchar();
*/
						/* Add P for current node */
/*
						fre1[m40*4]+=fre0[m4*4]*swt;
*/
						ival2=fre0[m4*4]*nd[m4]/nd[m40]*swt;
						fre1[m40*4]+=ival2;
						/* Add Wt for current node */
						sol1[m40*4]+=swt;
						}
					}
				}
			}
/*
printf("Pr %ld %ld %ld %ld  %e %e",m1,m2,m3,m4,ival,pr[m4]);getchar();
*/
		}
	}
/**/
/**/
/**/
/* Calc New  Residuals */
#pragma omp parallel for shared(sol0,sol1,fre1) private(m1) firstprivate(mgi,mgp,mgn) schedule(static)
/*
*/
for (m1=mgp[mgi]*4;m1<(mgp[mgi]+mgn[mgi])*4;m1++)
	{
	if(sol1[m1])
		{
/*
if(fre1[m1]){printf("RESID %ld %e %e %ld",m1,fre1[m1],sol1[m1],bondm[m1]);getchar();}
*/
		fre1[m1]/=sol1[m1];
		sol1[m1]=0;
		}
	sol0[m1]=0;
	}
}
/* End Vx,Vy,Vz,P residuals "Restriction" from finer to coarser level */



