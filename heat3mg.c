/* Left side, Right side or Err of Lagrangian Thermal Conductivity equation calc */
/* Explicit form  */
/* dT/dt = -1/Cp/ro*[dQx/dX + dQy/dY + dQz/Dz + dH/dt] */
/* Implicit form */
/* dT/dt + 1/Cp/ro*[dQx/dX + dQy/dYx + dQz/Dz] =  1/Cp/ro*(dH/dt] */
double heaterr(long int m1, long int m2, long int m3, int ynerr, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2 - node X,Y number */
/* ynerr - Calc mode: 0-Right, 1-Left, 2-Err */
{
/* Counters */
int n1;
long int v[7],p[8],wn1[100];
/* Buffer */
double mpb,ival=0,ival1=0,ival2=0,epsval,dpdx,dpdy,dpdz,x,y,z,dx,dy,dz;
/* Distances */
double xkf,ykf,zkf,rocpkf,vxyz1[100];
/**/
/**/
/**/
/* Distances Calc */
xkf=2.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
ykf=2.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
zkf=2.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/* T-Nodes num */
/*     2 6   */
/*     |/    */
/* 1---3---5 */
/*    /|     */
/*   0 4     */
v[3]=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
v[0]=v[3]-mgxy[mgi];
v[1]=v[3]-mgy[mgi];
v[2]=v[3]-1;
v[4]=v[3]+1;
v[5]=v[3]+mgy[mgi];
v[6]=v[3]+mgxy[mgi];
/**/
/**/
/**/
/* (Ro*Cp) calc */
rocpkf=ro[v[3]]*cp[v[3]];
/*
if(rocpkf<2e+6) rocpkf=2e+6;
*/
/**/
/**/
/**/
/* Mode of Calculation */
switch(ynerr)
	{
	/* RIght Side Calculation ---------------------- */
	case 0:
	/**/
	/* Lagrangian Thermal conductivity Equation in Simple Direct form use  */
	/* dT/dt = K/Cp/ro*(d2T/dX2+d2T/dY2)+1/Cp/ro*(dT/dX*dK/dX+dT/dY*dK/dY) + dH/dt/Cp/ro */
	/**/
	/* Heat Sources dH/dt/Cp/ro */
	/* Radioactive heat dH/dt, Wt/m3 */
	ival+=ht[v[3]];
	/**/
	/**/
	/**/
	/* Adiabatic, Shear heating terms koefficients */
	if(adiabyn==2 || frictyn)
		{
		/* P-SIGXX_EPSXX-Nodes num */
		/*    4----6 */
		/*   /|   /| */
		/*  / 5--/-7 */
		/* 0----2    */
		/* |    |    */
		/* 1----3    */
		p[0]=v[3];
		p[1]=p[0]+1;
		p[2]=p[0]+mgy[mgi];
		p[3]=p[2]+1;
		p[4]=p[0]+mgxy[mgi];
		p[5]=p[4]+1;
		p[6]=p[4]+mgy[mgi];
		p[7]=p[6]+1;
		/* Relative distances */
		x=(mggx[mgi][m1]-mggx[mgi][m1-1])/2.0*xkf;
		y=(mggy[mgi][m2]-mggy[mgi][m2-1])/2.0*ykf;
		z=(mggz[mgi][m3]-mggz[mgi][m3-1])/2.0*zkf;
		}
	/**/
	/**/
	/**/
	/* Vx, Vy in current node calc   */
	if(adiabyn) vxyzcalc1(mggx[mgi][m1],mggy[mgi][m2],mggz[mgi][m3],vxyz1,wn1);
	/**/
	/**/
	/**/
	/* Correct adiabatic term alphV*T*(vX*dP/dX+vY*dP/dY)/Cp/ro */
	if(adiabyn==2 && ro[v[3]]>2e+3 && mggy[mgi][m2]<ysize-stp100)
		{
		/* P-Nodes num */
		/*    4----6 */
		/*   /|   /| */
		/*  / 5--/-7 */
		/* 0----2    */
		/* |    |    */
		/* 1----3    */
		/* P derivatives */
		dpdx=xkf*((1.0-y)*(1.0-z)*(pr[p[2]]-pr[p[0]])+y*(1.0-z)*(pr[p[3]]-pr[p[1]])+(1.0-y)*z*(pr[p[6]]-pr[p[4]])+y*z*(pr[p[7]]-pr[p[5]]));
		dpdy=ykf*((1.0-x)*(1.0-z)*(pr[p[1]]-pr[p[0]])+x*(1.0-z)*(pr[p[3]]-pr[p[2]])+(1.0-x)*z*(pr[p[5]]-pr[p[4]])+x*z*(pr[p[7]]-pr[p[6]]));
		dpdz=zkf*((1.0-x)*(1.0-y)*(pr[p[4]]-pr[p[0]])+x*(1.0-y)*(pr[p[6]]-pr[p[2]])+(1.0-x)*y*(pr[p[5]]-pr[p[1]])+x*y*(pr[p[7]]-pr[p[3]]));
		ival+=et[v[3]]*tk0[v[3]]*(vxyz1[0]*dpdx+vxyz1[1]*dpdy+vxyz1[2]*dpdz);
		/* Adiabate computing */
		if(1==0 && timesum<(3.15576e+7*1e+3))
			{
			mpb =(1.0-x)*(1.0-y)*(1.0-z)*pr[p[0]];
			mpb+=(1.0-x)*(    y)*(1.0-z)*pr[p[1]];
			mpb+=(    x)*(1.0-y)*(1.0-z)*pr[p[2]];
			mpb+=(    x)*(    y)*(1.0-z)*pr[p[3]];
			mpb+=(1.0-x)*(1.0-y)*(    z)*pr[p[4]];
			mpb+=(1.0-x)*(    y)*(    z)*pr[p[5]];
			mpb+=(    x)*(1.0-y)*(    z)*pr[p[6]];
			mpb+=(    x)*(    y)*(    z)*pr[p[7]];
			ival+=et[v[3]]*tk0[v[3]]*mpb/(3.15576e+7*1e+3);
			}
/*
printf("%ld %ld   %e %e %e %e   %e %e %e\n",m1,m2,e,n,dpdx,dpdy,et[v[3]],ival1,ival2); getchar();
*/
		}
	/* Simplified adiabatic term alphV*T/Cp*(gY*vY+gX*vX+gZ*vZ) */
	if(adiabyn==1 && ro[v[3]]>2e+3 && mggy[mgi][m2]<ysize-stp100)
		{
		ival+=et[v[3]]*tk0[v[3]]*ro[v[3]]*(vxyz1[0]*GXKOEF+vxyz1[1]*GYKOEF+vxyz1[2]*GZKOEF);
		}
	/**/
	/**/
	/**/
	/* Viscouse friction dH/dt=NU*EPSik*dVi/dXk */
	/* Viscouse friction dH/dt=2NU*EPSik*EPSik=SIGik*EPSik, where SIGik=2NU*EPSik, EPSik=1/2(dVxi/dXk+dVxk/dXi) */
	if(frictyn && ro[v[3]]>2e+3 && mggy[mgi][m2]<ysize-stp100)
		{
		/* SIGxx*EPSxx+SIGyy*EPSyy+SIGzz*EPSzz
		/* SIGxx-SIGyy-SIGzz-Nodes num */
		/*    4----6 */
		/*   /|   /| */
		/*  / 5--/-7 */
		/* 0----2    */
		/* |    |    */
		/* 1----3    */
		epsval =(1.0-x)*(1.0-y)*(1.0-z)*(sxx[p[0]]*exx[p[0]]+syy[p[0]]*eyy[p[0]]+szz[p[0]]*ezz[p[0]]);
		epsval+=(1.0-x)*(    y)*(1.0-z)*(sxx[p[1]]*exx[p[1]]+syy[p[1]]*eyy[p[1]]+szz[p[1]]*ezz[p[1]]);
		epsval+=(    x)*(1.0-y)*(1.0-z)*(sxx[p[2]]*exx[p[2]]+syy[p[2]]*eyy[p[2]]+szz[p[2]]*ezz[p[2]]);
		epsval+=(    x)*(    y)*(1.0-z)*(sxx[p[3]]*exx[p[3]]+syy[p[3]]*eyy[p[3]]+szz[p[3]]*ezz[p[3]]);
		epsval+=(1.0-x)*(1.0-y)*(    z)*(sxx[p[4]]*exx[p[4]]+syy[p[4]]*eyy[p[4]]+szz[p[4]]*ezz[p[4]]);
		epsval+=(1.0-x)*(    y)*(    z)*(sxx[p[5]]*exx[p[5]]+syy[p[5]]*eyy[p[5]]+szz[p[5]]*ezz[p[5]]);
		epsval+=(    x)*(1.0-y)*(    z)*(sxx[p[6]]*exx[p[6]]+syy[p[6]]*eyy[p[6]]+szz[p[6]]*ezz[p[6]]);
		epsval+=(    x)*(    y)*(    z)*(sxx[p[7]]*exx[p[7]]+syy[p[7]]*eyy[p[7]]+szz[p[7]]*ezz[p[7]]);
		/**/
		/* 2*(SIGxy*EPSxy+SIGxz*EPSxz+SIGyz*EPSyz)
		/* SIGxy-SIGxz-SIGyz-T-Nodes num */
		/*     SIGxz2 SIGxy3  */
		/*          |/        */
		/* SIGyz1--T3--SIGyz3 */
		/*         /|         */
		/*   SIGxy0 SIGxz3    */
		/* SIGxy*EPSxy */
		epsval+=2.0*((1.0-z)*(sxy[v[0]]*exy[v[0]])+z*(sxy[v[3]]*exy[v[3]]));
		/* SIGxz*EPSxz */
		epsval+=2.0*((1.0-y)*(sxz[v[2]]*exz[v[2]])+y*(sxz[v[3]]*exz[v[3]]));
		/* SIGyz*EPSyz */
		epsval+=2.0*((1.0-x)*(syz[v[1]]*eyz[v[1]])+x*(syz[v[3]]*eyz[v[3]]));
		/**/
		ival+=epsval;
		/**/
/*
printf("A %ld %ld   %e\n",m1,m2,epsval); getchar();
printf("%ld %ld %ld %e %e %e %e",m1,m2,v[3],tk0[v[3]],kt[v[3]],cp[v[3]],ro[v[3]]); getchar();
*/
		}
	/**/
	/**/
	/**/
	/* Save Right Part for next calculations */
	tk1[v[3]]=ival;
	/**/
	/**/
	/**/
	/* Rate of i ternal heating return */
	return ival/rocpkf;
	/**/
	/**/
	/**/
	/* Left Side Calculation ---------------------- */
	case 1:
	/* Implicit form */
	/* dT/dt + 1/Cp/ro*[dQx/dX + dQy/dY + dQz/dZ] =  1/Cp/ro*(dH/dt] */
	/**/
	/* Add Right part  1/Cp/ro*(dH/dt) */
	if(mgi==0) ui1[0]=tk1[v[3]]+tk0[v[3]]/timestep*rocpkf; else ui1[0]=0;
	/**/
	/* Add Left part */
	/* dT/dt */
	ui1[1]=rocpkf/timestep;
	un1[1]=v[3];
	/* T-k-Nodes num */
	/*     2 6   */
	/*     |/    */
	/* 1---3---5 */
	/*    /|     */
	/*   0 4     */
	/* + [dQx/dX] where Qx=-k*dT/dX */
	/* dQx/dX=(Qx35-Qx13)/(dX135)=((k13)*(T3-T1)/(dX13) - (k35)*(T5-T3)/(dX35))/(dX135) */
	ival=(kt[v[1]]+kt[v[3]])/2.0*xkf/(mggx[mgi][m1]-mggx[mgi][m1-1]);
	ui1[2]=-ival; un1[2]=v[1];
	ui1[1]+=ival;
	ival=(kt[v[3]]+kt[v[5]])/2.0*xkf/(mggx[mgi][m1+1]-mggx[mgi][m1]);
	ui1[1]+=ival;
	ui1[3]=-ival; un1[3]=v[5];
	/* + [dQy/dY] where Qy=-k*dT/dY */
	/* dQy/dY=(Qy34-Qy23)/(dY234)=((k23)*(T3-T2)/(dY23) - (k34)*(T4-T3)/(dY34))/(dY234) */
	ival=(kt[v[2]]+kt[v[3]])/2.0*ykf/(mggy[mgi][m2]-mggy[mgi][m2-1]);
	ui1[4]=-ival; un1[4]=v[2];
	ui1[1]+=ival;
	ival=(kt[v[3]]+kt[v[4]])/2.0*ykf/(mggy[mgi][m2+1]-mggy[mgi][m2]);
	ui1[1]+=ival;
	ui1[5]=-ival; un1[5]=v[4];
	/* + [dQz/dz] where Qz=-k*dT/dZ */
	/* dQz/dZ=(Qz36-Qz03)/(dZ036)=((k03)*(T3-T0)/(dZ03) - (k36)*(T6-T3)/(dX36))/(dX036) */
	ival=(kt[v[0]]+kt[v[3]])/2.0*zkf/(mggz[mgi][m3]-mggz[mgi][m3-1]);
	ui1[6]=-ival; un1[6]=v[0];
	ui1[1]+=ival;
	ival=(kt[v[3]]+kt[v[6]])/2.0*zkf/(mggz[mgi][m3+1]-mggz[mgi][m3]);
	ui1[1]+=ival;
	ui1[7]=-ival; un1[7]=v[6];
	/* Total number of positions */
	un1[0]=7;
	/**/
/*
if(mgi>1){printf("%d  %ld %ld %ld   %ld   %e %e %e %e ",mgi,m1,m2,m3,v[3],kt[v[3]],ro[v[3]],cp[v[3]],rocpkf); getchar();}
if(mgi>0){for (n1=0;n1<=13;n1++){printf("%d  %ld %ld %ld  %d %ld %e",mgi,m1,m2,m3,n1,un[n1],ui[n1]); getchar();}}
*/
	return 0;
	/**/
	/**/
	/**/
	/* Error Calculation ---------------------- */
	case 2:
	/**/
	/* Left part Add */
	/* dT/dt */
	ival1=(tk[v[3]]-tk0[v[3]])/timestep;
	/**/
	/* T-k-Nodes num */
	/*     2 6   */
	/*     |/    */
	/* 1---3---5 */
	/*    /|     */
	/*   0 4     */
	/* + 1/Cp/ro*[dQx/dX] where Qx=-k*dT/dX */
	/* dQx/dX=(Qx35-Qx13)/(dX135)=((k13)*(T3-T1)/(dX13) - (k35)*(T5-T3)/(dX35))/(dX135) */
	ival=(kt[v[1]]+kt[v[3]])/2.0*xkf/(mggx[mgi][m1]-mggx[mgi][m1-1])/rocpkf;
	ival1+=(tk[v[3]]-tk[v[1]])*ival;
	ival=(kt[v[3]]+kt[v[5]])/2.0*xkf/(mggx[mgi][m1+1]-mggx[mgi][m1])/rocpkf;
	ival1-=(tk[v[5]]-tk[v[3]])*ival;
	/* + 1/Cp/ro*[dQy/dY] where Qy=-k*dT/dY */
	/* dQy/dY=(Qy34-Qy23)/(dY234)=((k23)*(T3-T2)/(dY23) - (k34)*(T4-T3)/(dY34))/(dY234) */
	ival=(kt[v[2]]+kt[v[3]])/2.0*ykf/(mggy[mgi][m2]-mggy[mgi][m2-1])/rocpkf;
	ival1+=(tk[v[3]]-tk[v[2]])*ival;
	ival=(kt[v[3]]+kt[v[4]])/2.0*ykf/(mggy[mgi][m2+1]-mggy[mgi][m2])/rocpkf;
	ival1-=(tk[v[4]]-tk[v[3]])*ival;
	/* + 1/Cp/ro*[dQz/dz] where Qz=-k*dT/dZ */
	/* dQz/dZ=(Qz36-Qz03)/(dZ036)=((k03)*(T3-T0)/(dZ03) - (k36)*(T6-T3)/(dX36))/(dX036) */
	ival=(kt[v[0]]+kt[v[3]])/2.0*zkf/(mggz[mgi][m3]-mggz[mgi][m3-1])/rocpkf;
	ival1+=(tk[v[3]]-tk[v[0]])*ival;
	ival=(kt[v[3]]+kt[v[6]])/2.0*zkf/(mggz[mgi][m3+1]-mggz[mgi][m3])/rocpkf;
	ival1-=(tk[v[6]]-tk[v[3]])*ival;
	/**/
	/* Right part Add */
	ival1-=tk1[v[3]]/rocpkf;
	return ival1;
	}
return 0;
}
/* Left side, Right side or Err of Thermal Conductivity equation calc */


/* Calculation of new T after time step */
void titerate(int m0)
/* m0 - Circle number */
{
/* Counters */
long int m2prn=0,m1prn=0;
long int m1,m2,m3,m4,mcmax,mcmax1,m1min,m1max,m2min,m2max;
int n1,n2,mgi;
/**/
/* Val buffer */
double ival,maxtkrate,tktimestep=1e+30,tktimesum,tktimesum0,celdx,celdy,swt1;
double eps1[100],wi1[100],ival2;
long int wn1[100];
long int un1[MAXPOS],pos0cur1;
double ui1[MAXPOS];
/* Err koef */
double heatsum=0,heatnum=0,bondsum=0,bondnum=0,mkt,mtk,rocpkf,dtk,dtk1,mintk,maxtk;
/**/
/**/
/**/
/* Right par of equation T step limit definition */
/* Save Old temperatures */
/* Save TK */
maxtk=-1e+30;
mintk=1e+30;
for (m1=0;m1<nodenum;m1++)
	{
	tk0[m1]=tk[m1];
	maxtk=MAXV(maxtk,tk[m1]);
	mintk=MINV(mintk,tk[m1]);
/*
if(tk[m1]>1640.0){printf("%ld %e %e   %e %e %e %e %e",m1,tk[m1],pr[m1],cp[m1],kt[m1],ro[m1],ht[m1],et[m1]);getchar();}
printf("%ld %e %e   %e %e %e %e",m1,tk[m1],pr[m1],cp[m1],kt[m1],ro[m1],ht[m1]);getchar();
*/
	}
printf("TEMPERATURE MIN = %e MAX = %e \n",mintk,maxtk);
/*
return ;
*/
/**/
mgi=0;
maxtkrate=0;
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Pos in tk[], kt[], cp[], etc. */
	mcmax=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
/*
if(tk[mcmax]>1640.0){printf("%ld %ld %ld %ld %e %e   %e %e %e %e %e",m1,m2,m3,mcmax,tk[mcmax],pr[mcmax],cp[mcmax],kt[mcmax],ro[mcmax],ht[mcmax],et[mcmax]);getchar();}
*/
	/**/	
        /* Cur Line Num in bondm[] */
	mcmax1=mcmax+nodenum8;
        if(bondm[mcmax1]==0)
                {
                /* Max dT/dt definition */
                ival=heaterr(m1,m2,m3,0,mgi,eps1,un1,ui1);
		ival=ABSV(ival);
                maxtkrate=MAXV(maxtkrate,ival);
/*
printf("%ld %ld %ld %e",m1,m2,m3,ival); getchar();
*/
                }
	}
/**/
/**/
/*
printf("\n maxtkrate = %e maxtkstep = %e !!!\n",maxtkrate,maxtkstep); getchar();
*/
/**/
/* Time step calc */
if (maxtkrate)
	{
	tktimestep=maxtkstep/maxtkrate;
	if (printmod) printf("\n !!! MAX VALID TIME STEP FOR TEMPERATURE %e YEAR !!!\n",tktimestep/3.15576e+7);
	}
/*
getchar();
*/
/* Check timestep */
if (timestep==0 || tktimestep==0) return;
/**/
/* Thermal Solution Cycle ------------------------------------------ */
/* Save TK */
for (m1=0;m1<nodenum;m1++)
	{
	tk2[m1]=tk[m1];
	}
/**/
/**/
/**/
tktimesum=0;
tktimesum0=timestep;
do
{
/* Timestep for temperature selection */
timestep=MINV(tktimesum0,tktimestep);
if (timestep>(tktimesum0-tktimesum)) timestep=tktimesum0-tktimesum;
/**/
/* Save TK */
for (m1=0;m1<nodenum;m1++)
	{
	tk0[m1]=tk[m1];
	}
/**/
/* Add Matrix by TK Equations */
if (printmod) printf("Adding Matrix By TK ...\n");
/* TK - Node Cycle */
leftnum=0;
rightnum=0;
/* Reload T OLD Results to sol0[] */
for (m1=0;m1<nodenum;m1++)
	{
	sol0[m1]=tk[m1];
	}
/**/
/**/
/**/
/* Add Matrix by T Equations */
/* Node  Cycle */
for(mgi=0;mgi<=multinumt;mgi++)
{
#pragma omp parallel for shared(bondm) private(m1,m2,m3,mcmax,mcmax1,pos0cur1,eps1,un1,ui1) firstprivate(mgi,multinum,mgp,mgxy,mgy) schedule(static)
/*
*/
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Cur Pos in sol0[] */
	mcmax=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
        /* Cur Line Num in bondm[] */
	mcmax1=mcmax+nodenum8;
	/* Cur position in val1[] */
	pos0cur1=mcmax*7;
	/**/
/*
if(mgi>0 && m3>0 ){printf("%ld %ld %ld   %ld  %ld %d",m1,m2,m3,mcmax,mcmax1,bondm[mcmax1]);getchar();}
printf("TEMPER %ld",pos0cur);getchar();
*/
	/* Add T-Equations --------------------------------------------------- */
        if(bondm[mcmax1]==0)
                {
                /* Add Heat Eq */
		heaterr(m1,m2,m3,1,mgi,eps1,un1,ui1);
		gausmat2(1,mcmax,pos0cur1,un1,ui1);
                }
	else
		{
		/* Add T-Boundary */
		tbonderr(mcmax,0,eps1,un1,ui1);
		/* Plume BC deactivation */
		if(timesum>timebeg && m2==mgy[mgi]-1)
			{
			if(timesum>=timeend)
				{
				ui1[0]=tempend;
				}
			else
				{
				ui1[0]=ui1[0]+(tempend-ui1[0])/(timeend-timebeg)*(timesum-timebeg);
				}
                	}
		gausmat2(1,mcmax,pos0cur1,un1,ui1);
		}
	}
}
/* End Add Matrix by T Equations */
mcmax=mgp[multinum]+mgn[multinum];
pos0cur1=mcmax*7;
mcmax1=mcmax;
if(printmod) printf("mcmax = %ld mcmax1 = %ld Pos0cur1 = %ld \n",mcmax,mcmax1,pos0cur1);
/**/
/**/
/**/
/* T ITERATION CYCLE =================================== */
do
{
/* Mulrigrid V cycle */
if(0==0)
for(n1=0;n1<multicyct;n1++)
	{
	printf("Iter%d/%d: ",n1+1,multicyct);
	/* Multigrid cylcle */
	mgi=0;
	gausmatt(multittt[0][mgi],mgp[mgi]+mgn[mgi]-1,mgp[mgi]);
	if(printmod) printf("Simple%d/%d ",mgi,multittt[0][mgi]);
	if(printmod) printf("\n");
	}
mgi=0;
/*
gausmatt(4,mgp[mgi]+mgn[mgi]-1,mgp[mgi]);
*/
/* End Mulrigrid V cycle */
/**/
/**/
/**/
/* Reload T Results */
/* T max-min definition */
mintk=1e+30;maxtk=-1e+30;
/* Reload T NEW Results from sol0[] */
for (m1=0;m1<nodenum;m1++)
	{
	tk[m1]=sol0[m1];
	mintk=MINV(mintk,tk[m1]);
	maxtk=MAXV(maxtk,tk[m1]);
	}
/* End Reload T Results */
/**/	
/**/	
/**/	
/* Check Error for T equations */
/* Err koef */
bondsum=0;bondnum=0;
heatsum=0;heatnum=0;
/* Node  Cycle */
mgi=0;
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Cur Pos in sol0[] */
	mcmax=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
        /* Cur Line Num in bondm[] */
	mcmax1=mcmax+nodenum8;
	/**/
	/**/
	/* Check Err T-Equations --------------------------------------------------- */
        if(bondm[mcmax1]==0)
                {
                /* Check Err Heat Eq */
		ival=heaterr(m1,m2,m3,2,mgi,eps1,un1,ui1);
/*
{printf("%ld %ld %ld   %ld  %ld %d   %e %e   %e %e",m1,m2,m3,mcmax,mcmax1,bondm[mcmax1],tk[mcmax],ival,maxtk,timestep);getchar();}
*/
		ival*=timestep;
		heatsum+=ival*ival;
		heatnum+=1.0;
                }
	else
		{
		/* Add vY-Boundary */
		ival=tbonderr(mcmax,1,eps1,un1,ui1);
		bondsum+=ival*ival;
		bondnum+=1.0;
		}
	}
heatsum=pow(heatsum/heatnum,0.5);
bondsum=pow(bondsum/bondnum,0.5);
/* End Check Error for T equations */
/**/
/**/
/**/
/* Print Results */
m1prn=m1prn+1;
m2prn=m2prn+1;
if (printmod && printmod<=m1prn)
	{
	printf("\n KRUG %d Cycle %ld \n ",m0+1,m2prn);
	printf("THERMAL STEP = %e YEARS    THERMAL TIME = %e YEARS \n",timestep/3.15576e+7,tktimesum/3.15576e+7);
	printf("TEMPERATURE    : min = %e max = %e \n",mintk,maxtk);
	printf("TERMAL EQUATION: num = %e err = %e \n",heatnum,heatsum);
	printf("BOUND T        : num = %e err = %e \n",bondnum,bondsum);
	printf("\n");
	m1prn=0;
/*
getchar();
*/
	}
/**/
/**/
/**/
}
while(heatsum>HEATMIN);
/* End T ITERATION CYCLE =================================== */
/**/
/**/
/**/
/* Add thermal step */
tktimesum+=timestep;
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	printf("\n KRUG %2d \n",m0+1);
	printf("THERMAL STEP = %e YEARS    THERMAL TIME = %e YEARS \n",timestep/3.15576e+7,tktimesum/3.15576e+7);
	printf("TEMPERATURE    : min = %e max = %e \n",mintk,maxtk);
	printf("TERMAL EQUATION: num = %e err = %e \n",heatnum,heatsum);
	printf("BOUND T        : num = %e err = %e \n",bondnum,bondsum);
	printf("\n");
/*
getchar();
*/
	}
/**/
}
while(tktimesum<tktimesum0);
/* Restore timestep */
timestep=tktimesum0;
/* End Thermal Solution Cycle ------------------------------------------ */
/**/
/**/
/**/
/* Recalc temperature for markers */
/* Clear nodes wt, save increment */
#pragma omp parallel for shared(sol0,sol1) private(m1) firstprivate(nodenum) schedule(static)
/*
*/
for (m1=0;m1<nodenum;m1++)
	{
	sol0[m1]=sol1[m1]=0;
/*
if(tk[m1]>1640.0){printf("%ld %e %e   %e %e %e %e %e",m1,tk[m1],pr[m1],cp[m1],kt[m1],ro[m1],ht[m1],et[m1]);getchar();}
*/
	}
/**/
/* Reset T for fluid markers */
#pragma omp parallel for shared(markx,marky,markz,markt,markk,tk,tk0,tk1,tk2) private(m3,eps1,wn1,wi1) firstprivate(marknum,xsize,ysize,zsize,gx,gy,gz) schedule(static)
/*
*/
for (m3=0;m3<marknum;m3++) 
/* Check markers out of grid */
if (markt[m3]>=50 && markt[m3]<100 && markx[m3]>0 && marky[m3]>0 && markz[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && (double)(markz[m3])<zsize)
	{
	/* Interpolate temperature */
	tkcalc1((double)(markx[m3]),(double)(marky[m3]),(double)(markz[m3]),eps1,wn1,wi1);
	/**/
	/* Reset marker temperature for newly coming markers */
	markk[m3]=(float)(eps1[2]);
	}
/**/
/* Recalc marker T + Diffusion, Add  nodes wt */
#pragma omp parallel for shared(markx,marky,markz,markt,markk,sol0,sol1,tk,tk0,tk1,tk2,cp,kt,ro) private(n1,ival,ival2,m3,m4,eps1,wn1,wi1,mkt,mtk,rocpkf,dtk,dtk1) firstprivate(heatdif,marknum,xsize,ysize,zsize,zdeep,tdeep,timestep,gx,gy,gz) schedule(static)
/*
*/
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if (markt[m3]<50 && markx[m3]>0 && marky[m3]>0 && markz[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && (double)(markz[m3])<zsize)
	{
	/* Interpolate temperature */
	tkcalc1((double)(markx[m3]),(double)(marky[m3]),(double)(markz[m3]),eps1,wn1,wi1);
	/**/
	/* Reset marker temperature for newly coming markers and for markers at the boundary */
/*
	if(markk[m3]<=0 || wn[0]<=0 || wn[0]>=xnumx-2 || wn[1]<=0 || wn[1]>=ynumy-2 || wn[2]<=0 || wn[2]>=znumz-2) 
*/
	if(markk[m3]<=0) 
		{
		markk[m3]=(float)(eps1[3]);
		/* Marker coming from the depth */
		if(marky[m3]>zdeep && markk[m3]<tdeep) markk[m3]=tdeep;
/*
if(markk[m3]>1700.0){printf("%ld %d %e %e %e  %e",m3,markt[m3],markx[m3],marky[m3],markz[m3],markk[m3]);getchar();}
		if(markk[m3]<tdeep) markk[m3]=tdeep;
*/
		}
	/**/
/*
printf("DIFFUSION1 %ld %e",m3,markk[m3]); getchar();
*/
	/* Numerical diffusion add to markers -------------*/
	if(heatdif)
		{
		/**/
		/* Interpolate of k, ro, Cp from 8-Nodes using interpolation weights */
		mkt=rocpkf=0;
		mtk=tk[wn1[3]];
		for (n1=0;n1<8;n1++)
			{
			/* Current node num, wt */
			m4=wn1[3+n1];
			ival=wi1[3+n1];
			mkt+=kt[m4]*ival;
			rocpkf+=ro[m4]*cp[m4]*ival;
			mtk=MAXV(mtk,tk[m4]);
/*
printf("CP %ld %e  %ld %e %e %e %e %e %e",m3,markk[m3],m4,mkt,mcp,mro,kt[m4],cp[m4],ro[m4]); getchar();
*/
			}
		/**/
		/* Correction for the boundary region */
		if((wn1[0]<=0 || wn1[0]>=xnumx-2 || wn1[1]<=0 || wn1[1]>=ynumy-2 || wn1[2]<=0 || wn1[2]>=znumz-2) && (markk[m3]+eps1[2]-eps1[3])>mtk) markk[m3]=mtk-eps1[2]+eps1[3];
		/**/
		/* Calc difference in marker temperature */
		dtk=eps1[3]-(double)(markk[m3]);
		/**/
		/* Calc, check diffusion changes in marker temperature: dT/dt=k/ro/cp*(d2T/dX2+d2T/dY2) */
		/* dT=k/ro/cp*(Tgrid-2Tmark+Tgrid)*(1/(xstep)^2+1/(ystep)^2)*dt */
/*
		if(rocpkf<2e+6) rocpkf=2e+6;
*/
		dtk1=-heatdif*mkt/rocpkf*2.0*(1.0/wi1[0]/wi1[0]+1.0/wi1[1]/wi1[1]+1.0/wi1[2]/wi1[2])*timestep;
		if(dtk1<-150.0) dtk1=-150.0;
		dtk1=dtk*(1.0-exp(dtk1));
		/**/
 	      	/* Diffuse Marker Temperature */
		markk[m3]+=(float)(dtk1);
/*
printf("DIFFUSION %ld %e  %ld %e %e %e %e %e",m3,markk[m3],m4,mkt,mcp,mro,dtk,dtk1); getchar();
*/
		/**/
		/* Wt for 8 nodes add */
		for (n1=0;n1<8;n1++)
			{
			/* Current node num, wt */
			m4=wn1[3+n1];
			ival=wi1[3+n1];
			/**/
			/* Add Node wt, T */
/*
			ival*=rocpkf;
*/
			ival2=dtk1*ival;
			sol1[m4]+=ival2;
			sol0[m4]+=ival;
			}
		}
	/* End Numerical diffusion add to markers -------------*/
	/**/
	/**/
	/**/
	/* Change marker temperature after solution */
	markk[m3]+=(float)(eps1[2]-eps1[3]);
	/**/
/*
if(markk[m3]>1700.0){printf("%ld %d %e %e %e  %e %e %e",m3,markt[m3],markx[m3],marky[m3],markz[m3],markk[m3],eps[2],eps[3]);getchar();}
*/
	/**/
	/**/
/*
printf("DIFFUSION2 %ld %e",m3,markk[m3]); getchar();
*/
	}
/**/
/* Numerical antidiffusion add to markers -------------*/
if(heatdif)
	{
	/* Recalc changes in nodes T */
#pragma omp parallel for shared(sol0,sol1) private(m1) firstprivate(nodenum) schedule(static)
/*
*/
	for (m1=0;m1<nodenum;m1++)
		{
/*
printf("A %ld %e %e %e",m1,tk0[m1],sol0[m1],sol1[m1]); getchar();
*/
		if(sol0[m1]) 
			{
			/* Averaged Changes in temperature due to smoothing */
			sol1[m1]/=sol0[m1];
			}
/*
printf("B %ld %e %e %e",m1,tk0[m1],sol0[m1],sol1[m1]); getchar();
*/
		}
	/**/
	/* Recalc marker T + Antidiffusion) */
#pragma omp parallel for shared(markx,marky,markz,markt,markk,sol1,tk,tk0,tk1,tk2) private(n1,ival,ival2,m3,m4,eps1,wn1,wi1) firstprivate(marknum,xsize,ysize,zsize,gx,gy,gz) schedule(static)
/*
*/
	for (m3=0;m3<=marknum;m3++) 
	/* Check markers out of grid */
	if (markt[m3]<50 && markx[m3]>0 && marky[m3]>0 && markz[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && (double)(markz[m3])<zsize)
		{
		/* Interpolate temperature */
		tkcalc1((double)(markx[m3]),(double)(marky[m3]),(double)(markz[m3]),eps1,wn1,wi1);
		/**/
/*
printf("C %ld %e",m3,markk[m3]); getchar();
*/
		/* Wt for nodes calc, add */
		/* Wt for 8 nodes add */
		ival2=0;
		for (n1=0;n1<8;n1++)
			{
			/* Current node num, wt */
			m4=wn1[3+n1];
			ival=wi1[3+n1];
			/**/
       		 	/* Antidiffuse Marker Temperature */
			ival2-=sol1[m4]*ival;
/*
printf("D %ld %ld %e %e ",m3,m4,markk[m3],sol1[m4]); getchar();
*/
			/**/
			}
       		 /* Antidiffuse Marker Temperature */
		markk[m3]+=(float)(ival2);
/*
printf("E %ld %e",m3,markk[m3]); getchar();
*/
		}
	}
/* End Numerical antidiffusion add to markers -------------*/
/**/
/* End T interpolation ------------------------ */
/**/
/**/
/**/
}
/* Calculation of new T after time step */


/* Left side or Err for T Boundary Condition Equation */ 
/* Ai=CONST+KOEF*An */
double tbonderr(long int mcmax, int ynerr, double *eps1, long int *un1, double *ui1)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftb=0;
int n1;
/* Cur Line Num in bondm[] */
long int mcmax1=mcmax+nodenum8;
/**/
/**/
/**/
/* Error Calc */
if (ynerr==1)
	{
	/* Add Const */
	leftb=sol0[mcmax]-bondv1[bondm[mcmax1]][0];
	/* Add Koef */
	if(bondn1[bondm[mcmax1]]) 
		{
		leftb-=bondv1[bondm[mcmax1]][1]*sol0[bondn1[bondm[mcmax1]]-1];
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
ui1[0]=bondv1[bondm[mcmax1]][0];
un1[1]=mcmax;
ui1[1]=1.0;
/* Add T PAR1,PAR2,PAR3 */
if(bondn1[bondm[mcmax1]]) 
	{
	un1[0]+=1;
	un1[un1[0]]=bondn1[bondm[mcmax1]]-1;
	ui1[un1[0]]=-bondv1[bondm[mcmax1]][1];
	}
/*
printf("Bo %ld %e %ld %e \n",mcmax1,bondv1[mcmax1][0],bondn1[mcmax1],bondv1[mcmax1][1]);getchar();
printf("Bo %ld %e %e \n",mcmax,sol0[mcmax],leftb);getchar();
for(n1=0;n1<3;n1++)printf("%e %d \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for T Boundary Condition Equation */ 



/* tk[] Correct for boundary conditions */
void tkrecalc()
{
/* Counters */
long int m1,m2,m3,mgi,mcmax,mcmax1;
int n1;
/**/
/**/
/**/
if (printmod) printf("Correcting Temperature by TK Boundaries ...");
/* Save TK */
for (m1=0;m1<nodenum;m1++)
	{
	tk1[m1]=tk[m1];
	}
/**/
/**/
/**/
/* T-Node  Cycle */
mgi=0;
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	/* Cur Pos in sol0[] */
	mcmax=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
        /* Cur Line Num in bondm[] */
	mcmax1=mcmax+nodenum8;
	/**/
	/**/
	/* Check Err T-Equations --------------------------------------------------- */
        if(bondm[mcmax1]>0)
                {
		/* Add Const */
		tk[mcmax]=bondv1[bondm[mcmax1]][0];
		/* Add Koef */
		if(bondn1[bondm[mcmax1]]) 
			{
			tk[mcmax]+=bondv1[bondm[mcmax1]][1]*tk1[bondn1[bondm[mcmax1]]-1];
			}
                }
/*
printf("%ld %ld %ld   %ld %ld %e   %d %e",m1,m2,m3,mcmax,mcmax1,tk[mcmax],bondm[mcmax1],bondv1[bondm[mcmax1]][0]);getchar();
printf("%e %d",bondv1[bondm[mcmax1]][1],bondn1[bondm[mcmax1]]-1);getchar();
printf("%ld %ld %ld   %ld %ld %e",m1,m2,m3,mcmax,mcmax1,tk[mcmax]);getchar();
*/
	}
/**/
if (printmod) printf("OK!\n");
/* End Boundary Solution Cycle ------------------------------------------ */
/**/
}
/* End tk[] Correct for boundary conditions */




/* SOLVE T MATRIX BY ITERATIVE METHOD */
int gausmatt(int am, long int mcmax, long int mcmin)
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val1[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val1[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* mcmax - current line in num0[] */
/* mcmin - Number of iteration cycles */
{
/* Counters */
long int m1,m2,m3,m4;
/* Val Buffer */
double ival,ival1,ival2,ivalsum;
/**/
/**/
/**/
/* Seudel iteration */
/* Major line cycle */
for (m3=0;m3<am;m3++)
{
/**/
#pragma omp parallel for shared(fre1,num0,pos0,sol0,lin0,val1) private(m1,m2,ival,ival1) firstprivate(m3,mcmin,mcmax) schedule(static)
for (m1=mcmin;m1<=mcmax;m1++)
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
	sol0[m1]-=ival/ival1*t0koef;
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
	for (m2=0;m2<num0[m1];m2++) 
		{
		/* Recalc increment */
		ival-=val1[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]];
		}
	fre0[m1]=ival;
	}
return 0;
}
/* End SOLVE T MATRIX BY ITERATIVE METHOD */





