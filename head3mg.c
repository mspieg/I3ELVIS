/* MACROS --------------------------------------------- */
/* Min, Max, Abs */
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
/**/
/* Main Sizes */
#define MAXMAT  2600000000 /* General Matrix Size */
#define MAXXLN    521 /* Max Node num for Grid in X-direction */
#define MAXYLN    281 /* Max Node num for Grid in Y-direction */
#define MAXZLN    521 /* Max Node num for Grid in Z-direction */
#define MAXTMR    200 /* Max markers types */
#define MAXPOS    200 /* Max Pos on Line num for wi[], wn[] buffers */
#define MAXFLN   5000 /* Max Output file names Num */
#define MAXNOD   80000000  /* Max Total nodes num for all grids */
#define MAXNNN   40000000  /* Max Total nodes num for finest grid */
#define MAXPAR  400000000 /* Max Total par num for all grids */
#define MAXPPP  320000000 /* Max vx,vy,vz,p par num for all grids */
#define MAXMRK  330000000 /* Max total marker Num */
#define MAXBON   10000000  /* Max Total bondary condition Equations */
/* End MACROS --------------------------------------------- */

/* ARRAYS --------------------------------------------- */
/* Processing+Service Arrays for solutions */
/* val0[] - matrix contents */
/* fre0[] - free member for lines */
/* bufv[] - buffer for matrix organisation */
/* lin0[] - line numbers for matrix contents */
/* num0[] - pos numbers for line in  val0[], pos0[] */
/* pos0[] - first pos numbers for line in  val0[], pos0[] */
/* sol0[],sol1[] - solution buffers */
/* un[],ui[] - Line Num,Koef buffer */
double sol0[MAXPPP],sol1[MAXPPP];
double val1[MAXMAT],fre1[MAXPPP],fre0[MAXPPP],bufv[MAXPPP];
int lin0[MAXMAT],num0[MAXPPP];
long int pos0[MAXPPP];
int bufn[MAXPPP]; 
long int un[MAXPOS];
double ui[MAXPOS];
long int leftnum,rightnum;
int mgn[100],mgp[100],mgs[100],mgx[100],mgy[100],mgz[100],mgxy[100];
double mggx[10][MAXXLN],mggy[10][MAXYLN],mggz[10][MAXZLN];
/**/
/* Processing+Service Arrays for solutions */
/* mat[], sol[] - intermediat LSQ arrays */
/* wn[],wi[] - Line Num,Koef buffer */
long int wn[MAXPOS];
double wi[MAXPOS];
double mat[100],sol[10];
/**/
/**/
/**/
/* Nodes information */
/* nu[], dh[], ss[] - Reological Eq Koef */
/* ro[] - densyty, kg/m^3 */
/* tk[] - themperature K */
/* cp[] - Heat capacity, J/kg */
/* kt[] - Thermal conductivity koef, Wt/m/K */
/* et[] - Thermal expansivity, 1/K */
/* ht[] - Heat Sources, Wt/kg */
/* vx[],vy[],vz[] - Vx,Vy,Vz in cur cycle,  m/sek */
/* vx0[],vy0[],vz0[] - Vx,Vy,Vz in last cycle, m/sek */
/* gx[], gy[], gz[] - coordinates of grid surfaces */
/* exx[],eyy[],ezz[],exy[],exz[],eyz[] - deviatoric strain rates */
/* sxx[],syy[],szz[],sxy[],sxz[],syz[] - deviatoric stresses */
/* bondm[] - bondary position in bondn[],bondv[] for cur par (0) - no bondary */
/* bondn[] - PAR1+1 num in bond equat: CURPAR=CONST+KOEF1*PAR1 */
/* bondv[] - CONST,KOEF1 val in bond equat: CURPAR=CONST+KOEF1*PAR1 */
double gx[MAXXLN],gy[MAXYLN],gz[MAXZLN];
float nu[MAXNOD],ro[MAXNOD],nd[MAXNOD];
float nxx[MAXNOD],nxy[MAXNOD],nxz[MAXNOD],nyz[MAXNOD];
float cp[MAXNOD],kt[MAXNOD];
float ht[MAXNOD],et[MAXNOD],rr[MAXNOD];
double pr[MAXNOD];
float tk[MAXNOD], tk0[MAXNNN],tk1[MAXNNN],tk2[MAXNNN];
float vx[MAXNOD],vy[MAXNOD],vz[MAXNOD];
float exx[MAXNNN],eyy[MAXNNN],ezz[MAXNNN],exy[MAXNNN],exz[MAXNNN],eyz[MAXNNN];
float sxx[MAXNNN],syy[MAXNNN],szz[MAXNNN],sxy[MAXNNN],sxz[MAXNNN],syz[MAXNNN];
float bondv1[MAXBON][2];
int bondn1[MAXBON];
int bondm[MAXPAR];
float td[351][351][15][3];
/**/
/**/
/**/
/* Markers and Rock information */
/* Markers information */
/* markx[], marky[], markz[] - X,Y,Z of markers */
/* markk[] - temperature in K for marker */
/* mareii[] - second strain rate invariant for marker, 1/s */
/* markt[] - rock type of markers */
/* Information for different rock types */
/* marknu[], markdh[], markdv[], markss[] markmm[] -  Koef in ductile rheology Eq */
/* markll[], marka0[], marka1[], markb0[] markb1[], marke0[], marke1[] - Koef in brittle rheology Eq */
/* markn0[], markn1[], marks0[], marks1[] - viscosity and stress limits for individual rock types */
/* markro[], markaa[], markbb[], markcp[], markkt[], markkf[], markkv[], markht[] - ro,aro,bro,cp,kt, Tkt, Pkt, ht */
float markx[MAXMRK],marky[MAXMRK],markz[MAXMRK];
float markk[MAXMRK],markv[MAXMRK];
float markd[MAXMRK],markw[MAXMRK],markto[MAXMRK];
float marke[MAXMRK],markex[MAXMRK],markft[MAXMRK],markdt[MAXMRK],markrr[MAXMRK],markrn[MAXMRK];
char markt[MAXMRK];
double marknu[MAXTMR],markdh[MAXTMR],markdv[MAXTMR],markss[MAXTMR],markmm[MAXTMR];
double markll[MAXTMR],marka0[MAXTMR],marka1[MAXTMR],markb0[MAXTMR],markb1[MAXTMR],marke0[MAXTMR],marke1[MAXTMR];
double markro[MAXTMR],markaa[MAXTMR],markbb[MAXTMR],markcp[MAXTMR],markkt[MAXTMR],markht[MAXTMR],markkf[MAXTMR],markkp[MAXTMR];
double markn0[MAXTMR],markn1[MAXTMR],marks0[MAXTMR],marks1[MAXTMR];
/**/
/**/
/**/
/* Service, Connection Buffers */
/* vxyz[] - Cur Vx,Vy,Vz */
/* eps[] - Cur Eps tenzors */
/* errbuf[] - Error buffer */
/* nunu[] - NUik buffer */
/* sa[] - input string buffer */
/* fl1in[] - input data file name */
/* fl1itp - input data file type */
/* fl1otp - output data file type */
/* fl0out[] - output data All file names */
/* fl0otp[] - output data files types */
/* fl0stp[],fl0cyc[] - output data file X,T,time steps cyc0max for Out files */
/* fl1out[] - output data Cur file name */
/* *fl - Streem for File Input/Output */
double vxyz[20],errbuf[20],nunu[6][4],eps[100];
char sa[250];
char fl1in[60],fl1out[60];
int fl1itp,fl1otp;
char fl0out[MAXFLN][60];
int fl0otp[MAXFLN];
double fl0stp[MAXFLN][20];
int fl0cyc[MAXFLN][5];
FILE *fl,*fl1;
/**/
/* End ARRAYS --------------------------------------------- */


/* VARIABLES -------------------------------------------------- */
/**/
/* Grid Parameters */
/* xnumx,ynumy,znumz, - num of nodes in grid for X,Y,Z directions */
/* xnumx1,ynumy1,znumz1, - num of cells in grid for X,Y,Z directions */
/* xynumxy,xynumxy1 - num of nodes,cells in one xy slide */
/* mnumx,mnumy,mnumz - num of markers in one cell for X,Y,Z directions */
/* xsize,ysize,zsize - size of grid in X,Y,Z directions, m */
/* GXKOEF,GYKOEF,GZKOEF - Gravitation g for X,Y,Z directions, m/sek^2 */
/* pinit - Pressure in standard pressure Cell, Pa */
/* pxinit,pyinit,pzinit - address of the standard Pressure Cell */
/* xstpx,ystpy,zstpz - X,Y,Z steps for grid */
/* kfx,kfy,kfz,kfxx,kfxy,kfyy,kfxz,kfyz,kfzz - Koef for numeric differenciation */
/* cellnum - total Cell Num */
/* nodenum - total Node Num */
/* marknum - total Marker Num */
/* rocknum - rock types num */
/* bondnum - bondary condition equation num */
/* nodenum3,nodenum4 etc - nodenum*3,nodenum*4 */
/* spheryn - spherical gravity Y(1)/N(0) */
long int xnumx,ynumy,znumz;
long int xnumx1,ynumy1,znumz1;
long int xynumxy,xynumxy1;
long int mnumx,mnumy,mnumz;
double xsize,ysize,zsize;
double GXKOEF,GYKOEF,GZKOEF;
double pinit;
long int pxinit,pyinit,pzinit;
double xstpx,ystpy,zstpz;
double kfx,kfy,kfz,kfxx,kfxy,kfyy,kfxz,kfyz,kfzz;
long int cellnum,nodenum;
long int marknum,bondnum;
int rocknum,spheryn;
long int nodenum2,nodenum3,nodenum4,nodenum8;
clock_t clockbeg, clockend;
/**/
/**/
/**/
/* Service parameters */
/* loadmod - load from data file(1) or set initial conditions(0) */
/* printmod - print information on the monitor Y(1)/N(0) */
/* fl0num - number of otput file Names */
/* pos0cur - Cur pos Counter in val0[] lin0[] */
/* movemod - solve Stoks+Continuity equations Y(1)/N(0) */
/* tempmod - solve Heat Transfer equation Y(1)/N(0) */
/* markmod - move markers Y(1-simple,2-Runge-Kutta4)/N(0) */
/* gridmod - recalc grid parameters Y(1)/N(0) */
/* densimod - mode of  density calculation: 0-constant, 1-PT-dependent, 2-TDbase 3-PT-dependent+WaterTDbase */
int cyc0max,movemod,tempmod,markmod,gridmod;
long int pos0cur=0;
int printmod,fl0num,markmod,movemod,densimod,filesjob;
double gridcur,gridtot,stp100;
/**/
/**/
/**/
/* Errosion/Sedimentation parameters */
/* erosmod - errosion/sedimentation Y(1)/N(0)  */
/* eroslev - errosion level, m */
/* eroscon - constant erosion rate, m/s */
/* eroskoe - increment of erosion rate with elevation, 1/s */
/* sedilev - sedimentation level, m */
/* sedicon - constant sedimentation rate, m/s */
/* sedikoe - increment of sedimentation rate with burial, 1/s */
/* sedimnum - Num of cycles of sedimentation */
/* sedimcyc - Num of cycles of sedimentation for One Layer */
/* waterlev - water/air boundary, m */
/* slopemax - Max slope y/x ratio */
/* basalty - basalt melting depth, m */
/* dehydrmin, dehydrmax - serpentine dehydration min, max depth, m */
int erosmod;
double eroslev,eroscon,eroskoe;
double sedilev,sedicon,sedikoe;
double waterlev,slopemax;
int sedimnum=0,sedimcyc=3;
/**/
/**/
/**/
/* Motion parameters */
/* cyc0max - Num of circles for cur calculation */
/* maxxystep - max Distance change for one time step, m */
/* maxtkstep - max Themperature for one time step, K */
/* maxtmstep - max Time for one time step, sek */
/* timestep - time step, sek */
/* timesum - time from start, sek */
/* outgrid - marker move out of grid Y(0)/N(1) Orthogonal Only (2) */
int cyc0max,outgrid=0;
double maxxystep,maxtkstep,maxtmstep,timestep,timesum,timedir;
/**/
/**/
/**/
/* Termodynamic database parameters <loadjxxx.c> */
double tkmin,pbmin,tkstp,pbstp;
int tknum,pbnum;
double tkpor,zmpor,vyfluid,vymelt,dmwamin,tdeep,dtdeep,drdeep,vdeep,zdeep,nudeep,dxwater,dywater,dzwater,maxpmelt,maxwater;
double dserp,yserp,tserp,lambmlt,ydike,dydike,nusseltno,circdepth,circtemp,healrate,meltmin,xvolc,xtkplut,teclmin,teclmax;
/**/
/**/
/**/
/* General iteration parameters in <gaus.c> */
/* ckoef - Cur Koef for Zeydel method of iteration */
double ckoef;
/**/
/**/
/**/
/* V Iteration parameters in <move.c> */
/* cyc1max - V Max Cycle num for iteration cycle */
/* oldbeg,oldend,oldstp - Init,Final,Step wt for V old solution */
/* rimbeg,rimend,rimstp - Init,Final,Step wt for P recalc after contin Eq */
/* DIVVMIN,STOKSMIN - Min Valid absolut Err val for Contin,Stokes Eq */
/* EPSERR - Min Valid absolut Err val for EPSik, in Reological Eq */
/* EPSMIN - Min Valid relativ Err val for EPSik, in Reological Eq */
/* nukoef - Avrige Nu for Pressure optimization */
/* dprmax - max Value for One time pressure change */
/* dvvmax - max Value for One time Vx,Vy,Vz change */
/* cv0koef,cv1koef,cv2koef - 1-st koef,2-nd koef,ERR svitch for Zeydel method, Cur,V,T Koef for Zeydel method of iteration */
/* dvbeg,dvend - Min-Max limits for Vx,Vy difference in all grid */
/* nubeg,nuend -  Min-Max limits of viscozity */
int cyc1max,multinum,multimax,multicyc,multinnn[5][50];
double multinnv[4][50];
double oldbeg,oldend,oldstp;
double rimbeg,rimend,rimstp;
double DIVVMIN,STOKSMIN,DIVVMAX,STOKSMAX,EPSERR,EPSMIN;
double p0koef,p1koef,p2koef,v0koef,v1koef,vxkoef,vykoef,vzkoef,nukoef,numult;
int maxcyc,multstp;
double cv0koef,cv1koef,cv2koef;
double dvbeg,dvend;
double nubeg,nuend;
double prbeg,prend,prlow,prhigh;
int viscmod;
double stoksum=100.0,stoksum1=100.0,stoksum2=100.0,stoknum=100.0,stoksum0=-1e+30,contsum0=-1e+30;
/**/
/**/
/**/
/* T Iteration parameters in <heat.c> */
/* cyc2max - T Max Cycle num for iteration cycle */
/* oldbegt,oldendt,oldstpt - Init,Final,Step wt for T old solution */
/* HEATMIN - Min Valid absolut Err val for Heat Equation */
/* ctkoef - Cur,V,T Koef for Zeydel method of iteration */
/* dtbeg,dtend - Min-Max limits for Max Memb in T Eq */
/* frictyn - Viscouse friction heat Y(1)/N(0) */
/* adiabyn - adiabatic heat calculation: N(0)/Y(1) */
/* heatdif - numerical heat diffusion koef */
int cyc2max,multinumt,multicyct,multittt[2][50];
double oldbegt,oldendt,oldstpt;
double HEATMIN;
double ctkoef;
double dtbeg,dtend;
double heatdif;
int frictyn,adiabyn;
double t0koef,t1koef;
/* Thermal BC change */
double convend,timebeg,timeend,tempbeg,tempend;
/**/
/**/
/* End VARIABLES -------------------------------------------------- */



/* FUNCTONS PROTOTYPES ------------------------------------------- */
/**/
/* <load.c> FILE */
/* ffscanf() - Load single word from file without empty lines */
/* ffscanf1() - Load single word from fl1 without empty lines */
/* loadconf() - Load configuration from mode.t3c */
/* loader() - Load information from data file */
/* saver() - Save information to data file */
/* gridcheck() - Calc,Check parameters of Grid */
void ffscanf();
void ffscanf1();
int loadconf();
void loader();
void saver(int,int);
void gridcheck();
/**/
/* <init.c> FILE */
/* setinit() - Organise initial condition for processing */
void setinit();
/**/
/* <gaus.c> FILE */
/* gausmat() - Solve interm LSQ matrix by Gauss method */
/* matclear() - Clear interm LSQ matrix */
/* gausmat2() - Add system of equations */
int gausmat2(int, long int, long int, long int *, double *);
void gausmat(int, double *, double *);
void matclear(int, double *);
/**/
/* <move.c> FILE */
/* vpiterate() - General soubr for Vx,Vy,P Calc by Iterativ method */
/* xstokserr() -  Right part or  Err in Stokes Equat for cur node calc */
/* continerr() - New Pressure and Err in continuity Equat for cur Cell calc */
/* sxxcalc() etc. - Value or add EPS and SIG equations */
/* ronurestrict() - Finer->Coarcer level recalc for Ro, Nu */
/* vpprolong() - Vx,Vy,Vz,P "Prolongation" from coarser to finer level */
/* vprestrict() - Vx,Vy,Vz,P residuals "Restriction" from finer to coarser level */
void vpiterate(int, int);
double xstokserr(long int, long int, long int, int, int, double *, long int *, double *);
double ystokserr(long int, long int, long int, int, int, double *, long int *, double *);
double zstokserr(long int, long int, long int, int, int, double *, long int *, double *);
double continerr(long int, long int, long int, int, int, double *, long int *, double *);
double continadd(long int, long int, long int, double, int, double *, long int *, double *);
double sxxcalc(long int, long int, long int, double, int, double *, long int *, double *);
double syycalc(long int, long int, long int, double, int, double *, long int *, double *);
double szzcalc(long int, long int, long int, double, int, double *, long int *, double *);
double sxycalc(long int, long int, long int, double, int, double *, long int *, double *);
double sxzcalc(long int, long int, long int, double, int, double *, long int *, double *);
double syzcalc(long int, long int, long int, double, int, double *, long int *, double *);
void ronurestrict(int);
void vpprolong(int, int);
void vprestrict(int);
int gausmat4(int, long int, long int, int);
/**/
/* <heat.c> FILE */
/* titerate() -  Themperature recalc after time step */
/* tkrecalc() -  Themperature recalc after old + New solution */
/* heaterr() -  Left+Right part or  Err in Heat Equat for cur node calc */
/* tbonderr() -  Left+Right part or  Err in Temperature boundary Equat for cur node calc */
void titerate(int);
void tkrecalc();
double heaterr(long int, long int, long int, int, int, double *, long int *, double *);
double tbonderr(long int, int, double *, long int *, double *);
int gausmatt(int, long int, long int);
/**/
/* <mark.c> FILE */
/* movemark() - move markers by Runge-Kutta method */
/* vxyzcalc() - Vx,Vy,Vz calc for marker by interpolation */
/* ronurecalc() - recalc ro[],nu[] etc after new marker position */
/* tkcalc() - T K calc for marker by interpolation */
/* prcalc() - P Pa calc for marker by interpolation */
/* dencalc() - Ro kg/m^3 calc for marker */
/* viscalc() - Nu Pa.s calc for marker */
/* epscalc() - EPSij 1/s calc for marker by interpolation */
/* m1serch(), m2serch(), m3serch() - serch of nearest upper-left-frontal node of cell for  current location */
/* antigor() - Mantle transformation to Antigorite */
/* hydration2() - Hydration front progress recalc */
/* melting(), meltpart(), meltpart1() - melting of rocks account */
/* gridchange() - change of the grid */
/* meltextract() - melt extraction, crustal growth */
/* tdbasecalc() - compute TD properties at given P,T,Composition */
void movemark();
void ronurecalc();
void tkcalc1(double, double, double, double *, long int *, double *);
double viscalc(double, double, long int, int, double, double, double, double, long int *);
long int m1serch(double);
long int m2serch(double);
long int m3serch(double);
void meltextract();
void antigor(double, double, long int);
void vxyzcalc1(double, double, double, double *, long int *);
void epscalc1(double, double, double, int, double *, long int *);
long int m1serch1(double, long int *);
long int m2serch1(double, long int *);
long int m3serch1(double, long int *);
void prcalc1(double, double, double, double *, long int *);
double meltpart0(double, double, long int, int, double *);
void meltpart11(double, double, long int, double *);
void melting1(double, double, long int, double *);
double dencalc1(double, double, int, double *);
void tdbasecalc1(double, double, int, long int, double *);
double hydration2();
/**/
/* End FUNCTONS PROTOTYPES ------------------------------------------- */

