/* CLEAR mat(rm) */
void matclear(int rm, double *adr)
/* rm - rank of matrix */
/* adr - adres of first element of matrix */
{
/* String lenth */
/* int rm1=rm+1; */
/* Counters */
int m1,m2;
/**/
for (m1=0;m1<rm;m1++)
for (m2=0;m2<=rm;m2++)
	{
	*adr=0;
	adr++;
	}
}
/* End CLEAR mat0(rm) */



/* SOLVE mat(rm) BY USING GAUSS METHOD */
void gausmat(int rm, double *adr, double *adr1)
/* rm - rank of matrix */
/* adr - adres of first element of matrix */
/* adr1 - adres of first element of solutions */
{
/* Counters */
int m1,m2,m3;
/* Buffer for val */
double val1, val2;
/* String lenth */
int rm1=rm+1;
/**/
/* STEP 2: KOEF OPERATIONS */
for (m1=0;m1<rm-1;m1++)
if ((val1=*(adr+m1*rm1+m1)))
	for (m2=m1+1;m2<rm;m2++)
	if ((val2=*(adr+m2*rm1+m1)))
		for (m3=m1+1;m3<=rm;m3++)
		/* F[]  G[] */
		*(adr+m2*rm1+m3)=(*(adr+m2*rm1+m3))/val2-(*(adr+m1*rm1+m3))/val1;
/* End STEP 2: KOEF OPERATIONS */
/**/
/* STEP 3: SOLUTION CALC CHECK */
for (m1=rm-1;m1>=0;m1--)
	{
	/* Clear Sol */
	*(adr1+m1)=0;
	if ((val1=*(adr+m1*rm1+m1))!=0)
		{
		/* Calc Sol */
		val2=*(adr+m1*rm1+rm)/val1;
		*(adr1+m1)=val2;
		/* Use Sol */
		for (m2=0;m2<=m1-1;m2++)
		if ((val1=*(adr+m2*rm1+m1))!=0)
		*(adr+m2*rm1+rm)-=val1*val2;
		}
	}
/* End STEP 3: SOLUTION CALC CHECK */
}
/* End SOLVE m0() BY USING GAUSS METHOD */


/* ADD MATRIX */
int gausmat2(int am, long int mcmax, long int pos0cur1, long int *un1, double *ui1)
/* un[] - line koef numbers */
/* ui[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val1[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* pos0cur1 - position number in val0[], lin0[] */
/* mcmax - current line in num0[] */
/* mcmin - first line in num0[] */
{
/* Counters */
long int m1,m2,m3;
/* Val Buffer */
double ival,ival1;
/**/
/* Space Check */
if (mcmax>=MAXPAR)
	{
	printf("EXIT PROGRAM: Space out in fre1[] %ld",mcmax);
	exit(0);
	}
/**/
/**/
/**/
/* STEP 0: RELOAD KOEF FROM ui[] to val1[] */
/*
printf("GAUS %ld   %e ",mcmax,fre1[mcmax]); getchar();
*/
	/**/
	/* Free member reload from buffer */
	fre1[mcmax]=ui1[0];
/*
if(un[0]>1) {printf(" %ld %ld   %ld %e \n",pos0cur,mcmax,un[0],ui[0]); getchar();}
*/
	/**/
	/* Line koef reload to val1[] from buffer */
	pos0[mcmax]=pos0cur1;
	num0[mcmax]=0;
	/* Reload koefficients */
	for (m2=1;m2<=un1[0];m2++)
		{
		/* Save Cur Koef */
		lin0[pos0cur1]=un1[m2];
		val1[pos0cur1]=ui1[m2];
		pos0cur1++;
		num0[mcmax]++;
		}
/* Check Space */
if (pos0cur1>=MAXMAT) 
	{
	printf("EXIT PROGRAM: Space out in val1[] lin0[] %ld %ld",mcmax,pos0cur1);
	exit(0);
	}
/*
if (mcmax>4000) getchar();
*/
	return 0;
}
/* End STEP 0: RELOAD KOEF FROM ui[] to val1[] */


