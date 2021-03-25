/* Load information from configuration file mode.t3c ============== */
int loadconf()
{
/* Counter */
int m1,m2,n1,n2,fln3;
/**/
/**/
/**/
/* Open File file.t3c */
fl = fopen("file.t3c","rt");
ffscanf(); fln3=atoi(sa)-1;
fclose(fl);
/**/
/**/
/**/
/* Open File mode.t3c */
fl = fopen("mode.t3c","rt");
/**/
/* Data File name */
ffscanf();
for (m1=0;m1<50;m1++) fl1in[m1]=sa[m1];
ffscanf(); if(sa[0] == 'b') fl1itp=1;
/**/
/* Load first Results File names */
ffscanf();
fl0num=0;
while(sa[0]!='~')
	{
	/* Check file Counter */
	if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
	/**/
	/* Save results file name */
	for (m1=0;m1<50;m1++) fl0out[fl0num][m1]=sa[m1];
	/**/
	/* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep */
	ffscanf(); if(sa[0] == 'b') fl0otp[fl0num]=1;
	ffscanf();fl0cyc[fl0num][0]=atoi(sa);
	ffscanf();fl0stp[fl0num][0]=atof(sa);
	ffscanf();fl0stp[fl0num][1]=atof(sa);
	ffscanf();fl0stp[fl0num][2]=atof(sa)*3.15576e+7;
	ffscanf();fl0stp[fl0num][3]=atof(sa);
	ffscanf();fl0stp[fl0num][4]=atof(sa);
	ffscanf();fl0stp[fl0num][5]=atof(sa);
	ffscanf();fl0stp[fl0num][6]=atof(sa);
	ffscanf();fl0stp[fl0num][7]=atof(sa);
	ffscanf();fl0stp[fl0num][8]=atof(sa);
	ffscanf();fl0stp[fl0num][9]=atof(sa);
	ffscanf();fl0stp[fl0num][10]=atof(sa);
	ffscanf();fl0cyc[fl0num][1]=atoi(sa);
	ffscanf();fl0cyc[fl0num][2]=atoi(sa);
	ffscanf();fl0cyc[fl0num][3]=atoi(sa);
	/**/
	/* Incr File Counters */
	fl0num++;
	/**/
	/* Load Next Results File names */
	ffscanf();
	}
/**/
/* Service */
ffscanf();printmod=atoi(sa);
ffscanf();movemod=atoi(sa);
ffscanf();tempmod=atoi(sa);
ffscanf();markmod=atoi(sa);
ffscanf();gridmod=atoi(sa);
ffscanf();outgrid=atoi(sa);
ffscanf();densimod=atoi(sa);
ffscanf();stp100=atof(sa);
ffscanf();filesjob=atoi(sa);
if(stp100<=0) {printf("Invalid stp100 %e",stp100);exit(0);}
printf("Current stp100 %e",stp100);
/**/
/* Errosion/Sedimentation */
ffscanf();eroslev=atof(sa);
ffscanf();sedilev=atof(sa);
ffscanf();waterlev=atof(sa);
/**/
/* V Iteration */
ffscanf();cyc1max=atoi(sa);
ffscanf();DIVVMIN=atof(sa);
ffscanf();STOKSMIN=atof(sa);
ffscanf();DIVVMAX=atof(sa);
ffscanf();STOKSMAX=atof(sa);
ffscanf();multimax=multinum=atoi(sa);
ffscanf();multicyc=atoi(sa);
for(m1=0;m1<=1;m1++) 
for(m2=0;m2<=multinum;m2++) 
	{
	ffscanf();multinnn[m1][m2]=atoi(sa);
	}
ffscanf();p0koef=atof(sa);
ffscanf();p1koef=atof(sa);
ffscanf();p2koef=atof(sa);
ffscanf();v0koef=atof(sa);
ffscanf();v1koef=atof(sa);
ffscanf();nubeg=atof(sa);
ffscanf();nuend=atof(sa);
ffscanf();nukoef=atof(sa);
ffscanf();viscmod=atoi(sa);
ffscanf();spheryn=atoi(sa);
/**/
/* T Iteration */
ffscanf();cyc2max=atoi(sa);
ffscanf();HEATMIN=atof(sa);
ffscanf();multinumt=atoi(sa);
ffscanf();multicyct=atoi(sa);
for(m1=0;m1<=1;m1++) 
for(m2=0;m2<=multinumt;m2++) 
	{
	ffscanf();multittt[m1][m2]=atoi(sa);
	}
ffscanf();t0koef=atof(sa);
ffscanf();t1koef=atof(sa);
ffscanf();heatdif=atof(sa);
ffscanf();frictyn=atof(sa);
ffscanf();adiabyn=atof(sa);
/**/
/* Water */
ffscanf();tkpor=atof(sa);
ffscanf();zmpor=atof(sa);
ffscanf();vyfluid=atof(sa);
ffscanf();vymelt=atof(sa);
ffscanf();dmwamin=atof(sa);
ffscanf();tdeep=atof(sa);
ffscanf();dtdeep=atof(sa);
ffscanf();drdeep=atof(sa);
ffscanf();zdeep=atof(sa);
ffscanf();vdeep=atof(sa);
ffscanf();nudeep=atof(sa);
ffscanf();dxwater=atof(sa);
ffscanf();dywater=atof(sa);
ffscanf();dzwater=atof(sa);
ffscanf();dserp=atof(sa);
ffscanf();yserp=atof(sa);
ffscanf();tserp=atof(sa);
ffscanf();lambmlt=atof(sa);
ffscanf();ydike=atof(sa);
ffscanf();dydike=atof(sa);
ffscanf();nusseltno=atof(sa);
ffscanf();circdepth=atof(sa);
ffscanf();circtemp=atof(sa);
ffscanf();healrate=atof(sa);
ffscanf();meltmin=atof(sa);
ffscanf();maxwater=atof(sa);
ffscanf();maxpmelt=atof(sa);
ffscanf();teclmin=atof(sa);
ffscanf();teclmax=atof(sa);
/* Volcanic ratio, source temperature influence for plutonics */
ffscanf();xvolc=atof(sa);
ffscanf();xtkplut=atof(sa);
/* Thermal BC for the plume */
ffscanf();timebeg=atof(sa);timebeg*=3.15576e+7;
ffscanf();timeend=atof(sa);timeend*=3.15576e+7;
ffscanf();tempbeg=atof(sa);
/* Convergence condition end */
ffscanf();convend=atof(sa);convend*=3.15576e+7;
/*
printf("%e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin);getchar();
*/
/* Data File name change after number */
if(fln3>=0 && fln3<fl0num)
	{
	for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[fln3][n1];
	fl1itp=fl0otp[fln3];
	}
else
	{
	fln3=-1;
	}
/**/
fclose(fl);
/* End Load information from configuration file mode.t3c */
/**/
/* stop.yn file creation */
fl = fopen("stop.yn","wt");
fprintf(fl,"n \n");
fclose(fl);
/**/
return fln3;
/*
*/
/**/
/* Load thermodynamic database */
/* Dry peridotite */
/* RO - density */
fl = fopen("pdry_rho","rt");
ffscanf();
ffscanf(); tknum=atoi(sa);
ffscanf(); pbnum=atoi(sa);
ffscanf(); tkmin=atof(sa);
ffscanf(); pbmin=atof(sa);
ffscanf(); tkstp=atof(sa);
ffscanf(); pbstp=atof(sa);
ffscanf();
ffscanf();
/*
printf("%d %d %e %e %e %e",tknum,pbnum,tkmin,pbmin,tkstp,pbstp);getchar();
*/
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][0]=atof(sa)/1000.0;
/*
printf("%d %d %e",n1,n2,td[n1][n2][0][0]);getchar();
*/
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][0][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][0][1]-td[n2-1][n1][0][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PDRY %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
if (printmod) printf("pdry_hh OK \n");
/**/
/* Wet peridotite */
/* RO - density */
fl = fopen("pwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][1][1]-td[n2-1][n1][1][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("pwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][1][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("pwet_h2o OK \n");
/**/
/**/
/* Molten peridotite */
/* RO - density */
fl = fopen("pwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("pwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][2][1]-td[n2-1][n1][2][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("PMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("pwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][2][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("pwetmelt_h2o OK \n");
/**/
/**/
/* Wet Gabbro */
/* RO - density */
fl = fopen("gabwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("gabwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][3][1]-td[n2-1][n1][3][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("GWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("gabwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][3][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("gabwet_h2o OK \n");
/**/
/**/
/* Molten Gabbro */
/* RO - density */
fl = fopen("gabwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("gabwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][4][1]-td[n2-1][n1][4][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("GMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("gabwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][4][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("gabwetmelt_h2o OK \n");
/**/
/**/
/* Wet sediments */
/* RO - density */
fl = fopen("swet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("swet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][5][1]-td[n2-1][n1][5][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("SWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("swet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][5][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("swet_h2o OK \n");
/**/
/**/
/* Molten sediments */
/* RO - density */
fl = fopen("swetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("swetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][6][1]-td[n2-1][n1][6][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("SMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("swetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][6][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("swetmelt_h2o OK \n");
/**/
/**/
/* Wet Basalt */
/* RO - density */
fl = fopen("bwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("bwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][7][1]-td[n2-1][n1][7][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("bwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][7][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("bwet_h2o OK \n");
/**/
/**/
/* Molten Basalt */
/* RO - density */
fl = fopen("bwetmelt_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("bwetmelt_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][8][1]-td[n2-1][n1][8][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("bwetmelt_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
	}
fclose(fl);
if (printmod) printf("bwetmelt_h2o OK \n");
/**/
/**/
/* Dry Upper crust */
/* RO - density */
fl = fopen("ucdry_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][11][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("ucdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][11][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][11][1]-td[n2-1][n1][11][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
printf("ucdry_hh OK \n");
/**/
/* Wet Upper crust */
/* RO - density */
fl = fopen("ucwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("ucwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][12][1]-td[n2-1][n1][12][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("ucwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][12][2]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("ucwet_h2o OK \n");
/**/
/**/
/* Dry Lower crust */
/* RO - density */
fl = fopen("lcdry_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("lcdry_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][13][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][13][1]-td[n2-1][n1][13][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
printf("lcdry_hh OK \n");
/**/
/* Wet Lower crust */
/* RO - density */
fl = fopen("lcwet_rho","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][0]=atof(sa)/1000.0;
	}
fclose(fl);
/**/
/* H  - enthalpy */
fl = fopen("lcwet_hh","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][1]=atof(sa)/1000.0/4.1837;
/*
if (n2>0) 
{
ival=1000.0*4.1837*(td[n2][n1][14][1]-td[n2-1][n1][14][1])/tkstp;
if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
}
*/
	}
fclose(fl);
/**/
/* Wa - water contents, wt% */
fl = fopen("lcwet_h2o","rt");
for (n1=0;n1<9;n1++) ffscanf();
for (n1=0;n1<pbnum;n1++)
for (n2=0;n2<tknum;n2++)
	{
	ffscanf(); td[n2][n1][14][2]=MAXV(0,atof(sa));
	}
fclose(fl);
printf("lcwet_h2o OK \n");
/**/
/**/
/*
printf("%d %d %e %e %e %e ",n1,n2,td[n1-1][n2-1][0][0],td[n1-1][n2-1][0][1],td[n1-1][n2-1][0][2],td[n1-1][n2-1][0][3]);getchar();
*/
/**/
/* Load global thermodynamic database */
/*
printf("%d %d %e %e %e %e ",n1,n2,td[n1-1][n2-1][0][0],td[n1-1][n2-1][0][1],td[n1-1][n2-1][0][2],td[n1-1][n2-1][0][3]);getchar();
*/
/**/
/* Return counter */
/* stop.yn file creation */
return fln3;
}
/* Load information from configuration file mode.t3c ============== */



/* Load Information from data file ------------------------------- */
void loader()
/* bondv[] - bondary value */
/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
{
/* Counter */
int n1;
char nn1;
long int m1,m2,m3,m4;
long int mm1;
char szint,szlong,szfloat,szdouble,szcur;
float ival0;
double ival1;
/**/
/**/
/**/
/* Load Past Results from data file-------------------------------- */
if (printmod) printf("Load Past results from %s ...",fl1in);
/**/
/**/
/**/
/* Load in Binary Format ---------------------------- */
	{
	fl = fopen(fl1in,"rb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	/* Check sizes of variables */
	fread(&szcur,1,1,fl);
	if (szcur!=szint) {printf("Current INT size <%d> is different from given in file <%d> \n",szint,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szlong) {printf("Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szfloat) {printf("Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); exit(0);}
	fread(&szcur,1,1,fl);
	if (szcur!=szdouble) {printf("Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); exit(0);}
	/**/
	/* Grid Parameters */
	fread(&xnumx,szlong,1,fl);
	fread(&ynumy,szlong,1,fl);
	fread(&znumz,szlong,1,fl);
	fread(&mnumx,szlong,1,fl);
	fread(&mnumy,szlong,1,fl);
	fread(&mnumz,szlong,1,fl);
	fread(&xsize,szdouble,1,fl);
	fread(&ysize,szdouble,1,fl);
	fread(&zsize,szdouble,1,fl);
	fread(&pxinit,szlong,1,fl);
	fread(&pyinit,szlong,1,fl);
	fread(&pzinit,szlong,1,fl);
	fread(&pinit,szdouble,1,fl);
	fread(&GXKOEF,szdouble,1,fl);
	fread(&GYKOEF,szdouble,1,fl);
	fread(&GZKOEF,szdouble,1,fl);
	fread(&rocknum,szint,1,fl);
	fread(&bondnum,szlong,1,fl);
	fread(&marknum,szlong,1,fl);
	fread(&n1,szint,1,fl);
	fread(&timesum,szdouble,1,fl);timesum*=3.15576e+7;
	fread(&gridcur,szdouble,1,fl);
	fread(&gridtot,szdouble,1,fl);
	/**/
	/* Calc,Check Grid parameters */
	gridcheck();
	/**/
	/* Rock Types information */
	fread(markn0,szdouble,rocknum,fl);
	fread(markn1,szdouble,rocknum,fl);
	fread(marks0,szdouble,rocknum,fl);
	fread(marks1,szdouble,rocknum,fl);
	fread(marknu,szdouble,rocknum,fl);
	fread(markdh,szdouble,rocknum,fl);
	fread(markdv,szdouble,rocknum,fl);
	fread(markss,szdouble,rocknum,fl);
	fread(markmm,szdouble,rocknum,fl);
	fread(markll,szdouble,rocknum,fl);
	fread(marka0,szdouble,rocknum,fl);
	fread(marka1,szdouble,rocknum,fl);
	fread(markb0,szdouble,rocknum,fl);
	fread(markb1,szdouble,rocknum,fl);
	fread(marke0,szdouble,rocknum,fl);
	fread(marke1,szdouble,rocknum,fl);
	fread(markro,szdouble,rocknum,fl);
	fread(markbb,szdouble,rocknum,fl);
	fread(markaa,szdouble,rocknum,fl);
	fread(markcp,szdouble,rocknum,fl);
	fread(markkt,szdouble,rocknum,fl);
	fread(markkf,szdouble,rocknum,fl);
	fread(markkp,szdouble,rocknum,fl);
	fread(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,ro[],nu[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		fread(&ival0,szfloat,1,fl);pr[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);vz[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ro[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nu[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nxx[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nxy[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nxz[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);nyz[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);tk[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);cp[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);et[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);kt[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);ht[m1]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);rr[m1]=(double)(ival0);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		fread(&ival0,szfloat,1,fl);gx[m1]=(double)(ival0);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		fread(&ival0,szfloat,1,fl);gy[m2]=(double)(ival0);
		}
	for (m3=0;m3<znumz;m3++)
		{
		fread(&ival0,szfloat,1,fl);gz[m3]=(double)(ival0);
		}
	/**/
	/* Bondary Conditions Equations */
	for (m1=0;m1<bondnum;m1++)
		{
		fread(&m2,szlong,1,fl);
		fread(&m3,szlong,1,fl);bondm[m2]=m3;
		/* Check boundary array */
		if(m3>MAXBON) {printf("Space out in bondv[]"); exit(0);}
		fread(&ival0,szfloat,1,fl);bondv1[m3][0]=(double)(ival0);
		fread(&ival0,szfloat,1,fl);bondv1[m3][1]=(double)(ival0);
		fread(&m2,szlong,1,fl);bondn1[m3]=m2;
/*
printf("Bond %ld %ld %ld %e %ld %e",bondnum,m1,m3,bondv1[m3][0],bondn1[m3],bondv1[m3][1]);getchar();
*/
		}
/*
printf("Bond %ld %ld %ld ",bondnum,m3,m2);getchar();
*/
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
		{
		/* General information load */
		fread(&ival0,szfloat,1,fl);markx[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);marky[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markz[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markk[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markw[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markd[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markv[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);marke[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markex[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markft[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markrr[mm1]=ival0;
		fread(&ival0,szfloat,1,fl);markrn[mm1]=ival0;
		fread(&nn1,1,1,fl);markt[mm1]=(int)nn1;
/*
printf("MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
*/
		}
	}
/* Load in Binary Format ---------------------------- */
/*
for (m3=0;m3<marknum;m3++) 
if (markt[m3]==100)
	{
{printf("%ld %d  %e %e %e ",m3,markt[m3], markx[m3],marky[m3],markz[m3]);getchar();}
	}
*/
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
/**/
/**/
/**/
/* Calc,Check Grid parameters */
gridcheck();
/**/
}
/* Load Information from data file ------------------------------- */



/* Print Results to data file ----------------------------------- */
void saver(int f0, int n0)
/* n0 - circle number */
/* f0 - file number */
{
/* Counters */
int n1;
char nn1;
long int m1,m2,m3,m4;
long int mm1;
/* Buffers for XY */
double x,y,z;
char szint,szlong,szfloat,szdouble,szcur;
float ival0;
double ival1;
/**/
/**/
/**/
if (printmod) printf("Print %d circle results to %s...",n0+1,fl1out);
/**/
/**/
/**/
/* Save data in binary format ---------------------------- */
	{
	fl = fopen(fl1out,"wb");
	/**/
	/* Sizes of var definition */
	szint=sizeof(n1);
	szlong=sizeof(m1);
	szfloat=sizeof(ival0);
	szdouble=sizeof(ival1);
	fwrite(&szint,1,1,fl);
	fwrite(&szlong,1,1,fl);
	fwrite(&szfloat,1,1,fl);
	fwrite(&szdouble,1,1,fl);
	/**/
	/* Grid Parameters */
	fwrite(&xnumx,szlong,1,fl);
	fwrite(&ynumy,szlong,1,fl);
	fwrite(&znumz,szlong,1,fl);
	fwrite(&mnumx,szlong,1,fl);
	fwrite(&mnumy,szlong,1,fl);
	fwrite(&mnumz,szlong,1,fl);
	fwrite(&xsize,szdouble,1,fl);
	fwrite(&ysize,szdouble,1,fl);
	fwrite(&zsize,szdouble,1,fl);
	fwrite(&pxinit,szlong,1,fl);
	fwrite(&pyinit,szlong,1,fl);
	fwrite(&pzinit,szlong,1,fl);
	fwrite(&pinit,szdouble,1,fl);
	fwrite(&GXKOEF,szdouble,1,fl);
	fwrite(&GYKOEF,szdouble,1,fl);
	fwrite(&GZKOEF,szdouble,1,fl);
	fwrite(&rocknum,szint,1,fl);
	fwrite(&bondnum,szlong,1,fl);
	fwrite(&marknum,szlong,1,fl);
	fwrite(&n0,szint,1,fl);
	ival1=timesum/3.15576e+7;fwrite(&ival1,szdouble,1,fl);
	fwrite(&gridcur,szdouble,1,fl);
	fwrite(&gridtot,szdouble,1,fl);
	/**/
	/* Rock Types information */
	fwrite(markn0,szdouble,rocknum,fl);
	fwrite(markn1,szdouble,rocknum,fl);
	fwrite(marks0,szdouble,rocknum,fl);
	fwrite(marks1,szdouble,rocknum,fl);
	fwrite(marknu,szdouble,rocknum,fl);
	fwrite(markdh,szdouble,rocknum,fl);
	fwrite(markdv,szdouble,rocknum,fl);
	fwrite(markss,szdouble,rocknum,fl);
	fwrite(markmm,szdouble,rocknum,fl);
	fwrite(markll,szdouble,rocknum,fl);
	fwrite(marka0,szdouble,rocknum,fl);
	fwrite(marka1,szdouble,rocknum,fl);
	fwrite(markb0,szdouble,rocknum,fl);
	fwrite(markb1,szdouble,rocknum,fl);
	fwrite(marke0,szdouble,rocknum,fl);
	fwrite(marke1,szdouble,rocknum,fl);
	fwrite(markro,szdouble,rocknum,fl);
	fwrite(markbb,szdouble,rocknum,fl);
	fwrite(markaa,szdouble,rocknum,fl);
	fwrite(markcp,szdouble,rocknum,fl);
	fwrite(markkt,szdouble,rocknum,fl);
	fwrite(markkf,szdouble,rocknum,fl);
	fwrite(markkp,szdouble,rocknum,fl);
	fwrite(markht,szdouble,rocknum,fl);
	/**/
	/* Nodes information */
	/* Vx,Vy,ro[],nu[],tk[],cp[],kt[],ht[] */
	for (m1=0;m1<nodenum;m1++)
		{
		ival0=(float)(pr[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(vz[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ro[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nu[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nxx[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nxy[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nxz[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(nyz[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(tk[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(cp[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(et[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(kt[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(ht[m1]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(rr[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Gridlines positions */
	for (m1=0;m1<xnumx;m1++)
		{
		ival0=(float)(gx[m1]);fwrite(&ival0,szfloat,1,fl);
		}
	for (m2=0;m2<ynumy;m2++)
		{
		ival0=(float)(gy[m2]);fwrite(&ival0,szfloat,1,fl);
		}
	for (m3=0;m3<znumz;m3++)
		{
		ival0=(float)(gz[m3]);fwrite(&ival0,szfloat,1,fl);
		}
	/**/
	/* Bondary Conditions Equations */
m3=0;
	for (m1=0;m1<nodenum*10;m1++)
	if(bondm[m1])
		{
		fwrite(&m1,szlong,1,fl);
		m2=bondm[m1];fwrite(&m2,szlong,1,fl);
		ival0=(float)(bondv1[m2][0]);fwrite(&ival0,szfloat,1,fl);
		ival0=(float)(bondv1[m2][1]);fwrite(&ival0,szfloat,1,fl);
		m2=bondn1[m2];fwrite(&m2,szlong,1,fl);
m3++;
		}
printf("BOND %ld %ld ",bondnum,m3);
	/**/
	/* Markers X,Y,types */
	for (mm1=0;mm1<marknum;mm1++)
		{
		/* General information Save */
		ival0=markx[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=marky[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markz[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markk[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markw[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markd[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markv[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=marke[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markex[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markft[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markrr[mm1];fwrite(&ival0,szfloat,1,fl);
		ival0=markrn[mm1];fwrite(&ival0,szfloat,1,fl);
		nn1=(char)markt[mm1];fwrite(&nn1,1,1,fl);
/*
if(nn1==13){printf("MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();}
printf("MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
*/
		}
printf("MARK %ld ",marknum);
	}
/* Save data in binary format ---------------------------- */
/**/
/**/
/**/
fclose(fl);
if (printmod) printf("OK!\n");
/**/
/* file.t3c file creation */
fl = fopen("file.t3c","wt");
fprintf(fl,"%d \n",f0);
fclose(fl);
/**/
/**/
/*
for (mm1=0;mm1<marknum;mm1++)
{
printf("SAVE MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
}
*/
/* stop.yn file information read */
fl = fopen("stop.yn","rt");
/* Read String */
ffscanf();
/**/
/* Stop Y/N */
if (sa[0]=='y' || sa[0]=='Y')
	{
	fclose(fl);
	printf("PROGRAM TERMINATED FROM stop.yn \n");
	exit(0);
	}
/**/
/* Change printmod */
if (sa[0]>='0' && sa[0]<='9')
	{
	printmod=atoi(sa);
	}
fclose(fl);
}
/* Print Results to data file ----------------------------------- */



/* LOAD WITHOUT EMPTY LINES =================================== */
void ffscanf()
{
/* Counter */
int m1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	m1=fscanf(fl,"%s",sa);
	/* Check end of file */
	if (m1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(m1=strlen(sa)-1;m1>=0;m1--)
	if (*(sa+m1)<=32)
	*(sa+m1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES =================================== */


/* LOAD WITHOUT EMPTY LINES from fl1 =================================== */
void ffscanf1()
{
/* Counter */
int n1;
/**/
/* Read cycle */
do
	{
	/* Input string */
	n1=fscanf(fl1,"%s",sa);
	/* Check end of file */
	if (n1<1)
		{
		printf("\n Unexpected end of file\n");
		fclose(fl1);
		exit(0);
		}
	/* Delete last symbol <32 */
	for(n1=strlen(sa)-1;n1>=0;n1--)
	if (*(sa+n1)<=32)
	*(sa+n1)=0;
	else
	break;
	}
while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */


/* Calc,Check Parameters of Grid */
void gridcheck()
{
/* Counter */
int n1,n2;
/* Nodes Num */
nodenum=xnumx*ynumy*znumz;
if(nodenum>MAXNOD) {printf("Space out in vx[],vy[],vz[]"); exit(0);}
printf(" Node NUM %ld  ",nodenum);
/**/
/* Parameters Num */
if(nodenum*10>MAXPAR) {printf("Space out in sol0[],bondm[]"); exit(0);}
printf(" Parameter NUM %ld  ",nodenum*10);
/**/
/* Cells Num */
cellnum=(xnumx-1)*(ynumy-1)*(znumz-1);
/**/
/* Mark num */
if(marknum>MAXMRK) {printf("Space out in markx[]"); exit(0);}
printf(" Marker NUM %ld  ",marknum);
/**/
/* Rock types Num */
if (rocknum>MAXTMR){printf("Space out in marknu[]"); exit(0);}
/**/
/* Koef for processing */
xstpx=xsize/(double)(xnumx-1);
ystpy=ysize/(double)(ynumy-1);
zstpz=zsize/(double)(znumz-1);
kfx=1.0/xstpx;
kfy=1.0/ystpy;
kfz=1.0/zstpz;
kfxx=kfx*kfx;
kfyy=kfy*kfy;
kfxy=kfx*kfy;
kfxz=kfx*kfz;
kfyz=kfy*kfz;
kfzz=kfz*kfz;
/**/
/* Spec counters */
/* Node num in one xy slide */
xynumxy=xnumx*ynumy;
/* Cell num in one xy slide */
xynumxy1=(xnumx-1)*(ynumy-1);
/* Cell num in x,y,z directions */
xnumx1=xnumx-1;
ynumy1=ynumy-1;
znumz1=znumz-1;
/* Counters for sol[] */
nodenum2=nodenum*2;
nodenum3=nodenum*3;
nodenum4=nodenum*4;
nodenum8=nodenum*8;
/* Multigrid variables */
mgn[0]=nodenum;
mgp[0]=0;
mgs[0]=1;
mgx[0]=xnumx;
mgy[0]=ynumy;
mgz[0]=znumz;
mgxy[0]=xnumx*ynumy;
/*
printf("%ld %ld %ld %d",xnumx,ynumy,znumz,multinum);getchar();
*/
for(n1=1;n1<=multinum;n1++) 
	{
	mgp[n1]=mgp[n1-1]+mgn[n1-1];
	mgs[n1]=mgs[n1-1]*2;
	mgx[n1]=((mgx[n1-1]-5)/2)+5;
	mgy[n1]=((mgy[n1-1]-5)/2)+5;
	mgz[n1]=((mgz[n1-1]-5)/2)+5;
	mgxy[n1]=mgx[n1]*mgy[n1];
	mgn[n1]=mgx[n1]*mgy[n1]*mgz[n1];
	}
/* Multigrid nodes positions */
for(n1=0;n1<=multinum;n1++) 
	{
	/* Set external X grid nodes */
	mggx[n1][0]=gx[0];
	mggx[n1][1]=gx[1];
	mggx[n1][2]=gx[2];
	mggx[n1][mgx[n1]-1]=gx[xnumx-1];
	mggx[n1][mgx[n1]-2]=gx[xnumx-2];
	mggx[n1][mgx[n1]-3]=gx[xnumx-3];
	/* Set internal X grid nodes */
	for(n2=3;n2<mgx[n1]-3;n2++) 
		{
		mggx[n1][n2]=gx[(n2-2)*mgs[n1]+2];
		}
/*
for(n2=0;n2<mgx[n1];n2++)printf("GX %d %d %e %e\n",n1,n2,mggx[n1][n2],gx[n2]);getchar(); 
*/
	/* Set external Y grid nodes */
	mggy[n1][0]=gy[0];
	mggy[n1][1]=gy[1];
	mggy[n1][2]=gy[2];
	mggy[n1][mgy[n1]-1]=gy[ynumy-1];
	mggy[n1][mgy[n1]-2]=gy[ynumy-2];
	mggy[n1][mgy[n1]-3]=gy[ynumy-3];
	/* Set internal Y grid nodes */
	for(n2=3;n2<mgy[n1]-3;n2++) 
		{
		mggy[n1][n2]=gy[(n2-2)*mgs[n1]+2];
		}
/*
for(n2=0;n2<mgy[n1];n2++)printf("GY %d %d %e \n",n1,n2,mggy[n1][n2]);getchar(); 
*/
	/* Set external Z grid nodes */
	mggz[n1][0]=gz[0];
	mggz[n1][1]=gz[1];
	mggz[n1][2]=gz[2];
	mggz[n1][mgz[n1]-1]=gz[znumz-1];
	mggz[n1][mgz[n1]-2]=gz[znumz-2];
	mggz[n1][mgz[n1]-3]=gz[znumz-3];
	/* Set internal Z grid nodes */
	for(n2=3;n2<mgz[n1]-3;n2++) 
		{
		mggz[n1][n2]=gz[(n2-2)*mgs[n1]+2];
		}
/*
for(n2=0;n2<mgz[n1];n2++)printf("GZ %d %d %e \n",n1,n2,mggz[n1][n2]);getchar(); 
*/
	}
}
/* Calc,Check Parameters of Grid */
