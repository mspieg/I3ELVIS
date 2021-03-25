
/*
 
 Overview:

	Three different paraview writers are provided in this file
		i3vis_writer_vtr_ascii()
		i3vis_writer_vtr_binary_appended_()
		i3vis_writer_vtr_binary_appended_gz()
	Presently, only three nodal fields are outputted, namely
		viscosity, temperature and density

	This can be modified via the function calls 
		i3vis_WriteField_scalar_XXX(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
		i3vis_WriteField_scalar_XXX(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
		i3vis_WriteField_scalar_XXX(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
	where
		XXX = { ascii, binary_appended, binary_appended_gz }
 
	One should be careful to match the data type of the field with the definition of the
	variable used in head3mg.c
 
  The first time saver_paraview() is called, a PVD file is created with the 
	following naming convention; the 4-letter model name followed by a time stamp formatted like: 
		YYYY.MM.DD_HH.MM.SS
	For example
		s3ce_2012.02.19_17:22:13.pvd
		comp_s3ce_2012.02.19_17:22:13.pvd

	The time stamp is required to ensure a unique pvd file name is assigned each time ivis is restarted.
	[I couldn't find a variable which defined the total number of time steps completed by the simulation....]
 
 Usage:

	To utilise the functionality in this file, you need to do the following things
 
	1) Edit in3mg.c i3mg.c and insert the line at 
		#include "i3vis_output_paraview.c"
	underneath the other files included.
 
	2) Edit head3mg.c and insert the function prototype for saver_paraview().
	I.e., insert the line
		void saver_paraview(int f0, int n0);
	underneath all the other function prototypes.
 
	3) Include a function call to saver_paraview() where ever saver() is called in the ivis source files.
 OR
		Leave the source in original form and add this macro into in3mg.c and i3mg.c

 #define saver(a,b); \
 { \
	saver((a),(b)); \
	saver_paraview((a),(b)); \
 } \
 
	4) Figure out whether the machine is a big endian or little endian.
	To do this, compile and run checkEndian.c ON THE MACHINE you intend to run i3vis on.

	gcc -o checkEndian checkEndian.c
	
	Then modify the lines, marked via 
	// [B] //  
	below to indicate whether the target machine uses big or little endian integers.
 
	5) If you want to use zlib to compress your data files, either
	i) Re-compile your code with the addition of the cflag
			-D__HAVE_ZLIB__
	and indicate you wannt to link against libz.so, via
			-lz.
	For example:
		gcc -Wall -O2 -D__HAVE_ZLIB__ -o main main.c -lz

 ii) You don't have to include -D__HAVE_ZLIB__ statement in the compilation line,
 instead you can simply uncomment the line marked below via // [A] //.
 
 
	~DAM, Feb 16, 2012. (dave.mayhem23@gmail.com)

*/

#define _GNU_SOURCE

// [A] //
//#define __HAVE_ZLIB__

#include "time.h"
#include "sys/time.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>




#ifdef __HAVE_ZLIB__
	#include "zlib.h"
#endif

// [B] //  
//#define WORDSIZE_BIGENDIAN
#define WORDSIZE_LITTLEENDIAN

/* definitions */
typedef enum { CTYPE_INT=0, CTYPE_LONGINT, CTYPE_FLOAT, CTYPE_DOUBLE } CDataType;

/* protoypes */
void i3vis_WriteField_scalar_ascii(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field);
void i3vis_WriteHeaderField_scalar_binary_appended(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,int *offset,void *field);
void i3vis_WriteField_scalar_binary_appended(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field);
#ifdef __HAVE_ZLIB__
void i3vis_WriteHeaderField_scalar_binary_appended_gz(gzFile fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,int *offset,void *field);
void i3vis_WriteField_scalar_binary_appended_gz(gzFile fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field);
#endif



/* ==================================================================================================================== */
/* ascii data */
void ivis_writefield_nodalfields_a(FILE *fp,long int M,long int N,long int P)
{
	i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
	i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
}

/* binary data */
void ivis_writeheader_nodalfields_b(FILE *fp,long int M,long int N,long int P,int *byte_offset)
{
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"viscosity",  M,N,P,byte_offset,nu);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"temperature",M,N,P,byte_offset,tk);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"density",    M,N,P,byte_offset,ro);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"EII",        M,N,P,byte_offset,EII);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"vy", M,N,P,byte_offset,vy);
    i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"vx",    M,N,P,byte_offset,vx);
    i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"vz",    M,N,P,byte_offset,vz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xx",    M,N,P,byte_offset,exx);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain yy",    M,N,P,byte_offset,eyy);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain zz",    M,N,P,byte_offset,ezz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xy",    M,N,P,byte_offset,exy);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xz",    M,N,P,byte_offset,exz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain yz",    M,N,P,byte_offset,eyz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xx",    M,N,P,byte_offset,sxx);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress yy",    M,N,P,byte_offset,syy);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress zz",    M,N,P,byte_offset,szz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xy",    M,N,P,byte_offset,sxy);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xz",    M,N,P,byte_offset,sxz);
	//i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress yz",    M,N,P,byte_offset,syz);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_DOUBLE,"pressure",    M,N,P,byte_offset,pr);
}


void ivis_writefield_nodalfields_b(FILE *fp,long int M,long int N,long int P)
{
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"EII",        M,N,P,EII);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"vy", M,N,P,vy);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"vx",    M,N,P,vx);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"vz",    M,N,P,vz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xx",    M,N,P,exx);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain yy",    M,N,P,eyy);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain zz",    M,N,P,ezz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xy",    M,N,P,exy);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain xz",    M,N,P,exz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"strain yz",    M,N,P,eyz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xx",    M,N,P,sxx);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress yy",    M,N,P,syy);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress zz",    M,N,P,szz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xy",    M,N,P,sxy);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress xz",    M,N,P,sxz);
	//i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"stress yz",    M,N,P,syz);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_DOUBLE,"pressure",    M,N,P,pr);
}

/* binary zippped data */
#ifdef __HAVE_ZLIB__
void ivis_writeheader_nodalfields_bgz(gzFile fp,long int M,long int N,long int P,int *byte_offset)
{
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"viscosity",  M,N,P,byte_offset,nu);
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"temperature",M,N,P,byte_offset,tk);
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"density",    M,N,P,byte_offset,ro);	
}

void ivis_writefield_nodalfields_bgz(gzFile fp,long int M,long int N,long int P)
{
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
}
#endif
/* ==================================================================================================================== */

/* utilities */
void GetTime(double *time)
{
	struct timeval t1;
	double tt;
	
	gettimeofday(&t1,NULL);
	tt = (double)(t1.tv_sec*1000000 + t1.tv_usec);
	
	*time = tt / 1000000.0; /* convert from micro-sec to seconds */
}

void get_formatted_date_time(char date_time[])
{
	time_t      currTime;
	struct tm*  timeInfo;
	int         adjustedYear;
	int         adjustedMonth;
	char        user_name[200];
	
	currTime = time( NULL );
	timeInfo = localtime( &currTime );
	/* See man localtime() for why to adjust these */
	adjustedYear = 1900 + timeInfo->tm_year;
	adjustedMonth = 1 + timeInfo->tm_mon;
	
	/* Format; YYYY.MM.DD_HH.MM.SS */	

	sprintf( date_time, "%.4d.%.2d.%.2d_%.2d.%.2d.%.2d",  
					adjustedYear, adjustedMonth, timeInfo->tm_mday,
					timeInfo->tm_hour, timeInfo->tm_min, timeInfo->tm_sec );
}


/* PVD writer */
void ParaviewPVDOpen(const char pvdfilename[])
{
	FILE *fp = NULL;
	char filesc[100];
	sprintf(filesc,"%s",pvdfilename);
	////
	fp = fopen(filesc,"w");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	///
	
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef WORDSIZE_BIGENDIAN
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	
	fprintf(fp,"<Collection>\n");
	
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");
	fclose(fp);
}

void ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[])
{
	FILE *fp;
	char line[10000];
	int key_L;
	char key[] = "</Collection>";
	char *copy,*tmp;
	char filesc[100];
	sprintf(filesc,"%s",pvdfilename);
	
	fp = fopen(filesc,"r");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	/* reset to start of file */
	rewind(fp);
	
	copy = NULL;
	key_L = strlen( key );
	while( !feof(fp) ) {
		fgets( line, 10000-1, fp );
		if ( strncmp(key,line,key_L)!=0 ) {
			
			/* copy line */
			if (copy!=NULL) {
			  asprintf(&tmp,"%s",copy);
			  free(copy);
				asprintf(&copy,"%s%s",tmp,line);
				free(tmp);
			}
			else {
			  asprintf(&copy,"%s",line);
			}
		}
		else {
			break;
		}
	}
	fclose(fp);
	
	/* open new file - clobbering the old */
	fp = fopen(filesc,"w");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	
	/* write all copied chars */
	fprintf(fp,"%s",copy);
	
	/* write new data */
	if (DirectoryName) {
		fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"%s/%s\"/>\n",time, DirectoryName, datafile );
	} else {
		fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"./%s\"/>\n",time, datafile );
	}
	
	/* close tag */
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");
	
	fclose(fp);
	free(copy);
}


/* reorder ivis array */
void ivis_2_pv_permute_nodal_array_F32(const float field[],const long int nx,const long int ny,const long int nz,float *field1[])
{
	long int ivis_idx,i,j,k;
	long int pv_idx;
	float *tmp;
	
	tmp = malloc(sizeof(float)*nx*ny*nz);
	memset(tmp,0,sizeof(float)*nx*ny*nz);
	
	
	for( k=0; k<nz; k++ ) {
		for( j=0; j<ny; j++ ) {
			for( i=0; i<nx; i++ ) {
				ivis_idx = k * (nx*ny) + i*ny + j;
				pv_idx  =  i + j*nx + k*(nx * ny);
				
				tmp[pv_idx] = field[ivis_idx];
			}
		}
	}
	*field1 = tmp;
}

void ivis_2_pv_permute_nodal_array_F64(const double field[],const long int nx,const long int ny,const long int nz,double *field1[])
{
	long int ivis_idx,i,j,k;
	long int pv_idx;
	double *tmp;
	
	tmp = malloc(sizeof(double)*nx*ny*nz);
	memset(tmp,0,sizeof(double)*nx*ny*nz);
	
	
	for( k=0; k<nz; k++ ) {
		for( j=0; j<ny; j++ ) {
			for( i=0; i<nx; i++ ) {
				ivis_idx = k * (nx*ny) + i*ny + j;
				pv_idx  =  i + j*nx + k*(nx * ny);
				
				tmp[pv_idx] = field[ivis_idx];
			}
		}
	}
	*field1 = tmp;
}

void i3vis_WriteField_scalar_ascii(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	long int i,j,k;
	long int cnt,m4;
	
	switch (type) {
		
		case CTYPE_INT:
		{
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_INT. \n",__FUNCTION__,__LINE__);
			exit(EXIT_FAILURE);
		}
			break;
		
		case CTYPE_FLOAT:
		{
			float *field,*field1;
			
			field = (float*)_field;
			ivis_2_pv_permute_nodal_array_F32(field,nx,ny,nz,&field1);
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\">\n",field_name);

			for( k=0; k<nz; k++ ) {
				for( j=0; j<ny; j++ ) {
					for( i=0; i<nx; i++ ) {
						m4 = i + j*nx + k * (nx*ny);
						
						fprintf( fp, "%1.13e ", field1[m4] );
			}}}	fprintf( fp, "\n");
			fprintf(fp,"        </DataArray>\n");
			free(field1);
		}
			break;

		case CTYPE_DOUBLE: 
		{
			double *field,*field1;
			
			field = (double*)_field;
			ivis_2_pv_permute_nodal_array_F64(field,nx,ny,nz,&field1);
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n",field_name);
			
			for( k=0; k<nz; k++ ) {
				for( j=0; j<ny; j++ ) {
					for( i=0; i<nx; i++ ) {
						m4 = i + j*nx + k * (nx*ny);
						
						fprintf( fp, "%1.13e ", field1[m4] );
					}}}	fprintf( fp, "\n");
			fprintf(fp,"        </DataArray>\n");
			free(field1);
		}
			break;

		default:
			break;
	}
}
void i3vis_WriteField_scalar_binary_appended(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	long int len,L;
	
	L = nx * ny * nz;
	
	switch (type) {
			
		case CTYPE_INT:
		{
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_INT. \n",__FUNCTION__,__LINE__);
			exit(EXIT_FAILURE);
		}
			break;
			
		case CTYPE_FLOAT:
		{
			float *field,*field1;
			
			field = (float*)_field;
			ivis_2_pv_permute_nodal_array_F32(field,nx,ny,nz,&field1);
			len = L * sizeof(float);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field1,sizeof(float),L,fp);
			free(field1);
		}
			break;
			
		case CTYPE_DOUBLE:
		{
			double *field,*field1;
			
			field = (double*)_field;
			ivis_2_pv_permute_nodal_array_F64(field,nx,ny,nz,&field1);
			len = L * sizeof(double);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field1,sizeof(double),L,fp);
			free(field1);
		}
			break;
			
		default:
			break;
	}
	
}

#ifdef __HAVE_ZLIB__
void i3vis_WriteField_scalar_binary_appended_gz(gzFile fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	int len,L;
	
	L = nx * ny * nz;
	
	switch (type) {
			
		case CTYPE_INT:
		{
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_INT. \n",__FUNCTION__,__LINE__);
			exit(EXIT_FAILURE);
		}
			break;
			
		case CTYPE_LONGINT:
		{
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_LONGINT. \n",__FUNCTION__,__LINE__);
			exit(EXIT_FAILURE);
		}
			break;
			
		case CTYPE_FLOAT:
		{
			float *field,*field1;
			
			field = (float*)_field;
			ivis_2_pv_permute_nodal_array_F32(field,nx,ny,nz,&field1);
			len = L * sizeof(float);
			
			gzwrite(fp,&len,sizeof(int));
			gzwrite(fp,field1,sizeof(float)*L);
			free(field1);
		}
			break;
			
		case CTYPE_DOUBLE:
		{
			double *field,*field1;
			
			field = (double*)_field;
			ivis_2_pv_permute_nodal_array_F64(field,nx,ny,nz,&field1);
			len = L * sizeof(double);
			
			gzwrite(fp,&len,sizeof(int));
			gzwrite(fp,field1,sizeof(double)*L);
			free(field1);
		}
			break;
			
		default:
			break;
	}
	
}
#endif

/* check memory is < 2 GB in file - known bug with appended data in paraview */
void check_offset_size(int offset)
{
	float memory = (float)offset;
	
	if (memory > 4.0e9) {
		printf("ERROR(%s,%d): \n",__FUNCTION__,__LINE__);
		printf("!!! i3vis_writer_vtr_binary_appended(): |-> Max memory offset = %1.2e GBytes !!!\n", memory*1.0e-9 );
		printf("!!! WARNING: Paraview binary appended data files only supported IF memory offsets are less than 2 GB !!!\n");
		exit(1);
	}
}

void i3vis_WriteHeaderField_scalar_binary_appended(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,int *offset,void *field)
{
	long int i,j,k;
	long int cnt,m4;
	
	switch (type) {
		
		case CTYPE_INT:
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Int32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			fprintf(fp,"        </DataArray>\n");
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(int);
			check_offset_size(*offset);
			break;

		case CTYPE_LONGINT:
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Int64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			fprintf(fp,"        </DataArray>\n");
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(long int);
			check_offset_size(*offset);
			break;
			
		case CTYPE_FLOAT:
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			fprintf(fp,"        </DataArray>\n");
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(float);
			check_offset_size(*offset);
			break;

		case CTYPE_DOUBLE:
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			fprintf(fp,"        </DataArray>\n");
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(double);
			check_offset_size(*offset);
			break;

		default:
			break;
	}
}

#ifdef __HAVE_ZLIB__
void i3vis_WriteHeaderField_scalar_binary_appended_gz(gzFile fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,int *offset,void *field)
{
	
	switch (type) {
			
		case CTYPE_INT:
		{
			gzprintf(fp,"        <DataArray Name=\"%s\" type=\"Int32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			gzprintf(fp,"        </DataArray>\n");
			
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(int);
			check_offset_size(*offset);
		}
			break;
			
		case CTYPE_LONGINT:
		{
			gzprintf(fp,"        <DataArray Name=\"%s\" type=\"Int64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			gzprintf(fp,"        </DataArray>\n");
			
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(long int);
			check_offset_size(*offset);
		}
			break;
			
		case CTYPE_FLOAT:
			gzprintf(fp,"        <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			gzprintf(fp,"        </DataArray>\n");
			
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(float);
			check_offset_size(*offset);
			break;
			
		case CTYPE_DOUBLE:
			gzprintf(fp,"        <DataArray Name=\"%s\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",field_name, *offset);
			gzprintf(fp,"        </DataArray>\n");
			
			(*offset) += sizeof(int) + nx*ny*nz * sizeof(double);
			check_offset_size(*offset);
			break;
			
		default:
			break;
	}
}
#endif

void i3vis_writer_vtr_ascii(const char filename[],long int M,long int N,long int P)
{
	FILE *fp;
	long int i,j,k;
	char filesc[100];
	sprintf(filesc,"%s",filename);
	
	////
	fp = NULL;
	fp = fopen(filesc,"w");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
	
  fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	fprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	
	/* grid coords */
	fprintf(fp,"      <Coordinates>\n");
	fprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( i=0; i<M; i++ ) {
		fprintf( fp, "%1.13e ", gx[i] );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( j=0; j<N; j++ ) {
		fprintf( fp, "%1.13e ", -gy[j] );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( k=0; k<P; k++ ) {
		fprintf( fp, "%1.13e ", gz[k] );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"      </Coordinates>\n");
	
	
	/* grid data */
	fprintf(fp,"      <PointData>\n");
	
	ivis_writefield_nodalfields_a(fp,M,N,P);
	/*
	 i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"viscosity",M,N,P,nu);
	 i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	 i3vis_WriteField_scalar_ascii(fp,CTYPE_FLOAT,"density",M,N,P,ro);	
	 */
	fprintf(fp,"      </PointData>\n");
	
	fprintf(fp,"    </Piece>\n");
  fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"</VTKFile>\n");
	
	fclose(fp);
}

void i3vis_writer_vtr_binary_appended(const char filename[],long int M,long int N,long int P)
{
	FILE *fp;
	long int i,j,k;
	int byte_offset,len;
	char filesc[100];
	sprintf(filesc,"%s",filename);
	
	////
	fp = NULL;
	fp = fopen(filesc,"wb");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
  fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	fprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	
	byte_offset = 0;
	
	/* grid coords */
	fprintf(fp,"      <Coordinates>\n");
	fprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + M * sizeof(double);
	check_offset_size(byte_offset);

	fprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + N * sizeof(double);
	check_offset_size(byte_offset);
	
	fprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"      </Coordinates>\n");
	byte_offset += sizeof(int) + P * sizeof(double);
	check_offset_size(byte_offset);
	
	
	/* grid data */
	fprintf(fp,"      <PointData>\n");
	ivis_writeheader_nodalfields_b(fp,M,N,P,&byte_offset);
	/*
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"viscosity",  M,N,P,&byte_offset,nu);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"temperature",M,N,P,&byte_offset,tk);
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_FLOAT,"density",    M,N,P,&byte_offset,ro);
	*/
	fprintf(fp,"      </PointData>\n");
	
	fprintf(fp,"    </Piece>\n");
  fprintf(fp,"  </RectilinearGrid>\n");
	
	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	
	/* appended data */
	/* coords */
	fprintf(fp,"_");

	len = M*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gx, sizeof(double),M,fp);

	len = N*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gy, sizeof(double),N,fp);

	len = P*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gz, sizeof(double),P,fp);
	
	/* fields */
	ivis_writefield_nodalfields_b(fp,M,N,P);
	/*
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	i3vis_WriteField_scalar_binary_appended(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
	 */
	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");
	
	
	fclose(fp);
}

#ifdef __HAVE_ZLIB__
void i3vis_writer_vtr_binary_appended_gz(const char filename[],long int M,long int N,long int P)
{
	gzFile fp;
	long int i,j,k;
	int byte_offset,len;
	int errnum,ierr;
	const char *info;
	
	////
	if ((fp = gzopen(filename,"wb"))==NULL)  {
		info = gzerror(fp,&errnum);
		printf("ERROR(%s,%d): Cannot open gzip file %s\n",__FUNCTION__,__LINE__,filename);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	gzprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	gzprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
  gzprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	gzprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",M-1,N-1,P-1);
	
	byte_offset = 0;
	
	/* grid coords */
	gzprintf(fp,"      <Coordinates>\n");
	gzprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + M * sizeof(double);
	check_offset_size(byte_offset);
	
	gzprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + N * sizeof(double);
	check_offset_size(byte_offset);
	
	gzprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	gzprintf(fp,"      </Coordinates>\n");
	byte_offset += sizeof(int) + P * sizeof(double);
	check_offset_size(byte_offset);
	
	
	/* grid data */
	gzprintf(fp,"      <PointData>\n");
	ivis_writeheader_nodalfields_bgz(fp,M,N,P,&byte_offset);
	/*
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"viscosity",  M,N,P,&byte_offset,nu);
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"temperature",M,N,P,&byte_offset,tk);
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"density",    M,N,P,&byte_offset,ro);
	*/
	gzprintf(fp,"      </PointData>\n");
	
	gzprintf(fp,"    </Piece>\n");
  gzprintf(fp,"  </RectilinearGrid>\n");
	
	gzprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	
	/* appended data */
	/* coords */
	gzprintf(fp,"_");
	
	len = M*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gx,sizeof(double)*M);
	
	len = N*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gy,sizeof(double)*N);
	
	len = P*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gz,sizeof(double)*P);

	/* fields */
	ivis_writefield_nodalfields_bgz(fp,M,N,P);
	/*
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"viscosity",  M,N,P,nu);
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"temperature",M,N,P,tk);
	i3vis_WriteField_scalar_binary_appended_gz(fp,CTYPE_FLOAT,"density",    M,N,P,ro);
	*/
	
	gzprintf(fp,"  </AppendedData>\n");
	gzprintf(fp,"</VTKFile>\n");
	
	ierr = gzclose(fp);
	info = gzerror(fp,&errnum);
}
#endif


void I3VISView_nodalfields_paraview(void)
{
	char *name;
	int i,L,location;
	double t0,t1;
	static char *pvdfilename;
	static int beenhere = 0;
	double time;
	char jobname[10];
	char filesc[100];
	char filesc2[100];
	
	snprintf(jobname,5,"%s",fl1out); /* the trailing null character counts */

	/* fetch model time from ivis data */
	time = timestep/3.15576e+7;
	
	
	if (beenhere==0) {
		char time_date[500];
		
		/* make a UNIQUE file name based on the current time */
		printf("JOB NAME : %s \n", jobname );
		
		get_formatted_date_time(time_date);
		asprintf( &pvdfilename, "%s_%s.pvd", jobname,time_date );
		
		sprintf(filesc,"%s",pvdfilename);
		printf("  Creating PVD file : %s \n", filesc);
		ParaviewPVDOpen(pvdfilename);
		beenhere = 1;
	}
	
	/* make a file name */
	/*
	asprintf( &name, "%s.vtr", jobname );
	printf("  Creating file : %s \n", name );
	*/
	asprintf(&name,"%s",fl1out);
	/* replace prn with vtr */
	L = strlen(name);
	for (i=0; i<L; i++) {
		if (name[i]=='.') {
			location = i;
			break;
		}
	}
	sprintf(&name[location+1],"vtr");/* +1 skips the . in the name */
	sprintf(filesc2,"%s",name);
	printf("  Creating file : %s \n", filesc2);
	
/*	
  // Don't use the ascii writer unless it's absolutely necessary //
  // Writing is slow compared to binary, and reading ascii into paraview is much slower than binary //
  GetTime(&t0);
	i3vis_writer_vtr_ascii(name,xnumx,ynumy,znumz);
	ParaviewPVDAppend(pvdfilename,time,name,NULL);
	GetTime(&t1);
  printf("  <<nodal fields>> VTR writer (ascii) : %1.2e (sec) \n", t1-t0);
*/	

#ifndef __HAVE_ZLIB__

	GetTime(&t0);
	i3vis_writer_vtr_binary_appended(name,xnumx,ynumy,znumz);
	ParaviewPVDAppend(pvdfilename,time,filesc2,NULL);
	GetTime(&t1);
  printf("  <<nodal fields>> VTR writer (binary_appended) : %1.2e (sec) \n", t1-t0);

#endif
#ifdef __HAVE_ZLIB__
	{
		char *name2;
		
		asprintf(&name2,"%s.gz",name);
		printf("  Writing as : %s \n", name2 );
		GetTime(&t0);
		i3vis_writer_vtr_binary_appended_gz(name2,xnumx,ynumy,znumz);
		ParaviewPVDAppend(pvdfilename,time,filesc2,NULL);
		GetTime(&t1);
		printf("  <<nodal fields>> VTR writer (binary_appended_gzip) : %1.2e (sec) \n", t1-t0);
		free(name2);
	}
#endif
	
	free(name);
	//free(pvdfilename);
}

/* =============== composition ================================ */
void i3vis_noreorder_WriteField_scalar_ascii(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	long int i,j,k;
	long int m4;
	
	switch (type) {
			
		case CTYPE_INT:
		{
			int *field;
			
			field = (int*)_field;
			
			fprintf(fp,"        <DataArray Name=\"%s\" type=\"Int32\" NumberOfComponents=\"1\" format=\"ascii\">\n",field_name);
			
			for( k=0; k<nz; k++ ) {
				for( j=0; j<ny; j++ ) {
					for( i=0; i<nx; i++ ) {
						m4 = i + j*nx + k * (nx*ny);
						
						fprintf( fp, "%d ", field[m4] );
					}}}	fprintf( fp, "\n");
			fprintf(fp,"        </DataArray>\n");
		}
			break;
			
		case CTYPE_LONGINT:
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_LONGINT. \n",__FUNCTION__,__LINE__);
			break;

		case CTYPE_FLOAT:
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_FLOAT. \n",__FUNCTION__,__LINE__);
			break;
			
		case CTYPE_DOUBLE:
			printf("ERROR(%s,%d): Unsupported format provided CTYPE_DOUBLE. \n",__FUNCTION__,__LINE__);
			break;

		default:
			break;
	}
}

void i3vis_noreorder_WriteField_scalar_binary_appended(FILE *fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	long int len,L;
	
	L = nx * ny * nz;
	
	switch (type) {
			
		case CTYPE_INT:
		{
			int *field;
			
			field = (int*)_field;
			len = L * sizeof(int);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field,sizeof(int),L,fp);
		}
			break;
			
		case CTYPE_LONGINT:
		{
			long int *field;
			
			field = (long int*)_field;
			len = L * sizeof(long int);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field,sizeof(long int),L,fp);
		}
			break;

		case CTYPE_FLOAT:
		{
			float *field;
			
			field = (float*)_field;
			len = L * sizeof(float);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field,sizeof(float),L,fp);
		}
			break;
			
		case CTYPE_DOUBLE:
		{
			double *field;
			
			field = (double*)_field;
			len = L * sizeof(double);
			fwrite(&len,sizeof(int),1,fp);
			fwrite(field,sizeof(double),L,fp);
		}
			break;
			
		default:
			break;
	}
	
}

#ifdef __HAVE_ZLIB__
void i3vis_noreorder_WriteField_scalar_binary_appended_gz(gzFile fp,CDataType type,const char field_name[],const long int nx,const long int ny,const long int nz,void *_field)
{
	int len,L;
	
	L = nx * ny * nz;
	
	switch (type) {
		case CTYPE_INT:
			len = L * sizeof(int);
			break;
			
		case CTYPE_LONGINT:
			len = L * sizeof(long int);
			break;
			
		case CTYPE_FLOAT:
			len = L * sizeof(float);
			break;
			
		case CTYPE_DOUBLE:
			len = L * sizeof(double);
			break;
			
		default:
			break;
	}

	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,_field,len);
}
#endif

void i3vis_generate_compo_map(double Lx,double Ly,double Lz,long int mx,long int my,long int mz,
															long int np,float xp[],float yp[],float zp[],char rocktype[],
															int *_compo[])
{
	int *compo;
	long int *closest_marker;
	long int p,c,cm_p,ncells,owner;
	long int ii,jj,kk,cidx,cidx_n;
	double dx,dy,dz,cx,cy,cz,sep1,sep2;
	
	ncells = mx * my * mz;
	compo          = malloc( sizeof(int)*ncells );
	closest_marker = malloc( sizeof(long int)*ncells );

	/* init */
	dx = (Lx)/(double)mx;
	dy = (Ly)/(double)my;
	dz = (Lz)/(double)mz;
	
	for (c=0; c<ncells; c++) {
		compo[c] = -1;
	}
	memset(closest_marker,0,sizeof(long int)*ncells);
	
	/* compute distances */
	for (p=0; p<np; p++) {
		if (rocktype[p] < 0) { continue; }
		if (rocktype[p] > 50) { continue; }

		if (xp[p] < 0.0) { continue; }
		if (yp[p] < 0.0) { continue; }
		if (zp[p] < 0.0) { continue; }
		if (xp[p] > Lx) { continue; }
		if (yp[p] > Ly) { continue; }
		if (zp[p] > Lz) { continue; }
		
		ii = xp[p]/dx;
		jj = yp[p]/dy;
		kk = zp[p]/dz;
		if (ii==mx) { ii--; }
		if (jj==my) { jj--; }
		if (kk==mz) { kk--; }
		cidx = ii + jj * mx + kk * mx * my;
		
		owner = closest_marker[cidx];
		
		if (compo[cidx]!=-1) {
			/* claim */
			closest_marker[cidx] = p;	
			compo[cidx] = rocktype[p];
		} else {
			cm_p = closest_marker[cidx];
			
			/* centroid of cell */
			cx = ii*dx + 0.5*dx;
			cy = jj*dy + 0.5*dy;
			cz = kk*dz + 0.5*dz;
			
			sep1 = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
			sep2 = (xp[p]-cx)*(xp[p]-cx) + (yp[p]-cy)*(yp[p]-cy);
			if (sep2 < sep1) {
				/* new marker is closer */
				closest_marker[cidx] = p;	
				compo[cidx] = rocktype[p];
			}
		}
		
	}

	
	/* check and update using neighbours */
	for (kk=0; kk<mz; kk++) {
		for (jj=0; jj<my; jj++) {
			for (ii=0; ii<mx; ii++) {
				double sep,maxsep = 1.0e32;
				long int closest_p = -1;
				
				cidx = ii + jj*mx + kk*mx*my;
				
				cx = (ii  )*dx + 0.5*dx;
				cy = (jj  )*dy + 0.5*dy;
				cz = (kk  )*dz + 0.5*dz;
				
				if (compo[cidx]!=-1) { continue; }
				
				// check neighbours //
				if (ii-1>=0) {
					cidx_n = (ii-1) + (jj  )*mx + (kk)*mx*my;
					
					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}
				}
				if (ii+1<mx) {
					cidx_n = (ii+1) + (jj  )*mx + (kk)*mx*my;
					
					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}
				}

				if (jj-1>=0) {
					cidx_n = (ii  ) + (jj-1)*mx + (kk)*mx*my;

					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}
				}
				if (jj+1<my) {
					cidx_n = (ii  ) + (jj+1)*mx + (kk)*mx*my;

					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}
				}

				if (kk-1>=0) {
					cidx_n = (ii  ) + (jj  )*mx + (kk-1)*mx*my;

					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}					
				}
				if (kk+1<mz) {
					cidx_n = (ii  ) + (jj  )*mx + (kk+1)*mx*my;

					if (compo[cidx_n]!=-1) { 
						cm_p = closest_marker[cidx_n];
						
						sep = (xp[cm_p]-cx)*(xp[cm_p]-cx) + (yp[cm_p]-cy)*(yp[cm_p]-cy);
						if (sep<maxsep) {
							maxsep = sep;
							closest_p = cm_p;
							compo[cidx] = rocktype[closest_p];
						}
					}					
				}
				/*
				if (compo[cidx]==-1) {
					printf("  cell %ld has empty von-Neumann neighbours \n", cidx);
				}
				*/
				
				/*  NaN (??) */
				/* if(compo[cidx]==0){
					compo[cidx]=(float)(0.0/0.0);
				}*/
			}
		}
	}
	
	{
		long int nempty = 0;
		
		for (c=0; c<ncells; c++) {
			if (compo[c] == -1) {
				nempty++;
			}
		}
		if (nempty!=0) {
			printf("  detected %ld empty cells, implying its von-Neumann neighbours were also empty\n",nempty );
		}
	}
	
	
	free(closest_marker);
	
	*_compo = compo;
}

void i3vis_compo_writer_vtr_ascii(const char filename[],long int MX,long int MY,long int MZ)
{
	FILE *fp;
	long int i,j,k;
	long int M,N,P;
	double dx,dy,dz;
	int *composition;
	
	M = MX + 1;
	N = MY + 1;
	P = MZ + 1;
	
	//dx = (xsize)/(double)MX;
	//dy = (ysize)/(double)MY;
	//dz = (zsize)/(double)MZ;

	dx = xstpx;
        dy = ystpy;
        dz = zstpz;

	MX=xnumx-1;
	MY=ynumy-1;
	MZ=znumz-1;

	
	/* generate crude composition map */
	i3vis_generate_compo_map(xsize,ysize,zsize,MX,MY,MZ, marknum,markx,marky,markz,markt, &composition);
	
	char filesc[100];
	sprintf(filesc,"%s",filename);
	////
	fp = NULL;
	fp = fopen(filesc,"w");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
	
  fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	fprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	
	/* grid coords */
	fprintf(fp,"      <Coordinates>\n");
	fprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( i=0; i<M; i++ ) {
		fprintf( fp, "%1.13e ", i*dx );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( j=0; j<N; j++ ) {
		fprintf( fp, "%1.13e ", j*dy );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for( k=0; k<P; k++ ) {
		fprintf( fp, "%1.13e ", k*dz );
	} fprintf( fp, "\n" );
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"      </Coordinates>\n");
	
	
	/* grid data - nodal */
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"      </PointData>\n");
	
	
	/* grid data - cells */
	fprintf(fp,"      <CellData>\n");
	i3vis_noreorder_WriteField_scalar_ascii(fp,CTYPE_INT,"composition",MX,MY,MZ,composition);
	fprintf(fp,"      </CellData>\n");
	
	
	
	
	fprintf(fp,"    </Piece>\n");
  fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"</VTKFile>\n");
	
	fclose(fp);
	free(composition);
}

void i3vis_compo_writer_vtr_binary_appended(const char filename[],long int MX,long int MY,long int MZ)
{
	FILE *fp;
	int byte_offset,len;
	long int i,j,k;
	long int M,N,P;
	double dx,dy,dz;
	int *composition;
	double *gx,*gy,*gz;
	
	M = MX + 1;
	N = MY + 1;
	P = MZ + 1;


	//dx = xstpx;
        //dy = ystpy;
        //dz = zstpz;

        //MX=xnumx-1;
        //MY=ynumy-1;
        //MZ=znumz-1;

	//M = MX + 1;
        //N = MY + 1;
        //P = MZ + 1;


	
	gx = malloc( sizeof(double)*M );
	gy = malloc( sizeof(double)*N );
	gz = malloc( sizeof(double)*P );
	
	dx = (xsize)/(double)MX;
	dy = (ysize)/(double)MY;
	dz = (zsize)/(double)MZ;

	for (i=0; i<M; i++) { gx[i] = i * dx; }
	for (j=0; j<N; j++) { gy[j] = j * dy; }
	for (k=0; k<P; k++) { gz[k] = k * dz; }
	
	
	/* generate crude composition map */
	i3vis_generate_compo_map(xsize,ysize,zsize,MX,MY,MZ, marknum,markx,marky,markz,markt, &composition);
	
	char filesc[100];
	sprintf(filesc,"%s",filename);
	
	////
	fp = NULL;
	fp = fopen(filesc,"wb");
	if (fp == NULL) {
		printf("ERROR(%s,%d): Cannot open file %s\n",__FUNCTION__,__LINE__,filesc);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
  fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	fprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	
	byte_offset = 0;
	
	/* grid coords */
	fprintf(fp,"      <Coordinates>\n");
	fprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + M * sizeof(double);
	check_offset_size(byte_offset);
	
	fprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + N * sizeof(double);
	check_offset_size(byte_offset);
	
	fprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	fprintf(fp,"        </DataArray>\n");
	fprintf(fp,"      </Coordinates>\n");
	byte_offset += sizeof(int) + P * sizeof(double);
	check_offset_size(byte_offset);
	
	
	/* grid data - nodal */
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"      </PointData>\n");

	/* grid data - cells */
	fprintf(fp,"      <CellData>\n");
	i3vis_WriteHeaderField_scalar_binary_appended(fp,CTYPE_INT,"composition",  MX,MY,MZ,&byte_offset,composition);
	fprintf(fp,"      </CellData>\n");
	
	fprintf(fp,"    </Piece>\n");
  fprintf(fp,"  </RectilinearGrid>\n");
	
	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	
	/* appended data */
	/* coords */
	fprintf(fp,"_");
	
	len = M*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gx, sizeof(double),M,fp);
	
	len = N*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gy, sizeof(double),N,fp);
	
	len = P*sizeof(double);
	fwrite(&len,sizeof(int),1,fp);
	fwrite(gz, sizeof(double),P,fp);
	
	/* fields */
	i3vis_noreorder_WriteField_scalar_binary_appended(fp,CTYPE_INT,"composition",  MX,MY,MZ,composition);
	
	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");
	
	
	fclose(fp);
	free(composition);
	free(gx);
	free(gy);
	free(gz);
}

#ifdef __HAVE_ZLIB__
void i3vis_compo_writer_vtr_binary_appended_gz(const char filename[],long int MX,long int MY,long int MZ)
{
	gzFile fp;
	int errnum,ierr;
	const char *info;	
	int byte_offset,len;
	long int i,j,k;
	long int M,N,P;
	double dx,dy,dz;
	int *composition;
	double *gx,*gy,*gz;
	
	M = MX + 1;
	N = MY + 1;
	P = MZ + 1;
	
	gx = malloc( sizeof(double)*M );
	gy = malloc( sizeof(double)*N );
	gz = malloc( sizeof(double)*P );
	
	dx = (xsize)/(double)MX;
	dy = (ysize)/(double)MY;
	dz = (zsize)/(double)MZ;
	
	for (i=0; i<M; i++) { gx[i] = i * dx; }
	for (j=0; j<N; j++) { gy[j] = j * dy; }
	for (k=0; k<P; k++) { gz[k] = k * dz; }
	
	
	/* generate crude composition map */
	i3vis_generate_compo_map(xsize,ysize,zsize,MX,MY,MZ, marknum,markx,marky,markz,markt, &composition);
	
	////
	if ((fp = gzopen(filename,"wb"))==NULL)  {
		info = gzerror(fp,&errnum);
		printf("ERROR(%s,%d): Cannot open gzip file %s\n",__FUNCTION__,__LINE__,filename);
		exit(EXIT_FAILURE);
	}
	///
	
#ifdef WORDSIZE_LITTLEENDIAN	
	gzprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
#ifdef WORDSIZE_BIGENDIAN	
	gzprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#endif
  gzprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	gzprintf(fp,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">\n",MX,MY,MZ);
	
	byte_offset = 0;
	
	/* grid coords */
	gzprintf(fp,"      <Coordinates>\n");
	gzprintf(fp,"        <DataArray Name=\"xcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + M * sizeof(double);
	check_offset_size(byte_offset);
	
	gzprintf(fp,"        <DataArray Name=\"ycoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	byte_offset += sizeof(int) + N * sizeof(double);
	check_offset_size(byte_offset);
	
	gzprintf(fp,"        <DataArray Name=\"zcoor\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\">\n",byte_offset);
	gzprintf(fp,"        </DataArray>\n");
	gzprintf(fp,"      </Coordinates>\n");
	byte_offset += sizeof(int) + P * sizeof(double);
	check_offset_size(byte_offset);
	
	
	/* grid data - nodal */
	gzprintf(fp,"      <PointData>\n");
	gzprintf(fp,"      </PointData>\n");
	
	/* grid data - cells */
	gzprintf(fp,"      <CellData>\n");
	i3vis_WriteHeaderField_scalar_binary_appended_gz(fp,CTYPE_INT,"composition",  MX,MY,MZ,&byte_offset,composition);
	gzprintf(fp,"      </CellData>\n");
	
	gzprintf(fp,"    </Piece>\n");
  gzprintf(fp,"  </RectilinearGrid>\n");
	
	gzprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	
	/* appended data */
	/* coords */
	gzprintf(fp,"_");
	
	len = M*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gx,sizeof(double)*M);
	
	len = N*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gy,sizeof(double)*N);
	
	len = P*sizeof(double);
	gzwrite(fp,&len,sizeof(int));
	gzwrite(fp,gz,sizeof(double)*P);
	
	/* fields */
	i3vis_noreorder_WriteField_scalar_binary_appended_gz(fp,CTYPE_INT,"composition",  MX,MY,MZ,composition);
	
	gzprintf(fp,"  </AppendedData>\n");
	gzprintf(fp,"</VTKFile>\n");
	
	
	ierr = gzclose(fp);
	info = gzerror(fp,&errnum);

	free(composition);
	free(gx);
	free(gy);
	free(gz);
}
#endif

void I3VISView_compofield_paraview(void)
{
	char *name;
	int i,L,location;
	double t0,t1;
	static char *pvdfilename;
	static int beenhere = 0;
	double time;
	long int compo_grid_mx,compo_grid_my,compo_grid_mz;
	char jobname[10];
	char filesc[100];
	char filesc2[100];

	snprintf(jobname,5,"%s",fl1out); /* the trailing null character counts */
	
	/* fetch model time from ivis data */
	time = timestep/3.15576e+7;
	
	compo_grid_mx = 1.5*(xnumx-1);
	compo_grid_my = 1.5*(ynumy-1);
	compo_grid_mz = 1.5*(znumz-1);
	
	/*changed*/
	/*compo_grid_mx = round(1.75*(xnumx-1) );
        compo_grid_my = round(1.75*(ynumy-1) );
        compo_grid_mz = round(1.75*(znumz-1) );*/



	
	if (beenhere==0) {
		char time_date[500];
		
		/* make a UNIQUE file name based on the current time */
		printf("JOB NAME : %s \n", jobname );
		
		get_formatted_date_time(time_date);
		asprintf( &pvdfilename, "comp_%s_%s.pvd", jobname,time_date );
		
		sprintf(filesc,"%s",pvdfilename);
		printf("  Creating PVD file : %s \n", filesc );
		ParaviewPVDOpen(pvdfilename);
		beenhere = 1;
	}
	
	/* make a file name */
	/*
	asprintf( &name, "comp_%s.vtr", jobname );
	printf("  Creating file : %s \n", name );
	*/
	asprintf(&name,"comp_%s",fl1out);
	/* replace prn with vtr */
	L = strlen(name);
	for (i=0; i<L; i++) {
		if (name[i]=='.') {
			location = i;
			break;
		}
	}
	sprintf(&name[location+1],"vtr");/* +1 skips the . in the name */
	sprintf(filesc2,"%s",name);
	printf("  Creating file : %s \n", filesc2 );
	
	
	/*	
	 // Don't use the ascii writer unless it's absolutely necessary //
	 // Writing is slow compared to binary, and reading ascii into paraview is much slower than binary //
	GetTime(&t0);
	
	// choose a mesh for the comp field //
	i3vis_compo_writer_vtr_ascii(name,compo_grid_mx,compo_grid_my,compo_grid_mz);

	ParaviewPVDAppend(pvdfilename,time,name,NULL);
	GetTime(&t1);
	printf("  <<cell composition>> VTR writer (ascii) : %1.2e (sec) \n", t1-t0);
*/
	
//////	
#ifndef __HAVE_ZLIB__
	GetTime(&t0);
	// choose a mesh for the comp field //
	i3vis_compo_writer_vtr_binary_appended(name,compo_grid_mx,compo_grid_my,compo_grid_mz);
	ParaviewPVDAppend(pvdfilename,time,filesc2,NULL);
	GetTime(&t1);
  printf("  <<cell composition>> VTR writer (binary_appended) : %1.2e (sec) \n", t1-t0);
#endif
#ifdef __HAVE_ZLIB__
	{
		char *name2;
		
		asprintf(&name2,"%s.gz",name);
		printf("  Writing as : %s \n", name2 );
		GetTime(&t0);
		// choose a mesh for the comp field //
		i3vis_compo_writer_vtr_binary_appended_gz(name2,xnumx-1,ynumy-1,znumz-1);
		ParaviewPVDAppend(filesc,time,filesc2,NULL);
		GetTime(&t1);
		printf("  <<cell composition>> VTR writer (binary_appended_gzip) : %1.2e (sec) \n", t1-t0);
		free(name2);
	}
#endif
	
	free(name);
//	free(pvdfilename);
}


/*
Wrapper function to call paraview output routines  
*/

/* n0 - circle number */
/* f0 - file number */
void saver_paraview(int f0, int n0)
{

	I3VISView_nodalfields_paraview();
	I3VISView_compofield_paraview();
}
