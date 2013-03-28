/*    This file is part of VecVis
 *
 *    VecVis is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    VecVis is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>
 *==============================================================================
 */
/*File: FlwReader.c
 *
 * Author: David Warne (david.warne@qut.edu.au)
 * 
 * Summary: Implementation of Reader function for vector fields
 */
 
#include "FlwReader.h"
 
 /* V_ReadHeader(): Reader header information of custom vectorField format
  */
 unsigned char V_ReadHeader(const char* filename,
                   int* total_verts,
                   int* total_elements, 
                   int* numBoundaries,
                   int* tsteps,
                   float* bounds,
                   int* numScalars,
                   fpos_t* vertpos, 
                   fpos_t* elementpos, 
                   fpos_t* fluxpos,
                   fpos_t* scalarpos,
                   fpos_t* boundarypos)
{
  	FILE* fp;
    char buffer[255]; 
    char c = 'a';
 
    /* open the file for writing */
    if (!(fp = fopen(filename,"r")))
    {     
        return 1;
    }
    
    c = fgetc(fp);
    while(c!=EOF)
    {
        /* # indicates a tag */
        if (c == '#')
        {
            int i=0;
            while (c != ' ' && c!= '\n')
            {
                buffer[i] = c;
                c = fgetc(fp);
                i++;
            }
            buffer[i] = '\0'; /* null terminate the string */
            
            /*check the tag and act accordingly*/
            if (!strcmp(buffer,"#NUMVERTS"))
            {
                fscanf(fp,"%d",total_verts);
            }
            else if (!strcmp(buffer,"#NUMPOLYS"))
            {
                fscanf(fp,"%d",total_elements);
            }
            else if (!strcmp(buffer,"#NUMTIMESTEPS"))
            {
                fscanf(fp,"%d",tsteps);
            }
            else if (!strcmp(buffer,"#NUMBOUNDARIES"))
            {
                fscanf(fp,"%d",numBoundaries);
            }
            else if (!strcmp(buffer,"#NUMSCALARS"))
            {
            	fscanf(fp,"%d",numScalars);
            }
            else if(!strcmp(buffer,"#VERTS"))
            {
               fgetpos(fp, vertpos); 
            }
            else if (!strcmp(buffer,"#POLYS"))
            {
                 fgetpos(fp, elementpos);
            }
            else if(!strcmp(buffer,"#FLUXES"))
            {
                fgetpos(fp, fluxpos); 
            }
            else if(!strcmp(buffer,"#XMAX"))
            {
                fscanf(fp,"%f",&(bounds[0]));
            }
            else if(!strcmp(buffer,"#XMIN"))
            {
                fscanf(fp,"%f",&(bounds[1]));
            }
            else if(!strcmp(buffer,"#YMAX"))
            {
                fscanf(fp,"%f",&(bounds[2]));
            }
            else if(!strcmp(buffer,"#YMIN"))
            {
                fscanf(fp,"%f",&(bounds[3]));
            }
            else if(!strcmp(buffer,"#SCALARS"))
            {
                fgetpos(fp,scalarpos);
            }
            else if(!strcmp(buffer,"#BOUNDARIES"))
            {
            
                fgetpos(fp,boundarypos);
            }
        }
        c = fgetc(fp);
    }
    fclose(fp);
    return 0;
} 
 
 /*V_ReadData(): Reads data into data structures in which memeory has been allocated
  *              (X[i],Y[i]) is the position vector of vertex i. Vx[i][t] Vy[i][t] is the 
  *              flow field sampled at vector i at time t
  */
 unsigned char V_ReadData(const char* filename,
                   int total_verts,
                   int total_elements,
                   int numBoundaries, 
                   int tsteps,
                   int numScalars,
                   float* X,
                   float* Y,
                   float** Vx,
                   float** Vy,
                   int** elements,
                   int* numVerts,
				   int* numVertsB,
                   int** boundaries,
                   float*** scalars,
                   fpos_t* vertpos, 
                   fpos_t* elementpos, 
                   fpos_t* fluxpos,
                   fpos_t* scalarpos,
                   fpos_t* boundarypos)
{
	FILE* fp;
    char c = 'a';
    double tmp1,tmp2,tmp3,tmp4;
    double tmpc;
    int i,j,t;
    /*open the file*/
    if (!(fp=fopen(filename,"r")))
    {
        fprintf(fp,"ERROR: [%s]\n",filename);
        return 1;
    }
 
    /* read vertex data*/
    fsetpos(fp,vertpos);
     
    for(i=0;i<total_verts;i++)
    {       
        fscanf(fp,"%f %f",&(X[i]),&(Y[i]));
    }
    
    fsetpos(fp,elementpos);
    for(i=0;i<total_elements;i++)
    {
    	fscanf(fp,"%d",&(numVerts[i]));
        for(j=0;j<numVerts[i];j++)
        {
            fscanf(fp,"%d",&(elements[i][j]));
		
        }
    }
    
    fsetpos(fp,fluxpos);
    /* for now just read the volumetric fluxes */ 
    
    for (t=0;t<tsteps;t++)
    {
        for(i=0; i<total_verts;i++)
        {
        
            fscanf(fp,"%lf %lf %lf %lf",&tmp1,&tmp2,&tmp3,&tmp4);
            Vx[t][i] = (float)tmp1;
            Vy[t][i] = (float)tmp2;
        }
    }
    
    fsetpos(fp,scalarpos);
    
    for (t=0;t<tsteps;t++)
	{
        for (i=0;i<total_verts;i++)
        {
        	for (j=0;j<numScalars;j++)
        	{
            	fscanf(fp,"%lf",&tmpc);
            	scalars[j][t][i] = (float)tmpc;
            }
        }
    }
    
    fsetpos(fp,boundarypos);
    for (i=0;i<numBoundaries;i++)
	{
    	fscanf(fp,"%d",&(numVertsB[i]));
        for (j=0;j<numVertsB[i];j++)
		{
            fscanf(fp,"%d",&(boundaries[i][j]));            
        }
    }
    return 1; 
}
