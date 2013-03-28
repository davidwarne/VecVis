/*   VecVis :  A texture based vector field visualisation tool
 *   Copyright (C) 2011  David J. Warne, Joe Young
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>
 *==============================================================================
 */
/*==============================================================================
 * File: VecVis.c
 *
 * Authors: David J. Warne (david.warne@qut.edu.au)
 * 			Joe Young	(j.young@qut.edu.au)
 *
 * Date Created: 31/10/2009
 * Last Modified: 29/08/2011
 * 
 * Description: Vector field Visualisation tool for both steady and unsteady 2d 
 *				flow fields. This tool implements the Image Based Flow Vizualisation
 * 				algorithm. Allow the user to inject dye (up to 10 injection points).
 *==============================================================================
 */

/* Standard Headers */
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* OpenGL and GLUT */
#include "GL/glut.h"


#include "Flw/FlwReader.h"
#include "BitMap/BitMapWriter.h"
#include "BitMap/BitMapReader.h"

#define VERSION 0.9
#define PI 3.141592653589793



/*==============================================================================
 * 		GLOBALS
 *==============================================================================
 */
/* noise testure size*/
#ifndef NSIZE
#define NSIZE 512
#endif
/*number of noise frames*/
#ifndef NOISE
#define NOISE 64
#endif

#ifndef NSPOT
#define NSPOT 1
#endif
/* default view port dimensions */
int 			HSIZE = 512; 				/* image size */
int      		WSIZE = 512; 				/* image size*/


 
/* texture properties */
float    		nspot = NSPOT; 					/* size of the noise spots*/
float    		T;							/* texture coordinate upper limit*/
GLuint  	    tex[NOISE+2]; 				/* noise and work texture*/
GLuint			jittertex[1];				/* jitter grid texture*/
GLuint 			ImageTex[1]; 				/* background image*/
GLuint			DyeTex[1];					/* Texture for dye overlay*/
GLfloat 		jitterGrid[NSIZE][NSIZE][4];/* jitter grid texture pattern */
int 			phase[NSIZE][NSIZE];		/* temporal noise pattern*/
GLfloat 		noise[NSIZE][NSIZE][2];     /* spatial noise pattern*/
GLfloat			dyeImage[NSIZE][NSIZE][4]; 	/* dye texture*/
float*			TX;
float*			TY;

/* frame counter*/
int     	    frame = 0;

/* time used for calculating fps*/
float	        seconds = 0.0;
unsigned int	tframe = 0; 
unsigned int 	numPolys = 0;
/* number of timesteps*/
int            	total_steps=101;

/* just keeping it the same for now*/
float          	curtime = 0.0;

/* commandline args  */
char*		    filename;  					/* input .flw data file*/
char* 			backgroundImagefile; 		/* .bmp file to show under flow*/
float           dt = 0.05; 					/* fraction of a timestep perframe*/
float 			oldDt = 0.0; 				/* previous timestep setting*/
float           scale = 1.0;           		/* vector field scaling factor*/
float           alpha = 0.1;				/* opacity factor*/
float			dyeAlpha = 0.0;				/* opacity of Dye*/
int 			numScalars = 0;				/* number of scalar data associated with each vertex*/
int 			CurScalar = 0;				/* Currently selected Scalar*/


/* data structures */
float 			bounds[4];		/*  x and y region bounds*/
int* 			numVerts; 					/* array of verts per elements*/
int*			numVertsB;					/* array of verts per boundary*/
float* 			X;							/* array of x-coordinates*/
float* 			Y;							/* array of y-coordinates*/
float** 		Vx;							/* x-velocity component*/
float** 		Vy;							/* y-velocity component*/
float** 		unVx;						/* un-nomalised x-velocity data*/
float** 		unVy;						/* un-nomalised y-velocity data*/
float** 		nVx;						/* normalised x-velocity data*/
float** 		nVy;						/* normalised y-velocity data*/
int				total_Verts;				/* number of vertices in mesh*/
float***		Scalars;					/* array of scalar data sets*/
int 			numElements;				/* number of elements in mesh*/
int 			numBoundaries;				/* number of boundaries in mesh*/
int** 			elements;					/* array of indices pointing to vertices in each element*/
int** 			Boundaries;					/* array of indices pointing to vertices in each boundary*/
float 			xmin, xmax, ymin, ymax;		
float			o_xMax,o_yMax,o_xMin,o_yMin;/*  x and y region bounds*/
float 			xrange_inv;					/* inverse of length of x domain*/
float			yrange_inv; 				/* inverse of length of y domain*/
float			xrange;						/* length of x domain*/
float			yrange;						/* length of y domain*/
float*			PX;							/* integral of X*/
float*			PY;							/* integral of Y*/
float 			Vmax = 0.0;					/* max spatial step*/

/* geometry of circle */
#ifndef CIRCLE_RES
#define	CIRCLE_RES 50
#endif
float			dtheta = (2.0*PI)/((float)(CIRCLE_RES-1));
float			cost[CIRCLE_RES];
float			sint[CIRCLE_RES];
float			circle[CIRCLE_RES];

/* location of dye injection*/
float 			curx = 0.0;
float 			cury = 0.0;

/* other properties of dye*/
#ifndef MAX_DYE
#define	MAX_DYE 10
#endif

#define CIRCLE_STYLE 0
#define CROSS_STYLE 1
#define SQUARE_STYLE 2
#define LINE_STYLE 3
#define NUMSTYLES 4
float 			ds[MAX_DYE];				/* size of dye spot*/ 
int				numDyeIn = 1;				/* number of dye spots in the view*/
int				curDye = 0;					/* currently selected dye spot*/
float 			dyeX[MAX_DYE];				/* dye x-coordinates*/
float			dyeY[MAX_DYE];				/* dye y-coordinates*/
float			dyePhi[MAX_DYE];			/* dye rotation angle*/
unsigned int	dyePulseRate[MAX_DYE];		/* rate of dye injection*/
unsigned short	dyeStyle[MAX_DYE];			/* reder style of dye points*/

/* color of dye*/
float			dyeCol[MAX_DYE][3] = {{1.0,0.0,0.0},
												{0.0,1.0,0.0},
												{0.0,0.0,1.0},
												{1.0,1.0,0.0},
												{1.0,0.0,1.0},
												{0.0,1.0,1.0},
												{0.5,0.5,1.0},
												{0.0,0.5,0.5},
												{0.2,0.3,0.4},
												{0.4,0.2,0.1}}; 
/* our colour maps*/
#ifndef MAXCONTROLPOINTS
#define	MAXCONTROLPOINTS 20
#endif
#ifndef MAXSCALARS
#define MAXSCALARS 10
#endif
#ifndef MAXVERTS
#define MAXVERTS 6
#endif

int				ControlPointNums[MAXSCALARS]; /* transfer function control points*/
float			ScalarMins[MAXSCALARS];		  /* ranges of the scalars*/
float			ScalarMaxs[MAXSCALARS];		 
float			ColorMaps[MAXSCALARS][MAXCONTROLPOINTS][5];	/* Transfer functions structures*/
float			inv_Diffs[MAXSCALARS][MAXCONTROLPOINTS]; /* difference in control points*/
char			ScalarNames[MAXSCALARS][100];
float***		colors;
/* anntotions*/
char* 			title;
char* 			dir;
char * 			annotationsFile;
char* 			tunits; /* label for the time units*/
char 			sunits[100]; /* label for spatial units*/

/* mapping world space size of a pixel*/
float 			xh;
float 			yh;

/* the the time step in real time units */
float        	real_time_step = 0.0;
#ifndef MAX_OBSERVATION_POINTS
#define MAX_OBSERVATION_POINTS 10
#endif
#define MAX_OBSERVATION_LINES 2
int				curObs[MAX_OBSERVATION_LINES] = {0,0};
int 			curLine = 0;
int 			numObsIn[MAX_OBSERVATION_LINES] = {1,1};
float			obsX[MAX_OBSERVATION_LINES][MAX_OBSERVATION_POINTS];
float			obsY[MAX_OBSERVATION_LINES][MAX_OBSERVATION_POINTS];
unsigned int*	obsVerts[MAX_OBSERVATION_LINES];
unsigned char*	obsVertsSeg[MAX_OBSERVATION_LINES];
unsigned int	numObsVerts[MAX_OBSERVATION_LINES] = {0,0};
float 			lineCol[MAX_OBSERVATION_LINES][4] = {{1.0,0.0,0.0,1.0},{0.0,0.0,1.0,1.0}};

float 			fluxdata[100];
float 			timedata[100];
#define INTERACT_DYE 0
#define INTERACT_OBS 1
int 			mouseSelectType = INTERACT_DYE;

/* Run-time options*/
unsigned char	dye = 0;
unsigned char	drawmesh = 0;
unsigned char	showfps = 0;
unsigned char	outputframes = 0;
unsigned char	jitter = 0;
unsigned char	overlay = 0;
unsigned char	norm = 0;	
unsigned char	colormap = 0;
unsigned char	drawPlot = 0;	
unsigned char	obs = 0;		
	
/*==============================================================================
 * 		FUNCTION PROTOTYPES
 *==============================================================================
 */

/* IBFV Functions*/
float 	F(int);
void 	InitNoiseTextures(float);
void 	InitJitterTexture(float);
void 	InitDyeTexture(float);
void 	InitDye(void);
void 	Euler(float, float , float , float , float* , float* );
void 	Advect(float ,float ,int ,int,float,float, float*,float*, float** ,float** ,float*, float*,float*,float*);
void 	DrawBackground(void);
void 	DrawJitter(void);
void 	DrawDomain(void);
void 	Inject(void);
void 	InjectDye(void);
void 	DrawColorOverlay(float,float,int,int);

/*GLUT callbacks*/
void 	InitGL(void);
void 	Reshape(int, int);
void 	Timer(int);
void 	Timer2(int);
void 	Display(void);
void 	KeyPressed(unsigned char, int, int);
void 	ButtonClick(int, int, int, int);

/* Annotations*/
void 	InterpolateColor(float*, float, float[][5], int, float[]);
void 	BuildColorMaps(void);
void 	DrawText(const char*, void*);
void 	ColorBar(float [][5], int, char*, float, float, float);
void 	AddSpatialExtent(void);
void 	AddTitle(float, float, char*, char*);
void 	AddLabel(float, float, char*, unsigned char);
void 	DrawPlot(float* ,float* ,float*,char*,float*, float*,int ,char*,float*);
void 	PlotProfile(float,float,float,float,char*,float,float,int,int,unsigned int*,unsigned char*, unsigned int,float*,int);
/*Data Import Functions*/
int 	LoadFlow(void);
void 	LoadAnnotations(char*);
int 	ParseCommandline(int, char**);
void 	LoadImage(char *);

/*Output functions*/
int 	OutputFrame(void);
void 	PrintState(void);
float 	getTotalFlux(float*,float*,float*,float*,float*,int);
float 	getTotalFlow(float,float,int,int,unsigned int*,unsigned char*, unsigned int,int);
/*Main and Help Functions*/
void 	PrintHelp(void);
void 	DisplayCommands(void);
int 	main(int, char**);

/*==============================================================================
 * 		IBFV ALGORITHM FUNCTIONS
 *==============================================================================
 */
 
/* F(): Periodic funtion used to generate the random noise image
 *

 * Parameters:
 *		t - integer in [0, 255]
 */
float F(int t)
{
	return (t>127) ? 1 : 0;
}
 
/* InitNoiseTextures(): make the noise images
 *
 * Parameters: 
 *		a - opacity
 */	
void InitNoiseTextures(float a)
{
	unsigned int i,j,k;
	register int t;
	/* make spatial noise */
	for (i=0;i<NSIZE;i++)
	{
		for(j=0;j<NSIZE;j++)
		{
			phase[i][j] = rand()%256;
		}
	}
    
    /* make temporal noise */
	for(k=0;k<NOISE+1;k++)
	{
		t = k*255/NOISE;
		     
		for(i=0;i<NSIZE;i++)
		{
			for(j=0;j<NSIZE;j++)
			{
				t=(t+phase[i][j])%255;
             
				noise[i][j][0] = F(t);
				noise[i][j][1] = a;            
			}            
		}
		
		/* set up texture parameters */
		glBindTexture(GL_TEXTURE_2D,tex[k]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);   
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D,0,2,NSIZE,NSIZE,0,GL_LUMINANCE_ALPHA,GL_FLOAT,noise);
	}
    
    /* work texture */
    for(i=0;i<NSIZE;i++)
	{
		for(j=0;j<NSIZE;j++)
		{    
			noise[i][j][0] = F(0);
			noise[i][j][1] = 1.0;
		}
	}
	
	/* set up texture parameters */
	glBindTexture(GL_TEXTURE_2D,tex[NOISE+1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);   
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D,0,2,NSIZE,NSIZE,0,GL_LUMINANCE_ALPHA,GL_FLOAT,noise);
}

/* InitDyeTexture(): generates a  texture for to use in for the dye advection
 *					 overlay.
 * Parameters:
 *		a - opacity
 */
void InitDyeTexture(float a)
{
	unsigned int i,j;
	 for(i=0;i<NSIZE;i++)
	{
		for(j=0;j<NSIZE;j++)
		{    
			dyeImage[i][j][0] = 0.0;
			dyeImage[i][j][1] = 0.0;
			dyeImage[i][j][2] = 0.0;
			dyeImage[i][j][3] = a;
		}
	}
	
	/* set up texture parameters*/
	glBindTexture(GL_TEXTURE_2D,DyeTex[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);   
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D,0,4,NSIZE,NSIZE,0,GL_RGBA,GL_FLOAT,dyeImage);
}

/* InitJitterTexture(): generates a jitterGrid texture for to use in for the jitter
 *						gird overlay.
 * Parameters:
 *		a - opacity
 */
void InitJitterTexture(float a)
{
	unsigned int i,j,step;
	step = 16;
	for (i=0;i<NSIZE;i++)
	{
		for (j=0;j<NSIZE;j++)
		{
			jitterGrid[i][j][0] = 0.0;
			jitterGrid[i][j][1] = 0.0;
			jitterGrid[i][j][2] = 0.0;
			jitterGrid[i][j][3] = 0.0;
		}
	}
	/* make random perturbation of jitter grid points */
	for (i=0;i<NSIZE;i++)
	{
		for(j=0;j<NSIZE;j++)
		{
			phase[i][j] = rand()%256;
		}
	}
	
	/* draw jitter grid to texture */
	for (i=step;i<NSIZE;i+=step)
	{
		for (j=step;j<NSIZE;j+=step)
		{
			int t=(int)F(phase[i][j]%255)*5.0;
			jitterGrid[i+t][j+t][0] = 1.0;
			jitterGrid[i+t][j+t][1] = 0.0;
			jitterGrid[i+t][j+t][2] = 0.0;
			jitterGrid[i+t][j+t][3] = 1.0;
			
			jitterGrid[i+t+1][j+t][0] = 1.0;
			jitterGrid[i+t+1][j+t][1] = 0.0;
			jitterGrid[i+t+1][j+t][2] = 0.0;
			jitterGrid[i+t+1][j+t][3] = 1.0;
			
			jitterGrid[i+t+1][j+t+1][0] = 1.0;
			jitterGrid[i+t+1][j+t+1][1] = 0.0;
			jitterGrid[i+t+1][j+t+1][2] = 0.0;
			jitterGrid[i+t+1][j+t+1][3] = 1.0;
			
			jitterGrid[i+t-1][j+t+1][0] = 1.0;
			jitterGrid[i+t-1][j+t+1][1] = 0.0;
			jitterGrid[i+t-1][j+t+1][2] = 0.0;
			jitterGrid[i+t-1][j+t+1][3] = 1.0;
			
			jitterGrid[i+t-1][j+t-1][0] = 1.0;
			jitterGrid[i+t-1][j+t-1][1] = 0.0;
			jitterGrid[i+t-1][j+t-1][2] = 0.0;
			jitterGrid[i+t-1][j+t-1][3] = 1.0;
			
			jitterGrid[i+t+1][j+t-1][0] = 1.0;
			jitterGrid[i+t+1][j+t-1][1] = 0.0;
			jitterGrid[i+t+1][j+t-1][2] = 0.0;
			jitterGrid[i+t+1][j+t-1][3] = 1.0;
			
			jitterGrid[i+t][j+t+1][0] = 1.0;
			jitterGrid[i+t][j+t+1][1] = 0.0;
			jitterGrid[i+t][j+t+1][2] = 0.0;
			jitterGrid[i+t][j+t+1][3] = 1.0;
			
			jitterGrid[i+t-1][j+t][0] = 1.0;
			jitterGrid[i+t-1][j+t][1] = 0.0;
			jitterGrid[i+t-1][j+t][2] = 0.0;
			jitterGrid[i+t-1][j+t][3] = 1.0;
			
			jitterGrid[i+t][j+t-1][0] = 1.0;
			jitterGrid[i+t][j+t-1][1] = 0.0;
			jitterGrid[i+t][j+t-1][2] = 0.0;
			jitterGrid[i+t][j+t-1][3] = 1.0;
		}
	}
	
	/* set up texture parameters */
	glBindTexture(GL_TEXTURE_2D,jittertex[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);   
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D,0,4,NSIZE,NSIZE,0,GL_RGBA,GL_FLOAT,jitterGrid);
}

void InitDye()
{
	unsigned int i;
	for (i=0;i<MAX_DYE;i++)
	{
		dyePhi[i] = 0.0;
		ds[i] = 0.05;
		ds[i] /= xrange_inv;
		dyePulseRate[i] = 1;
		dyeStyle[i] = CIRCLE_STYLE;
	}

	
	circle[0] = 0.0;
   	sint[0] = 0.0;
   	cost[0] = 1.0;
   	
   	/* to draw circle dye */
   	for (i=1;i<CIRCLE_RES;i++)
   	{
   		circle[i] = circle[i-1]+dtheta;
   		sint[i] = sin(circle[i]);
   		cost[i] = cos(circle[i]);
   	}
}

/* Euler(): integrates the vector field based on the loaded data
 *           uses first order forward Euler method. assumes t as 1
 *
 *	Parameters:
 *		x - x coordinate of vector (vx,vy)
 *		y - y coordinate of vector (vx, vy)
 *		vx - x vector component
 *		vy - y vector component 
 *		px - address to write integral result in x
 *		py - address to write integral result in y
 */
void Euler(float x, float y, float vx, float vy, float* px, float* py)
{
	float dx, dy;
    
	dx = vx*scale;
	dy = vy*scale;
   
	*px = x + dx;
	*py = y + dy;
}
 

/* Advect(): Advect the mesh vertex
 *
 *	Parameters:
 *		i1 - index of previous timestep
 *		i2 - index of next timestep
 *		s1,s2 - interpolation weights
 */
void Advect(float s1,float s2,int i1,int i2,float n,float scale, float* px,float* py, float** vx,float** vy,float* x, float* y,float* tx,float* ty)
{
	unsigned int len,i;
	len = n;
	/* integrate the field at each vertex to warp the mesh */
	for(i=0;i<len;i++)
	{
		px[i] = x[i] + (vx[i1][i]*s2+vx[i2][i]*s1)*scale;
	}
	for (i=0;i<len;i++)
	{
		py[i] = y[i] + (vy[i1][i]*s2+vy[i2][i]*s1)*scale;
	}
	/* compute texture coordinates of previous time step*/		
	for (i=0;i<len;i++)
	{
		tx[i] = (x[i] - xmin)*xrange_inv;
	}
	for (i=0;i<len;i++)
	{
		ty[i] = (y[i] - ymin)*yrange_inv;
	}
}


/* DrawWarpedMesh(): Draw the mesh using the integrated positions 
 */
void DrawWarpedMesh()
{
	register int tmp;
	unsigned int i,j;
	/* Draw each warped elements */	 
	for(i=0;i<numElements;i++)
	{
		glBegin(GL_POLYGON);
			for (j=0;j<numVerts[i];j++)
			{
				tmp = elements[i][j]; 
				glTexCoord2f(TX[tmp],TY[tmp]); 
				glVertex2f(PX[tmp],PY[tmp]);    
			}    
		glEnd();
	}
}

/* DawMesh(): Draws the mesh using the original positions
 */
void DrawMesh()
{
	register int tmp;
	unsigned int i,j;
	/* Draw each element */ 		 
	for(i=0;i<numElements;i++)
	{
		glBegin(GL_POLYGON);
			for (j=0;j<numVerts[i];j++)
			{
				tmp = elements[i][j];
				glTexCoord2f(TX[tmp],TY[tmp]); 
				glVertex2f(X[tmp],Y[tmp]);        
			}    
		glEnd();
	}
}
 
/* DrawBackground(): This Draws black background or a user loaded bitmap image 
 */
void DrawBackground()
{
	if(overlay)
	{
		glDisable(GL_BLEND);
		glEnable(GL_TEXTURE_2D);
       	glBindTexture(GL_TEXTURE_2D,ImageTex[0]);
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,0.0);  
			glVertex2f(o_xMin,o_yMin);
			
			glTexCoord2f(1.0,0.0);    
			glVertex2f(o_xMax,o_yMin);
			
			glTexCoord2f(1.0,1.0);      
			glVertex2f(o_xMax,o_yMax);
			
			glTexCoord2f(0.0,1.0);    
			glVertex2f(o_xMin,o_yMax);		
		glEnd();
       	glEnable(GL_BLEND);
	}
	else
	{
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_BLEND);
		glColor4f(0.0,0.0,0.0,0.0);
		glBegin(GL_QUADS);
			glVertex2f(o_xMin,o_yMin);
			glVertex2f(o_xMax,o_yMin);
			glVertex2f(o_xMax,o_yMax);
			glVertex2f(o_xMin,o_yMax);		
		glEnd();
		glEnable(GL_BLEND);
		glEnable(GL_TEXTURE_2D);
	}
}

/* DrawJitter(): Overlays the jitter grid
 */
void DrawJitter()
{
	glEnable(GL_BLEND);
	glBindTexture(GL_TEXTURE_2D,jittertex[0]);
	glBegin(GL_QUADS);
		glTexCoord2f(0.0,0.0);  
		glVertex2f(o_xMin-0.05,o_yMin-0.05);
		
		glTexCoord2f(1.0,0.0);    
		glVertex2f(o_xMax+0.05,o_yMin-0.05);
		
		glTexCoord2f(1.0,1.0);      
		glVertex2f(o_xMax+0.05,o_yMax+0.05);
		
		glTexCoord2f(0.0,1.0);    
		glVertex2f(o_xMin-0.05,o_yMax+0.05);		
	glEnd();
	glDisable(GL_BLEND);
}
 
/* DrawDomain(): adds a black overlay for the area outside of the domain, 
 *				 TODO: may not be required
 */
void DrawDomain()
{
	register int tmp;
	unsigned int i,j;
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	for (i=0;i<numBoundaries;i++)
	{
		glBegin(GL_POLYGON);
	        for (j=0;j<numVertsB[i];j++)
	        {
	        	tmp = Boundaries[i][j];
	            glColor4f(0.0,0.0,0.0,1.0); 
	            glVertex2f(X[tmp],Y[tmp]);
	        }
		glEnd();
	}
	glDisable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);
}
 
/*Inject(): Inject Noise
 */
void Inject()
{
	glEnable(GL_BLEND);
	/* next noise texture to use */
	glBindTexture(GL_TEXTURE_2D,tex[frame%NOISE]);

	/* noise injection */
	glBegin(GL_QUADS);
		glTexCoord2f(0.0,0.0);  
		glVertex2f(o_xMin,o_yMin);
		
		glTexCoord2f(T,0.0);    
		glVertex2f(o_xMax,o_yMin);
		
		glTexCoord2f(T,T);      
		glVertex2f(o_xMax,o_yMax);
		
		glTexCoord2f(0.0,T);    
		glVertex2f(o_xMin,o_yMax);
	glEnd();
    
	glDisable(GL_BLEND);  
}

/*InjectDye(): Inject Dye.
 */
void InjectDye()
{
	unsigned int i,t;
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/* dye injection if enabled */
    if (dye)
    {
    	for (i=0;i<numDyeIn;i++)
    	{
    		if (!(frame%dyePulseRate[i]))
    		{
    			
    			glDisable(GL_TEXTURE_2D);
    			glPushMatrix();
    			
    			/* set up local coordinte transformations */
    			glTranslatef(dyeX[i],dyeY[i],0.0);
    			glRotatef(dyePhi[i],0.0,0.0,1.0);
    			
    			/* set the color */
    			glColor4f(dyeCol[i][0],dyeCol[i][1],dyeCol[i][2],1.0);
    			
    			/* draw dye according to the selected dye draw style */
    			switch (dyeStyle[i])
    			{
    				case CIRCLE_STYLE:
    				{
    					glBegin(GL_TRIANGLES);
    					for (t=0;t<CIRCLE_RES-1;t++)
    					{
    						glVertex2f(0.0,0.0);
    						glVertex2f(cost[t]*ds[i],sint[t]*ds[i]);
    						glVertex2f(cost[t+1]*ds[i],sint[t+1]*ds[i]);
    					}
    					/*glVertex2f(0.0,0.0);
    					glVertex2f(cost[CIRCLE_RES-1]*ds[i],sint[CIRCLE_RES-1]*ds[i]);
    					glVertex2f(cost[0]*ds[i],sint[0]*ds[i]);*/
    					glEnd();
    				}
    					break;
    				case CROSS_STYLE:
    				{
    					glLineWidth(5);
    					glBegin(GL_LINES);
    						glVertex2f(ds[i],0.0);
    						glVertex2f(-ds[i],0.0);
    						glVertex2f(0.0,ds[i]);
    						glVertex2f(0.0,-ds[i]);
    					glEnd();
    					glLineWidth(1);
    				}
    					break;
    				case SQUARE_STYLE:
    				{
    					glBegin(GL_QUADS);
    						glVertex2f(-ds[i],-ds[i]);
    						glVertex2f(-ds[i],ds[i]);
    						glVertex2f(ds[i],ds[i]);
    						glVertex2f(ds[i],-ds[i]);
    					glEnd();
    				}
    					break;
    				case LINE_STYLE:
					{
						glLineWidth(8);
						glBegin(GL_LINES);
							glVertex2f(-ds[i],0.0);
							glVertex2f(ds[i],0.0);
						glEnd();
						glLineWidth(1);
					}
						break;
    			}
    			glPopMatrix();
    		}
    	
		}
	}
    
	
    glEnable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);    
}

/* DrawColorOverlay(): Overlay colormapped mesh, based on currently selected scalar
 *	
 *	Parameters:
 *		i1 - index of previous timestep
 *		i2 - index of next timestep
 *		s1,s2 - interpolation weights
 */
void DrawColorOverlay(float s1,float s2,int i1,int i2)
{
	register float val;
	register int tmp;
	float r,g,b,a,tr1,tg1,tb1,ta1,tr2,tg2,tb2,ta2;
    unsigned int ind0,ind1,ind2,ind3;
    unsigned int i,j;
    
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
   
    ind0 = CurScalar*4;
    ind1 = ++ind0;
    ind2 = ++ind0;
    ind3 = ++ind0;
    ind0 -=3;
    
    /* draw the mesh */
    for(i=0;i<numElements;i++)
    {
    	glBegin(GL_POLYGON);
    		for (j=0;j<numVerts[i];j++)
			{
				tmp = elements[i][j];
				
				tr1 = colors[tmp][i1][ind0];
				tg1 = colors[tmp][i1][ind1];
				tb1 = colors[tmp][i1][ind2];
				ta1 = colors[tmp][i1][ind3];
				tr2 = colors[tmp][i2][ind0];
				tg2 = colors[tmp][i2][ind1];
				tb2 = colors[tmp][i2][ind2];
				ta2 = colors[tmp][i2][ind3];
				
				r = tr1*s2 + tr2*s1;
				g = tg1*s2 + tg2*s1	;
				b = tb1*s2 + tb2*s1	;
				a = ta1*s2 + ta2*s1	;
				
				glColor4f(r,g,b,a);
	   			glVertex2f(X[tmp],Y[tmp]);
			}    
		glEnd();
    }    
    
    glDisable(GL_BLEND);
    glEnable(GL_TEXTURE_2D);
}

/* DrawWireFrame(): Draws the mesh skeleton.
 */
void DrawWireFrame()
{
	register int tmp;
    unsigned int i,j;
    
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
   
    for(i=0;i<numElements;i++)
    {
       	glBegin(GL_LINE_LOOP);
       	for (j=0;j<numVerts[i];j++)
       	{
       		tmp = elements[i][j];
           	glColor4f(0.0,0.0,1.0,1.0);
       		glVertex2f(X[tmp],Y[tmp]);
                
    	}    
	   	glEnd();
    }
 
    for(i=0;i<numElements;i++)
    {
    	glBegin(GL_POINTS);
    	for (j=0;j<numVerts[i];j++)
    	{
    		tmp = elements[i][j];
       		glColor4f(0.0,1.0,0.0,1.0);
       		glVertex2f(X[tmp],Y[tmp]);
                
    	}    
       	glEnd();

    }       
    glDisable(GL_BLEND);
}

void CollectObs(int lines)
{
	unsigned int i,j,k; 
	float avgX,avgY;
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_TEXTURE_2D);
	for (i=0;i<numObsIn[lines];i++)
    {
    	glPushMatrix();
    	glColor4fv(lineCol[lines]);
    	/* set up local coordinte transformations */
    	glTranslatef(obsX[lines][i],obsY[lines][i],0.0);
		glLineWidth(5);
    	glBegin(GL_LINES);
    		glVertex2f(ds[0],0.0);
    		glVertex2f(-ds[0],0.0);
    		glVertex2f(0.0,ds[0]);
    		glVertex2f(0.0,-ds[0]);
    	glEnd();
    	glLineWidth(1);
    	glPopMatrix();
    }
    
    avgX = (obsX[lines][0]+obsX[lines][1])*0.5;
	avgY = (obsY[lines][0]+obsY[lines][1])*0.5;
	
	glPushMatrix();	
    if (numObsIn[lines] >=2)
    {	
    	numObsVerts[lines] = 0;
    	for (k=0;k<numObsIn[lines]-1;k++)
    	{
    		float tmp2,tmp3;
			float tmpMin;
			int tmpMax, tmpMinInd;
			float m2,m3,m;
    		avgX = (obsX[lines][k]+obsX[lines][k+1])*0.5;
			avgY = (obsY[lines][k]+obsY[lines][k+1])*0.5;
    		glPushMatrix();
    		glColor4fv(lineCol[lines]);
    		/* set up local coordinte transformations */
    		glTranslatef(avgX,avgY,0.0);
			glLineWidth(5);
    		glBegin(GL_LINES);
    			glVertex2f(obsX[lines][k+1]-avgX,obsY[lines][k+1]-avgY);		
    			glVertex2f(obsX[lines][k]-avgX,obsY[lines][k]-avgY);
    		glEnd();
    		glLineWidth(1);
    		glPopMatrix();
    		tmpMax = 0.0;
			tmpMin = 1000000.0;
			
			m3 = sqrt((obsX[lines][k] - avgX)*(obsX[lines][k] - avgX) + (obsY[lines][k]-avgY)*(obsY[lines][k]-avgY)); 
			m = sqrt((obsY[lines][k+1]-avgY)*(obsY[lines][k+1]-avgY) +(-obsX[lines][k+1]+avgX)*(-obsX[lines][k+1]+avgX) );
			
			/* get intersecting elements */
			for (i=0;i<numElements;i++)
			{
				unsigned char intersects = 0;
				m2 = sqrt((X[elements[i][0]] - avgX)*(X[elements[i][0]] - avgX) + (Y[elements[i][0]]-avgY)*(Y[elements[i][0]]-avgY));
				tmp2 = ((X[elements[i][0]] - avgX)/m2)*((obsY[lines][k+1]-avgY)/m) + ((Y[elements[i][0]]-avgY)/m2)*((-obsX[lines][k+1]+avgX)/m);
		
				tmpMin = fabs(tmp2);
				tmpMinInd = 0;
				for (j=1;j<numVerts[i];j++)
				{
					m2 = sqrt((X[elements[i][j]] - avgX)*(X[elements[i][j]] - avgX) + (Y[elements[i][j]]-avgY)*(Y[elements[i][j]]-avgY));
					
					tmp3 = ((X[elements[i][j]] - avgX)/m2)*((obsY[lines][k+1]-avgY)/m) + ((Y[elements[i][j]]-avgY)/m2)*((-obsX[lines][k+1]+avgX)/m);
			
					intersects = intersects || (tmp2*tmp3 < 0.0 && m2 < m3);
					if (fabs(tmp3) < tmpMinInd)
					{
						tmpMin = fabs(tmp3);
						tmpMinInd = j;
					}
				}
		
				if (intersects)
				{
					unsigned char already_in = 0;
					unsigned int ii;
					for (ii=0;ii<numObsVerts[lines];ii++)
					{
						already_in = already_in || obsVerts[lines][ii] == elements[i][tmpMinInd];
					}
					if (!already_in)
					{	
						obsVerts[lines][numObsVerts[lines]] = elements[i][tmpMinInd];
						obsVertsSeg[lines][numObsVerts[lines]] = k;
						numObsVerts[lines]++;	
					}
				}
			}
		
			glBegin(GL_POINTS);
			for (i=0;i<numObsVerts[lines];i++)
			{	
				glColor4f(1.0,0.0,0.0,1.0);
				glVertex2f(X[obsVerts[lines][i]],Y[obsVerts[lines][i]]);
				
			}
			glEnd();
			for (i=0;i<numObsIn[lines];i++)
    		{
    			char label[11];
    			sprintf(label,"%u",i);
    			AddLabel(obsX[lines][i],obsY[lines][i],label,i == curObs[lines]);
    			
    		}
    	}
    }
    glPopMatrix();	
    
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_BLEND); 
}


 
/*==============================================================================
 *		GLUT CALLBACKS
 *==============================================================================
 */ 
 
/* InitGL(): Initialise OpenGL state
 */
void InitGL()
{
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.0,0.0,0.0,0.0); 
	glEnable(GL_TEXTURE_2D); 
	glGenTextures(NOISE+1,tex);
	glGenTextures(1,jittertex);
	glGenTextures(1,ImageTex); 
	glGenTextures(1,DyeTex);
}
 
/* Reshape(): Reshape callback function -handle the resize
 */
void Reshape(int w, int h)
{
	if ( h == 0) 
	{
		h = 1;
	}
	
	WSIZE = w;
	HSIZE = h;
	
	T = WSIZE/(0.1*NSIZE);
	xh = ((float)(xmax - xmin))/((float)WSIZE);
    yh = ((float)(ymax - ymin))/((float)HSIZE);
    
    glViewport(0,0,w,h);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	/* set up the view transformation */
	gluOrtho2D(xmin,xmax,ymin,ymax);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
 
/* Timer(): timer callback function - refresh screen 60 times a second
 */
void Timer(int x)
{
	frame++;
	tframe++;
	glutPostRedisplay();
	curtime = (curtime <= (float) total_steps) ? curtime + dt : 0.0;
	if(outputframes)
    {
    	OutputFrame();
    }
	glutTimerFunc(16.6666667,Timer,x);
}

/*Timer2(): timer callback function - display performance 1 time a second
 */
void Timer2(int x)
{
	seconds++;
	if (showfps)
	{
		fprintf(stdout,"Frames: %d, Secs: %f FPS: %d Polys/sec: %d\n",frame,seconds,tframe,numPolys);
		tframe ^= tframe;
		numPolys ^= numPolys;   
	}
	glutTimerFunc(1000.0,Timer2,x);
}

/* display(): display callback function - render the frames 
 *
 * TODO: Clean up this function... too much has been added, break up into functions 
 */
void Display()
{ 
	/* for temporal linear interpolation */   
	int tt1 = (int)curtime;
	float s1 = (curtime-(float)tt1);
	float s2 = 1.0-s1;
	int i1 = tt1%total_steps;
	int i2 = (tt1+1)%total_steps;
    unsigned int i;
    
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
    glEnable(GL_TEXTURE_2D);
    
    /*--------------------------------------------------------------------------
     *run the next iteration of IBFV
     *--------------------------------------------------------------------------
     */
    
    Advect( s1, s2, i1, i2, total_Verts, scale, PX, PY, Vx, Vy, X, Y, TX, TY); 
    glBindTexture(GL_TEXTURE_2D,tex[NOISE]);
    DrawWarpedMesh();numPolys+=numElements;
    Inject(); 
	glBindTexture(GL_TEXTURE_2D,tex[NOISE]);    
	glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGB,0,0,WSIZE,HSIZE,0);
    
    /*--------------------------------------------------------------------------
     * Run IBFV for dye advection
     *--------------------------------------------------------------------------    
	 */
	/* advection already computed.. so we dont need to redo this */
	glBindTexture(GL_TEXTURE_2D,DyeTex[0]);
	DrawWarpedMesh();numPolys+=numElements;
	InjectDye();
	if(jitter)
    {
    	DrawJitter(); /*TODO: fix bug here... seems to screw up dye advection*/
    }
	glBindTexture(GL_TEXTURE_2D,DyeTex[0]);    
	glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,0,0,WSIZE,HSIZE,0);
    
    /* -------------------------------------------------------------------------
     * Final render
     *--------------------------------------------------------------------------
     */
    glClear(GL_COLOR_BUFFER_BIT);
    
    DrawBackground(); /*either black or image*/
		
	glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);
	
	/* overlay work texture on the background */
	glBindTexture(GL_TEXTURE_2D,tex[NOISE]);
	DrawMesh();numPolys+=numElements;
	
	/* overlay dye if required */	
	if(dye)
	{
		glBindTexture(GL_TEXTURE_2D,DyeTex[0]);
		DrawMesh();numPolys+=numElements;
	}
	
	/*--------------------------------------------------------------------------
	 * Draw overlays
	 *--------------------------------------------------------------------------
	 */
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_BLEND);
	/* overlay color maps */
	if (colormap)
	{
		DrawColorOverlay(s1,s2,i1,i2);numPolys+=numElements;	
    }
	
	/* overlay wire frame */
    if (drawmesh)
    {
    	DrawWireFrame();numPolys+=numElements;
    }
    
    /* overlay domain */ 
    DrawDomain();
    
	/* draw observation points */ 
	if (obs)
	{
		for(i=0;i<MAX_OBSERVATION_LINES;i++)
	  		CollectObs(i);
	}
	
	/* overlay text annotations */
    glDisable(GL_TEXTURE_2D);
    glLoadIdentity();
    
    if (drawPlot)
	{
		float col[4] = {1.0,1.0,0.0,1.0};
		float pos[2];
		float size[2];
		float lims[4];
		for (i=1;i<100;i++)
		{
    		fluxdata[i-1] = fluxdata[i];
    	}
    	for (i=1;i<100;i++)
    	{
    		timedata[i-1] = timedata[i];
		}
		for(i=0;i<MAX_OBSERVATION_LINES;i++)
		{
    		PlotProfile(xmin+20*xh,ymax-500*yh - 170*i*yh,300,150,ScalarNames[CurScalar],
    					s1,s2,i1,i2,obsVerts[i],obsVertsSeg[i],numObsVerts[i],lineCol[i],i);
		
		/*getTotalFlow(s1,s2,i1,i2,obsVerts[i],obsVertsSeg[i],numObsVerts[i],i);*/
		}
		
		/*timedata[99] = real_time_step*curtime;
		fluxdata[99] =  getTotalFlow(s1,s2,i1,i2,obsVerts[0],obsVertsSeg[0],numObsVerts[0],0);
		pos[0] = xmin+20*xh;
		pos[1] = ymax-500*yh - 170*2*yh;
		size[0] = 300;
		size[1] = 150;
		lims[0] = timedata[0];
		lims[1] = timedata[99];
		lims[2] = -10000;
		lims[3] = 10000;
		*/
		/*DrawPlot( pos,size,lims,"flux",timedata, fluxdata,100,tunits,col);*/
	}
    
    glColor4f(1.0,1.0,1.0,1.0);
    AddTitle(xmin+10*xh,ymax-30*yh,title,tunits);
    AddSpatialExtent();
    
    if (dye)
    {
    	/* label dye injection points and high light selected point */
    	for (i=0;i<numDyeIn;i++)
    	{
    		char label[11];
    		sprintf(label,"%u",i);
    		AddLabel(dyeX[i],dyeY[i],label,i == curDye);
    	}
    }
    
    /* overlay colorbar */
   	if (colormap)
    {
    	ColorBar(ColorMaps[CurScalar],ControlPointNums[CurScalar], ScalarNames[CurScalar],xmin+20*xh,ymax-300*yh,200.0);
    }
    
	glutSwapBuffers();    
}

/* KeyPressed(): handle the ASCII key strokes
 *
 * Parameters:
 *		key - ACSII code of key stroke
 *		x,y - mouse location in window relative coordinates
 */
void KeyPressed(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'd': /* toggle dye on and off */
			dye  = !dye;
			PrintState();
			break;
		case '-': /* move to next dye point */
			curDye = (curDye+1)%numDyeIn;
			PrintState();
			break;
		case '+': /* add a new dye point */
			numDyeIn++;
			curDye = (curDye+1)%numDyeIn;
			PrintState();
			break;
		case 't': /* increase temporal speed */
			dt += 0.005; 
			PrintState();
			break;
		case 'y': /* decrease temporal speed*/
			dt -= 0.005; 
			PrintState();			
			break; 
		case 's': /* toggele temporal pause on and off*/
		{
			register float tmp = dt;
			dt = oldDt;
			oldDt = tmp;
			PrintState();
			break;
		}
		case 'g': /* increase dye size parameter*/
			ds[curDye] +=0.01/xrange_inv;
			PrintState();
			break;
		case 'h': /* decrease dye size parameter*/
			ds[curDye] -= 0.01/xrange_inv;  
			PrintState();
			break;
		case 'a':  /* increase noise texture opacity*/
			alpha += 0.01;
			InitNoiseTextures(alpha);
			PrintState();
			break;
		case 'z': /* decrease noise texture opacity*/
			alpha -= 0.01;
			InitNoiseTextures(alpha); 
			PrintState();
			break;
		case 'c': /* toggle color overlays on and off*/
			colormap = !colormap;
			PrintState();
			break;    
		case 'p': /* move to next scalar to overlay*/
			CurScalar=(CurScalar+1)%numScalars;
            PrintState();
            break;
		case 'l': /* move to previous scalar to overlay*/
			CurScalar=(CurScalar-1)%numScalars;
			PrintState();
			break;
		case 'm': /* toggle wireframe overlay on and off*/
			drawmesh = !drawmesh;
			PrintState();
			break;
		case 'f': /* toggle performance printout on and off*/
			showfps = !showfps;
			break;
		case '?': /* print run-time command list*/
			DisplayCommands();	
			break;
		case 'O': /* toggle frame output on and off*/
			outputframes = !outputframes;
			PrintState();
			break;
		/* zooming*/
		case 'n': /* zoom out*/
			xmax+=0.1*xrange;
			xmin-=0.1*xrange;
			ymin-=0.1*yrange;
			ymax+=0.1*yrange;
			xrange_inv = 1.0/(xmax-xmin);
			yrange_inv = 1.0/(ymax-ymin);
			xrange = 1.0/xrange_inv;
			yrange = 1.0/yrange_inv;
			Reshape(WSIZE,HSIZE);
			break;
		case 'b': /* zoom in*/
			xmax-=0.1*xrange;
			xmin+=0.1*xrange;
			ymax-=0.1*yrange;
			ymin+=0.1*yrange;
			xrange_inv = 1.0/(xmax-xmin);
			yrange_inv = 1.0/(ymax-ymin);
			xrange = 1.0/xrange_inv;
			yrange = 1.0/yrange_inv;
			Reshape(WSIZE,HSIZE);
			break;
		/* panning*/
		case '8': /* pan up*/
			ymin+=0.1*yrange;
			ymax+=0.1*yrange;
			yrange_inv = 1.0/(ymax-ymin);
			yrange = 1.0/yrange_inv;
			Reshape(WSIZE,HSIZE);
			break;
		case '6': /* pan right*/
			xmin+=0.1*xrange;
			xmax+=0.1*xrange;
			xrange_inv = 1.0/(xmax-xmin);
			xrange = 1.0/xrange_inv;
			Reshape(WSIZE,HSIZE);
			break;
		case '2': /* pan down*/
			ymin-=0.1*yrange;
			ymax-=0.1*yrange;
			yrange_inv = 1.0/(ymax-ymin);
			yrange = 1.0/yrange_inv;
			Reshape(WSIZE,HSIZE);
			break;
		case '4': /* pan left*/
			xmin-=0.1*xrange;
			xmax-=0.1*xrange;
			xrange_inv = 1.0/(xmax-xmin);
			xrange = 1.0/xrange_inv;
			Reshape(WSIZE,HSIZE);	
			break;
		case ']': /* increase vector scale*/
			scale*=1.5;
			PrintState();
			break;
		case '[': /* decrease vector scale*/
			scale/=1.5;
			PrintState();
			break;
		case 'j': /* ovelay jitter grid*/
			jitter = !jitter;
			PrintState();
			break;
		case 'N': /* switch between normalised and un-normalise vectors*/
		{
			norm = !norm;
			if (norm)
			{
				Vx = nVx;
				Vy = nVy;
			}
			else
			{
				Vx = unVx;
				Vy = unVy;
			}
			PrintState();
			break;
		}
		case 'I': /* show image background*/
			overlay = !overlay;
			PrintState();
			break;
		case 'A':  /* increase dye alpha*/
			dyeAlpha += 0.01;
			//InitDyeTexture(dyeAlpha);
			InitNoiseTextures(dyeAlpha);
			PrintState();
			break;
		case 'Z':  /* decrease dye alpha*/
			dyeAlpha -= 0.01;
			InitNoiseTextures(dyeAlpha);
			//InitDyeTexture(dyeAlpha);
			PrintState(); 
			break;
		case ',': /* rotate dye counter clock-wise*/
			dyePhi[curDye] -= 5.0;
			PrintState(); 
			break;
		case '.': /* rotate dye clock-wise*/
			dyePhi[curDye] += 5.0;
			PrintState(); 
			break;
		case '{': /* switch to previous dye style*/
			dyeStyle[curDye] = (dyeStyle[curDye]-1)%NUMSTYLES;
			PrintState(); 
			break;
		case '}': /* switch to next dye style*/
			dyeStyle[curDye] = (dyeStyle[curDye]+1)%NUMSTYLES;
			PrintState(); 
			break;
		case '(': /* decrease dye pulses per sec*/
			dyePulseRate[curDye] += 10;
			PrintState(); 
			break;
		case ')': /* increase dye pulses per sec*/
			dyePulseRate[curDye] -= (dyePulseRate[curDye] == 1) ? 1 : 10; 
			PrintState(); 
			break;
		case '7':
			mouseSelectType = INTERACT_DYE;
			break;
		case '9':
			mouseSelectType = INTERACT_OBS;
			break;
		case '=':
			numObsIn[curLine]++;
			curObs[curLine] = (curObs[curLine]+1)%numObsIn[curLine];
			PrintState();
			break;
		case 'P':
			drawPlot = !drawPlot;
			PrintState();
			break;
		case 'o':
			obs = !obs;
			PrintState();
			break;
		case 'L':
			curLine ^= 1;
	}
	  
}



/* ButtonClick(): handle the mouse input
 *
 * Parameters:
 *		button - the button clicked, either left middle or right
 *		state  - state of button either up or down
 *		x,y    - mouse location in window relative coordinates
 */
void ButtonClick(int button, int state, int x, int y)
{
	unsigned int i;
	switch(mouseSelectType)
	{
		case INTERACT_DYE:
		{
			if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
			{
				/* place dye injection point */
				curx = (((float)x)/WSIZE)*(xmax - xmin) +xmin;
				cury =  -(((float)y)/HSIZE)*(ymax-ymin)+ymax;
				dyeX[curDye] = curx;
				dyeY[curDye] = cury;
			}
			else if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
			{
				/* select a placed dye injection point */
				curx = (((float)x)/WSIZE)*(xmax - xmin) +xmin;
				cury =  -(((float)y)/HSIZE)*(ymax-ymin)+ymax;
				for (i=0;i<numDyeIn;i++)
				{
					if(fabs(dyeX[i]-curx)/fabs(xmax-xmin) < 0.01 && fabs(dyeY[i]-cury)/fabs(ymax-ymin) < 0.01)
					{
						curDye = i;
					}
				}
			}
			
		}

			break;
		case INTERACT_OBS:
		{
			if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
			{
				/* place observation point */
				curx = (((float)x)/WSIZE)*(xmax - xmin) +xmin;
				cury =  -(((float)y)/HSIZE)*(ymax-ymin)+ymax;
				obsX[curLine][curObs[curLine]] = curx;
				obsY[curLine][curObs[curLine]] = cury;
			}
			else if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
			{
				/* select a placed observation point */
				curx = (((float)x)/WSIZE)*(xmax - xmin) +xmin;
				cury =  -(((float)y)/HSIZE)*(ymax-ymin)+ymax;
				for (i=0;i<numObsIn[curLine];i++)
				{
					if(fabs(obsX[curLine][i]-curx)/fabs(xmax-xmin) < 0.01 && fabs(obsY[curLine][i]-cury)/fabs(ymax-ymin) < 0.01)
					{
						curObs[curLine] = i;
					}
				}
			}
		}
			break;
	}
}
 
/*==============================================================================
 *		ANNOTATIONS
 *==============================================================================
 */

/* InterpolateColor(): Finds the color of a data point based on a given color map, uses linear interpolation
 *
 * Parameters:
 *		Color 		 - output color RGBA
 *		val 		 - scalar value to map to a color
 *		TransferFunc - Color map to use
 *		numCntlPts	 - number of control points the color map has
 *		dist_inv     _ the inverse of the distance between each control point
 * 
 * TODO: add handling of HSV color model
 */
void InterpolateColor(float* Color, float val, float TransferFunc[][5], int numCntrlPts, float dist_inv[])
{
	/* find the interval containing the Value */
    float s1,s2;
	
	/* if its outside the bounds then we just set it to max or min*/
	if (val <= TransferFunc[0][0])
	{
		Color[0] = TransferFunc[0][1];
		Color[1] = TransferFunc[0][2];
		Color[2] = TransferFunc[0][3];
		Color[3] = TransferFunc[0][4];
	}
	else if (val >= TransferFunc[numCntrlPts-1][0])
	{
		Color[0] = TransferFunc[numCntrlPts-1][1];
		Color[1] = TransferFunc[numCntrlPts-1][2];
		Color[2] = TransferFunc[numCntrlPts-1][3];
		Color[3] = TransferFunc[numCntrlPts-1][4];
	}
	else /* find what interval its in */
	{   
		int i = 0; 
		while (val > TransferFunc[i][0])
		{
			i++;
		}
        
		/* linear interpolation of RGBA values*/
		s1 = (val - TransferFunc[i-1][0])*dist_inv[i];
        s2 = 1.0-s1;
        
		Color[0] = TransferFunc[i-1][1]*s2 + TransferFunc[i][1]*s1;
		Color[1] = TransferFunc[i-1][2]*s2 + TransferFunc[i][2]*s1;
		Color[2] = TransferFunc[i-1][3]*s2 + TransferFunc[i][3]*s1;
		Color[3] = TransferFunc[i-1][4]*s2 + TransferFunc[i][4]*s1;
	}
}
 
/* BuildColorMaps(): sets up color maps for the scalar properties, uses default 
 *					 color map of blue-> white -> red.
 */
void BuildColorMaps()
{
	unsigned int i,j,t;
	/* get the range of each scalar attribute */
 	for (i=0;i<numScalars;i++)
 	{
 		ScalarMins[i] = Scalars[i][0][0];
 		ScalarMaxs[i] = Scalars[i][0][0];
 		for (t=0; t<total_steps;t++)
 		{
 			for (j=0;j<total_Verts;j++)
 			{
 			
 				if (ScalarMins[i] > Scalars[i][t][j])
 				{
 					ScalarMins[i] = Scalars[i][t][j];			
 				}
 				else if (ScalarMaxs[i] < Scalars[i][t][j])
 				{
 					ScalarMaxs[i] = Scalars[i][t][j];
 				}
 			}
 		}
 		
 		/*set all color maps to default blue -> white -> red */
 		ControlPointNums[i] = 3;
 		
 		/* min -> blue */
 		ColorMaps[i][0][0] = ScalarMins[i];
 		ColorMaps[i][0][1] = 0.0;
 		ColorMaps[i][0][2] = 0.0;
 		ColorMaps[i][0][3] = 1.0;
 		ColorMaps[i][0][4] = 0.5;
 		
 		/* midpoint -> white */
 		ColorMaps[i][1][0] = (ScalarMins[i]+ScalarMaxs[i])*0.5;
 		ColorMaps[i][1][1] = 1.0;
 		ColorMaps[i][1][2] = 1.0;
 		ColorMaps[i][1][3] = 1.0;
 		ColorMaps[i][1][4] = 0.5;
 		
 		/* max -> red */
 		ColorMaps[i][2][0] = ScalarMaxs[i];
 		ColorMaps[i][2][1] = 1.0;
 		ColorMaps[i][2][2] = 0.0;
 		ColorMaps[i][2][3] = 0.0;
 		ColorMaps[i][2][4] = 0.5;
 		
 		fprintf(stdout,"Scalar Ranges: (%f -> %f)\n",ScalarMins[i],ScalarMaxs[i]);
 	}
}

void DrawPlot(float* pos,float* size,float* lims,char* title,float* xdata, float* ydata,int dim,char* units,float* col)
{
	char buffer[100];
	float x1,x2,y1,y2;
	float xdata_range,ydata_range;
	float xdata_val, ydata_val;
	float xplot, yplot;
	unsigned int i;
	x1 = pos[0];
	x2 = pos[0]+size[0]*xh;
	y1 = pos[1];
	y2 = pos[1]+size[1]*yh;
	xdata_range = lims[1] - lims[0];
	ydata_range = lims[3] - lims[2];
	
	glColor4f(1.0,1.0,1.0,1.0);
	glRasterPos2f(x1,y2+4*yh);
	DrawText(title,GLUT_BITMAP_HELVETICA_18);
    
    /*draw plot background*/
	glColor4f(0.1,0.1,0.1,1.0);
	glBegin(GL_QUADS); 
		glVertex2f(x1,y1);   
		glVertex2f(x1,y2);
		glVertex2f(x2,y2);
		glVertex2f(x2,y1);
	glEnd();
	
	/* draw grid lines */
	glColor4f(0.3,0.3,0.3,1.0);
	glBegin(GL_LINES);
		for (i=0;i<=10;i++)
		{
			glVertex2f(x1,y1+((float)i)*0.1*(y2-y1));
			glVertex2f(x2,y1+((float)i)*0.1*(y2-y1));
		}
	glEnd();
	
	/* draw border */
	glColor4f(1.0,1.0,1.0,1.0);
	glBegin(GL_LINE_LOOP); 
		glVertex2f(x1,y1);   
		glVertex2f(x1,y2);
		glVertex2f(x2,y2);
		glVertex2f(x2,y1);
	glEnd();
	
	/* Draw the plot */
	glBegin(GL_LINES);
	for (i=0;i<dim-1;i++)
	{
		ydata_val = ydata[i];
	 	xdata_val = xdata[i];
		xplot = x1 + ((xdata_val-lims[0])/xdata_range)*(x2-x1);
		yplot = y1 +((ydata_val-lims[2])/ydata_range)*(y2-y1);		
		glColor4fv(col);
		glVertex2f(xplot,yplot);
		ydata_val = ydata[i+1];
	 	xdata_val = xdata[i+1];
		xplot = x1 + ((xdata_val-lims[0])/xdata_range)*(x2-x1);
		yplot = y1 +((ydata_val-lims[2])/ydata_range)*(y2-y1);		
		glColor4fv(col);
		glVertex2f(xplot,yplot);
	}
	glEnd();
	
	
	glColor4f(1.0,1.0,1.0,1.0);
	sprintf(buffer,"%f (%s)",0.0,units);
	glRasterPos2f(x1,y1-19*yh);
	DrawText(buffer,GLUT_BITMAP_HELVETICA_10);
    sprintf(buffer,"%f",xdata_range);
    glRasterPos2f(x2,y1-19*yh);
    DrawText(buffer,GLUT_BITMAP_HELVETICA_10);
	
	/* draw axis ticks */
	for (i=0;i<10;i++)
	{
		sprintf(buffer,"%0.02f",lims[2] + ((float)i)*0.1*ydata_range);
		glRasterPos2f(x2+10*xh,y1+((float)i)*0.1*(y2-y1));
		DrawText(buffer,GLUT_BITMAP_HELVETICA_10);
	}
}

float getTotalFlow(float s1,float s2,int i1,int i2,unsigned int* selectedVerts,unsigned char* Segment, unsigned int numSelectedVerts,int line)
{
	float* ydata;
	float* xdata;
	float* vydata;
	float* vxdata;
	float* depth;
	float* dist;
	unsigned int total;
	unsigned char curSeg;
	float dx_,dy_;
	unsigned int i,j;
	float flux;
	if (numSelectedVerts <= 0)
	{
		return 0;
	}
	ydata = (float*)malloc(numSelectedVerts*sizeof(float));
	xdata = (float*)malloc(numSelectedVerts*sizeof(float));
	vydata = (float*)malloc(numSelectedVerts*sizeof(float));
	vxdata = (float*)malloc(numSelectedVerts*sizeof(float));
	depth = (float*)malloc(numSelectedVerts*sizeof(float));
	dist = (float*)malloc(numSelectedVerts*sizeof(float));
	total = 0;
	curSeg = 0;

	while (total < numSelectedVerts)
	{
		depth[total] = Scalars[1][i1][selectedVerts[total]]*s2 + Scalars[1][i2][selectedVerts[total]]*s1;
		xdata[total] = X[selectedVerts[total]];
		ydata[total] = Y[selectedVerts[total]];		
		vxdata[total] = Vx[i1][selectedVerts[total]]*s2 + Vx[i2][selectedVerts[total]]*s1;
		vydata[total] = Vy[i1][selectedVerts[total]]*s2 + Vy[i2][selectedVerts[total]]*s1;		
		dx_ = (X[selectedVerts[total]]-obsX[line][curSeg]);
		dy_ = (Y[selectedVerts[total]]-obsY[line][curSeg]);
		
		for (i = curSeg; i != 0 ;i--)
		{
			dx_ += obsX[line][curSeg] - obsX[line][curSeg-1];
			dy_ += obsY[line][curSeg] - obsY[line][curSeg-1];
		}
		dist[total] = sqrt(dx_*dx_+ dy_*dy_);
		if (curSeg != Segment[total])
		{
			curSeg = Segment[total];
		}
		total++;		
	}
	
	for (i=0;i<numSelectedVerts;i++)
	{
		float min = dist[i];
		unsigned int min_i = i;
		for (j=i+1;j<numSelectedVerts;j++)
		{
			if (dist[j] < min)
			{
				min = dist[j];
				min_i = j;
			}
		}
		dist[min_i] = dist[i];
		dist[i] = min;
		
		min = xdata[min_i];
		xdata[min_i] = xdata[i];
		xdata[i] = min;
		
		min = ydata[min_i];
		ydata[min_i] = ydata[i];
		ydata[i] = min;
		
		min = vxdata[min_i];
		vxdata[min_i] = vxdata[i];
		if (isnan(min))
		{
			vxdata[i] = 0;
		}
		else
		{
			vxdata[i] =min;
		}
		min = vydata[min_i];
		vydata[min_i] = vydata[i];
		if (isnan(min))
		{
			vydata[i] = 0;
		}
		else
		{
			vydata[i] = min;
		}
		min = depth[min_i];
		depth[min_i] = depth[i];
		depth[i] = min;
	}
	
	flux = getTotalFlux(vxdata,vydata,depth,xdata,ydata,numSelectedVerts);
	free(depth);
	free(dist);
	free(xdata);
	free(ydata);
	free(vxdata);
	free(vydata);
	return line, flux;
}

void PlotProfile(float xpos,float ypos,float xsize,float ysize,char* title,float s1,float s2,int i1,int i2,unsigned int* selectedVerts,unsigned char* Segment, unsigned int numSelectedVerts,float* col,int line)
{
	char buffer[100];
	float x1,x2,y1,y2;
	float pos[2];
	float size[2];
	float lims[4];
	float * ydata;
	float* xdata;
	unsigned int total;
	unsigned char curSeg;
	float dx_,dy_;
	float m2;
	float val;
    unsigned int i,j;
	pos[0] = xpos;
	pos[1] = ypos;
	size[0] = xsize;
	size[1] = ysize;
	lims[2] = ColorMaps[CurScalar][0][0];
	lims[3] = ColorMaps[CurScalar][ControlPointNums[CurScalar]-1][0];
	
	if (numSelectedVerts <= 0)
	{
		return;
	}
	ydata = (float*)malloc(numSelectedVerts*sizeof(float));
	xdata = (float*)malloc(numSelectedVerts*sizeof(float));
	
	total = 0;
	curSeg = 0;

	while (total < numSelectedVerts)
	{
		ydata[total] = Scalars[CurScalar][i1][selectedVerts[total]]*s2 + Scalars[CurScalar][i2][selectedVerts[total]]*s1;
		dx_ = (X[selectedVerts[total]]-obsX[line][curSeg]);
		dy_ = (Y[selectedVerts[total]]-obsY[line][curSeg]);
		
		for (i = curSeg; i != 0 ;i--)
		{
			dx_ += obsX[line][curSeg] - obsX[line][curSeg-1];
			dy_ += obsY[line][curSeg] - obsY[line][curSeg-1];
		}
		xdata[total] = sqrt(dx_*dx_+ dy_*dy_);
		if (curSeg != Segment[total]){
			curSeg = Segment[total];
		}
		total++;	
		
	}
	/* quick hack sort ()...*/
	for (i=0;i<numSelectedVerts;i++)
	{
		float min = xdata[i];
		unsigned int min_i = i;
		for (j=i+1;j<numSelectedVerts;j++)
		{
			if (xdata[j] < min)
			{
				min = xdata[j];
				min_i = j;
			}
		}
		xdata[min_i] = xdata[i];
		xdata[i] = min;
		min = ydata[min_i];
		ydata[min_i] = ydata[i];
		ydata[i] = min;
	}
	
	
	lims[0] = xdata[0];
	lims[1] = xdata[numSelectedVerts-1];
	sprintf(buffer,"Profile of %s",title);
    
	DrawPlot(pos,size,lims,buffer,xdata,ydata,numSelectedVerts,sunits,col);
	free(ydata);
	free(xdata);
}

/* DrawText(): Draws the specified message.
 *
 *	Parameters:
 *		Title - text to print
 *		font - font to render text as
 */
void DrawText(const char* Title, void* font)
{
	unsigned int i;
	/* write using bitmap chars */ 
	for (i=0;Title[i]!='\0';i++)
	{
		glutBitmapCharacter(font, Title[i]);
	}
}

/* ColorBar(): A color bar for the given color map at the screen position,
 *             xpos and ypos in world coordinates.
 *
 * Parameters:
 *		TransferFunction - Color map to use
 *		numControlPoints - number of control points the color map has
 *		title            - heading of color bar (typically name of scalar and units)
 *		xpos,ypos        - position of lower left corner of the color bar in 
 *						   world coordinates
 *		size 			 - height of color bar
 */
void ColorBar(float TransferFunction[][5],int numControlPoints, char* title, float xpos,float ypos,float size)
{  
	float dy;  /* y-step per controlpoint*/
	float minT; /* min-val*/
	float maxT; /* max-val*/
	float lenT_inv; /* scalar range*/
	unsigned int i;
	char tick[25];
	
	dy = size/((float)numControlPoints-1.0);
	minT = TransferFunction[0][0];
	maxT = TransferFunction[numControlPoints-1][0];
	lenT_inv = 1.0/(maxT - minT);
	
	/* draw black background */
	glColor3f(0.0,0.0,0.0);
	glBegin(GL_QUADS);
		glVertex2f(xpos-10*xh,ypos-3*yh);
		glVertex2f(xpos-10*xh,ypos+((numControlPoints-1)*dy+3)*yh);
		glVertex2f(xpos+20*xh,ypos+((numControlPoints-1)*dy+3)*yh);
		glVertex2f(xpos+20*xh,ypos-3*yh);
	glEnd();

	/* draw the colorbar (two thick lines together) TODO: change to draw quads */
	glLineWidth(10);
	glBegin(GL_LINE_STRIP);
	    for (i = 0;i<numControlPoints;i++)
	    {
	        dy = (size)*(TransferFunction[i][0] - minT)*lenT_inv;
	        glColor3f(TransferFunction[i][1],TransferFunction[i][2],TransferFunction[i][3]);
	        glVertex2f(xpos,ypos+dy*yh);
	    }
	glEnd();
	glBegin(GL_LINE_STRIP);
		for (i = 0;i<numControlPoints;i++)
		{
			dy = (size)*(TransferFunction[i][0] - minT)*lenT_inv;
			glColor3f(TransferFunction[i][1],TransferFunction[i][2],TransferFunction[i][3]);
			glVertex2f(xpos+10*xh,ypos+dy*yh);
		}
	glEnd();
	glLineWidth(1);
    
	/* annotate the colorbar */
	glColor3f(1.0,1.0,1.0);
	
	/* title */
	glRasterPos2f(xpos-5*xh,ypos-20*yh);
	DrawText(title,GLUT_BITMAP_HELVETICA_18);
	
	/* draw control point ticks */
	for (i=0; i<numControlPoints;i++)
	{		
		sprintf(tick,"    %0.2f",TransferFunction[i][0]);
		dy = (size)*(TransferFunction[i][0] - minT)*lenT_inv;
		glRasterPos2f(xpos+10*xh,ypos+(dy-5)*yh);
		DrawText(tick,GLUT_BITMAP_HELVETICA_12);	
	}
}
 
/* AddSpatialExtent(): adds the x and  y axis on the visualisation
 */
void AddSpatialExtent()
{
	unsigned int i;
	char tick[25];
	char toPrint[255];
	
	glColor3f(1.0,1.0,1.0);
    
	glBegin(GL_LINES);
		/* X and Y axes */
		glVertex2f(xmin+3*xh,ymin+3*yh);
		glVertex2f(xmax-3*xh,ymin+3*yh);
		
		glVertex2f(xmax-3*xh, ymin+3*yh);
		glVertex2f(xmax-3*xh,ymax-3*yh);


		/* tick marks */
        for (i=0;i<=10;i++)
		{
			glVertex2f(xmin+i*(xmax - xmin - 6*xh)*0.1 +3*xh,ymin+3*yh);
			glVertex2f(xmin+i*(xmax - xmin - 6*xh)*0.1 +3*xh,ymin+6*yh);
		}
		
		for (i=0;i<=10;i++)
		{
			glVertex2f(xmax-6*xh,ymin+i*(ymax - ymin - 6*yh)*0.1+3*yh);
			glVertex2f(xmax-3*xh,ymin+i*(ymax - ymin - 6*yh)*0.1+3*yh);
		}
	glEnd();

	/* label ticks */
	for (i=0;i<10;i++)
	{
		sprintf(tick,"%0.2f",xmin+i*(xmax - xmin)*0.1);
		glRasterPos2f(xmin+i*(xmax - xmin - 6*xh)*0.1 ,ymin+6*yh);
		DrawText(tick,GLUT_BITMAP_HELVETICA_12);
	}

	for (i=1;i<=10;i++)
	{
		sprintf(tick,"%0.2f",ymin+i*(ymax - ymin)*0.1);
		glRasterPos2f(xmax-50*xh,ymin+i*(ymax - ymin - 6*yh)*0.1-5*yh);
		DrawText(tick,GLUT_BITMAP_HELVETICA_12);
	}
	
	/* mark the units */
	glRasterPos2f(xmin+3*xh,ymin+40*yh);
	sprintf(toPrint,"(%s)",sunits);
	DrawText(toPrint,GLUT_BITMAP_HELVETICA_18);
}
 
/* AddTitle(): Adds a title at the given  x y location, where x and y are in 
 *             world coordinates
 *
 * Parameters:
 *		x,y   - position of Title in world coordinates
 *		title - text for title
 *		units - temporal time units
 */
void AddTitle(float x,float y,char* title,char* units)
{ 
	char toPrint[255];
	glColor3f(1.0,1.0,1.0);
	glRasterPos2f(x,y);
	sprintf(toPrint,"%s         Time: %.02f (%s)",title,real_time_step*curtime,units);
	DrawText(toPrint,GLUT_BITMAP_TIMES_ROMAN_24);	
}

/* AddLabel: adds a smaller text lable, used to annotate other objects
 *
 *	Parameters:
 *		x,y - position of label in world coordinates
 *		label - text for label
 *		highlight - indicates to bold label or not
 */
void AddLabel(float x, float y, char* label, unsigned char highlight)
{
	
	if (highlight)
	{
		glColor3f(1.0,1.0,1.0);
		glRasterPos2f(x,y);		
		DrawText(label,GLUT_BITMAP_HELVETICA_18);	
	}
	else
	{
		glColor3f(0.0,0.0,0.0);
		glRasterPos2f(x,y);
		DrawText(label,GLUT_BITMAP_HELVETICA_10);	
	}
}

/*==============================================================================
 *		DATA CALULATION FUNCTIONS
 *==============================================================================
 */

float getTotalFlux(float* vx,float* vy,float* depth,float* x,float* y,int numpoints)
{
	float flux;
	float m;
	float n[2];
	float v[2];
	float d;
	float area;
	unsigned int i;
	flux = 0.0;
	for (i=1;i<numpoints;i++)
	{
		/* get normal vector */
		n[0] = y[i] - y[i-1];
		n[1] = -(x[i] - x[i-1]);
		m = sqrt(n[0]*n[0] + n[1]*n[1]);
		n[0] /= m;
		n[1] /= m;
		/* for now we just average the flow vectors */
		v[0] = (vx[i]+vx[i-1])*0.5;
		v[1] = (vy[i]+vy[i-1])*0.5;
		d = v[0]*n[0]+v[1]*n[1];
		/* now compute the area */
		area = m*(depth[i] + depth[i-1])*0.5;
		flux += area*d;	
	}
	
	return flux;
}

/*==============================================================================
 *		DATA IMPORT FUNCTIONS
 *==============================================================================
 */

/* LoadFlow(): Read data files and load into memory
 *
 */
int LoadFlow()
{
	/*size parameters*/
	int NV, NE, NT;
	/*loop counters*/
	unsigned int i,j,t;
	/* locations in file of data */
	fpos_t startVerts, startElements,startfluxes,startscalars, Boundarypos;
    
	/* read data header*/
	if(V_ReadHeader(filename,&NV,&NE,&numBoundaries,&NT,bounds,&numScalars,&startVerts,&startElements,&startfluxes,&startscalars,&Boundarypos))
	{
		fprintf(stderr,"Error: could not read [%s]\n",filename);
		return 1;
	}
    
	/* allocate memory to data structures */
	X = (float*)malloc(NV*sizeof(float));
	Y = (float*)malloc(NV*sizeof(float));
	PX = (float*)malloc(NV*sizeof(float));
	PY = (float*)malloc(NV*sizeof(float));
	numVerts = (int*)malloc(NE*sizeof(int));
	numVertsB = (int*)malloc(numBoundaries*sizeof(int));
	
	TX = (float*)malloc(NV*sizeof(float));
	TY = (float*)malloc(NV*sizeof(float));
	
	unVx = (float**)malloc(NT*sizeof(float*));
	for (i=0;i<NT;i++)
	{
		unVx[i] = (float*)malloc(NV*sizeof(float));
	}
    
	unVy = (float**)malloc(NT*sizeof(float*));
	for (i=0;i<NT;i++)
	{
		unVy[i] = (float*)malloc(NV*sizeof(float));
	}
	
	nVx = (float**)malloc(NT*sizeof(float*));
	for (i=0;i<NT;i++)
	{ 
		nVx[i] = (float*)malloc(NV*sizeof(float));
	}
    
	nVy = (float**)malloc(NT*sizeof(float*));
	for (i=0;i<NT;i++)
	{
		nVy[i] = (float*)malloc(NV*sizeof(float));
	}
    
	Scalars = (float***)malloc(MAXSCALARS*sizeof(float**));
	for (i=0;i<MAXSCALARS;i++)
	{
		Scalars[i] = (float**)malloc(NT*sizeof(float*));
		for (j=0;j<NT;j++)
		{
			Scalars[i][j] = (float*)malloc(NV*sizeof(float));
		}
	}

    elements = (int**)malloc(NE*sizeof(int*));
	for(i=0;i<NE;i++)
	{
		elements[i] = (int*)malloc(MAXVERTS*sizeof(int));
	}
	
	
    
	Boundaries = (int**)malloc(numBoundaries*sizeof(int*));
	for (i=0;i<numBoundaries;i++)
	{
		Boundaries[i] = (int*)malloc(MAXVERTS*sizeof(int)); 
	}
    
    for (i=0;i<MAX_OBSERVATION_LINES;i++)
    {
    	obsVerts[i] = (unsigned int*)malloc((NV/10)*sizeof(unsigned int));
    }
    for (i=0;i<MAX_OBSERVATION_LINES;i++)
    {
    	obsVertsSeg[i] = (unsigned char*)malloc((NV/10)*sizeof(unsigned char));
    }
    for (i=0;i<100;i++)
    {
    	fluxdata[i] =0.0;
    }
    for (i=0;i<100;i++)
    {
    	timedata[i] =0.0;
    }
    /* read the data */
	V_ReadData(filename,NV, NE,numBoundaries, NT,numScalars, X,Y,unVx,unVy,elements,numVerts,numVertsB,Boundaries,Scalars,&startVerts,&startElements,&startfluxes,&startscalars,&Boundarypos);
	
	/* get the data bounds */
	for (t=0;t<NT;t++)
	{
		register float r;
		for (i=0;i<NV;i++)
		{
			r = sqrt(unVx[t][i]*unVx[t][i] + unVy[t][i]*unVy[t][i]);
			if (r != 0.0)
			{
           		nVx[t][i] = unVx[t][i]/r;
				nVy[t][i] = unVy[t][i]/r;
			}
			else 
			{
				nVx[t][i] = unVx[t][i];
				nVy[t][i] = unVy[t][i];
				/*NOTE: uncomment below to make zeros flow disappear*/
				/*nVx[t][i] = unVx[t][i]/r;
				nVy[t][i] = unVx[t][i]/r;
				unVx[t][i] = unVx[t][i]/r;
				unVy[t][i] = unVx[t][i]/r;
				*/
			}
		
			if (r>Vmax)
			{
				Vmax = r;
			} 
		}
	}
	
	colors = (float***)malloc(NV*sizeof(float**));
	for (i=0;i<NV;i++)
	{
		colors[i] = (float**)malloc(NT*sizeof(float*));
		for (j=0;j<NT;j++)
			colors[i][j] = (float*)malloc(4*numScalars*sizeof(float));
	}
	
    
	xmax = bounds[0]; 
	xmin = bounds[1]; 
	ymax = bounds[2]; 
	ymin = bounds[3];

	/* we divide by these alot, so here we reduce division required */
	xrange_inv = 1.0/(xmax-xmin);
	yrange_inv = 1.0/(ymax-ymin);

    /* get the size of a pixel */
    xh = ((float)(xmax - xmin))/((float)WSIZE);
    yh = ((float)(ymax - ymin))/((float)HSIZE);
    
	total_steps = NT; 
	total_Verts = NV;
	numElements = NE;
	 
	if (norm)
	{
		Vx = nVx;
		Vy = nVy;
	}
	else
	{
		Vx = unVx;
		Vy = unVy;
	}
    
    /* size of dye square */
	
   	o_xMax = xmax;
   	o_yMax = ymax;
   	o_xMin = xmin;
   	o_yMin = ymin;
   	
	fprintf(stdout,"Vertices: %d\n Elements: %d\n Timesteps: %d\n\n",NV,NE,NT);
	fprintf(stdout,"Scalars: %d\n Boundaries: %d\n",numScalars,numBoundaries);
	fprintf(stdout,"Bounds: x (%f -> %f), y (%f -> %f)\n",xmin,xmax,ymin,ymax);
	fprintf(stdout,"Max velocity: %f\n",Vmax);

	return 0;
}

/* LoadAnnotations(): loads a user provided .annotations file. The file provides
 *					  information like colourmaps, titles, units etc.
 */ 
void LoadAnnotations(char* annotationsFile)
{
	FILE* fp;
	char c = 'a';
	int k=0;		
	int tmp;
	unsigned int i,j,t,s;	
	float* tmpcol;
 	register float val;
 	
	fp = fopen(annotationsFile,"r");
 	
	if (fp != NULL)
	{
 		title = (char*)malloc(255*sizeof(char));
 		tunits = (char*)malloc(100*sizeof(char));
		c = fgetc(fp);
		k=0;		
		while (c != '\n')
		{
			title[k] = c;
			c = fgetc(fp);
			k++;
		}
		title[k] = '\0';
 		
		c = fgetc(fp);
		k=0;						
		while (c != '\n')
		{
			tunits[k] = c;
			c = fgetc(fp);
			k++;
		}
		tunits[k] = '\0';
 		
 		fscanf(fp,"%f",&real_time_step);
		c = fgetc(fp);
		c = fgetc(fp);
		k=0;						
		while (c != '\n')
		{
			sunits[k] = c;
			c = fgetc(fp);
			k++;
		}
		sunits[k] = '\0';
 		
 		for (i=0;i<numScalars;i++)
 		{
			c = fgetc(fp);
			k=0;						
			while (c != '\n')
			{
				ScalarNames[i][k] = c;
				c = fgetc(fp);
				k++;
			}
			ScalarNames[i][k] = '\0'; 
						
			fscanf(fp,"%d",&tmp);
 			if (tmp > 1)
 			{
 				ControlPointNums[i] = tmp;
 				for (j=0;j<ControlPointNums[i];j++)
 				{
 					fscanf(fp,"%f %f %f %f %f",ColorMaps[i][j],ColorMaps[i][j]+1,ColorMaps[i][j]+2,ColorMaps[i][j]+3,ColorMaps[i][j]+4);
 				}
 			}
			c = fgetc(fp);
 		}
 	}
 	
 	/* compute color distances */
 	for (i=0;i<numScalars;i++)
 	{
 		inv_Diffs[i][0] = 0.0;
 		for (j=1;j<ControlPointNums[i];j++)
 		{
 			inv_Diffs[i][j] = 1.0/(ColorMaps[i][j][0] - ColorMaps[i][j-1][0]);
 		}
 	}
 	
 	
 	fprintf(stdout,"computing colors\n");
 	for (t=0;t<total_steps;t++)
 	{
 		for (j=0;j<total_Verts;j++)
 		{
 			for (s=0;s<numScalars;s++)
 			{
 				tmpcol = colors[j][t] + 4*s; 
 				val = Scalars[s][t][j];
 				InterpolateColor(tmpcol,val,ColorMaps[s],ControlPointNums[s],inv_Diffs[s]);	
 			}
 		}
 	}
 	fprintf(stdout,"finished\n");

}

/* ParseCommandline(): read commandline args
 */
int ParseCommandline(int argc, char** argv)
{
	unsigned int i;
	for (i=1; i<argc;i++)
	{
		if (!strcmp(argv[i],"-f"))
		{
			filename = argv[++i];
		}
		else if(!strcmp(argv[i],"-d"))
		{
			dt = atof(argv[++i]);
		}
		else if(!strcmp(argv[i],"-s"))
		{
			scale = atof(argv[++i]);
		}
		else if(!strcmp(argv[i],"-a"))
		{
			alpha = atof(argv[++i]);
		}
		else if(!strcmp(argv[i],"-sp"))
		{
			nspot = atof(argv[++i]);
			T = WSIZE/(nspot*NSIZE);
		}
		else if (!strcmp(argv[i],"-st"))
		{
			curtime = atof(argv[++i]);
		}
		else if(!strcmp(argv[i],"-t"))
		{
			title=argv[++i];
		}
		else if (!strcmp(argv[i],"-dt"))
		{
			real_time_step = atof(argv[++i]);
		}
		else if (!strcmp(argv[i],"-tu"))
		{
			tunits = argv[++i];
		}
		else if (!strcmp(argv[i],"-col"))
		{
			CurScalar = atoi(argv[++i]);
			colormap = 1;
		}
		else if (!strcmp(argv[i],"-dye"))
		{
			dye = 1;
		}
		else if (!strcmp(argv[i],"-norm"))
		{
			norm = 1;
		}
		else if (!strcmp(argv[i],"-an"))
		{
			annotationsFile = argv[++i];
		}
		else if (!strcmp(argv[i],"-v"))
		{
			fprintf(stdout,"VecVis Version: %2.2f\n",VERSION);
			return 1;
		}
		else if (!strcmp(argv[i],"-h"))
		{
			PrintHelp();
			return 1;
		}
		else if (!strcmp(argv[i],"-dim"))
		{
			WSIZE = atoi(argv[++i]);
			HSIZE = atoi(argv[++i]);
			T = WSIZE/(nspot*NSIZE);
			xh = ((float)(xmax - xmin))/((float)WSIZE);
    		yh = ((float)(ymax - ymin))/((float)HSIZE);
		}
		else if (!strcmp(argv[i],"-image"))
		{
			backgroundImagefile = argv[++i];
		}
		else
		{
			fprintf(stdout,"Unrecognised option %s",argv[i]);
		}           
	}

	if (dir == NULL)
	{
		dir = "./";
	}	

	return 0;
}

/* LoadImage(): Loads a backgroung image bitmap.
 *
 */
void LoadImage(char* fileName)
{
	unsigned char* imageData;
	unsigned int i;
	if (fileName != NULL)
	{
		BMPImage bmpImage;
		BMPFILE* bmpfile = CreateBMPFILE((char* )fileName);
		BMP_OpenBitMap(bmpfile,"rb");
		BMP_ReadHeaders(bmpfile);
		imageData = (unsigned char*)malloc(bmpfile->bmpInfoHeader->width*bmpfile->bmpInfoHeader->height*3*sizeof(unsigned char));
	
		BMP_ReadImageData(bmpfile);
		for (i = 0;i<(bmpfile->bmpInfoHeader->width*bmpfile->bmpInfoHeader->height);i++)
		{
			imageData[i*3] = bmpfile->imageData[i*3+2];
			imageData[i*3+1] = bmpfile->imageData[i*3+1];
			imageData[i*3+2] = bmpfile->imageData[i*3];
		}
	
	
		glBindTexture(GL_TEXTURE_2D,ImageTex[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);   
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,bmpfile->bmpInfoHeader->width,bmpfile->bmpInfoHeader->height,0,GL_RGB,GL_UNSIGNED_BYTE,imageData);
	}
}

/*==============================================================================
 *			OUTPUT FUNCTIONS
 *==============================================================================
 */
int OutputFrame()
{
	float r,g,b;
	char name[20];
	unsigned int i,j;
	BMPImage* image;
	GLfloat* pixels;
	
	image = (BMPImage*)malloc(sizeof(BMPImage));
	image->width = (unsigned long int)WSIZE;
	image->height = (unsigned long int)HSIZE;
	image->RGB = (unsigned char**)malloc(HSIZE*sizeof(unsigned char*));
	
	for (i=0;i<HSIZE;i++)
	{
		image->RGB[i] = (unsigned char*)malloc(WSIZE*3*sizeof(unsigned char));
	}
	pixels = (GLfloat*)malloc(WSIZE*HSIZE*3*sizeof(GLfloat));
	glReadPixels(0,0,WSIZE,HSIZE,GL_RGB,GL_FLOAT,(GLvoid*)pixels);
	
	fprintf(stdout,"WritingFrame %d\n",frame);
	for (j=0;j<HSIZE;j++)
	{
		for (i=0; i<WSIZE;i++)
		{
			r = pixels[j*WSIZE*3+i*3];
			g = pixels[j*WSIZE*3+i*3+1];
			b = pixels[j*WSIZE*3+i*3+2];
			
			image->RGB[j][i*3] = (unsigned char)(r*255.0);
			image->RGB[j][i*3+1] = (unsigned char)(g*255.0);
			image->RGB[j][i*3+2] = (unsigned char)(b*255.0);
		}
	}
	sprintf(name,"frame%016d.bmp",frame);
	printf("%s\n",name);
	WriteBMP(name,image);
	
	/*clean up*/
	for (i=0;i<HSIZE;i++)
	{
		free(image->RGB[i]);
	}
	free(image->RGB);
	free(pixels);
	free(image);
	return 0;
}

void PrintState()
{
	fprintf(stdout,"-----------------------------------------------------------\n");
	fprintf(stdout,"Dye on: %d\n",(int)dye);
	fprintf(stdout,"Selected Dye: %d\n",curDye);
	fprintf(stdout,"Number of Dye sources: %d\n",numDyeIn);
	fprintf(stdout,"Change in time: %f\n",dt);
	fprintf(stdout,"Dye Radii: %f\n",ds[curDye]);
	fprintf(stdout,"Alpha: %f\n",alpha);
	fprintf(stdout,"Color on: %d\n",(int)colormap);		
	fprintf(stdout,"Color overlay %d selected\n",CurScalar);
	fprintf(stdout,"Mesh on: %d\n",(int)drawmesh);
	fprintf(stdout,"Vector Scale: %f\n",scale);
	fprintf(stdout,"Jitter Grid on: %d\n",(int)jitter);
	fprintf(stdout,"Normalised: %d\n",(int)norm);
	fprintf(stdout,"BackGround Image: %d\n",(int)overlay);
	fprintf(stdout,"Dye Alpha: %f\n",dyeAlpha);
	fprintf(stdout,"Dye Style: %d\n",dyeStyle[curDye]);
	fprintf(stdout,"Dye Pulse Rate (pulse/sec): %f\n",60.0/((float)dyePulseRate[curDye]));
	fprintf(stdout,"Dye Rotation: %f\n",dyePhi[curDye]);
	fprintf(stdout,"-----------------------------------------------------------\n");
}

/*==============================================================================
 *			MAIN AND HELP FUNCTIONS
 *==============================================================================
 */

void PrintHelp()
{
	fprintf(stdout,"\nVecVis Commandline Option List:\n\n");
	fprintf(stdout," -f string\t\tIndicates the .flw file to load (required).\n");
	fprintf(stdout," -dim int int\t\tSet the window size (default 512 x 512).\n");
	fprintf(stdout," -d float\t\tSet the fraction of a timestep taken per frame (default = 0.05).\n");
	fprintf(stdout," -s float\t\tSet the scale factor to be applied to the field (default = 1).\n");
	fprintf(stdout," -a float\t\tSet the alpha value for the blending (default = 0.1).\n");
	fprintf(stdout," -sp float\t\tSet the spot size of the noise images (default = 1).\n");
	fprintf(stdout," -st float\t\tSet the start time for the animation (default = 0.0).\n");
	fprintf(stdout," -t string\t\tSet a Title to display.\n");
	fprintf(stdout," -dt float\t\tThe real world time interval for each timestep.\n");
	fprintf(stdout," -tu string\t\tThe real world units for the timestep.\n");
	fprintf(stdout," -col int\t\tEnable colormapping for the given scalar attribute.\n");
	fprintf(stdout," -dye\t\t\tEnable Dye injection (best in combination with -a 0.01).\n");
	fprintf(stdout," -norm\t\t\tNormalise the vector field.\n");
	fprintf(stdout," -an string\t\tIndicates a .ans file to load, this will override any relevant options.\n");
	fprintf(stdout," -v\t\t\tPrint version number.\n");
	fprintf(stdout," -image string\t\tSpecify Background image (only .bmp supported in version %0.02f).\n",VERSION);
	fprintf(stdout," -h\t\t\tPrint this help menu.\n");
	
}

void DisplayCommands()
{
	fprintf(stdout,"\n VecVis Commands:\n\n");
	
	fprintf(stdout,"\nDye Injection:\n");
	fprintf(stdout," 7\tSet mouse interaction to dye (default).\n");
	fprintf(stdout," d\tToggle dye injection.\n");
	fprintf(stdout," +\tAdd Dye Injection Point\n");
	fprintf(stdout," -\tMove to Next Dye Injection Point\n");
	fprintf(stdout," A Z\tincrease/decrease dye alpha value.\n");
	fprintf(stdout," g h \tIncrease, Decrease dye size.\n");
	fprintf(stdout," , .\trotate dye counter-clockwise/clockwise.\n");
	fprintf(stdout," { }\tcycle through dye styles.\n");
	fprintf(stdout," ( )\tdecrease/increase dye pulse rate.\n");
	
	fprintf(stdout,"\nObservation point Interaction:\n");
	fprintf(stdout," 9\tSet mouse interaction to observation points.\n");
	fprintf(stdout," o\tToggle Observation points.\n");
	fprintf(stdout," =\tAdd an Observation point.\n");
	fprintf(stdout," P\tDraw Plot.\n");
	fprintf(stdout," L\tToggle between profiles 1 and 2\n");
	
	fprintf(stdout,"\nGeneral Interaction:\n");
	fprintf(stdout," 8 2\t Pan in y.\n");
	fprintf(stdout," 6 4\t Pan in x.\n");
	fprintf(stdout," b n\tzoon in, out.\n");
	fprintf(stdout," t, y\tIncrease, Decrease timestep.\n");
	fprintf(stdout," s\tPause/continue.\n");
	
	fprintf(stdout,"\nIBFV Parameter Interaction:\n");
	fprintf(stdout," a z\tIncrease, Decrease alpha value.\n");
	fprintf(stdout," N\tToggle flow Normalisation\n");
	fprintf(stdout," [ ]\t Decrease/Increase vector scaling.\n");
	
	fprintf(stdout,"\nOverlays:\n");
	fprintf(stdout," c\tToggle colour maps.\n");
	fprintf(stdout," p, l\tCycle forward and backward through scalar attributes to show.\n");
	fprintf(stdout," m\tToggle mesh draw.\n");
	fprintf(stdout," j\tToggle jittered grid/\n");
	fprintf(stdout," I\tToggle background image.\n");
	
	fprintf(stdout,"\nOutput:\n");
	fprintf(stdout," f\tToggle FPS print out.\n");
	fprintf(stdout," O\tOutput Movie Frames.\n");
	fprintf(stdout," ?\tShow command list.\n");
}

void printGPL()
{
	fprintf(stdout,"\nVecVis Copyright (C) 2011  David J. Warne, Joe Young\n");
    fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY;\n");
    fprintf(stdout,"This is free software, and you are welcome to redistribute under certain conditions;\n\n");
}

/* main(): Progam starting point
 */
int main(int argc, char** argv)
{
	char heading[100]; 
	sprintf(heading,"VecVis Version: %0.2f (press ? for help)",VERSION);
	printGPL();
	
	
	
	/* parse the commandline, check for errors */
	if(ParseCommandline(argc, argv))
	{
		return 0;
	}
    T = WSIZE/(nspot*NSIZE);
    fprintf(stdout,"Loading %s ...\n",filename);
	/* load data    */
	if(LoadFlow())
	{
		exit(1);
	}
	BuildColorMaps();
	LoadAnnotations(annotationsFile);
	
	
	/* Initialise GLUT */
	fprintf(stdout,"Initialising GLUT...");
	glutInit(&argc, argv);
	/* Initialise the GL context and window */
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA  | GLUT_DEPTH); 
	glutInitWindowSize(WSIZE, HSIZE);
	glutCreateWindow(heading);
	InitGL();
	
	/*Load background image*/
	LoadImage(backgroundImagefile);
	/* texture Initialisation */
  	InitNoiseTextures(alpha);  
	InitJitterTexture(alpha);
	InitDyeTexture(dyeAlpha);
	
	/* initialise dye injection */
	InitDye();
	fprintf(stdout,"done\n");
	
	/* register glut callbacks */
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	/* main timer at 60 times per second*/
	glutTimerFunc(16.66667,Timer,1);
	glutKeyboardFunc(KeyPressed);
	glutMouseFunc(ButtonClick);
	glutTimerFunc(1000,Timer2,1);
	
	fprintf(stdout,"Starting...\n");
	/* now run the main loop */
	glutMainLoop();
	return 0;
}
 
