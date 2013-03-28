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
/*File: FlwReaders.h
 *
 * Author: David Warne (david.warne@qut.edu.au)
 * 
 * Summary: Specification of Reader function for vector fields
 */
 
 #ifndef READERS_H
 #define READERS_H
 
 #include<stdio.h>
 #include<stdlib.h>
 #include<string.h>
 
 /* Header info reader */
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
                   fpos_t* boudarypos); 
                   
 /*Data Readers*/                  
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
                   fpos_t* boudarypos);
                   
 #endif
