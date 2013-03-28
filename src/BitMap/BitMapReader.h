
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
/* File: BitMapReader.h
 * Author: David Warne (david.warne@qut.edu.au)
 * Date Created: 18/04/2010
 *
 * Summary: Definitions of BitMap Reader Functions.
 *
 */


#ifndef BITMAPREADER_H
#define BITMAPREADER_H


#include "BitMapFile.h"

/*Core Reader Functions*/
ERROR BMP_ReadHeaders(BMPFILE*);
ERROR BMP_ReadColourPalette(BMPFILE*);
ERROR BMP_ReadImageData(BMPFILE*);

/*Higher order functions for reading, more user frendly*/
BMPImage ReadBMP(char*);

#endif
