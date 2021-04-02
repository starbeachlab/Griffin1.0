/////////////////////////////////////////////////////////////////////////
//!									 
//  This file is part of the GRIFFIN program				 
//									 
//  Copyright (C) 2010 by Rene Staritzbichler		      	 
//  rene.staritzbichler@biophys.mpg.de			       	 
//									 
//  GRIFFIN is free software; you can redistribute it and/or modify	 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.				 
//									 
//  GRIFFIN is distributed in the hope that it will be useful,	 
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	 
//  GNU General Public License for more details.			 
//!									 
//!									 
//!									 
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace math
{

const float Pi( acos( -1));
const float HalfPi( .5 * acos( -1));
const float OneAndAHalfPi( 1.5 * acos( -1));
const float TwoPi( 2.0 * acos( -1));
const float GradToRad( acos( -1) / 180.0);
const float RadToGrad( 180.0 / acos( -1));

} // end namespace $




#endif /* CONSTANTS_H_ */
