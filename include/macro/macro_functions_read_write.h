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


#ifndef MACRO_FUNCTIONS_READ_WRITE_H_
#define MACRO_FUNCTIONS_READ_WRITE_H_

// these macro definitions allow to control output according to DEBUG or STANDARD message level
// this implementation allows to
// switching to another message level requires recompilation
// ===>  DO NOT write stuff like:  StandardWrite( arg++)  it can cause a float increment!  <===

#ifdef DEBUG
#define DebugWriteNoFlush( ARGUMENT) std::cout << ARGUMENT << " "
#define DebugWrite( ARGUMENT) std::cout << ARGUMENT << std::endl
#define DebugPrint( STREAM, ARGUMENT) STREAM << ARGUMENT << std::endl;
#else
#define DebugWriteNoFlush( ARGUMENT)
#define DebugWrite( ARGUMENT)
#define DebugPrint( STREAM, ARGUMENT)
#endif

#if defined STANDARD or defined DEBUG
#define StandardWriteNoFlush(ARGUMENT) std::cout << ARGUMENT << " "
#define StandardWrite(ARGUMENT) std::cout << ARGUMENT << std::endl
#define StandardPrint( STREAM, ARGUMENT) STREAM << ARGUMENT << std::endl;
#else
#define StandardWriteNoFlush( ARGUMENT)
#define StandardWrite( ARGUMENT)
#define StandardPrint( STREAM, ARGUMENT)
#endif


#endif /* MACRO_FUNCTIONS_READ_WRITE_H_ */
