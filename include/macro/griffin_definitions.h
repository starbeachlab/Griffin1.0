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


#ifndef MACRO_GRIFFIN_H_
#define MACRO_GRIFFIN_H_


#define STANDARD 		            // enables standard output
// #define DEBUG                    // loads of output

// #define MONITORING               // enable memory profiling
// #define TIMERS                   // enable time profiling of major functions (for a complete profile use -pg flag for gprof)

// #define VMD_OUTPUT               // vmd output for visualization of forces

// #define FORCE_COORDINATES        // print force files with coordinates (molecule_force_grid.cc) 
// #define FORCE_INDICES            // print force files with indices // DON'T USE WITH MPI_PARALLEL or frequency flag, check for interpolation!
// #define SPLIT_FORCES		        // OUTDATED print individual force contributions (makes most sense with #define FORCE_COORDINATES)or frequency flag
// #define NO_INTERACTION_FORCES    // ignore interaction forces (vdw+coulomb)
// #define NO_VDW		            // ignore vdw forces // still there?
// #define NO_COULOMB		        // ignore coulomb forces // still there?
// #define ANALYZE_BURIAL



////  DO NOT EDIT !!  ////
#ifdef MPI_PARALLEL
#define FULL_PARALLEL               
#else
#define DAEMON
#endif


#include "griffin_global.h"

#endif /* MACRO_GRIFFIN_H_ */
