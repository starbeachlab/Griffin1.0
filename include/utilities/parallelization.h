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


#ifndef UTIL_PARALLIZATION_H
#define UTIL_PARALLIZATION_H

#include <cassert>

#include "../math/double_functions.h"
#include "../external/std_functions.h"
#include "../math/vector.h"
#include "../math/vector3N.h"
#include "../math/double_functions.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


namespace util
{

  bool ParallelizeProcess( math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA, const size_t &NR_PROCESSES, const size_t &PROCESS_ID);

  std::vector< std::vector< size_t> > PossibleDevisors( math::Vector &MIN, math::Vector &MAX, const math::Vector &DELTA);

  std::vector< std::vector< size_t> > MatchingDevisors( const std::vector< std::vector< size_t> > &SUB_BLOCKS, const size_t &NR_PROCESSES);

  std::vector< size_t> SelectMostBalancedDevisors( const std::vector< std::vector< size_t> > &SUB_BLOCKS);

  void AdjustLimits( math::Vector &MIN, math::Vector &MAX, const std::vector< size_t> &SUB_BLOCK, const size_t &PROCESS_ID);

  std::vector< size_t> NrToIDs( const std::vector< size_t> &LIMITS, const size_t &NR);

} // namespace util


#endif
