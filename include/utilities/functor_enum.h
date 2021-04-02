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


#ifndef FUNCTOR_ENUM_H
#define FUNCTOR_ENUM_H

namespace util
{
	enum FunctorEnum
	{
		e_Function,
		e_VoidGridPoint,
		e_SurfGridPoint,
		e_ConstantGridPoint,
		e_InteractionGridPoint,
		e_TypeMappedGridPoint,
		e_RecursiveTypeMappedGridPoint,
		e_SumGridPoint,
		e_TransientInteractionGridPoint,
		e_PotentialElectrostaticForce,
		e_PotentialAttractiveVanDerWaalsForce,
		e_PotentialRepulsiveVanDerWaalsForce,
		e_UNDEFINED_Functor = 88888888
	};
} // end namespace util

#endif

