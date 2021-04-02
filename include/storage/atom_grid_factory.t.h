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


#ifndef ATOM_GRID_FACTORY_T_H_
#define ATOM_GRID_FACTORY_T_H_

#include "factory.h"
#include "../molecule/atom.h"
#include "void_gridpoint.t.h"


//namespace store
//{
//
//	namespace factory
//	{
//
//
//		template< >
//		inline
//		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
//		BuildNeutralObject< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >()
//		{
//			DebugWrite( __PRETTY_FUNCTION__);
//		    return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
//		    (
//		    		new VoidGridPoint< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >()
//		    );
//		}
//
//	} // end namespace factory
//}// end namespace store


namespace store
{
    namespace factory
    {

		template<>
		inline
		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        BuildNeutralObject< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >()
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
            (
                    new VoidGridPoint< mol::Atom, math::Vector3N >()
            );
        }
    } // end namespace factory
} // end namespace store




#endif /* ATOM_GRID_FACTORY_T_H_ */
