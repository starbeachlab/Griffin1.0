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


#ifndef INTERACTION_GRID_POINT_H
#define INTERACTION_GRID_POINT_H

#include "../math/function.t.h"
#include "../math/vector3N.h"
#include "../phys/potential_force_container.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace store
{

    class InteractionGridPoint
      : public math::Function< mol::Atom, math::Vector3N >
    {
     protected:
         phys::PotentialForceContainer                     m_Forces;
         int                                               m_NSP;       //!< nearest surface point id

     public:

         InteractionGridPoint
         (
                 const phys::PotentialForceContainer &FORCES = phys::PotentialForceContainer(),
                 const int &NSP = std::numeric_limits< int>::max()
         )
         : m_Forces( FORCES),
           m_NSP( NSP)
         {}

         InteractionGridPoint
         (
                 const int &NSP 
         )
         : m_Forces(),
           m_NSP( NSP)
         {}

        InteractionGridPoint( const InteractionGridPoint & IGP)
          : m_Forces( IGP.m_Forces),
            m_NSP( IGP.m_NSP)
        {}

        virtual ~InteractionGridPoint(){}

        virtual InteractionGridPoint *Clone() const
        { return new InteractionGridPoint( *this);}

        virtual math::Vector3N & operator()( const mol::Atom &ATOM) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return math::Function< mol::Atom, math::Vector3N >::s_Tmp = m_Forces( ATOM);
        }

        virtual const int &GetNSP() const
        {
            return m_NSP;
        }

	virtual void SetNSP( const int &NSP)
	{
	    m_NSP = NSP;
	}


        virtual float Energy( const mol::Atom &ATOM) const
        {
        	return m_Forces.Energy( ATOM);
        }

//        virtual
//        InteractionGridPoint *
//        MolTypeFilter( const std::string &MOL_TYPE)
//        {
//	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
//	    return new InteractionGridPoint( *this);
//        }

        phys::PotentialForceContainer &
        GetSetPotentialForcesContainer()
        {
        	return m_Forces;
        }

        virtual std::istream &Read( std::istream &STREAM)
        {
        	DebugWrite( __PRETTY_FUNCTION__);
			std::string str;

			STREAM >> m_Forces >> str >> m_NSP;
			if( str != "nsp:")
			{
				std::cout << "===> should be nsp: <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
				std::cout << "nsp: " << m_NSP << std::endl;
				exit( -1);
			}
			return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
			STREAM << util::EnumHandler< util::FunctorEnum>().String( util::e_InteractionGridPoint) << " ";
			STREAM.flush();
			STREAM << m_Forces;
			STREAM << "nsp: " << m_NSP << std::endl;
			return STREAM;
        }


        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_InteractionGridPoint;
        }

    }; // end class InteractionGridPoint

} // end namespace store


#endif

