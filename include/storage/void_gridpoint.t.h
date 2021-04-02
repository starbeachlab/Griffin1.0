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


#ifndef VOID_GRIDPOINT_H
#define VOID_GRIDPOINT_H

#include "../math/function.t.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"
#include "factory.h"



namespace store
{


    template< typename t_INPUT, typename t_RETURN>
    class VoidGridPoint
    : public math::Function< t_INPUT, t_RETURN>
    {
    private:
		static const t_RETURN s_Zero;
    public:
        VoidGridPoint()
        : math::Function< t_INPUT, t_RETURN>()
        {}

        VoidGridPoint( const VoidGridPoint & GPOINT)
        : math::Function< t_INPUT, t_RETURN>( GPOINT)
        {}

        virtual ~VoidGridPoint(){}

        virtual VoidGridPoint *Clone() const
        { return new VoidGridPoint( *this);}


        virtual t_RETURN & operator()( const t_INPUT &DATA) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return (t_RETURN &) s_Zero;
        }


        virtual std::ostream &Write( std::ostream &STREAM) const
        {
	    STREAM << util::EnumHandler< util::FunctorEnum>().String( util::e_VoidGridPoint) << " ";
            return STREAM;
        }

        virtual float Energy( const t_INPUT &DATA) const
        {
        	return 0.0;
        }

        virtual
        VoidGridPoint *
        MolTypeFilter( const std::string &MOL_TYPE)
        {
	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
	    return new VoidGridPoint( *this);
        }

        virtual std::istream &Read( std::istream &STREAM)
        {
        	DebugWrite( __PRETTY_FUNCTION__ );
            return STREAM;
        }

        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_VoidGridPoint;
        }


    }; // end class VoidGridPoint


    template< typename t_INPUT, typename t_RETURN>
    const t_RETURN VoidGridPoint< t_INPUT, t_RETURN>::s_Zero = factory::BuildNeutralObject< t_RETURN>();

} // end namespace store



#endif
