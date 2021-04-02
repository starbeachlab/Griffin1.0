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


#ifndef FUNCTION_H
#define FUNCTION_H

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../utilities/functor_enum.h"
#include "../utilities/enum_handler.t.h"

namespace math
{



    template< typename t_INPUT, typename t_RETURN>
    class Function
    : public StreamOperator
    {
    protected:
    	static t_RETURN   			   s_Tmp;
    public:
        virtual t_RETURN & operator()( const t_INPUT &DATA) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return s_Tmp = t_RETURN();
        }

        virtual t_INPUT dF() const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return t_INPUT();
        }

        virtual const int &GetNSP() const /// TODO: derive grid point parent and keep function class clean
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return g_IntMax;
        }

        virtual float Energy( const t_INPUT &DATA) const  // TODO: take out of function class !
        {
            DebugWrite( __PRETTY_FUNCTION__);
        	return g_FloatMax;
        }

//        virtual
//        Function *
//        MolTypeFilter( const std::string &MOL_TYPE)
//        {
//	    StandardWrite( "====> WHY?? " <<__PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
//	    return new Function( *this);
//        }

        virtual Function *Clone() const
        {
            return new Function( *this);
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            return STREAM;
        }

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_Function;
        }


    }; // end class Function

    template< typename t_INPUT, typename t_RETURN>
    t_RETURN Function< t_INPUT, t_RETURN>::s_Tmp = t_RETURN();

} // end namespace math

#endif
