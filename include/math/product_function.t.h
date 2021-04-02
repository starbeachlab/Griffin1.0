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


#ifndef ProductFunction_T_H_
#define ProductFunction_T_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"

#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace math
{
    template< typename t_INPUT, typename t_RETURN, typename t_SCALAR>
    class ProductFunction
    : public Function< t_INPUT, t_RETURN>
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        t_SCALAR m_Scalar;
        boost::shared_ptr< Function< t_INPUT, t_RETURN> >  m_Function;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ProductFunction()
        : Function< t_INPUT, t_RETURN>(),
        m_Scalar(),
        m_Function()
        {}

        //! construct from data
        ProductFunction( const t_SCALAR &DATA, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION)
        : Function< t_INPUT, t_RETURN>(),
        m_Scalar( DATA),
        m_Function( FUNCTION)
        {}

        //! copy constructor
        ProductFunction( const ProductFunction &ORIGINAL)
        : Function< t_INPUT, t_RETURN>( ORIGINAL),
        m_Scalar( ORIGINAL.m_Scalar),
        m_Function( ORIGINAL.m_Function)
        {}

        //! virtual destructor
        virtual ~ProductFunction(){}

        //! virtual copy constructor
        virtual ProductFunction *Clone() const{ return new ProductFunction( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        t_RETURN operator()( const t_INPUT &INPUT) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << ": " << m_Scalar << ": " << INPUT << ": " << m_Function->GetClassName());
            return m_Scalar * m_Function->operator()( INPUT);
        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            STREAM << m_Scalar << "  " << m_Function << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class ProductFunction
} // end namespace math




#endif /* ProductFunction_T_H_ */
