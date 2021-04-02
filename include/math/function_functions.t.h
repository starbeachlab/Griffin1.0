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


#ifndef FUNCTION_FUNCTIONS_T_H
#define FUNCTION_FUNCTIONS_T_H

#include "sum_function.t.h"
#include "product_function.t.h"

#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace math
{

    template< typename t_INPUT, typename t_RETURN>
    inline
    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > operator + ( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        return boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> >( new SumFunction< t_INPUT, t_RETURN>( FIRST_FUNCTION, SECOND_FUNCTION));
    }

    template< typename t_INPUT, typename t_RETURN>
    inline
    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > operator + ( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        return ( SUM_FUNCTION->operator += ( FIRST_FUNCTION));
    }

    template< typename t_INPUT, typename t_RETURN>
    inline
    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> >& operator + ( boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        return ( SUM_FUNCTION->operator += ( SECOND_FUNCTION));
    }

    template< typename t_INPUT, typename t_RETURN>
    boost::shared_ptr< Function< t_INPUT, t_RETURN> > operator += ( boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        if( FIRST_FUNCTION.use_count() > 0)
        {
            FIRST_FUNCTION =  boost::shared_ptr< Function< t_INPUT, t_RETURN> >( new SumFunction< t_INPUT, t_RETURN>( FIRST_FUNCTION, SECOND_FUNCTION));
        }
        else
        {
            FIRST_FUNCTION = SECOND_FUNCTION;
        }
        return FIRST_FUNCTION;
    }

//    template< typename t_INPUT, typename t_RETURN>
//    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > operator + ( const boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
//    {
//        return boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> >( new SumFunction< t_INPUT, t_RETURN>( FIRST_FUNCTION, SECOND_FUNCTION));
//    }
//
//    template< typename t_INPUT, typename t_RETURN>
//    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > operator + ( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
//    {
//        return ( SECOND_FUNCTION->operator += ( FIRST_FUNCTION));
//    }
//
//    template< typename t_INPUT, typename t_RETURN>
//    boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > operator + ( const boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
//    {
//        return ( SUM_FUNCTION->operator += ( SECOND_FUNCTION));
//    }
//

    template< typename t_INPUT, typename t_RETURN>
    boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN, float> > operator * ( const float &FACTOR, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION)
    {
        DebugWrite( __PRETTY_FUNCTION__);

#if defined DEBUG or defined STANDARD
        if( FUNCTION.use_count() == 0)
        {
            std::cout << "===> no object behind pointer!" << std::endl;
        }
        else
        {
            std::cout << "object: " << FUNCTION->GetClassName() << std::endl;
        }
#endif

        return boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN, float> >( new ProductFunction< t_INPUT, t_RETURN, float>( FACTOR, FUNCTION));
    }

    template< typename t_INPUT, typename t_RETURN>
    boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN, float> > operator * ( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION, const float &FACTOR)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        return boost::shared_ptr< ProductFunction< t_INPUT, t_RETURN, float> >( new ProductFunction< t_INPUT, t_RETURN, float>( FACTOR, FUNCTION));
    }

    template< typename t_INPUT, typename t_RETURN>
    boost::shared_ptr< Function< t_INPUT, t_RETURN> > operator *= ( boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION, const float &FACTOR)
    {
        DebugWrite( __PRETTY_FUNCTION__);
        FUNCTION =  boost::shared_ptr< Function< t_INPUT, t_RETURN> >( new ProductFunction< t_INPUT, t_RETURN, float>( FACTOR, FUNCTION));
        return FUNCTION;
    }


}

#endif
