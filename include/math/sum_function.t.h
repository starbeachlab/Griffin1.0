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


#ifndef SUM_FUNCTION_H
#define SUM_FUNCTION_H

#include "function.t.h"
#include "../storage/shared_pointer_vector.t.h"

namespace math
{
    template< typename t_INPUT, typename t_RETURN>
    class SumFunction
      : public Function< t_INPUT, t_RETURN>
    {
    protected:
        store::ShPtrVec< Function< t_INPUT, t_RETURN> >  m_Function;

    public:
        SumFunction( )
        : m_Function()
        {}

        SumFunction( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FIRST_FUNCTION, const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &SECOND_FUNCTION)
        : m_Function(2)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            m_Function( 0) = FIRST_FUNCTION;
            m_Function( 1) = SECOND_FUNCTION;
        }


        SumFunction( const SumFunction &FCT)
        : m_Function( FCT.m_Function)
        {}

        virtual SumFunction *Clone() const{ return new SumFunction( *this);}

        boost::shared_ptr< SumFunction< t_INPUT, t_RETURN> > &
        operator +=( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            m_Function.insert( m_Function.end(), FUNCTION);
            return *this;
        }


        SumFunction< t_INPUT, t_RETURN> &
        operator +=( const SumFunction< t_INPUT, t_RETURN> &FUNCTION)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            m_Function.insert( m_Function.end(), FUNCTION.begin(), FUNCTION.end());
            return *this;
        }

        SumFunction< t_INPUT, t_RETURN> &
        operator = ( const boost::shared_ptr< Function< t_INPUT, t_RETURN> > &FUNCTION)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            m_Function.clear();
            m_Function.push_back( FUNCTION);
            return *this;
        }

        SumFunction< t_INPUT, t_RETURN> &
        operator = ( const SumFunction< t_INPUT, t_RETURN> & SUM_FUNCTION)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            m_Function = SUM_FUNCTION.m_Function;
            return *this;
        }

        t_RETURN
        operator()( const t_INPUT &INPUT) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << std::endl << "sum function size: " << m_Function.size());
            t_RETURN result;
            for( typename std::vector< boost::shared_ptr< Function< t_INPUT, t_RETURN> > >::const_iterator itr = m_Function.begin(); itr != m_Function.end(); ++itr)
            {
#ifdef DEBUG
            if( itr->use_count() == 0)
            {
                std::cout << "===> no object behind pointer!" << std::endl;
            }
            else
            {
                std::cout << "object behind pointer in sum function: " << ( *itr)->GetClassName() << std::endl;
            }
#endif
                result += ( *itr)->operator()( INPUT);
            }
            return result;
        }



        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            std::cout << GetClassName() << std::endl;
            for( typename std::vector< boost::shared_ptr< Function< t_INPUT, t_RETURN> > >::const_iterator itr = m_Function.begin(); itr != m_Function.end(); ++itr)
            {
                ( *itr)->Write( STREAM);
            }
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }




    }; // end class SumFunction
} // end namespace math

#endif /*SUM_FUNCTION_H*/
