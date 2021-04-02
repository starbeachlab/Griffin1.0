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


#ifndef INTEGRATOR_T_H_
#define INTEGRATOR_T_H_

#include <cstdio>

#include "function.t.h"
#include "iteration_vector.t.h"

namespace math
{
    template< typename t_DATA, typename t_RETURN>
    class Integrator
    : public Function< Function< t_DATA, t_RETURN>, t_RETURN>
    {
    protected:
        boost::shared_ptr< IterationVector< size_t>    >                              m_Iterate;  //!< the n-dim loop
        boost::shared_ptr< Function< std::vector< size_t>, t_DATA> >           m_Transformation;  //!< coordinate system in which to integrate in

    public:
        Integrator()
        : m_Iterate(), m_Transformation()
        { std::cerr << "calling default constructor makes no sense"; exit( -1);}

        Integrator( const boost::shared_ptr< IterationVector< size_t> > &ITER, const boost::shared_ptr< Function< std::vector< size_t>, t_DATA> > &TRANS)
        : m_Iterate( ITER),
        m_Transformation( TRANS)
        {}

        Integrator( const Integrator &INTEGRATOR)
        : m_Iterate( INTEGRATOR.m_Iterate),
        m_Transformation( INTEGRATOR.m_Transformation)
        {}

        virtual ~Integrator(){}

        virtual Integrator *Clone() const{ return new Integrator( *this);}

        virtual t_RETURN operator()( const Function< t_DATA, t_RETURN> &FUNCTION) const
        {
            t_RETURN result;
            while( m_Iterate->Iterate())
            {
                result += FUNCTION( ( *m_Transformation)( m_Iterate->GetValues()));
            }
            return result;
        }


        virtual t_RETURN operator()() const
        {
            t_RETURN result;
            while( m_Iterate->Iterate())
            {
                result += ( *m_Transformation)( m_Iterate->GetValues());
            }
            return result;
        }

    }; // end class Integrator
} // end namespace math

#endif /* INTEGRATOR_T_H_ */
