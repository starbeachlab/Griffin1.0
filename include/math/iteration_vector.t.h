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


#ifndef ITERATION_VECTOR_H_
#define ITERATION_VECTOR_H_

#include <vector>

namespace math
{
    template< typename t_TYPE>
    class IterationVector
    {
    protected:
        std::vector< t_TYPE>       m_Min;
        std::vector< t_TYPE>       m_Max;
        std::vector< t_TYPE>       m_Delta;
        std::vector< t_TYPE>       m_CurrentPoint;
//        size_t                     m_NrShells;
        size_t                       m_NrIncrements;

    public:
        IterationVector( const std::vector< t_TYPE> &MIN, const std::vector< t_TYPE> &MAX, const std::vector< t_TYPE> &DELTA)
        : m_Min( MIN),
        m_Max( MAX),
        m_Delta( DELTA),
        m_CurrentPoint( MIN),
        m_NrIncrements( 0)
        {}


        IterationVector( const IterationVector &ITER)
        : m_Min( ITER.m_Min),
        m_Max( ITER.m_Max),
        m_Delta( ITER.m_Delta),
        m_CurrentPoint( ITER.m_CurrentPoint),
        m_NrIncrements( ITER.m_NrIncrements)
        {}


        virtual ~IterationVector(){}


        virtual IterationVector *Clone() const
        { return new IterationVector( *this);}


        const std::vector< t_TYPE> &GetValues() const
        { return m_CurrentPoint;}


        size_t GetNrIncrements() const
        { return m_NrIncrements;}

        bool Iterate(){ return Increment();}

        bool Increment()
        {
            typename std::vector< t_TYPE>::const_reverse_iterator
				min_itr( m_Min.rbegin()),
				max_itr( m_Max.rbegin()),
				delta_itr( m_Delta.rbegin());
            typename std::vector< t_TYPE>::reverse_iterator
				current_itr( m_CurrentPoint.rbegin());

            while( current_itr != m_CurrentPoint.rend())
            {
                if( *current_itr + *delta_itr < *max_itr)
                {
                    *current_itr += *delta_itr;
                    ++m_NrIncrements;
                    return true;
                }
                else
                {
                    *current_itr = *min_itr;
                    ++current_itr;
                    ++max_itr;
                    ++min_itr;
                }
            }
            return false;
        } // end Increment()

    }; // end class IterationVector

} // end namespace math

#endif /* ITERATION_VECTOR_H_ */
