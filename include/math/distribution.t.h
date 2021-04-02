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


#ifndef DISTRIBUTION_T_H_
#define DISTRIBUTION_T_H_

//#include "../readwrite/stream_operator.h"
//#include "../string/io_string_functions.h"

#include "function.t.h"

namespace math
{
    template< typename t_INPUT>
    class Distribution
    : public Function< t_INPUT, t_INPUT>
    {
    protected:
        ///////////////
        //    std::vector< size_t>     //
        ///////////////
        t_INPUT               m_Min;
        t_INPUT                  m_Delta;
        size_t                  m_NrBins;
        std::vector< t_INPUT> m_Data;
        t_INPUT                m_BelowLowestBin;
        t_INPUT                  m_AboveHighestBin;
        t_INPUT                  m_InverseDelta;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! construct from std::vector< size_t>
        Distribution( const t_INPUT &MIN, const t_INPUT &DELTA, const size_t &NR_BINS)
        : Function< t_INPUT, t_INPUT>(),
        m_Min( MIN),
        m_Delta( DELTA),
        m_NrBins( NR_BINS),
        m_Data( NR_BINS),
        m_BelowLowestBin( 0),
        m_AboveHighestBin( 0),
        m_InverseDelta( 1.0 / DELTA)
        {}

        //! copy constructor
        Distribution( const Distribution &ORIGINAL)
        : Function< t_INPUT, t_INPUT>( ORIGINAL),
        m_Min( ORIGINAL.m_Min),
        m_Delta( ORIGINAL.m_Delta),
        m_NrBins( ORIGINAL.m_NrBins),
        m_Data( ORIGINAL.m_Data),
        m_BelowLowestBin( ORIGINAL.m_BelowLowestBin),
        m_AboveHighestBin( ORIGINAL.m_AboveHighestBin),
        m_InverseDelta( ORIGINAL.m_InverseDelta)
        {}

        //! virtual destructor
        virtual ~Distribution(){}

        //! virtual copy constructor
        virtual Distribution *Clone() const{ return new Distribution( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        t_INPUT & operator() ( const t_INPUT &INPUT) const
        {
            size_t id( ValueToBin( INPUT));
            assert( id != std::numeric_limits< size_t>::max() && id != -std::numeric_limits< size_t>::max());
            return (t_INPUT &) m_Data[ id];
        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const t_INPUT &GetMinimum() const
        {
            return m_Min;
        }

        const t_INPUT &GetDelta() const
        {
            return m_Delta;
        }

        const size_t &GetNrBins() const
        {
            return m_NrBins;
        }

        const std::vector< t_INPUT> &GetData() const
        {
            return m_Data;
        }

        const t_INPUT &GetInverseDelta() const
        {
            return m_InverseDelta;
        }

        const t_INPUT &GetAboveHighestBin() const
        {
            return m_AboveHighestBin;
        }

        const t_INPUT &BelowLowestBin() const
        {
            return m_BelowLowestBin;
        }


        void InsertValue( const t_INPUT &POSITION, const t_INPUT &VALUE)
        {
            size_t id( ValueToBin( POSITION));
            if( id == -std::numeric_limits< size_t>::max())
            {
                m_BelowLowestBin += VALUE;
            }
            else if( id == std::numeric_limits< size_t>::max())
            {
                m_AboveHighestBin += VALUE;
            }
            else
            {
                m_Data[ id] += VALUE;
            }
        }


        size_t ValueToBin( const t_INPUT &VALUE) const
        {
            if( VALUE < m_Min)
            {
                return -std::numeric_limits< size_t>::max();
            }
            else if( VALUE > m_Min + m_NrBins * m_Delta)
            {
                return std::numeric_limits< size_t>::max();
            }
            return size_t( ( VALUE - m_Min) * m_InverseDelta);
        }

        t_INPUT Sum() const
        {
            t_INPUT sum( 0);
            for( typename std::vector< t_INPUT>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                sum += *itr;
            }
            return sum;
        }

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
//            std::cout << "events_below_bins: " << m_BelowLowestBin << std::endl;
            t_INPUT lower( m_Min), upper( m_Min + m_Delta);
            STREAM << "#sum of values found below lowest bin: " << m_BelowLowestBin << std::endl;
            for( typename std::vector< t_INPUT>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM << lower + 0.5 * m_Delta << "  " << *itr << std::endl;
                lower += m_Delta;
                upper += m_Delta;
            }
            STREAM << "#sum of values found above hightest bin: " << m_AboveHighestBin << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Distribution
} // end namespace math




#endif /* Distribution_T_H_ */
