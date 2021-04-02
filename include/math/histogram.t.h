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


#ifndef HISTOGRAM_T_H_
#define HISTOGRAM_T_H_

//#include "../readwrite/stream_operator.h"
//#include "../string/io_string_functions.h"

#include "function.t.h"
#include <iterator>

namespace math
{
    template< typename t_INPUT>
    class Histogram
    : public Function< t_INPUT, size_t>
    {
    protected:
        ///////////////
        //    std::vector< size_t>     //
        ///////////////
        t_INPUT              m_Min;
        t_INPUT                 m_Delta;
        size_t                 m_NrBins;
        std::vector< size_t> m_Data;
        size_t               m_BelowLowestBin;
        size_t                 m_AboveHighestBin;
        t_INPUT                 m_InverseDelta;
	std::vector< t_INPUT>  m_Above;
	std::vector< t_INPUT>  m_Below;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! construct from std::vector< size_t>
        Histogram( const t_INPUT &MIN, const t_INPUT &DELTA, const size_t &NR_BINS)
        : Function< t_INPUT, size_t>(),
        m_Min( MIN),
        m_Delta( DELTA),
        m_NrBins( NR_BINS),
        m_Data( NR_BINS),
        m_BelowLowestBin( 0),
        m_AboveHighestBin( 0),
	    m_InverseDelta( 1.0 / DELTA),
	    m_Above(),
	    m_Below()
        {}

        //! copy constructor
        Histogram( const Histogram &ORIGINAL)
        : Function< t_INPUT, size_t>( ORIGINAL),
        m_Min( ORIGINAL.m_Min),
        m_Delta( ORIGINAL.m_Delta),
        m_NrBins( ORIGINAL.m_NrBins),
        m_Data( ORIGINAL.m_Data),
        m_BelowLowestBin( ORIGINAL.m_BelowLowestBin),
        m_AboveHighestBin( ORIGINAL.m_AboveHighestBin),
	    m_InverseDelta( ORIGINAL.m_InverseDelta),
	    m_Above( ORIGINAL.m_Above),
	    m_Below( ORIGINAL.m_Below)
        {}

        //! virtual destructor
        virtual ~Histogram(){}

        //! virtual copy constructor
        virtual Histogram *Clone() const{ return new Histogram( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        size_t & operator() ( const t_INPUT &INPUT) const
        {
            size_t id( ValueToBin( INPUT));
            assert( id != std::numeric_limits< size_t>::max() && id != -std::numeric_limits< size_t>::max());
            return (size_t &) m_Data[ id];
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

        const std::vector< size_t> &GetData() const
        {
            return m_Data;
        }

        const t_INPUT &GetInverseDelta() const
        {
            return m_InverseDelta;
        }

        const size_t &GetAboveHighestBin() const
        {
            return m_AboveHighestBin;
        }

        const size_t &BelowLowestBin() const
        {
            return m_BelowLowestBin;
        }


        void InsertValue( const t_INPUT &VALUE)
        {
            int id( ValueToBin( VALUE));
	    DebugWrite(__FUNCTION__ <<   "id: " << id);
            if( id == -std::numeric_limits< int>::max())
            {
		DebugWrite(  " belowcount++");
                ++m_BelowLowestBin;
		m_Below.push_back( VALUE);
            }
            else if( id == std::numeric_limits< int>::max())
            {
                ++m_AboveHighestBin;
		DebugWrite(  " abovecount++");
		m_Above.push_back( VALUE);
            }
            else
            {
                ++m_Data[ id];
		DebugWrite(  " insidecount = " << m_Data[id] );
            }
        }

        void InsertValue( const t_INPUT &VALUE, const int &INCREMENT)
        {
            int id( ValueToBin( VALUE));
	    DebugWrite( __FUNCTION__ <<  "id: " << id);
            if( id == -std::numeric_limits< int>::max())
            {
		DebugWrite(  " belowcount++" );
                m_BelowLowestBin += INCREMENT;
		m_Below.push_back( VALUE);
            }
            else if( id == std::numeric_limits< int>::max())
            {
		DebugWrite(  " abovecount++" );
                m_AboveHighestBin += INCREMENT;
		m_Above.push_back( VALUE);
            }
            else
            {
                m_Data[ id] += INCREMENT;
		DebugWrite(  " insidecount = " << m_Data[id] );
            }
        }

        int ValueToBin( const t_INPUT &VALUE) const
        {
	    DebugWrite(  __FUNCTION__ << "(" << VALUE << ") min: " << m_Min);
            if( VALUE < m_Min)
            {
		DebugWrite(  " below min! ");
                return -std::numeric_limits< int>::max();
            }
            else if( VALUE > m_Min + m_NrBins * m_Delta)
            {
		DebugWrite(  " above max (" << m_Min + m_NrBins * m_Delta << ") ");
                return std::numeric_limits< int>::max();
            }
	    std::cout.flush();
            return int( ( VALUE - m_Min) * m_InverseDelta);
        }

        size_t Sum() const
        {
            size_t sum( 0);
            for( std::vector< size_t>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                sum += *itr;
            }
            return sum;
        }

        void Clear()
        {
        	std::fill( m_Data.begin(), m_Data.end(), 0);
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
            for( std::vector< size_t>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM << lower + 0.5 * m_Delta << "  " << *itr << std::endl;
                lower += m_Delta;
                upper += m_Delta;
            }
            STREAM << "#sum of values found above hightest bin: " << m_AboveHighestBin << std::endl << "# ";
	    std::copy( m_Above.begin(), m_Above.end(), std::ostream_iterator< t_INPUT>( STREAM, " "));
	    STREAM << std::endl << "# ";
	    std::copy( m_Below.begin(), m_Below.end(), std::ostream_iterator< t_INPUT>( STREAM, " "));
	    STREAM << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Histogram
} // end namespace math




#endif /* HISTOGRAM_T_H_ */
