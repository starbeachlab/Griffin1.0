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


#ifndef POSITION_GRID_H
#define POSITION_GRID_H

#include "../math/vector3N.h"
#include "../math/function.t.h"
#include "vector3N.t.h"

namespace store
{
    template< typename t_RETURN>
    class PositionGrid
    {
    protected:
        Vector3N< size_t>           m_NrBins;
        math::Vector3N              m_Min;
        math::Vector3N              m_Max;
        math::Vector3N              m_Delta;

    public:
        PositionGrid()
        : m_NrBins(),
        m_Min(),
        m_Max(),
        m_Delta()
        {
            DebugWrite( "PositionGrid default constructor");
        }

        PositionGrid( const math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA, const store::Vector3N< size_t> &NR_BINS)
        : m_NrBins( NR_BINS),
        m_Min( MIN),
        m_Max( MAX),
        m_Delta( DELTA)
        {
            DebugWrite( "PositionGrid initialized min: " << m_Min << " max: " << m_Max << " delta: " << m_Delta << " nrbins: " << m_NrBins
                << "\nMIN: " << MIN << " MAX: " << MAX << " DELTA: " << DELTA << " nrbins: " << NR_BINS);
        }

        PositionGrid( const PositionGrid &POSITION_GRID)
        : m_NrBins( POSITION_GRID.m_NrBins),
        m_Min( POSITION_GRID.m_Min),
        m_Max( POSITION_GRID.m_Max),
        m_Delta( POSITION_GRID.m_Delta)
        {
        	std::cout << __FUNCTION__ << " copy constructor!" << std::endl;
        	DebugWrite(  "PositionGrid min: " << m_Min << " max: " << m_Max << " delta: " << m_Delta);
        }

        virtual ~PositionGrid(){}

        virtual PositionGrid *Clone() const = 0;

        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION) const
        {
        	std::cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " should not be called! needs to be overwritten in derived classes!" << std::endl;
        	exit( -1);
        	return t_RETURN();
        }

        virtual t_RETURN GetGridPoint( const store::Vector3N< int> &INDICES) const
        {
        	std::cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " should not be called! needs to be overwritten in derived classes!" << std::endl;
        	exit( -1);
        	return t_RETURN();
        }

        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE) const = 0;

#ifdef FORCE_INDICES
        virtual
        t_RETURN
        GetGridPoint
	    ( 
		const math::Vector3N &POSITION, 
		const std::string &MOL_TYPE, 
#ifdef CPPIO
	std::ostream &STREAM
#else
	FILE *STREAM
#endif
	    ) const = 0;
#endif

        virtual t_RETURN GetGridPoint( const store::Vector3N< int> &INDICES, const std::string &MOL_TYPE) const = 0;

        virtual
        util::FunctorEnum
        GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE) = 0;

        virtual const store::Vector3N< size_t> &GetNrBins() const{ return m_NrBins;}

        virtual int TotalNrBins() const
        {
        	return (int) m_NrBins[0] * m_NrBins[1] * m_NrBins[2];
        }


        virtual math::Vector3N PositionFromIDs( const store::Vector3N< size_t> &IDS) const
        {
            assert( IDS.size() == 3);
            math::Vector3N pos;
            std::vector< float>::iterator pos_itr( pos.begin());
            std::vector< float>::const_iterator min_itr( m_Min.begin()), delta_itr( m_Delta.begin());
            for( std::vector< size_t>::const_iterator itr( IDS.begin()); itr != IDS.end(); ++itr, ++pos_itr, ++min_itr, ++delta_itr)
            {
                *pos_itr = *min_itr + ( *itr + 0.5) * ( *delta_itr);
            }
            return pos;
        }

        // may return unreasonable values, check when applying!!
        virtual store::Vector3N< int> IDsFromPosition( const math::Vector3N &POS) const
        {
        	DebugWrite( __PRETTY_FUNCTION__);
            store::Vector3N< int> ids;
	    //	    std::cout << __FUNCTION__ << std::endl;
//        assert( POS( 0) >= m_Min( 0) && POS( 0) <= m_Max( 0) &&
//        POS( 1) >= m_Min( 1) && POS( 1) <= m_Max( 1) &&
//        POS( 2) >= m_Max( 2) && POS( 2) <= m_Max( 2));

            ids[ 0] = int( ( POS[ 0] - m_Min[ 0]) / m_Delta[ 0]); // time compare with bracket operator
            ids[ 1] = int( ( POS[ 1] - m_Min[ 1]) / m_Delta[ 1]);
            ids[ 2] = int( ( POS[ 2] - m_Min[ 2]) / m_Delta[ 2]);

	    //	    std::cout << __FUNCTION__ << std::endl;
            return ids;
        }

        virtual int IDFromIDs( const store::Vector3N< int> &IDS) const
        {
			if
			(
					IDS[0] < 0
					|| IDS[1] < 0
					|| IDS[2] < 0
					|| IDS[0] >= (int) m_NrBins[0]
					|| IDS[1] >= (int) m_NrBins[1]
					|| IDS[2] >= (int) m_NrBins[2]
			)
			{
//				std::cout << "===> out of limits in IDFromIDs " << std::endl;
//				std::cout << IDS[0] << " max: " << m_NrBins[0] << std::endl;
//				std::cout << IDS[1] << " max: " << m_NrBins[1] << std::endl;
//				std::cout << IDS[2] << " max: " << m_NrBins[2] << std::endl;
				//std::cout << "return " << g_IntMax << std::endl;
				return g_IntMax;
			}
	    return IDS[ 0] * m_NrBins[1] * m_NrBins[2] + IDS[1] * m_NrBins[2] + IDS[2];
		
//	    return IDS[ 0]  +  IDS[1] * m_NrBins[0]  +  IDS[2] * m_NrBins[0] * m_NrBins[1];
        }

        virtual int IDFromPosition( const math::Vector3N &POS) const
        {
        	return IDFromIDs( IDsFromPosition( POS));
        }

        virtual store::Vector3N< size_t> NeighborWithLowestIndices( const math::Vector3N &POSITION) const
        {
        	DebugWrite( __PRETTY_FUNCTION__);
            store::Vector3N< size_t> ids;

            ids[ 0] = size_t( ( POSITION[ 0] - m_Min[ 0] - 0.5 * m_Delta[ 0]) / m_Delta[ 0]); // time compare with bracket operator
            ids[ 1] = size_t( ( POSITION[ 1] - m_Min[ 1] - 0.5 * m_Delta[ 1]) / m_Delta[ 1]);
            ids[ 2] = size_t( ( POSITION[ 2] - m_Min[ 2] - 0.5 * m_Delta[ 2]) / m_Delta[ 2]);

            return ids;
        }

        virtual store::Vector3N< size_t> NeighborWithLowestIndices( const math::Vector3N &POSITION, const store::Vector3N< size_t> &NEAREST_NEIGHBOR) const
        {
            store::Vector3N< size_t> zero( NEAREST_NEIGHBOR);
            for( size_t i = 0; i < 3; ++i)
            {
                if
                (
                        m_Min[ i] + ( NEAREST_NEIGHBOR[ i] + 0.5) * m_Delta[ i] > POSITION[ i]
                        || NEAREST_NEIGHBOR[ i] == m_NrBins[ i] - 1
                )
                {
                    --zero[ i];
                }
            }
            return zero;
        }

        virtual void SetGridPoint( const math::Vector3N &POSITION, const t_RETURN &DATA)
        {
        	std::cout << __PRETTY_FUNCTION__ << "\n ===> should not be called!" << std::endl;
        	exit( -1);
	}
//        virtual void SetGridPoint( const math::Vector3N &POSITION, const t_RETURN &DATA){}
//        { SetGridPoint( IDsFromPosition( POSITION), DATA);}

        virtual void SetGridPoint( const int &I, const int &J, const int &K, const t_RETURN &DATA)
        {
        	std::cout << __PRETTY_FUNCTION__ << "\n ===> should not be called!" << std::endl;
        	exit( -1);
	}

        virtual void SetGridPoint( const store::Vector3N< int> &INDICES, const t_RETURN &DATA)
        {
        	std::cout << __PRETTY_FUNCTION__ << "\n ===> should not be called!" << std::endl;
        	exit( -1);
	}

        virtual t_RETURN &GetSetGridPoint(  const store::Vector3N< int> &INDICES)
        {
        	std::cout << __PRETTY_FUNCTION__ << "\n ===> should not be called!" << std::endl;
        	exit( -1);
        }

        virtual const math::Vector3N &GetMax() const{ return m_Max;}

        virtual const math::Vector3N &GetMin() const{ return m_Min;}

        virtual const math::Vector3N &GetDelta() const{ return m_Delta;}

        virtual void SetMax( const math::Vector3N &MAX){ m_Max = MAX;}

        virtual void SetMin( const math::Vector3N &MIN){ m_Min = MIN;}

        virtual void SetDelta( const math::Vector3N &DELTA){ m_Delta = DELTA;}

//        virtual math::Vector3N &GetNrBins() const{ return IndexGrid< t_RETURN>::GetNrBins();}


        virtual void InitializeDataStructure(){ std::cout << "InitializeDataStructure from PositionGrid should not be called" << std::endl; exit( -1);}

    protected:

        virtual void AdjustLimits( const float &PRECISION)
        {
            std::cout << "...adjust limits ..." << m_Max << std::endl;
            for( size_t i = 0; i < 3; ++i)
            {
                float d( ( m_Max[i] - m_Min[i]) / m_Delta[i]);
                float rest( d - size_t( d));
                if( rest > PRECISION && rest < 0.5)
                {
                    m_Max[i] = m_Min[i] + size_t( d) * m_Delta[i];
                    std::cout << "maximum for " << i << " decreased to " << m_Max[i] << "!" << std::endl;
                }
                else if( rest != 0)
                {
                    m_Max[i] = m_Min[i] + ( size_t( d) + 1) * m_Delta[i];
                    std::cout << "maximum for " << i << " increased to " << m_Max[i] << "!" << std::endl;
                }
            }
            std::cout << "... limit adjustment done: " << m_Max << std::endl;
        }


    public:

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << "bins: " << m_NrBins << std::endl;
            STREAM << "min: " << m_Min << std::endl;
            STREAM << "max: " << m_Max << std::endl;
            STREAM << "delta: " << m_Delta << std::endl;
            return STREAM;
        }

//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
//        {
//            return STREAM;
//        }

        virtual std::istream &Read( std::istream &STREAM)
        {
        	DebugWrite( __PRETTY_FUNCTION__);
            std::string str;

        	STREAM >> str;
        	if( str != "bins:")
        	{
        		std::cout << "===> not the correct header (" << str << ") for: " << __PRETTY_FUNCTION__ << std::endl;
        	}
        	STREAM >> m_NrBins; // [0] >> m_NrBins[1] >> m_NrBins[2];

        	STREAM >> str;
            if( str != "min:")
            {
            	std::cout << "===> not the correct \'min:\' string-id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
            }
            STREAM >> m_Min;//[0] >> m_Min[1] >> m_Min[2];
            STREAM >> str;
            if( str != "max:")
            {
            	std::cout << "===> not the correct \'max:\' string-id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
            }
            STREAM >> m_Max;//[0] >> m_Max[1] >> m_Max[2];
            STREAM >> str;
            if( str != "delta:")
            {
            	std::cout << "===> not the correct \'delta:\'string-id <" << str << "> in " << __PRETTY_FUNCTION__ << std::endl;
            }
            STREAM >> m_Delta;//[0] >> m_Delta[1] >> m_Delta[2];
            return STREAM;
        }


        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }


    }; // end class Grid

} // end namespace store

#endif //POSITION_GRID_H
