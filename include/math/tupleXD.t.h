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


#ifndef TUPLE_XD_H
#define TUPLE_XD_H

#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>

#include "../string/io_string_functions.h"
#include "../readwrite/stream_operator.h"

namespace math
{
  ////////////////////////////////////////////////////////////////////////////
  //! A tuple of flexible dimension and sizes.
  //! dim = 1: vector, dim = 2: matrix, dim = 3: tensor, ...
    //! all used vectors, except for data have to have size = dim
  ////////////////////////////////////////////////////////////////////////////

    template< typename t_DATA, int DIM>
    class TupleXD
    : public StreamOperator
    {
    protected:
        std::vector< int> 	 m_NrElements;
        std::vector< t_DATA> m_Data;
        std::vector< float>  m_Min;
        std::vector< float>  m_Delta;
        t_DATA               m_Default; //!< value to return for points outside of the grid
    
    public:
        TupleXD()
        : m_NrElements( DIM, 0),
          m_Data(),
          m_Min( DIM, 0),
          m_Delta( DIM, 0),
          m_Default()
          {}

        TupleXD( const std::vector< int> &SIZES, const t_DATA &DEFAULT = 0)
        : m_NrElements( SIZES),
          m_Data( CalcTotalSize( SIZES)),
          m_Min( DIM, 0),
          m_Delta( DIM, 0),
          m_Default( DEFAULT)
          {
            assert( m_NrElements.size() == DIM);
          }

          TupleXD( const float &DELTA, const std::vector< float> &MINIMUM, const std::vector< float> &MAXIMUM, const t_DATA &DEFAULT = 0)
            : m_NrElements( CalcNrElements( MINIMUM, MAXIMUM, DELTA)),
          m_Data(),
          m_Min( MINIMUM),
          m_Delta( DIM, DELTA), 
          m_Default( DEFAULT)
          {
        	  m_Data = std::vector< t_DATA>( CalcTotalSize( m_NrElements), t_DATA( 0));
          }

        TupleXD( const std::vector< size_t> &SIZES, const std::vector< float> &MINIMUM, const std::vector< float> &DELTA, const t_DATA &DEFAULT = 0)
        : m_NrElements( SIZES),
          m_Data( CalcTotalSize( SIZES)),
          m_Min( MINIMUM),
          m_Delta( DELTA), 
          m_Default( DEFAULT)
          {
            assert( m_Min.size() == DIM &&  m_NrElements.size() == DIM && m_Delta.size() == DIM);
          }

        TupleXD( const std::vector< size_t> &SIZES, const std::vector< float> &MINIMUM, const std::vector< float> &DELTA, const std::vector< t_DATA> &DATA, const t_DATA &DEFAULT = 0)
        : m_NrElements( SIZES),
          m_Data( DATA),
          m_Min( MINIMUM),
          m_Delta( DELTA), 
          m_Default( DEFAULT)
          {
            assert( m_Data.size() == CalcTotalSize( m_NrElements) && m_Min.size() == DIM &&  m_NrElements.size() == DIM && m_Delta.size() == DIM);
          }


        virtual ~TupleXD(){}

        virtual TupleXD *Clone() const
        {
            return new TupleXD();
        }

        virtual
        size_t GetSize() const
        { return m_Data.size();}


        virtual
        std::vector< t_DATA>&
        Data()
        {
            return m_Data;
        }


        virtual
        const std::vector< t_DATA>&
        GetData() const
        {
            return m_Data;
        }

        virtual
        void
        SetData( const std::vector< t_DATA> &DATA)
        {
            assert( size_t( CalcTotalSize()) == DATA.size());
            m_Data = DATA;
        }

        virtual
        const std::vector< int> &
        GetNrElements() const
        {
            return m_NrElements;
        }

        virtual
        void
        SetNrElements( const std::vector< int> &SIZES)
        {
            assert( SIZES.size() == DIM);
            m_NrElements = SIZES;
        }

        virtual
        const std::vector< float> &
        GetMinimum() const
        {
            return m_Min;
        }

        virtual
        void
        SetMinimum( const std::vector< float> &MIN)
        {
            assert( MIN.size() == DIM);
            m_Min = MIN;
        }

        virtual
        const std::vector< float> &
        GetDelta() const
        {
            return m_Delta;
        }

        virtual
        void
        SetDelta( const std::vector< float> &DELTA)
        {
            assert( DELTA.size() == DIM);
            m_Delta = DELTA;
        }

        virtual
        const t_DATA&
        GetDefault() const
        {
        	return m_Default;
        }

        virtual
        t_DATA &
        Value( const int &ID)
        {
            if( ID < 0 || ID >= int( m_Data.size()))
            {
                return m_Default;
            }
            return m_Data[ ID];
        }

        virtual
        t_DATA &
        operator()( const int &ID)
        {
            return Value( ID);
        }

        virtual
        t_DATA& operator()( const std::vector< int> &IDS)
        {
            if( !IsInsideGrid( IDS))
            {
                return m_Default;
            }
            return m_Data[ CalcID( IDS)];
        }

        virtual
        const t_DATA&
        operator()( const std::vector< int> &IDS) const
        {
            if( !IsInsideGrid( IDS))
            {
                return m_Default;
            }
            return m_Data[ CalcID( IDS)];
        }

        virtual
        t_DATA &
        operator()( const std::vector< float> &POINT)
        {
            assert( POINT.size() == DIM);
            if( !IsInsideGrid( POINT))
            {
                return m_Default;
            }
            return m_Data[ CalcID( CalcIndices( POINT))];
        }

        virtual
        const t_DATA &
        operator()( const std::vector< float> &POINT) const
        {
            if( !IsInsideGrid( POINT))
            {
                  return m_Default;
            }
            return m_Data[ CalcID( CalcIndices( POINT))];
        }

        virtual
        bool 
        IsInsideGrid(  const std::vector< float> &POINT) const
        {
            assert( POINT.size() == DIM);
            std::vector< int>::const_iterator size_itr = m_NrElements.begin();

            for( std::vector< float>::const_iterator p_itr = POINT.begin(), min_itr = m_Min.begin(), delta_itr = m_Delta.begin(); p_itr != POINT.end() && min_itr != m_Min.end(); ++p_itr, ++min_itr, ++size_itr, ++delta_itr)
            {
                if( *p_itr < *min_itr || *p_itr >= *min_itr + *size_itr * *delta_itr - 1e-6)
                {
                    return false;
                }
            }
            return true;
        }

        virtual
        bool 
        IsInsideGrid(  const std::vector< int> &INDICES) const
        {
            assert( INDICES.size() == DIM);
            for( std::vector< int>::const_iterator itr = INDICES.begin(), size_itr = m_NrElements.begin(); itr != INDICES.end(); ++itr, ++size_itr)
            {
                if( *itr < 0 || *itr >= *size_itr)
                {
                  return false;
                }
            }
            return true;
        }


        virtual
        typename std::vector< t_DATA>::iterator
        Iterator( const size_t &ID)
        {
            return m_Data.begin() + ID;
        }


        virtual
        typename std::vector< t_DATA>::iterator
        Iterator( const std::vector< int> &IDS)
        {
            return m_Data.begin() + CalcID( IDS);
        }


        virtual
        typename std::vector< t_DATA>::iterator
        Iterator( const std::vector< float> &POINT)
        {
            return m_Data.begin() + CalcID( CalcIndices( POINT));
        }

        virtual
        TupleXD 
        SubMatrix( const typename std::vector< t_DATA>::iterator & CENTER_POINT, const std::vector< size_t> & NR_SHELLS)
        {
            TupleXD result;
            return result;
        }

        virtual
        TupleXD 
        SubMatrix( const std::vector< size_t> &IDS_MINIMUM, const std::vector< size_t> &SIZES)
        {
            TupleXD result;
            return result;
        }

        virtual
        int
        CalcTotalSize( const std::vector< int> &SIZES) const
        {
            int size( 1);
            for( std::vector< int>::const_iterator itr = SIZES.begin(); itr != SIZES.end(); ++itr)
            { size *= *itr;}
            return size;
        }

        virtual
        int
        CalcTotalSize() const
        {
            int size( 1);
            for( std::vector< int>::const_iterator itr = m_NrElements.begin(); itr != m_NrElements.end(); ++itr)
            { size *= *itr;}
            return size;
        }


        virtual
        std::vector< float>
        CalcMax() const
        {
            std::vector< float> max( DIM);
            std::vector< float>::iterator
                max_itr = max.begin();
            std::vector< float>::const_iterator
                min_itr = m_Min.begin(),
                delta_itr = m_Delta.begin();
            std::vector< int>::const_iterator
                size_itr = m_NrElements.begin();
            for( ; min_itr != m_Min.end() && delta_itr != m_Delta.end(); ++min_itr, ++delta_itr, ++size_itr, ++max_itr)
            {
                *max_itr = *min_itr + *size_itr * *delta_itr;
            }
            return max;
        }

        virtual
        std::vector< float>
        CalcCenter() const
        {
            std::vector< float> center( DIM);
            std::vector< float>::iterator
                center_itr = center.begin();
            std::vector< float>::const_iterator
                min_itr = m_Min.begin(),
                delta_itr = m_Delta.begin();
            std::vector< int>::const_iterator
                size_itr = m_NrElements.begin();
            for( ; min_itr != m_Min.end() && delta_itr != m_Delta.end(); ++min_itr, ++delta_itr, ++size_itr, ++center_itr)
            {
                *center_itr = *min_itr + 0.5 * *size_itr * *delta_itr;
            }
            return center;
        }
          


        virtual
        int
        CalcID( const std::vector< int> &IDS) const
        {
            size_t dim( 1);
            size_t id( 0);
            for( std::vector< int>::const_reverse_iterator id_itr = IDS.rbegin(), size_itr = m_NrElements.rbegin(); id_itr != IDS.rend() && size_itr != m_NrElements.rend(); ++id_itr, ++size_itr)
            {
                id += *id_itr * dim;
                dim *= *size_itr;
            }
//            if( id >= m_Data.size())
//            {
//                std::cout << "===> " <<  __FUNCTION__ << ":  " << id << " should be smaller: " << m_Data.size() << std::endl << "ids: ";
//                std::copy( IDS.begin(), IDS.end(), std::ostream_iterator< int>( std::cout, " "));
//                std::cout << std::endl;
//                //                exit( -1);
//            }
            return id;
        }


        virtual
        std::vector< int>
        CalcIndicesFromID( const int &ID) const
        {
            std::vector< int>
                factors( m_NrElements.size()),
                indices( m_NrElements.size());
            int dim = 1;
            std::vector<int>::iterator itr = factors.begin();
            std::vector< int>::const_reverse_iterator size_itr = m_NrElements.rbegin();
            for( ; size_itr != m_NrElements.rend(); ++size_itr, ++itr)
            {
                *itr = dim;
                dim *= *size_itr;
            }
//            std::copy( factors.begin(), factors.end(), std::ostream_iterator< int>( std::cout, " "));
//            std::cout << std::endl;
            int id = ID;
            std::vector< int>::const_reverse_iterator fitr = factors.rbegin();
            for( std::vector< int>::iterator ind_itr = indices.begin(); fitr != ( std::vector< int>::const_reverse_iterator) factors.rend(); ++fitr, ++ind_itr)
            {
                *ind_itr = int( id / *fitr);
                id -= *ind_itr * *fitr;
//                std::cout << *ind_itr << " " << id << std::endl;
            }
            return indices;
        }

        virtual
        std::vector< int>
        CalcIndices( const std::vector< float> &POINT) const
        {
            assert( POINT.size() == DIM);
            std::vector< int>
                ids( DIM);
            std::vector< float>::const_iterator
                point_itr = POINT.begin(),
                min_itr = m_Min.begin(),
                delta_itr = m_Delta.begin();
            std::vector< int>::const_iterator
                size_itr = m_NrElements.begin();

            for( std::vector< int>::iterator id_itr = ids.begin(); id_itr != ids.end(); ++id_itr, ++point_itr, ++min_itr, ++delta_itr, ++size_itr)
            {
                *id_itr = int( ( *point_itr - *min_itr ) / *delta_itr);
#ifdef DEBUG
                if( *id_itr < 0 || *id_itr >= *size_itr)
                {
                    std::cout << "===> " << __FUNCTION__ << " " << *id_itr << "should be between 0 and " << *size_itr << std::endl;
                    std::copy( POINT.begin(), POINT.end(), std::ostream_iterator< float>( std::cout, " "));
                    std::cout << std::endl << std::endl;
                    std::cout << "min: " << std::endl;
                    std::copy( m_Min.begin(), m_Min.end(), std::ostream_iterator< float>( std::cout, " "));
                    std::cout << std::endl << std::endl;
                    std::vector< float> max = CalcMax();
                    std::cout << "max: " << std::endl;
                    std::copy( max.begin(), max.end(), std::ostream_iterator< float>( std::cout, " "));
                    std::cout << std::endl << std::endl;
                    std::cout << std::endl;
                    std::cout << "sizes: " << std::endl;
                    std::copy( m_NrElements.begin(), m_NrElements.end(), std::ostream_iterator< float>( std::cout, " "));
                    std::cout << std::endl << std::endl;
                }
#endif
            }
            return ids;
        }

        virtual
        std::vector< float>
        CalcPositionFromIndices( const std::vector< int> &INDICES) const
        {
            assert( INDICES.size() == DIM);
            std::vector< float> pos( DIM);
            std::vector< float>::const_iterator
                min_itr = m_Min.begin(),
                delta_itr = m_Delta.begin();
            std::vector< float>::iterator
                point_itr = pos.begin();

            for( std::vector< int>::const_iterator id_itr = INDICES.begin(); id_itr != INDICES.end(); ++id_itr, ++point_itr, ++min_itr, ++delta_itr)
            {
                *point_itr = *min_itr + ( *id_itr + 0.5) * *delta_itr;
            }
            return pos;
        }

        virtual
        std::vector< float>
        CalcPositionFromID( const int &ID) const
        {
        	std::vector< int>
				factors( DIM),
                indices( DIM);
            std::vector< float>
                pos( DIM);
            std::vector< float>::const_iterator
                min_itr = m_Min.begin(),
                delta_itr = m_Delta.begin();

            int dim = 1;
            std::vector<int>::iterator
                fac_itr = factors.begin();
            std::vector< int>::const_reverse_iterator
				size_itr = m_NrElements.rbegin();

            for( ; size_itr != m_NrElements.rend(); ++size_itr, ++fac_itr)
            {
                *fac_itr = dim;
                dim *= *size_itr;
            }
//            std::copy( factors.begin(), factors.end(), std::ostream_iterator< int>( std::cout, " "));
//            std::cout << std::endl;

            int id = ID, value;
            std::vector< int>::const_reverse_iterator fitr = factors.rbegin();
            for( std::vector< float>::iterator point_itr = pos.begin(); point_itr != pos.end(); ++fitr, ++point_itr, ++min_itr, ++delta_itr)
            {
                value = int( id / *fitr);
                id -= value * *fitr;
                *point_itr = *min_itr + ( value + 0.5) * *delta_itr;
            }
            return pos;
        }


        virtual 
        std::vector< int>
        ProjectOnGrid( const std::vector< int> &INDICES) const
        {
            assert( INDICES.size() == DIM);
            std::vector< int> indices( INDICES);
            std::vector< int>::iterator
                ind_itr = indices.begin();
            std::vector< int>::const_iterator
                size_itr = m_NrElements.begin();

            for( ; size_itr != m_NrElements.end(); ++ind_itr, ++size_itr)
            {
                if( *ind_itr < 0)
                {
                    *ind_itr = 0;
                }
                else if( *ind_itr >= *size_itr)
                {
                    *ind_itr = *size_itr - 1;
                }
            }
              return indices;
        }

        virtual
        std::vector< int>
        CalcNrElements
        (
                const std::vector< float> &MIN,
                const std::vector< float> &MAX,
               const float &DELTA
        ) const
        {
            assert( MIN.size() == DIM && MAX.size() == DIM);
            std::vector< int> sizes( DIM);
            std::vector< int>::iterator size_itr = sizes.begin();
            std::vector< float>::const_iterator
            min_itr = MIN.begin(),
            max_itr = MAX.begin();

            for( ; min_itr != MIN.end(); ++min_itr, ++max_itr, ++size_itr)
            {
                *size_itr = int( ( *max_itr - *min_itr) / DELTA);
            }
            return sizes;
        }



        virtual
        std::ostream &
        Write( std::ostream& STREAM) const
        {
            STREAM << "TupleXD" << std::endl;
            STREAM << "dim:   " << DIM << std::endl;
            STREAM << "bins:  ";
            for( std::vector< int>::const_iterator itr = m_NrElements.begin(); itr != m_NrElements.end(); ++itr)
            {
                STREAM << *itr << " ";
            }
            STREAM << std::endl;

            STREAM << "min:   ";
            for( std::vector< float>::const_iterator itr = m_Min.begin(); itr != m_Min.end(); ++itr)
            {
                STREAM << *itr << " ";
            }
            STREAM << std::endl;

            STREAM << "delta: ";
            for( std::vector< float>::const_iterator itr = m_Delta.begin(); itr != m_Delta.end(); ++itr)
            {
                STREAM << *itr << " ";
            }
            STREAM << std::endl;

            STREAM << "data: " << m_Data.size() << std::endl;
            int cc( 0);
            for( typename std::vector< t_DATA>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM << *itr << "  ";
                ++cc;
                int dim( 1);
//                for( std::vector<int>::const_reverse_iterator ii = m_NrElements.rbegin(); ii != m_NrElements.rend(); ++ii)
//                {
//                    dim *= *ii;
//                    //          std::cout << *ii << " " << dim << " " << std::endl;
//                    if( cc % dim == 0)
//                    {
//                        STREAM << std::endl;
//                    }
//                    //          else{ continue;}
//                }
                for( std::vector<int>::const_iterator ii = m_NrElements.begin(); ii != m_NrElements.end(); ++ii)
                {
                    dim *= *ii;
                    //          std::cout << *ii << " " << dim << " " << std::endl;
                    if( cc % dim == 0)
                    {
                        STREAM << std::endl;
                    }
                    //          else{ continue;}
                }
            }
            STREAM << "default: " << m_Default << std::endl;
            return STREAM;
        }

       virtual
        std::ostream &
        WriteGnuplot( std::ostream& STREAM) const
        {
            int cc( 0);
            for( typename std::vector< t_DATA>::const_iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM << *itr << std::endl;
                ++cc;
                int dim( 1);
                for( std::vector<int>::const_iterator ii = m_NrElements.begin(); ii != m_NrElements.end(); ++ii)
                {
                    dim *= *ii;
                    if( cc % dim == 0)
                    {
                        STREAM << std::endl;
                    }
                }
            }
            return STREAM;
        }

        virtual
        std::istream &
        Read( std::istream& STREAM)
        {
            std::string str;
            int size;


            STREAM >> str;
            assert( str == "TupleXD");
            STREAM >> str;
            assert( str == "dim:");
            STREAM >> size;
            assert( size == DIM);
            STREAM >> str;
            assert( str == "bins:");

            for( std::vector< int>::iterator itr = m_NrElements.begin(); itr != m_NrElements.end(); ++itr)
            {
                STREAM >> *itr;
            }

            STREAM >> str;
            assert( str == "min:");
            for( std::vector< float>::iterator itr = m_Min.begin(); itr != m_Min.end(); ++itr)
            {
                STREAM >> *itr;
            }

            STREAM >> str;
            assert( str == "delta:");
            for( std::vector< float>::iterator itr = m_Delta.begin(); itr != m_Delta.end(); ++itr)
            {
                STREAM >> *itr;
            }

            STREAM >> str;
            assert( str == "data:");
            STREAM >> size;
            m_Data = std::vector< t_DATA>( size);
            for( typename std::vector< t_DATA>::iterator itr = m_Data.begin(); itr != m_Data.end(); ++itr)
            {
                STREAM >> *itr;
            }

            STREAM >> str;
            assert( str == "default:");
            STREAM >> m_Default;

            return STREAM;
        }

        virtual
        std::ostream &
        WriteAsArray( std::ostream& STREAM, const char *SEPERATOR = " ") const
        {
        	std::copy( m_Data.begin(), m_Data.end(), std::ostream_iterator< t_DATA>( STREAM, SEPERATOR));
        	STREAM << std::endl;
        	return STREAM;
        }


  }; // end class TupleXD
} // end namespace math

#endif // TUPLE_XD_H
