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


#ifndef FLOATING_DISTANCE_MAP_H_
#define FLOATING_DISTANCE_MAP_H_

#include <map>

#include "../../include/storage/shared_pointer_vector.t.h"
#include "../../include/math/double_functions.h"
#include "../../include/math/vector_functions.h"


namespace util
{
    template< typename t_DATA>
    bool operator == ( boost::shared_ptr< t_DATA> &D1, boost::shared_ptr< t_DATA> &D2)
    { return *D1 == *D2;}


    template< typename t_DATA>
    class FloatingDistanceMap
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        float                   m_Cutoff;
        store::ShPtrVec< t_DATA> m_Data;
        store::ShPtrVec< t_DATA> m_Previous;
        std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > >  m_DistanceMap; // replace with multimap!!
        std::map< boost::shared_ptr< t_DATA>, float>  m_Translations;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        FloatingDistanceMap()
        : m_Cutoff(),
        m_Data(),
        m_Previous(),
        m_DistanceMap(),
        m_Translations()
        { exit( -8);}

        //! construct from store::ShPtrVec< t_DATA>
        FloatingDistanceMap( const store::ShPtrVec< t_DATA> &DATA, const float &CUTOFF)
        : m_Cutoff( CUTOFF),
        m_Data(),
        m_Previous( DATA),
        m_DistanceMap( CalculateDistanceMap( DATA)),
        m_Translations()
        {
            m_Data = DATA.SoftCopy();
//            m_Previous = DATA.HardCopy();  // better in constructor
        }

        //! copy constructor
        FloatingDistanceMap( const FloatingDistanceMap &ORIGINAL)
        : m_Cutoff( ORIGINAL.m_Cutoff),
        m_Data(),
        m_Previous( ORIGINAL.m_Previous),
        m_DistanceMap( ORIGINAL.m_DistanceMap),
        m_Translations( ORIGINAL.m_Translations)
        {
            m_Data = ORIGINAL.m_Data.SoftCopy();
        }

        //! virtual destructor
        virtual ~FloatingDistanceMap(){}

        //! virtual copy constructor
        virtual FloatingDistanceMap *Clone() const{ return new FloatingDistanceMap( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > >
        CalculateDistanceMap( const store::ShPtrVec< t_DATA> &DATA)
        {
            std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > > map;
            for( typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator itr( DATA.begin()); itr + 1 != DATA.end(); ++itr)
                for( typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator bitr( itr + 1); bitr != DATA.end(); ++bitr)
                {
                    map.insert( std::make_pair( math::Distance( ( *itr)->GetPosition(), ( *bitr)->GetPosition()), std::make_pair( *itr, *bitr)));
                }
            return map;
        }

        void CalculateTranslations()
        {
            for( typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator recent( m_Data.begin()), previous( m_Previous.begin()); recent != m_Data.end() && previous != m_Previous.end(); ++recent, ++previous)
            {
                m_Translations[ *recent] = math::Distance( ( *recent)->GetPosition(), ( *previous)->GetPosition());
            }
        }

        void AdjustMap()
        {
            std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > > tmp;
            for( typename std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > >::const_iterator itr( m_DistanceMap.begin()); itr != m_DistanceMap.end(); ++itr)
            {
                float dist( itr->first - m_Translations[ itr->second.first] - m_Translations[ itr->second.second]);
                if( dist < m_Cutoff)
                { dist = math::Distance( itr->second.first->GetPosition(), itr->second.second->GetPosition());}
                tmp[ dist] = itr->second;
            }
            m_DistanceMap = tmp;
        }

        void Update()
        {
            CalculateTranslations();
            AdjustMap();
            m_Previous = m_Data.HardCopy();
        }

        std::vector< std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > > GetElementsToBeRecalculated() const
        {
            std::vector< std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > > altered;
            for
            (
                    typename
                    std::map< float, std::pair< boost::shared_ptr< t_DATA>, boost::shared_ptr< t_DATA> > >::const_iterator itr( m_DistanceMap.begin());
                    itr != m_DistanceMap.end() && itr->first < m_Cutoff;
                    ++itr
            )
            {
                altered.push_back( itr->second);
            }
            return altered;
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
            return STREAM;
        }

    }; // end class FloatingDistanceMap
} // end namespace util




#endif /* FLOATING_DISTANCE_MAP_H_ */
