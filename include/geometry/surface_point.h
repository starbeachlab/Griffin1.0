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


#ifndef SURFACE_POINT_H
#define SURFACE_POINT_H

#include <deque>
#include <set>
#include <vector>
#include <cassert>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "../math/vector3N.h"

namespace geom
{

    //! SurfacePoint contains information about an point representing the surface
    class SurfacePoint
    : public math::Vector3N
    {
    protected:
        float m_Surf;                  //!< surface area represented by the SurfacePoint
//        ptr< Object>  m_PtrToObject;

    public:
        SurfacePoint()
        : math::Vector3N(),
        m_Surf()
        {/*std::cout << "SurfacePoint default constructor" << std::endl;*/}

        SurfacePoint( const math::Vector3N &POS, const float &SURF = 0.0)
        : math::Vector3N( POS),
        m_Surf( SURF)
        {/*std::cout << "SurfacePoint construct from data" << std::endl;*/}

        SurfacePoint( const SurfacePoint &SURF_POINT)
        : math::Vector3N( SURF_POINT),
        m_Surf( SURF_POINT.m_Surf)
        {/*std::cout << "SurfacePoint copy constructor" << std::endl;*/}

        SurfacePoint(  const float &X,  const float &Y, const float &Z)
        : math::Vector3N( X, Y, Z)
        {}


        virtual ~SurfacePoint(){ /*std::cout << "SurfacePoint destroyed" << std::endl;*/}

        virtual SurfacePoint *Clone() const{ /*std::cout << "SurfPoint Clone" << std::endl;*/ return new SurfacePoint( *this);}

        virtual math::Vector3N GetPosition() const{ return *this;}

        virtual void SetPosition( const math::Vector3N &POS){ *this = POS;}

        virtual const float &GetSurface() const{ return m_Surf;}

        virtual void SetSurface( const float &VALUE){ m_Surf = VALUE;}

        virtual std::ostream& Write( std::ostream &STREAM) const
        {
//            STREAM << SurfacePoint::GetClassName() << std::endl;
            STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
            STREAM << "position: " << GetPosition() << std::endl;
            STREAM << "surface: " << m_Surf << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
//            STREAM << "graphics " << MOL_ID << " color blue" <<  std::endl;
            STREAM << "graphics " << MOL_ID << " sphere {" << this->GetPosition()(0) << " " << this->GetPosition()(1) << " " << this->GetPosition()(2) << "} radius 0.1 resolution 3" << std::endl; // resolution 21" << std::endl;
            return STREAM;
        }

        virtual std::istream& Read( std::istream &STREAM)
        {
            std::string str;
            STREAM >> str;
            assert( str == GetClassName());
            math::Vector3N::Read( STREAM);
            STREAM >> m_Surf;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
//            std::cout << __FUNCTION__ << std::endl;
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }


    }; // end class SurfacePoint
} // end namespace geom

#define SurfaceContainer std::deque< boost::shared_ptr< SurfacePoint> >
#define SurfacePtr    boost::shared_ptr< SurfacePoint>
#define SurfaceIterator  std::deque< boost::shared_ptr< SurfacePoint> >::iterator


//class SurfaceMap  // should be replaced by a mysqlpp connector class
//  : public BinnedMap3D< float, std::set< size_t> >
//{
//  //  BinnedMap3D< float, std::set< size_t> > m_SurfMap; //!< a voxel map; its coordinates are represented as floats; each voxel has a list of neighboring surface points (as indices of the SurfaceContainer)
//
// public:
//
//  SurfaceMap()
//    : m_SurfMap()
//    {}
//
//    SurfaceMap( const std::vector< float> &MIN, const std::vector< float> &MAX, const std::vector< size_t> &NR_BINS)
//      : BinnedMap3D< float, std::set< size_t> >( MIN, MAX, NR_BINS)//,
//      //      m_SurfMap( ReadSurfaceContainer( SURF, NR_BINS))
//      {}
//
//      virtual ~SurfaceMap(){}
//
//      //! translating a surface file into a SurfaceMap, most time-consuming step
//      //! future step: move this into mysql database
//      virtual void
//    ReadSurfaceContainer( const SurfaceContainer &CONT)
//      {
//    for( size_t i = 0; i <= m_NrBins[0]; ++i)
//      for( size_t j = 0; j <= m_NrBins[1]; ++j)
//        for( size_t k = 0; k <= m_NrBins[2]; ++k)
//          {
//        std::vector< float> pos( 3);
//        // avoid for loop for speed
//        pos[0] = m_Minimum[0] + i * m_Delta[0];
//        pos[1] = m_Minimum[1] + j * m_Delta[1];
//        pos[2] = m_Minimum[2] + k * m_Delta[2];
//
//        size_t count( 0), shortest_id;
//        float shortest_dist( std::numeric_limits< float>::max());
//        for( std::deque< boost::shared_ptr< SurfacePoint> >::const_iterator itr = CONT.begin(); itr != CONT.end(); ++itr, ++count)
//          {
//            float dist( math::SquaredNorm( 3, pos, itr->GetPosition()));
//            if( dist < shortest_dist)
//              {
//            shortest_dist = dist;
//            shortest_id = count;
//              }
//          }
//        // more systematically???:
//        operator()( i, j, k).insert( shortest_id);
//        if( k > 0)
//          { operator()( i, j, k-1).insert( shortest_id);}
//        if( j > 0)
//          {
//            operator()( i, j-1, k).insert( shortest_id);
//            if( k > 0)
//              { operator()( i, j-1, k-1).insert( shortest_id);}
//          }
//        if( i > 0)
//          {
//            operator()( i-1, j, k).insert( shortest_id);
//            if( k > 0)
//              { operator()( i-1, j, k-1).insert( shortest_id);}
//            if( j > 0)
//              {
//            operator()( i-1, j-1, k).insert( shortest_id);
//            if( k > 0)
//              { operator()( i-1, j-1, k-1).insert( shortest_id);}
//              }
//          }
//          }
//
//      }
//
//      virtual std::ostream &
//    Write( std::ostream &STREAM) const
//      {
//    STREAM << "SurfaceMap" << std::endl;
//    STREAM << m_Minimum[0] << " " << m_Minimum[1] << " " << m_Minimum[2] << std::endl;
//    STREAM << m_Maximum[0] << " " << m_Maximum[1] << " " << m_Maximum[2] << std::endl;
//    STREAM << m_NrBins[0] << " " <<  m_NrBins[1] << " " << m_NrBins[2] << std::endl;
//    STREAM << m_Delta[0] << " " << m_Delta[1] << " " << m_Delta[2] << std::endl;
//    for( std::vector< std::vector< std::vector< std::set< size_t> > > >::const_iterator tensor_itr( m_Data.begin()); tensor_itr != m_Data.end(); ++tensor_itr)
//      for( std::vector< std::vector< std::set< size_t> > >::const_iterator matrix_itr( tensor_itr->begin()); matrix_itr != tensor_itr->end(); ++matrix_itr)
//        for( std::vector< std::set< size_t> >::const_iterator vector_itr( matrix_itr->begin()); vector_itr != matrix_itr->end(); ++vector_itr)
//          {
//        STREAM << vector_itr->size()<< "  ";
//        for( std::set<size_t>::const_itr list_itr( vector_itr->begin()); list_itr != vector_itr->end(); ++list_itr)
//          { STREAM << *list_itr << "  ";}
//          }
//    STREAM << std::endl;
//    return STREAM;
//      }
//
//
//
//      virtual std::istream&
//    Read( std::istream &STREAM) const
//      {
//    std::string str;
//    size_t nr;
//    STREAM >> str;
//    assert( str == "SurfaceMap");
//    STREAM >> m_Minimum[0] >> m_Minimum[1] >> m_Minimum[2];
//    STREAM >> m_Maximum[0] >> m_Maximum[1] >> m_Maximum[2];
//    STREAM >> m_NrBins[0] >>  m_NrBins[1] >> m_NrBins[2];
//    STREAM >> m_Delta[0] >> m_Delta[1] >> m_Delta[2];
//    m_Data = BuildMap( m_Minimum, m_Maximum, m_NrBins);
//    for( std::vector< std::vector< std::vector< std::set< size_t> > > >::const_iterator tensor_itr( m_Data.begin()); tensor_itr != m_Data.end(); ++tensor_itr)
//      for( std::vector< std::vector< std::set< size_t> > >::const_iterator matrix_itr( tensor_itr->begin()); matrix_itr != tensor_itr->end(); ++matrix_itr)
//        for( std::vector< std::set< size_t> >::const_iterator vector_itr( matrix_itr->begin()); vector_itr != matrix_itr->end(); ++vector_itr)
//          {
//        STREAM >> nr;
//        for( size_t i = 0; i < nr; ++i)
//          {
//            STREAM >> value;
//            vector_itr->insert( value);
//          }
//          }
//    return STREAM;
//      }
//};
//
//
//
//
//class SurfaceHandler
//{
//private:
//
//  CalculateNearestSurfaceNeighboursOfAtom( const LipidContainer &LIPIDS)
//  {
//
//    for( std::deque< boost::shared_ptr< LipidAtom> >::iterator itr = atoms.begin(); itr != atoms.end(); ++itr)
//      for( std::deque< boost::shared_ptr< SurfacePoint> >::iterator itr = surface.begin(); itr != surface.end(); ++itr)
//
//  }
//
//
//  std::map< size_t, std::vector< float> >
//  CalculateForceVectorForAtomsWithinSurface( const LipidContainer &LIPIDS)
//  {
//    CalculateNearestSurfaceNeighboursOfAtom
//  }
//
//}; // end class SurfaceManager
//
//




#endif
