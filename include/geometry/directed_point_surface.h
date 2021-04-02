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


#ifndef DIRECTED_POINT_SURFACE_H_
#define DIRECTED_POINT_SURFACE_H_

#include "../math/vector_functions.h"
#include "../storage/shared_pointer_vector.t.h"
#include "../storage/triplet.h"
#include "directed_surface_point.h"
#include "surfaceIF.h"
#include "object.h"


namespace geom
{

    class DirectedPointSurface
    : public SurfaceIF
    {
    protected:
        store::ShPtrVec< DirectedSurfacePoint>       m_Points;
        float                                     m_AreaSize;

    public:
        DirectedPointSurface()
        : m_Points(),
        m_AreaSize( 0.0)
        { /*std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " default constructor " << std::endl;*/}

        DirectedPointSurface( const std::vector< store::Triplet< float, math::Vector3N, math::Vector3N> > &DATA)
        : m_Points( DATA.size()),
        m_AreaSize( 0.0)
        {
            DebugWrite( __FUNCTION__);
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " data constructor " << std::endl;
            std::vector< store::Triplet< float, math::Vector3N, math::Vector3N> >::const_iterator data_itr( DATA.begin());
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator itr( m_Points.begin()); itr != m_Points.end() && data_itr != DATA.end(); ++itr, ++data_itr)
            {
                *itr = boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( data_itr->Third(), data_itr->Second(), data_itr->First()));
            }
        }

        DirectedPointSurface( const store::ShPtrVec< DirectedSurfacePoint> &POINTS, const float &AREA = 0.0)
        : m_Points( POINTS.size()),
        m_AreaSize( AREA)
        {
            DebugWrite( __FUNCTION__);
//            std::cout <<  "PointSurf construct from data" << std::endl;
            std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator orig_itr( POINTS.begin());
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++orig_itr)
            {
//                ( *this_itr)( new DirectedSurfacePoint( ( *orig_itr)->Clone())); // hard copy
                ( *this_itr) = boost::shared_ptr< DirectedSurfacePoint>( ( *orig_itr)->Clone()); // hard copy
            }
        }

//        DirectedPointSurface( const store::ShPtrVec< geom::DirectedPointSurfaceObject> &OBJECTS)
//        : m_Points(),
//        m_AreaSize()
//        {
//            for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::iterator itr( OBJECTS.begin()); itr != OBJECTS.end(); ++itr)
//            {
//                m_Points += ( *itr)->GetDirectedSurfacePoints();
//            }
//        }
//
//        DirectedPointSurface( const mol::SimpleMolecule< mol::SurfAtom> &MOL)
//        : m_Points(), m_Area()
//        {
//            for( std::vector< boost::shared_ptr< mol::SurfAtom> >::const_iterator itr = MOL.GetAtoms().begin(); itr != MOL.GetAtoms().end(); ++itr)
//            {
//                m_Points += ( *itr)->GetSurf()->GetSurface().GetData();
//            }
//        }

//        DirectedPointSurface( const store::ShPtrVec< Object> > & OBJECTS)
//        : m_Points()// GeometricalObjectsToSurface( OBJECTS))
//        {}

        DirectedPointSurface( const DirectedPointSurface &POINT_SURFACE)
        : m_Points( POINT_SURFACE.m_Points.size()),
        m_AreaSize( POINT_SURFACE.m_AreaSize)
        {
            DebugWrite( __FUNCTION__);
//            std::cout <<  "PointSurf copy constructor" << std::endl;
            std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator orig_itr( POINT_SURFACE.m_Points.begin());
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++orig_itr)
            {
//                ( *this_itr)( new DirectedSurfacePoint( ( *orig_itr)->Clone())); // hard copy
                ( *this_itr) = boost::shared_ptr< DirectedSurfacePoint>( ( *orig_itr)->Clone()); // hard copy
            }
        }


        virtual ~DirectedPointSurface(){}

        virtual DirectedPointSurface *Clone() const { return new DirectedPointSurface( *this);}

        virtual DirectedPointSurface operator * ( const float &FACTOR) const
        {
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " operator * " << std::endl;
            store::ShPtrVec< DirectedSurfacePoint> data( m_Points.size());
            std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator itr( data.begin());
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++itr)
            {
                boost::shared_ptr< DirectedSurfacePoint> tmp( ( *this_itr)->Clone());
                tmp->SetPosition( FACTOR * tmp->GetPosition());
                tmp->SetSurface( math::Square( FACTOR) * tmp->GetSurface());         /// <= SPHERE ONLY !!!!!!!!!!
                *itr = tmp;
            }
            return data;
        }

        virtual DirectedPointSurface operator + ( const math::Vector3N &SHIFT) const
        {
            DirectedPointSurface tmp( *this);
            tmp.Translate( SHIFT);
            return tmp;
        }

        virtual DirectedPointSurface &operator += ( const DirectedPointSurface &SURF)
        {
            m_Points += SURF.m_Points;
            return *this;
        }

        virtual DirectedPointSurface &PushBack( const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT)
        {
            m_Points.push_back( SURF_POINT);
            return *this;
        }

        virtual DirectedPointSurface &Translate( const math::Vector3N &SHIFT)
        {
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " translate " << std::endl;
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr)
            {
                ( *this_itr)->SetPosition( ( *this_itr)->GetPosition() + SHIFT);
            }
            return *this;
        }

        virtual const store::ShPtrVec< DirectedSurfacePoint> &GetData() const
        { return m_Points;}

        virtual store::ShPtrVec< DirectedSurfacePoint> &Data()
        { return m_Points;}

        virtual size_t GetNrPoints() const
        { return m_Points.size();}

        virtual float GetFreeSurface() const
        { return float();}

        virtual boost::shared_ptr< DirectedSurfacePoint> ClosestSurfacePoint( const math::Vector3N &POSITION) const
        {
            float
                distance,
                closest_distance( std::numeric_limits< float>::max());
            boost::shared_ptr< DirectedSurfacePoint>
                closest_surf_point;

            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                distance = math::Distance( POSITION, ( *itr)->GetPosition());
                if( distance < closest_distance)
                {
                    closest_distance = distance;
                    closest_surf_point = *itr;
                }
            }

            return closest_surf_point;
        }


        virtual boost::shared_ptr< std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> > >
        DistanceSortedSurfacePoints( const math::Vector3N &POSITION) const
        {
            boost::shared_ptr< std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> > >
                map( new std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> >());
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                map->insert( std::make_pair( math::Distance( POSITION, ( *itr)->GetPosition()), *itr));
            }
            return map;
        }


        virtual std::istream& Read( std::istream &STREAM)
        {
            std::string str;
            size_t nr;
            STREAM >> str;
            assert( str == GetClassName());
            STREAM >> nr;
            DirectedSurfacePoint point;
            for( size_t i = 0; i < nr; ++i)
            {
                STREAM >> point;
                m_Points.push_back( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint( point)));
            }
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }


        virtual std::ostream& Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName()  << std::endl;
            STREAM << m_Points.size() << std::endl;
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                STREAM << **itr;
            }
            STREAM << m_AreaSize << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
            for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                ( *itr)->WriteVmdCommands( STREAM, MOL_ID);
            }
            return STREAM;
        }


        virtual std::ostream &WriteAsPdb( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
          static size_t s_atom_count = 0;
          for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr, ++s_atom_count)
            {
              ( *itr)->WriteAsPdb( STREAM, MOL_ID, s_atom_count);
            }
            return STREAM;
        }


        virtual std::string GetClassName() const
        {
//            std::cout << __FUNCTION__ << std::endl;
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }
    };
}


#endif /* DIRECTED_POINT_SURFACE_H_ */
