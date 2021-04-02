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


#ifndef POINT_SURFACE_H_
#define POINT_SURFACE_H_

#include "../math/vector_functions.h"
#include "../storage/shared_pointer_vector.t.h"

#include "surface_point.h"
#include "surfaceIF.h"
#include "object.h"

namespace geom
{

    class PointSurface
    : public SurfaceIF
    {
    protected:
        store::ShPtrVec< SurfacePoint>       m_Points;
        float                                     m_AreaSize;

    public:
        PointSurface()
        : m_Points(),
        m_AreaSize( 0.0)
        { /*std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " default constructor " << std::endl;*/}

        PointSurface( const std::vector< std::pair< float, math::Vector3N> > &DATA)
        : m_Points( DATA.size()),
        m_AreaSize( 0.0)
        {
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " data constructor " << std::endl;
            std::vector< std::pair< float, math::Vector3N> >::const_iterator data_itr( DATA.begin());
            for( std::vector< boost::shared_ptr< SurfacePoint> >::iterator itr( m_Points.begin()); itr != m_Points.end() && data_itr != DATA.end(); ++itr, ++data_itr)
            {
                *itr = boost::shared_ptr< SurfacePoint>( new SurfacePoint( data_itr->second, data_itr->first));
            }
        }

        PointSurface( const store::ShPtrVec< SurfacePoint> &POINTS, const float &AREA = 0.0)
        : m_Points( POINTS.size()),
        m_AreaSize( AREA)
        {
//            std::cout <<  "PointSurf construct from data" << std::endl;
            std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator orig_itr( POINTS.begin());
            for( std::vector< boost::shared_ptr< SurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++orig_itr)
            {
//                ( *this_itr)( new SurfacePoint( ( *orig_itr)->Clone())); // hard copy
                ( *this_itr) = boost::shared_ptr< SurfacePoint>( ( *orig_itr)->Clone()); // hard copy
            }
        }

//        PointSurface( const store::ShPtrVec< geom::PointSurfaceObject> &OBJECTS)
//        : m_Points(),
//        m_AreaSize()
//        {
//            for( std::vector< boost::shared_ptr< geom::PointSurfaceObject> >::iterator itr( OBJECTS.begin()); itr != OBJECTS.end(); ++itr)
//            {
//                m_Points += ( *itr)->GetSurfacePoints();
//            }
//        }
//
//        PointSurface( const mol::SimpleMolecule< mol::SurfAtom> &MOL)
//        : m_Points(), m_Area()
//        {
//            for( std::vector< boost::shared_ptr< mol::SurfAtom> >::const_iterator itr = MOL.GetAtoms().begin(); itr != MOL.GetAtoms().end(); ++itr)
//            {
//                m_Points += ( *itr)->GetSurf()->GetSurface().GetData();
//            }
//        }

//        PointSurface( const store::ShPtrVec< Object> > & OBJECTS)
//        : m_Points()// GeometricalObjectsToSurface( OBJECTS))
//        {}

        PointSurface( const PointSurface &POINT_SURFACE)
        : m_Points( POINT_SURFACE.m_Points.size()),
        m_AreaSize( POINT_SURFACE.m_AreaSize)
        {
//            std::cout <<  "PointSurf copy constructor" << std::endl;
            std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator orig_itr( POINT_SURFACE.m_Points.begin());
            for( std::vector< boost::shared_ptr< SurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++orig_itr)
            {
//                ( *this_itr)( new SurfacePoint( ( *orig_itr)->Clone())); // hard copy
                ( *this_itr) = boost::shared_ptr< SurfacePoint>( ( *orig_itr)->Clone()); // hard copy
            }
        }


        virtual ~PointSurface(){}

        virtual PointSurface *Clone() const { return new PointSurface( *this);}

        virtual PointSurface operator * ( const float &FACTOR) const
        {
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " operator * " << std::endl;
            store::ShPtrVec< SurfacePoint> data( m_Points.size());
            std::vector< boost::shared_ptr< SurfacePoint> >::iterator itr( data.begin());
            for( std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr, ++itr)
            {
                boost::shared_ptr< SurfacePoint> tmp( ( *this_itr)->Clone());
                tmp->SetPosition( FACTOR * tmp->GetPosition());
                tmp->SetSurface( math::Square( FACTOR) * tmp->GetSurface());         /// <= SPHERE ONLY !!!!!!!!!!
                *itr = tmp;
            }
            return data;
        }

        virtual PointSurface operator + ( const math::Vector3N &SHIFT) const
        {
            PointSurface tmp( *this);
            tmp.Translate( SHIFT);
            return tmp;
        }

        virtual PointSurface &operator += ( const PointSurface &SURF)
        {
            m_Points += SURF.m_Points;
            return *this;
        }

        virtual PointSurface &PushBack( const boost::shared_ptr< SurfacePoint> &SURF_POINT)
        {
            m_Points.push_back( SURF_POINT);
            return *this;
        }

        virtual PointSurface &Translate( const math::Vector3N &SHIFT)
        {
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " translate " << std::endl;
            for( std::vector< boost::shared_ptr< SurfacePoint> >::iterator this_itr( m_Points.begin()); this_itr != m_Points.end(); ++this_itr)
            {
                ( *this_itr)->SetPosition( ( *this_itr)->GetPosition() + SHIFT);
            }
            return *this;
        }

        virtual const store::ShPtrVec< SurfacePoint> &GetData() const
        { return m_Points;}

        virtual store::ShPtrVec< SurfacePoint> &Data()
        { return m_Points;}

        virtual size_t GetNrPoints() const
        { return m_Points.size();}

        virtual float GetFreeSurface() const
        { return float();}


        virtual std::istream& Read( std::istream &STREAM)
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }


        virtual std::ostream& Write( std::ostream &STREAM) const
        {
//            STREAM << GetClassName() << std::endl;
            for( std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                STREAM << **itr;
            }
            STREAM << m_AreaSize << std::endl;
            return STREAM;
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
        {
            for( std::vector< boost::shared_ptr< SurfacePoint> >::const_iterator itr( m_Points.begin()); itr != m_Points.end(); ++itr)
            {
                ( *itr)->WriteVmdCommands( STREAM);
            }
            return STREAM;
        }

    };
}


#endif /* POINT_SURFACE_H_ */
