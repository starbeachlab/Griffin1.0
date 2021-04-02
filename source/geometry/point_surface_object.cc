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


#include "../../include/geometry/point_surface_object.h"

namespace geom
{

    PointSurfaceObject::PointSurfaceObject()
    : m_Surf(),
    m_Object()
    { /*std::cout << "PointSurfObj default constructor" << std::endl;*/}

    PointSurfaceObject::PointSurfaceObject( const PointSurfaceObject &OBJECT)
    : m_Surf( OBJECT.m_Surf),
    m_Object( OBJECT.m_Object)
    { /*std::cout << "PointSurfObj copy constructor" << std::endl;*/}

    PointSurfaceObject::PointSurfaceObject( const boost::shared_ptr< Object> &OBJECT, const PointSurface &SURF)
    : m_Surf( SURF),
    m_Object( OBJECT)
    { /*std::cout << "PointSurfObj construct from data" << std::endl;*/}

    PointSurfaceObject::~PointSurfaceObject(){}

    PointSurfaceObject *PointSurfaceObject::Clone() const
    { return new PointSurfaceObject( *this);}

    const math::Vector3N &PointSurfaceObject::GetPosition() const
    { return m_Object->GetPosition();}

    void PointSurfaceObject::SetPosition( const math::Vector3N &NEWPOS)
    { m_Object->SetPosition( NEWPOS);}

    bool PointSurfaceObject::IsPointWithin( const math::Vector3N &POS) const
    {
    //            return false;
        return m_Object->IsPointWithin( POS);
    }

    const float PointSurfaceObject::GetRadius() const
    { return m_Object->GetRadius();}

    const store::ShPtrVec< SurfacePoint> &PointSurfaceObject::GetSurfacePoints() const
    { return m_Surf.GetData();}

    store::ShPtrVec< SurfacePoint> &PointSurfaceObject::SurfacePoints()
    { return m_Surf.Data();}

    const PointSurface &PointSurfaceObject::GetSurface() const
    { return m_Surf;}

    void PointSurfaceObject::SetSurface( const PointSurface &SURF)
    { m_Surf = SURF;}

    const boost::shared_ptr< Object> &PointSurfaceObject::GetObject() const
    { return m_Object;}

    float PointSurfaceObject::GetTotalSurface() const
    { return m_Object->GetTotalSurface();}

    float PointSurfaceObject::GetFreeSurface() const
    { return m_Surf.GetFreeSurface();}

    std::string PointSurfaceObject::GetClassName() const
    { return mystr::GetClassName( __PRETTY_FUNCTION__);}

    /*
    *
    *         virtual const boost::shared_ptr< PointSurface> &CalculateSurface( const float &RESOLUTION) const;
    {
        boost::shared_ptr< math::IterationVector> iter( Object::GetIterationVector());
        boost::shared_ptr< Function< math::Vector3N, math::Vector3N> > trans( Object::GetCoordinateTransformation());
        m_Surf = boost::shared_ptr< PointSurface>( new PointSurface( math::Integrator< math::Vector3N, std::vector< math::Vector3N> >( iter, trans)()));
        return m_Surf;
    }
    */
    std::vector< float> PointSurfaceObject::GetSizes() const
    { return m_Object->GetSizes();}

    void PointSurfaceObject::PushBackSurfacePoint( const boost::shared_ptr< SurfacePoint> &POINT)
    { m_Surf.PushBack( POINT);}

    std::ostream& PointSurfaceObject::Write( std::ostream &STREAM) const
    {
        STREAM << *m_Object;
        STREAM << m_Surf;
        return STREAM;
    }

    std::ostream &PointSurfaceObject::WriteVmdCommands( std::ostream &STREAM) const
    {
        m_Object->WriteVmdCommands( STREAM);
        m_Surf.WriteVmdCommands( STREAM);
        return STREAM;
    }


    std::istream& PointSurfaceObject::Read( std::istream &STREAM)
    {
        // reading differen child classes by class name switch MISSING!!!!
        STREAM >> m_Surf;
        return STREAM;
    }


} // end namespace geom


/*
    PointSurfaceObject *PointSurfaceObject::Clone() const{ return new PointSurfaceObject( *this);}

//    const math::Vector3N &PointSurfaceObject::GetPosition() const{ return math::Vector3N( *this);}

//    void PointSurfaceObject::SetPosition( const math::Vector3N &NEWPOS){ *this = NEWPOS;} // PointSurfaceObject::Object::math::Vector3N::operator = ( NEWPOS);}


        const boost::shared_ptr< PointSurface> &PointSurfaceObject::GetSurface() const
        { return m_Surf;}

     const boost::shared_ptr< PointSurface> &PointSurfaceObject::CalculateSurface( const float &RESOLUTION) const
        {
            boost::shared_ptr< math::IterationVector< float> > iter( PointSurfaceObject::Object::GetIterationVector());
            boost::shared_ptr< math::Function< math::Vector3N, math::Vector3N> > trans( PointSurfaceObject::Object::GetCoordinateTransformation());
//            m_Surf = boost::shared_ptr< PointSurface>( new PointSurface( math::Integrator< math::Vector3N, std::vector< math::Vector3N> >( iter, trans)()));
            return m_Surf;
        }

     std::ostream& PointSurfaceObject::Write( std::ostream &STREAM) const
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }

     std::istream& PointSurfaceObject::Read( std::istream &STREAM)
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }


} // end namespace geom
*/
