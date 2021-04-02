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


#include "../../include/geometry/directed_point_surface_object.h"

namespace geom
{

    DirectedPointSurfaceObject::DirectedPointSurfaceObject()
    : m_Surf(),
    m_Object()
    { /*std::cout << "PointSurfObj default constructor" << std::endl;*/}

    DirectedPointSurfaceObject::DirectedPointSurfaceObject( const DirectedPointSurfaceObject &OBJECT)
    : m_Surf( OBJECT.m_Surf),
    m_Object( OBJECT.m_Object)
    { /*std::cout << "PointSurfObj copy constructor" << std::endl;*/}

    DirectedPointSurfaceObject::DirectedPointSurfaceObject( const boost::shared_ptr< Object> &OBJECT, const DirectedPointSurface &SURF)
    : m_Surf( SURF),
    m_Object( OBJECT)
    {
        DebugWrite( __FUNCTION__);
    }

    DirectedPointSurfaceObject::~DirectedPointSurfaceObject(){}

    DirectedPointSurfaceObject *DirectedPointSurfaceObject::Clone() const
    { return new DirectedPointSurfaceObject( *this);}

    const math::Vector3N &DirectedPointSurfaceObject::GetPosition() const
    { return m_Object->GetPosition();}

    void DirectedPointSurfaceObject::SetPosition( const math::Vector3N &NEWPOS)
    { m_Object->SetPosition( NEWPOS);}

    bool DirectedPointSurfaceObject::IsPointWithin( const math::Vector3N &POS) const
    {
    //            return false;
        return m_Object->IsPointWithin( POS);
    }

    float DirectedPointSurfaceObject::Distance( const math::Vector3N &POS) const
    {
        return m_Object->ClosestDistance( POS);
    }

    math::Vector3N DirectedPointSurfaceObject::ProjectionOnSurface( const math::Vector3N &POS) const
    {
        return m_Object->ProjectionOnSurface( POS);
    }

    const float DirectedPointSurfaceObject::GetRadius() const
    { return m_Object->GetRadius();}

    const store::ShPtrVec< DirectedSurfacePoint> &DirectedPointSurfaceObject::GetDirectedSurfacePoints() const
    { return m_Surf.GetData();}

    store::ShPtrVec< DirectedSurfacePoint> &DirectedPointSurfaceObject::DirectedSurfacePoints()
    { return m_Surf.Data();}

    const DirectedPointSurface &DirectedPointSurfaceObject::GetSurface() const
    { return m_Surf;}

    void DirectedPointSurfaceObject::SetSurface( const DirectedPointSurface &SURF)
    { m_Surf = SURF;}

    const boost::shared_ptr< Object> &DirectedPointSurfaceObject::GetObject() const
    { return m_Object;}

    float DirectedPointSurfaceObject::GetTotalSurface() const
    { return m_Object->GetTotalSurface();}

    float DirectedPointSurfaceObject::GetFreeSurface() const
    { return m_Surf.GetFreeSurface();}

    std::string DirectedPointSurfaceObject::GetClassName() const
    { return mystr::GetClassName( __PRETTY_FUNCTION__);}

    /*
    *
    *         virtual const boost::shared_ptr< DirectedPointSurface> &CalculateSurface( const float &RESOLUTION) const;
    {
        boost::shared_ptr< math::IterationVector> iter( Object::GetIterationVector());
        boost::shared_ptr< Function< math::Vector3N, math::Vector3N> > trans( Object::GetCoordinateTransformation());
        m_Surf = boost::shared_ptr< DirectedPointSurface>( new DirectedPointSurface( math::Integrator< math::Vector3N, std::vector< math::Vector3N> >( iter, trans)()));
        return m_Surf;
    }
    */
    std::vector< float> DirectedPointSurfaceObject::GetSizes() const
    { return m_Object->GetSizes();}

    std::vector< math::Vector3N> DirectedPointSurfaceObject::GetAxes() const
    { return m_Object->GetAxes();}


    void DirectedPointSurfaceObject::PushBackDirectedSurfacePoint( const boost::shared_ptr< DirectedSurfacePoint> &POINT)
    { m_Surf.PushBack( POINT);}

    std::ostream& DirectedPointSurfaceObject::Write( std::ostream &STREAM) const
    {
        STREAM << GetClassName() << std::endl;
        STREAM << "object: " << *m_Object << std::endl;
        STREAM << "surf: " << m_Surf << std::endl;
        return STREAM;
    }

    std::ostream &DirectedPointSurfaceObject::WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID) const
    {
        m_Object->WriteVmdCommands( STREAM, MOL_ID);
        m_Surf.WriteVmdCommands( STREAM, MOL_ID);
        return STREAM;
    }


    std::istream& DirectedPointSurfaceObject::Read( std::istream &STREAM)
    {
        // reading differen child classes by class name switch MISSING!!!!
        STREAM >> m_Surf;
        return STREAM;
    }

    size_t DirectedPointSurfaceObject::NumberSurfacePoints() const
    {
        return m_Surf.GetData().size();
    }



} // end namespace geom


/*
    DirectedPointSurfaceObject *DirectedPointSurfaceObject::Clone() const{ return new DirectedPointSurfaceObject( *this);}

//    const math::Vector3N &DirectedPointSurfaceObject::GetPosition() const{ return math::Vector3N( *this);}

//    void DirectedPointSurfaceObject::SetPosition( const math::Vector3N &NEWPOS){ *this = NEWPOS;} // DirectedPointSurfaceObject::Object::math::Vector3N::operator = ( NEWPOS);}


        const boost::shared_ptr< DirectedPointSurface> &DirectedPointSurfaceObject::GetSurface() const
        { return m_Surf;}

     const boost::shared_ptr< DirectedPointSurface> &DirectedPointSurfaceObject::CalculateSurface( const float &RESOLUTION) const
        {
            boost::shared_ptr< math::IterationVector< float> > iter( DirectedPointSurfaceObject::Object::GetIterationVector());
            boost::shared_ptr< math::Function< math::Vector3N, math::Vector3N> > trans( DirectedPointSurfaceObject::Object::GetCoordinateTransformation());
//            m_Surf = boost::shared_ptr< DirectedPointSurface>( new DirectedPointSurface( math::Integrator< math::Vector3N, std::vector< math::Vector3N> >( iter, trans)()));
            return m_Surf;
        }

     std::ostream& DirectedPointSurfaceObject::Write( std::ostream &STREAM) const
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }

     std::istream& DirectedPointSurfaceObject::Read( std::istream &STREAM)
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }


} // end namespace geom
*/
