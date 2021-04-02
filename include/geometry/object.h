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


#ifndef GEOM_OBJECT_H
#define GEOM_OBJECT_H

#include <vector>

#include "../math/vector3N.h"
#include "../math/function.t.h"
#include "../math/iteration_vector.t.h"
#include "../readwrite/stream_operator.h"
#include "../readwrite/class_name.h"
#include "../storage/limits3D.h"

namespace geom
{

    class Object
    : public StreamOperator//, public readwrite::ClassName
    {
    protected:
        math::Vector3N     m_Position;

    public:
        Object()
        : /*StreamOperator(), readwrite::ClassName(),*/ m_Position( 0.0, 0.0, 0.0)
        {}

        Object( const math::Vector3N & POS)
        : /*StreamOperator(), readwrite::ClassName(),*/ m_Position( POS)
        {}

        Object( const Object &OBJECT)
        : /*StreamOperator(), readwrite::ClassName(),*/ m_Position( OBJECT.m_Position)
        {}

        virtual ~Object(){}

        virtual Object *Clone() const = 0; //{ return new Object( *this);}

        virtual float MaxRadius() const{ std::cout << "Object::MaxRadius??" << std::endl; exit( -1); return float();}

        virtual const math::Vector3N& GetPosition() const{ return m_Position;}

        virtual void SetPosition( const math::Vector3N &POS){ m_Position = POS;}

        virtual float GetTotalSurface() const{ return std::numeric_limits< float>::max();}

        virtual store::Limits3D BoxLimits() const
        {
            return store::Limits3D();
        }

        virtual float GetRadius() const{ return std::numeric_limits< float>::max();}

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const{ return STREAM;}

        virtual bool IsPointWithin( const math::Vector3N &POS) const{ return false;}

        virtual float ClosestDistance( const math::Vector3N) const{ return std::numeric_limits< float>::max();}

        //! shortest distance between two objects - probably the slow but rock solid Monte Carlo stuff / maybe exploiting surface object?
        virtual float ClosestDistance( const Object &OBJECT) const{  return std::numeric_limits< float>::max();}

        virtual math::Vector3N ProjectionOnSurface( const math::Vector3N &POS) const
        { return math::Vector3N();}

//    virtual bool IsPointOnSurface( const math::Vector3N) const = 0;

        virtual boost::shared_ptr< math::Function< math::Vector3N, math::Vector3N> > GetCoordinateTransformation() const
        { return boost::shared_ptr< math::Function< math::Vector3N, math::Vector3N> >();}

        virtual boost::shared_ptr< math::IterationVector< float> > GetIterationVector() const
        { return boost::shared_ptr< math::IterationVector< float> >();}
//        {
//            static bool new_iterator( true);
//            static IterationVector iter; // build iteration vector only once and access it for all objects
//            if( new_iterator)
//            {
//                new_iterator = false;
//
//            }
//        }

        virtual std::vector< float> GetSizes() const{ std::cout << "Object::GetSizes()" << std::endl; return std::vector< float>();}
        virtual std::vector< math::Vector3N> GetAxes() const{ std::cout << "Object::GetAxes()" << std::endl; return std::vector< math::Vector3N>();}
        virtual math::Vector3N GetAxes( const size_t &ID) const{ std::cout << "Object::GetAxes( ID)" << std::endl; return math::Vector3N();}

    }; // end class Object

} // end namespace geom

#endif
