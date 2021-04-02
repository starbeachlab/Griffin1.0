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


#ifndef SURF_OBJECT_H_
#define SURF_OBJECT_H_

#include "object.h"
#include "point_surface.h"

#include "../math/integrator.t.h"

namespace geom
{
    class PointSurfaceObject
//    : public PointSurface
    {
    protected:
        PointSurface                m_Surf;
        boost::shared_ptr< Object>  m_Object;

    public:
        PointSurfaceObject();

        PointSurfaceObject( const PointSurfaceObject &OBJECT);

        PointSurfaceObject( const boost::shared_ptr< Object> &OBJECT, const PointSurface &SURF = PointSurface());

        virtual ~PointSurfaceObject();

        virtual PointSurfaceObject *Clone() const;

        virtual const math::Vector3N &GetPosition() const;

        virtual void SetPosition( const math::Vector3N &NEWPOS);

        virtual bool IsPointWithin( const math::Vector3N &POS) const;

        virtual const float GetRadius() const;

        virtual const store::ShPtrVec< SurfacePoint> &GetSurfacePoints() const;

        virtual store::ShPtrVec< SurfacePoint> &SurfacePoints();

        virtual const PointSurface &GetSurface() const;

        virtual void SetSurface( const PointSurface &SURF);

        virtual const boost::shared_ptr< Object> &GetObject() const;

        virtual float GetTotalSurface() const;

        virtual float GetFreeSurface() const;

        virtual std::string GetClassName() const;

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
        virtual std::vector< float> GetSizes() const;

        virtual void PushBackSurfacePoint( const boost::shared_ptr< SurfacePoint> &POINT);

        virtual std::ostream& Write( std::ostream &STREAM) const;

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const;


        virtual std::istream& Read( std::istream &STREAM);


    }; // end class PointSurfaceObject

    inline
    std::ostream &operator << ( std::ostream &STREAM, const PointSurfaceObject &SPHERE)
    { return SPHERE.Write( STREAM);}

} // end namespace geom

#endif /* SURF_OBJECT_H_ */
