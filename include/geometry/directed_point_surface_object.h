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


#ifndef DIRECTED_SURF_OBJECT_H_
#define DIRECTED_SURF_OBJECT_H_

#include "object.h"
#include "directed_point_surface.h"

#include "../math/integrator.t.h"

namespace geom
{
    class DirectedPointSurfaceObject
//    : public PointSurface
    {
    protected:
        DirectedPointSurface                     m_Surf;
        boost::shared_ptr< Object>               m_Object;

    public:
        DirectedPointSurfaceObject();

        DirectedPointSurfaceObject( const DirectedPointSurfaceObject &OBJECT);

        DirectedPointSurfaceObject( const boost::shared_ptr< Object> &OBJECT, const DirectedPointSurface &SURF = DirectedPointSurface());

        virtual ~DirectedPointSurfaceObject();

        virtual DirectedPointSurfaceObject *Clone() const;

        virtual const math::Vector3N &GetPosition() const;

        virtual void SetPosition( const math::Vector3N &NEWPOS);

        virtual bool IsPointWithin( const math::Vector3N &POS) const;

        virtual float Distance( const math::Vector3N &POS) const;

        virtual math::Vector3N ProjectionOnSurface( const math::Vector3N &POS) const;

        virtual const float GetRadius() const;

        virtual const store::ShPtrVec< DirectedSurfacePoint> &GetDirectedSurfacePoints() const;

        virtual store::ShPtrVec< DirectedSurfacePoint> &DirectedSurfacePoints();

        virtual size_t NumberSurfacePoints() const;

        virtual const DirectedPointSurface &GetSurface() const;

        virtual void SetSurface( const DirectedPointSurface &SURF);

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

        virtual std::vector< math::Vector3N> GetAxes() const;

        virtual void PushBackDirectedSurfacePoint( const boost::shared_ptr< DirectedSurfacePoint> &POINT);

        virtual std::ostream& Write( std::ostream &STREAM) const;

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const;


        virtual std::istream& Read( std::istream &STREAM);


    }; // end class DirectedPointSurfaceObject

    inline
    std::ostream &operator << ( std::ostream &STREAM, const DirectedPointSurfaceObject &SPHERE)
    { return SPHERE.Write( STREAM);}

} // end namespace geom

#endif /* SURF_OBJECT_H_ */
