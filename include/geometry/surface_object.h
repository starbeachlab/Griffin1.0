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



namespace geom
{
    class PointSurfaceObject
    : public Object
    {
    protected:
        boost::shared_ptr< PointSurface>  m_Surf;

    public:
        PointSurfaceObject()
        : Object(),
        m_Surf()
        {}

        PointSurfaceObject( const Object &OBJECT)
        : Object( OBJECT),
        m_Surf()
        {}

        PointSurfaceObject( const Object &OBJECT, const boost::shared_ptr< PointSurface> &SURF)
        : Object( OBJECT),
        m_Surf( SURF)
        {}

        PointSurfaceObject( const PointSurfaceObject &PSO)
        : Object( PSO),
        m_Surf( PSO.m_Surf)
        {}

        virtual PointSurfaceObject *Clone() const{ return new PointSurfaceObject( *this);}

        boost::shared_ptr< PointSurface> &GetSurface() const;

        boost::shared_ptr< PointSurface> &CalculateSurface( const float &RESOLUTION) const
        {
            boost::shared_ptr< math::IterationVector> iter( Object::GetIterationVector());
            boost::shared_ptr< Function< math::Vector3N, math::Vector3N> > trans( Object::GetCoordinateTransformation());
            m_Surf = boost::shared_ptr< PointSurface>( new PointSurface( math::Integrator< math::Vector3N, std::vector< math::Vector3N> >( iter, trans)()));
            return m_Surf;
        }

        virtual std::ostream& Write( std::ostream &STREAM) const
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }

        virtual std::istream& Read( std::istream &STREAM)
        {
//            STREAM << m_Position[0] << "  " << m_Position[1] << "  " << m_Position[2]<< "  " << m_NormalVector[0] << "  " << m_NormalVector[1] << "  " << m_NormalVector[2] << std::endl;
            return STREAM;
        }


    }; // end class
} // end namespace geom

#endif /* SURF_OBJECT_H_ */
