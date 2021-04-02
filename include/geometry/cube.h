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


#ifndef CUBE_H_
#define CUBE_H_

//#include <>

#include "object3D.h"


namespace geom
{

    class Cube
    : public Object3D
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        std::vector< float> m_Sizes;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Cube()
        : Object3D(),
        m_Sizes( 3)
        {}

        //! construct symmetric cube
        Cube( const float &SIZE)
        : Object3D(),
        m_Sizes( 3)
        {
            m_Sizes[ 0] = SIZE;
            m_Sizes[ 1] = SIZE;
            m_Sizes[ 2] = SIZE;
        }

        //! construct symmetric cube
        Cube( const math::Vector3N &POS, const float &SIZE)
        : Object3D( POS),
        m_Sizes( 3)
        {
            m_Sizes[ 0] = SIZE;
            m_Sizes[ 1] = SIZE;
            m_Sizes[ 2] = SIZE;
        }

        //! construct from size vector
        Cube( const std::vector< float> &SIZES)
        : Object3D(),
        m_Sizes( SIZES)
        { assert( m_Sizes.size() == 3);}

        Cube( const math::Vector3N &POS, const float &XSIZE, const float &YSIZE, const float &ZSIZE)
        : Object3D( POS),
        m_Sizes( 3)
        {
            DebugWrite( __FUNCTION__);
            m_Sizes[0] = XSIZE;
            m_Sizes[1] = YSIZE;
            m_Sizes[2] = ZSIZE;
        }

        //! copy constructor
        Cube( const Cube &ORIGINAL)
        : Object3D( ORIGINAL),
        m_Sizes( ORIGINAL.m_Sizes)
        {}

        //! virtual destructor
        virtual ~Cube(){}

        //! virtual copy constructor
        virtual Cube *Clone() const{ return new Cube( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual std::vector< float> GetSizes() const{ return m_Sizes;}

        virtual float GetTotalSurface() const
        {
            return
                2.0 * (
                        m_Sizes[0] * m_Sizes[1]
                        + m_Sizes[1] * m_Sizes[2]
                        + m_Sizes[0] * m_Sizes[2]
                );
        }

        virtual float MaxRadius() const
        { return std::sqrt( math::Square( m_Sizes[0]) + math::Square( m_Sizes[1]) + math::Square( m_Sizes[2]));}

        virtual store::Limits3D BoxLimits() const
        {
            return store::Limits3D
            (
                    m_Position[0] - m_Sizes[0],
                    m_Position[0] + m_Sizes[0],
                    m_Position[1] - m_Sizes[1],
                    m_Position[1] + m_Sizes[1],
                    m_Position[2] - m_Sizes[2],
                    m_Position[2] + m_Sizes[2]
            );
        }



        virtual bool IsPointWithin( const math::Vector3N &POS) const // ONLY AS LONG NO ROTATIONS ARE PERFORMED!!!!
        {
            math::Vector3N pos( POS - GetPosition());
            if( math::Unsign( pos(0)) <= 0.5 * m_Sizes[0] && math::Unsign( pos( 1)) <= 0.5 * m_Sizes[ 1] && math::Unsign( pos( 2)) <= 0.5 * m_Sizes[ 2])
            { return true;}
            return false;
        }

        virtual float ClosestDistance( const math::Vector3N &POS) const // TODO: QUATSCH!!!
        {
            math::Vector3N
                connector( POS - GetPosition());
            std::map< float, size_t>
                list;


            for( size_t i = 0; i < 3; ++i)
                for( size_t j = 0; j < 2; ++j)
                {
                    list.insert( std::make_pair( math::Angle( connector, pow( -1, j) * GetAxes( i)), i));
                }

            std::map< float, size_t>::const_iterator
                itr( list.begin());

            return connector.Length() - 0.5 * GetSizes()[ itr->second] / cos( itr->first);
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

        virtual std::string GetClassName() const
        { return mystr::GetClassName( __PRETTY_FUNCTION__);}


        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM, const size_t &MOL_ID = 0) const
        {
    //        STREAM << "mol new" << std::endl;
            STREAM << "graphics " << MOL_ID << " color red" <<  std::endl;

            std::vector< math::Vector3N> axes( GetAxes());

            std::vector< float> sizes( GetSizes());

            std::vector< math::Vector3N> points( 8, GetPosition());

            math::Vector3N
                a( 0.5 * sizes[0] * axes[0]),
                b( 0.5 * sizes[1] * axes[1]),
                c( 0.5 * sizes[2] * axes[2]);

            points[ 0] += -1.0 * a - b - c;
            points[ 1] +=        a - b - c;
            points[ 2] +=        a + b - c;
            points[ 3] += -1.0 * a + b - c;
            points[ 4] += -1.0 * a - b + c;
            points[ 5] +=        a - b + c;
            points[ 6] +=        a + b + c;
            points[ 7] += -1.0 * a + b + c;


            // TODO: systematic!
            // the bottom
            STREAM << "graphics " << MOL_ID << " triangle {" << points[0]( 0) << " " << points[0]( 1) << " " << points[0]( 2) << "} {" << points[1]( 0) << " " << points[1]( 1) << " " << points[1]( 2) << "} {" << points[2]( 0) << " " << points[2]( 1) << " " << points[2]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[0]( 0) << " " << points[0]( 1) << " " << points[0]( 2) << "} {" << points[2]( 0) << " " << points[2]( 1) << " " << points[2]( 2) << "} {" << points[3]( 0) << " " << points[3]( 1) << " " << points[3]( 2) << "} " << std::endl;
            // the top
            STREAM << "graphics " << MOL_ID << " triangle {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} {" << points[5]( 0) << " " << points[5]( 1) << " " << points[5]( 2) << "} {" << points[6]( 0) << " " << points[6]( 1) << " " << points[6]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} {" << points[6]( 0) << " " << points[6]( 1) << " " << points[6]( 2) << "} {" << points[7]( 0) << " " << points[7]( 1) << " " << points[7]( 2) << "} " << std::endl;
            // front side
            STREAM << "graphics " << MOL_ID << " triangle {" << points[0]( 0) << " " << points[0]( 1) << " " << points[0]( 2) << "} {" << points[1]( 0) << " " << points[1]( 1) << " " << points[1]( 2) << "} {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} {" << points[1]( 0) << " " << points[1]( 1) << " " << points[1]( 2) << "} {" << points[5]( 0) << " " << points[5]( 1) << " " << points[5]( 2) << "} " << std::endl;
            // right side
            STREAM << "graphics " << MOL_ID << " triangle {" << points[5]( 0) << " " << points[5]( 1) << " " << points[5]( 2) << "} {" << points[1]( 0) << " " << points[1]( 1) << " " << points[1]( 2) << "} {" << points[2]( 0) << " " << points[2]( 1) << " " << points[2]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[5]( 0) << " " << points[5]( 1) << " " << points[5]( 2) << "} {" << points[6]( 0) << " " << points[6]( 1) << " " << points[6]( 2) << "} {" << points[2]( 0) << " " << points[2]( 1) << " " << points[2]( 2) << "} " << std::endl;
            // back side
            STREAM << "graphics " << MOL_ID << " triangle {" << points[2]( 0) << " " << points[2]( 1) << " " << points[2]( 2) << "} {" << points[3]( 0) << " " << points[3]( 1) << " " << points[3]( 2) << "} {" << points[6]( 0) << " " << points[6]( 1) << " " << points[6]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[7]( 0) << " " << points[7]( 1) << " " << points[7]( 2) << "} {" << points[3]( 0) << " " << points[3]( 1) << " " << points[3]( 2) << "} {" << points[6]( 0) << " " << points[6]( 1) << " " << points[6]( 2) << "} " << std::endl;
            // left side
            STREAM << "graphics " << MOL_ID << " triangle {" << points[0]( 0) << " " << points[0]( 1) << " " << points[0]( 2) << "} {" << points[3]( 0) << " " << points[3]( 1) << " " << points[3]( 2) << "} {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} " << std::endl;
            STREAM << "graphics " << MOL_ID << " triangle {" << points[7]( 0) << " " << points[7]( 1) << " " << points[7]( 2) << "} {" << points[3]( 0) << " " << points[3]( 1) << " " << points[3]( 2) << "} {" << points[4]( 0) << " " << points[4]( 1) << " " << points[4]( 2) << "} " << std::endl;


            STREAM << "graphics " << MOL_ID << " sphere {" << GetPosition()(0) << " " << GetPosition()(1) << " " << GetPosition()(2) << "} radius 0.1 resolution 21" << std::endl;
            return STREAM;
        }

    }; // end class
} // end namespace $




#endif /* CUBE_H_ */
