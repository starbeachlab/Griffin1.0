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


#include "../../include/storage/limits3D.h"



namespace store
{


    //!  copy constructor
     Limits3D *Limits3D::Clone() const{ return new Limits3D( *this);}

    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////

     const float &Limits3D::GetXMin() const
    {
        return m_XMin;
    }

     void Limits3D::SetXMin( const float &XMAX)
    {
        m_XMin = XMAX;
    }

     const float &Limits3D::GetXMax() const
    {
        return m_XMax;
    }

     void Limits3D::SetXMax( const float &XMAX)
    {
        m_XMax = XMAX;
    }

     const float &Limits3D::GetYMin() const
    {
        return m_YMin;
    }

     void Limits3D::SetYMin( const float &YMIN)
    {
        m_YMin = YMIN;
    }

     const float &Limits3D::GetYMax() const
    {
        return m_YMax;
    }

     void Limits3D::SetYMax( const float &YMAX)
    {
        m_YMax = YMAX;
    }

     const float &Limits3D::GetZMin() const
    {
        return m_ZMin;
    }

     void Limits3D::SetZMin( const float &ZMIN)
    {
        m_ZMin = ZMIN;
    }

     const float &Limits3D::GetZMax() const
    {
        return m_ZMax;
    }

     void Limits3D::SetZMax( const float &ZMAX)
    {
        m_ZMax = ZMAX;
    }

     float Limits3D::XDelta() const
    {
        return m_XMax - m_XMin;
    }

     float Limits3D::YDelta() const
    {
        return m_YMax - m_YMin;
    }

     float Limits3D::ZDelta() const
    {
        return m_ZMax - m_ZMin;
    }

     float Limits3D::GetMin( const size_t ID) const
     {
         switch( ID)
         {
         case 0:
             return m_XMin;
         case 1:
             return m_YMin;
         case 2:
             return m_ZMin;
         default:
             StandardWrite( "WRONG ID; " << ID << " (should be 0..2)");
             exit( -1);
         }
         return float();
     }

     math::Vector3N Limits3D::GetMin() const
    {
         return math::Vector3N( m_XMin, m_YMin, m_ZMin);
    }


     float Limits3D::GetMax( const size_t ID) const
     {
         switch( ID)
         {
         case 0:
             return m_XMax;
         case 1:
             return m_YMax;
         case 2:
             return m_ZMax;
         default:
             StandardWrite( "WRONG ID; " << ID << " (should be 0..2)");
             exit( -1);
         }
         return float();
     }


     math::Vector3N Limits3D::GetMax() const
    {
         return math::Vector3N( m_XMax, m_YMax, m_ZMax);
    }


     float Limits3D::GetDelta( const size_t ID) const
     {
         switch( ID)
         {
         case 0:
             return m_XMax - m_XMin;
         case 1:
             return m_YMax - m_YMin;
         case 2:
             return m_ZMax - m_ZMin;
         default:
             StandardWrite( "WRONG ID; " << ID << " (should be 0..2)");
             exit( -1);
         }
         return float();
     }

     void Limits3D::SetMin( const size_t ID, const float &VALUE)
     {
         switch( ID)
         {
         case 0:
             m_XMin = VALUE;
             break;
         case 1:
             m_YMin = VALUE;
             break;
         case 2:
             m_ZMin = VALUE;
             break;
         default:
             StandardWrite( "WRONG ID; " << ID << " (should be 0..2)");
             exit( -1);
         }
     }

     void Limits3D::SetMax( const size_t ID, const float &VALUE)
     {
         switch( ID)
         {
         case 0:
             m_XMax = VALUE;
             break;
         case 1:
             m_YMax = VALUE;
             break;
         case 2:
             m_ZMax = VALUE;
             break;
         default:
             StandardWrite( "WRONG ID; " << ID << " (should be 0..2)");
             exit( -1);
         }
     }




     void Limits3D::AdjustToVoxelSize( const float &VOXELSIZE)
    {
        float fract, rest, offset;
        fract = XDelta() / VOXELSIZE;
        rest = fabs( fract - int( fract));
        if( rest > 1e-8)
        {
            offset = 0.5 * ( 1.0 - rest) * VOXELSIZE;
            m_XMin -= offset;
            m_XMax += offset;
        }

        fract = YDelta() / VOXELSIZE;
        rest = fabs( fract - int( fract));
        if( rest > 1e-8)
        {
            offset = 0.5 * ( 1.0 - rest) * VOXELSIZE;
            m_YMin -= offset;
            m_YMax += offset;
        }

        fract = ZDelta() / VOXELSIZE;
        rest = fabs( fract - int( fract));
        if( rest > 1e-8)
        {
            offset = 0.5 * ( 1.0 - rest) * VOXELSIZE;
            m_ZMin -= offset;
            m_ZMax += offset;
        }
    }

     bool Limits3D::AreLimitsInAgreementWithVoxelSize( const float &VOXELSIZE, const float &PRECISION) const
    {
        float fract;
        fract = XDelta() / VOXELSIZE;
        if( fabs( fract - int( fract + PRECISION)) > PRECISION)
        {
            DebugWrite( "deviation: " << fabs( fract - int( fract + PRECISION)) << " (" << fract << ", " << int( fract + PRECISION) << ")");
            return false;
        }

        fract = YDelta() / VOXELSIZE;
        if( fabs( fract - int( fract + PRECISION)) > PRECISION)
        {
            DebugWrite( "deviation: " << fabs( fract - int( fract + PRECISION)) << " (" << fract << ", " << int( fract + PRECISION) << ")");
            return false;
        }

        fract = ZDelta() / VOXELSIZE;
        if( fabs( fract - int( fract + PRECISION)) > PRECISION)
        {
            DebugWrite( "deviation: " << fabs( fract - int( fract + PRECISION)) << " (" << fract << ", " << int( fract + PRECISION) << ")");
            return false;
        }

        return true;
    }


     Limits3D &Limits3D::Merge( const Limits3D &LIMITS)
    {
        if( LIMITS.m_XMin < m_XMin)
        {
            m_XMin = LIMITS.m_XMin;
        }
        if( LIMITS.m_XMax > m_XMax)
        {
            m_XMax = LIMITS.m_XMax;
        }
        if( LIMITS.m_YMin < m_YMin)
        {
            m_YMin = LIMITS.m_YMin;
        }
        if( LIMITS.m_YMax > m_YMax)
        {
            m_YMax = LIMITS.m_YMax;
        }
        if( LIMITS.m_ZMin < m_ZMin)
        {
            m_ZMin = LIMITS.m_ZMin;
        }
        if( LIMITS.m_ZMax > m_ZMax)
        {
            m_ZMax = LIMITS.m_ZMax;
        }
        return *this;
    }


    bool Limits3D::IsWithin( const math::Vector3N &POSITION) const
    {
        if( POSITION( 0) < m_XMin || POSITION( 0) > m_XMax)
        {
            return false;
        }

        if( POSITION( 1) < m_YMin || POSITION( 1) > m_YMax)
        {
            return false;
        }

        if( POSITION( 2) < m_ZMin || POSITION( 2) > m_ZMax)
        {
            return false;
        }
        return true;
    }


    float Limits3D::Volume() const
    {
        return fabs( XDelta() * YDelta() * ZDelta());
    }

    bool Limits3D::AreDefined() const
    {
        return
            m_XMin != std::numeric_limits< float>::max()
            && m_XMax != -std::numeric_limits< float>::max()
            && m_YMin != std::numeric_limits< float>::max()
            && m_YMax != -std::numeric_limits< float>::max()
            && m_ZMin != std::numeric_limits< float>::max()
            && m_ZMax != -std::numeric_limits< float>::max();
    }

    /////////////////////////
    //      Read/Write     //
    /////////////////////////

     std::istream &Limits3D::Read( std::istream &STREAM)
    {
      //        std::string str;
        //        STREAM >> str;
      //        assert( str == GetClassName());
        STREAM >> m_XMin >> m_XMax;
        STREAM >> m_YMin >> m_YMax;
        STREAM >> m_ZMin >> m_ZMax;
        return STREAM;
    }

     std::ostream &Limits3D::Write( std::ostream &STREAM) const
    {
      //        STREAM << GetClassName() << std::endl;
        STREAM << m_XMin << "  " << m_XMax << std::endl;
        STREAM << m_YMin << "  " << m_YMax << std::endl;
        STREAM << m_ZMin << "  " << m_ZMax << std::endl;
        return STREAM;
    }

     std::string Limits3D::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace store



