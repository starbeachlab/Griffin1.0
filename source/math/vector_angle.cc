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


#include "../../include/math/vector_angle.h"

namespace math
{
    VectorAngle::VectorAngle( const float &ANGLE)
	: m_CosAngle(),
	  m_AngleThreshold( ANGLE)
#ifdef ANALYZE_BURIAL
	, m_Angles( 0.0, 5.0, 36), 
	  m_Count( 0)
#endif
    {
    	m_CosAngle = CalcFactor();
#ifdef ANALYZE_BURIAL
    	std::cout << __FUNCTION__;
        std::cout << ": threshold in degrees: " << RadiansToDegrees( ANGLE) << " rad: " <<  ANGLE << " cos: " << cos(  ANGLE) << std::endl;
#endif
    }

    const float &VectorAngle::GetCosThreshold() const
	{
    	return m_CosAngle;
	}

//    void VectorAngle::ResetLength( const float &LENGTH)
//    {
//    	StandardWrite( __PRETTY_FUNCTION__ << " length: " << LENGTH);
//    	m_Length = LENGTH;
//    	m_LengthCosAngle = CalcFactor();
//    }

    void VectorAngle::ResetAngle( const float &ANGLE)
    {
		DebugWrite( __PRETTY_FUNCTION__ << " angle: " << RadiansToDegrees( ANGLE));
		m_AngleThreshold = ANGLE;
		m_CosAngle = CalcFactor();
    }



    float & VectorAngle::operator()( const std::pair< math::Vector, math::Vector> &VECTORS) const
    {
        return Function< std::pair< math::Vector, math::Vector>, float>::s_Tmp = Angle( VECTORS.first, VECTORS.second);
    }


    bool VectorAngle::IsLargerThreshold( const std::pair< math::Vector, math::Vector> &VECTORS) const
    {
    	return IsAngleLargerThreshold( VECTORS.first, VECTORS.second);
    }


    bool VectorAngle::IsAngleLargerThreshold( const math::Vector &FIRST, const math::Vector &SECOND) const
    {

#ifdef ANALYZE_BURIAL
    	float angle = RadiansToDegrees( Angle( FIRST, SECOND));
    	m_Angles.InsertValue( angle);
//    	std::cout << __FUNCTION__ << ": angle: " << angle << "first length: " << FIRST.Length() << " second: " << SECOND.Length() << std::endl;
#endif


        if( FIRST * SECOND + 1e-4 < sqrt( FIRST.SquaredLength() * SECOND.SquaredLength()) * m_CosAngle) // sum 1e-4 to avoid rounding error for angle = 180 // todo: fix differently // float problem ?
        {

#ifdef ANALYZE_BURIAL
	    std::cout << __FUNCTION__ << ": angle: " << angle << " is larger than threshold" << std::endl;
#endif

            return true;
        }
        return false;
    }



    std::ostream &VectorAngle::Write( std::ostream &STREAM) const
    {
    	STREAM << "VectorAngle" << std::endl;
    	STREAM << RadiansToDegrees( m_AngleThreshold) << std::endl;
    	return STREAM;
    }


    std::istream &VectorAngle::Read( std::istream &STREAM)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
    	if( !STREAM)
    	{
    		std::cout << "===> dead stream in " << __PRETTY_FUNCTION__ << std::endl;
    		exit( -1);
    	}
    	std::string str;
    	STREAM >> str;
    	assert( str == "VectorAngle");
    	float value;
    	STREAM >> value;
    	m_AngleThreshold = DegreesToRadians( value);
    	m_CosAngle = CalcFactor();
    	return STREAM;
    }


    float VectorAngle::CalcFactor() const
    {
    	return  cos( m_AngleThreshold);
    }


    bool operator == ( const VectorAngle &A, const VectorAngle &B)
	{
		return A.GetCosThreshold() == B.GetCosThreshold();
	}


}
