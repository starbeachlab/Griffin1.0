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


#ifndef VECTOR_ANGLE_H_
#define VECTOR_ANGLE_H_

#include "function.t.h"
#include "vector.h"
#include "double_functions.h"
#include "vector_functions.h"

#ifdef ANALYZE_BURIAL
#include "histogram.t.h"
#include "double_functions.h"
#endif


namespace math
{

	class VectorAngle
	: public Function< std::pair< math::Vector, math::Vector>, float>
	{
	private:
		float m_CosAngle;
		float m_AngleThreshold;

#ifdef ANALYZE_BURIAL
        mutable math::Histogram< float>  m_Angles;
        int 				     m_Count;
#endif


    public:
        VectorAngle()
	    : m_CosAngle(),
	    m_AngleThreshold()
#ifdef ANALYZE_BURIAL
          , m_Angles( 0.0, 5.0, 36), m_Count( 0)
#endif
        {}

        ~VectorAngle()
        {
#ifdef ANALYZE_BURIAL
        	std::ofstream out( "hist/angle_histogram.txt");
        	out << m_Angles;
        	out.close();
        	out.clear();
        }

        void Reset()
        {
        	std::stringstream strstr;
        	strstr << m_Count;
        	std::string str;
			strstr >> str;
			std::ofstream out( ("hist/angle_histogram" + str + ".txt").c_str());
			out << m_Angles;
			out.close();
			out.clear();
			++m_Count;
			m_Angles.Clear();
#endif
        }

        VectorAngle( const float &THRESHOLD_GRAD);

        const float &GetCosThreshold() const;

//        void ResetLength( const float &LENGTH);

	void ResetAngle( const float &ANGLE_IN_DEGREES);

        float  & operator()( const std::pair< math::Vector, math::Vector> &VECTORS) const;

        bool IsLargerThreshold( const std::pair< math::Vector, math::Vector> &VECTORS) const;

        bool IsAngleLargerThreshold( const math::Vector &FIRST, const math::Vector &SECOND) const;

        std::ostream &Write( std::ostream &STREAM) const;

        std::istream &Read( std::istream &STREAM);

    private:
        float CalcFactor() const;

	}; // end class VectorAngle

	bool operator == ( const VectorAngle &A, const VectorAngle &B);


} // namespace math


#endif /* VECTOR_ANGLE_H_ */
