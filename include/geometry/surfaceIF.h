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


#ifndef SURFIF_H_
#define SURFIF_H_

#include "../readwrite/stream_operator.h"
#include "../readwrite/class_name.h"

namespace geom
{
    class SurfaceIF
    : public StreamOperator, readwrite::ClassName
    {
    protected:
        float     m_Size;   //!< the surface area

    public:

        virtual SurfaceIF *Clone() const = 0;

        float GetSize() const { return m_Size;}


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

#endif /* SURFIF_H_ */
