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


#ifndef CLASS_NAME_H
#define CLASS_NAME_H

#include "../string/io_string_functions.h"

namespace readwrite
{

    class ClassName
    {
    public:
        virtual ~ClassName(){}

        virtual std::string GetClassName() const
        { return mystr::GetClassName( __PRETTY_FUNCTION__);}

        virtual std::ostream & WriteClassName( std::ostream &STREAM) const
        {
            STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
            return STREAM;
        }

    }; // end class ClassName
} // end namespace readwrite

#endif
