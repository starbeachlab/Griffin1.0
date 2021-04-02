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


#include "../../include/utilities/class_id.h"

namespace util
{


    std::map< std::string, size_t> ClassID::s_TotalNumbers = std::map< std::string, size_t>();
    std::map< std::string, int> ClassID::s_ActualNumbers = std::map< std::string, int>();

    ClassID::ClassID()
    : m_ThisID()
    {
//        InsertNewAndIncrementCountForKnownClasses();
//        m_ThisID = s_TotalNumbers[ mystr::GetClassName( t_DATA()::__PRETTY_FUNCTION__)];
    }

    //ClassID::ClassID( const ClassID &ORIG)
    //: m_ThisID( ++s_TotalNumber)
    //{
    //    ++s_ActualNumber;
    //}

    ClassID::~ClassID()
    {
//        --s_ActualNumbers[].second;
//        assert( s_ActualNumber >= 0);
    }

    const size_t &ClassID::GetClassID() const{ return m_ThisID;}

//    int ClassID::GetActualNumber() const{ return s_ActualNumbers;}
} // end namespace util
