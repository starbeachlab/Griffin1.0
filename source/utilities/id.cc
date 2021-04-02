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


#include "../../include/utilities/id.h"

namespace util
{

    size_t ID::s_TotalNumber = 0;
    int ID::s_ActualNumber = 0;

    ID::ID()
    : m_ThisID( ++s_TotalNumber)
    {
        ++s_ActualNumber;
    }

    ID::ID( const ID &ORIG)
    : m_ThisID( ++s_TotalNumber)
    {
        ++s_ActualNumber;
    }

    ID::~ID()
    {
        --s_ActualNumber;
        assert( s_ActualNumber >= 0);
    }

    const size_t &ID::GetID() const{ return m_ThisID;}

    void ID::SetID( const size_t &NEW_ID){ m_ThisID = NEW_ID;}

    const int &ID::GetActualNumber() const{ return s_ActualNumber;}
} // end namespace util
