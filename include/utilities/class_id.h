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


#ifndef CLASS_ID_H_
#define CLASS_ID_H_

#include <cassert>
#include <map>
#include <string>

namespace util
{
    class ClassID
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        size_t        m_ThisID;
        static std::map< std::string, int> s_ActualNumbers;
        static std::map< std::string, size_t> s_TotalNumbers;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ClassID();


        //! copy constructor
        ClassID( const ClassID &ORIGINAL);

        //! destructor
        ~ClassID();

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const size_t &GetClassID() const;

        int GetActualNumber() const;

    }; // end class ClassID
} // end namespace util




#endif /* ClassID_H_ */
