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


#ifndef ELEMENT_DATA_MAP_H_
#define ELEMENT_DATA_MAP_H_

#include <map>
#include <iostream>
#include <fstream>

#include "element_data.h"
#include "../math/function.t.h"

namespace mol
{

    class ElementDataMap
    : public math::Function< std::string, boost::shared_ptr< ElementData> >
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        static const std::map< std::string, boost::shared_ptr< ElementData> > s_Map;
        static const boost::shared_ptr< ElementData>                          s_Empty;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ElementDataMap()
        : math::Function< std::string, boost::shared_ptr< ElementData> >()
        {}

        //! copy constructor
        ElementDataMap( const ElementDataMap &ORIGINAL)
        : math::Function< std::string, boost::shared_ptr< ElementData> >( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~ElementDataMap(){}

        //! virtual copy constructor
        virtual ElementDataMap *Clone() const{ return new ElementDataMap( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        virtual boost::shared_ptr< ElementData> & operator () ( const std::string &STRING) const;

        virtual bool IsValidElementType( const std::string &STRING) const{ return s_Map.find( STRING) != s_Map.end();}

    protected:
        static std::map< std::string, boost::shared_ptr< ElementData> > BuildMap();

    public:
        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

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

    }; // end class
} // end namespace $




#endif /* ELEMENT_DATA_MAP_H_ */
