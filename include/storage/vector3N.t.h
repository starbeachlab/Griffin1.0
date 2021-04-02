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


#ifndef STORE_VECTOR3N_T_H
#define STORE_VECTOR3N_T_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>


namespace store
{

    template< typename t_DATA>
    class Vector3N
    : public std::vector< t_DATA>
    {
    public:
        Vector3N()
        : std::vector<t_DATA>(3)
        {};

        Vector3N( const Vector3N &V)
        : std::vector< t_DATA>( V)
        {};

        Vector3N( const std::vector< t_DATA> &V)
        : std::vector< t_DATA>( V)
        { assert( V.size() == 3);};


//	//! constructor that enables implicit conversions (if data types are convertible)
//	template< typename t_OTHER>
//        Vector3N( const std::vector< t_OTHER> &V)
//        : std::vector< t_DATA>( V.size())
//        {
//	    typename std::vector< t_OTHER>::const_iterator other_itr = V.begin();
//	    typename std::vector< t_DATA>::iterator this_itr = this->begin();
//	    for( ; other_itr != V.end(); ++other_itr, ++this_itr)
//	    {
//		*this_itr = t_DATA( *other_itr);
//	    }
//	};

        Vector3N( const t_DATA &X, const t_DATA &Y, const t_DATA &Z)
        : std::vector< t_DATA>( 3)
        {
            ( *this)[0] = X;
            ( *this)[1] = Y;
            ( *this)[2] = Z;
        };

        ~Vector3N(){};


        t_DATA &operator() ( const size_t &ID)
        {
//            assert( ID < 3);
            return ( *this)[ ID];
        }

        const t_DATA &operator() ( const size_t &ID) const
        {
//            assert( ID < 3);
            return ( *this)[ ID];
        }

        std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << __PRETTY_FUNCTION__ << std::endl;
            STREAM << ( *this)[0] << "   " << ( *this)[1] << "   " << ( *this)[2] << std::endl;
            return STREAM;
        }


        std::istream &Read( std::istream &STREAM)
        {
            std::string str;
            STREAM >> str;
            STREAM >> ( *this)[0] >> ( *this)[1] >> ( *this)[2];
            return STREAM;
        }


    }; // end class Vector3N


} // end namespace store


#endif // STORE_VECTOR3N_T_H
