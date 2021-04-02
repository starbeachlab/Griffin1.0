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


#include "../../include/string/io_string_functions.h"

namespace mystr
{

    std::string GetFunctionName( const std::string &STRING)
    {
        std::vector< std::string> parts( SplitString( STRING));

        std::vector< std::string>::iterator itr = parts.begin(), first = itr;
        while( itr->find( ")") == std::string::npos)
        {
            if( itr->find( "<") != std::string::npos)
            { first = itr;}
            ++itr;
        }
        if( first == parts.begin())
        { first = itr;}

        std::string new_str;
        for( std::vector< std::string>::iterator tr = first; tr != itr + 1; ++tr)
        { new_str += *tr;}

//        std::cout << std::endl;
//        // replace template names if newstr contains '<' and '>'!!
//        size_t begin = new_str.find( "<");
//        if( begin != std::string::npos)
//        {
//            std::cout << STRING << std::endl;
//            size_t last = new_str.find( ">");
//            ++begin;
//            std::cout << "substr: <" << new_str.substr( begin, last - begin) << ">" << std::endl;
//            std::vector< std::string> templates( SplitString( new_str.substr( begin, last - begin), ","));
//            std::cout << templates;
//            std::vector< std::string> replace( templates.size());
//            std::vector< std::string>::iterator rtr = replace.begin();
//            for( std::vector< std::string>::const_iterator itr = templates.begin(); itr != templates.end(); ++itr, ++rtr)
//            {
//                begin = STRING.rfind( *itr);
////                std::cout << STRING.substr( begin, STRING.size() - begin) << std::endl;
//                begin = STRING.find( '=', begin);
//                begin = STRING.find_first_not_of( " ", begin+1);
//                if( itr + 1 != templates.end())
//                {
//                    last = STRING.rfind( *( itr + 1)) - 2;
//                }
//                else
//                {
//                    last = STRING.rfind( "]");
//                }
//                *rtr = STRING.substr( begin, last - begin);
//            }
////            replace.back() = replace.back().substr( 0, replace.back().size() - 1);
//            std::cout << "replace: " << replace;
//        }
//        // replace
        return new_str;
    }


    std::string GetClassName( const std::string &STRING)
    {
        std::string str( GetFunctionName( STRING));
//        std::cout << str << std::endl;
        str = str.substr( 0, str.find_first_of( "("));

        std::vector< std::string> vec( SplitString( str, "::"));

        std::string final;

        for( size_t i = 0; i < vec.size() - 1; ++i)
        {
            final += vec[i];
            if( i < vec.size() - 2)
            { final += "::";}
        }
        return final;
    }



} // end namespace mystr
