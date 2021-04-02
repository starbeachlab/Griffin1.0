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


#ifndef STRING_FUNCTIONS
#define STRING_FUNCTIONS

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>

#include <cctype>
#include <cassert>
#include <cmath>
#include <stdio.h>

#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


namespace mystr
{

    //! translate a char array into a string, removing eventual non alphanumerical endings
    void BufferToFileName( const char BUFFER[], std::string &STR, const int &SIZE);


    //! trims spaces from beginning and end of a copied string
    std::string TrimString( const std::string &STRING);


    //! test whether string is numerical (float, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    bool IsNumerical( const std::string &STRING);


    std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER = " ");

    std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER_A, const std::string &SPLITTER_B);

//    std::string &RemoveSubstr( std::string &STR, const std::string SUBSTR)
//    { exit( -1); return STR;}

    //! converts std::string to object of template class T1
    template< class T1>
    T1 ConvertStringToNumericalValue( const std::string &STRING)
    {
        // remove tabs, spaces, new lines, etc.
        const std::string string( TrimString( STRING));

        // if string is empty an undefined is returned
        if( string.empty())
        { return std::numeric_limits< T1>::max();}

        assert( IsNumerical( string));
        T1 new_t;
        std::stringstream iss( string);
        iss >> new_t;
        return new_t;
    }

    template<class T>
    std::string NumericalValueToString(const T& VALUE)
    {
        std::stringstream strm;
        strm << VALUE;
        return strm.str();
    }

//    template<class T>
//    std::string NumericalValueToString(const T& VALUE, const size_t & TOTAL, const size_t &PRECISION)  // TODO : not 100% reliable!
//    {
//    	std::stringstream strm;
//        strm.width( TOTAL);
//        T value =  T( int( VALUE * pow( 10, PRECISION) + 0.5)) / T( pow( 10, PRECISION)); // TODO: use bitwise operators rather than mathematical operations!!!
//        strm << value;
//        std::string str = strm.str();
////		std::cout << "now: " << str << std::endl;
//		if( str.find(".") == std::string::npos)
//		{
//			str += ".";
//		}
//		int dist =  str.size() - str.find(".") - 1;
//		for( size_t i = 0; i < std::max( 0, int( PRECISION) - dist); ++i)
//		{
//			str += "0";
//		}
//		return str;
//    }

    template<class T>
    std::string NumericalValueToString(const T& VALUE, const int & TOTAL, const int &PRECISION)
    {
    	std::string str( TOTAL, ' ');

    	char
			*nr = new char[4],
			*prec = new char[4],
			*char_ptr = new char[ TOTAL];

    	sprintf( prec, "%d", PRECISION);
    	sprintf( nr, "%d", TOTAL);

    	std::string
			format = "%";
    	format.append( nr);
    	format +=  ".";
    	format.append( prec);
    	format += "f";

    	//      std::cout << "format: " << format << std::endl;

    	int n = sprintf( char_ptr, format.c_str(), VALUE);

    	// the resulting char array might be longer than TOTAL
    	// adjust the precision in that case
    	if( n > TOTAL)
    	{
    		//	  std::cout << "long stuff" << std::endl;
    		int new_size = std::max( 0, PRECISION - ( n - TOTAL));
    		sprintf( prec, "%d", new_size);
    		format = "%";
    		format.append( nr);
    		format +=  ".";
    		format.append( prec);
    		format += "f";
    		//	  std::cout << "new format: " << format << std::endl;

    		sprintf( char_ptr, format.c_str(), VALUE);
    	}

    	// pass char array to correct sized string
    	std::copy( char_ptr, char_ptr + TOTAL, str.begin());

    	return str;
    }


//    template<class T>
//    std::string NumericalValueToString(const T& VALUE, const size_t & TOTAL, const size_t &PRECISION)  // TODO : not 100% reliable!
//    {
//    	std::stringstream strm;
//    	strm.width( TOTAL);
//    	T value = T( int( VALUE * pow( 10, PRECISION))) / T( pow( 10, PRECISION)); // TODO: use bitwise operators rather than mathematical operations!!!
//    	strm << value;
//    	return strm.str();
//    }

//    template<class T>
//    std::string NumericalValueToString(const T& VALUE, const size_t &PRECISION)
//    {
//        std::stringstream strm;
//        strm << std::setprecision( PRECISION);
////        T value( T( size_t( VALUE * pow( 10, PRECISION)) / pow( 10, PRECISION))
//        strm << VALUE;
//        return strm.str();
//    }

    //! transforms given string into lower case
    std::string &AllToLowerCase( std::string &STRING);

    //! returns a lower case copy of string
    std::string ToLower( const std::string &STRING);

    //! transforms given string to upper case
    std::string &AllToUpperCase( std::string &STRING);

    //! returns an upper case copy of string
    std::string ToUpper( const std::string &STRING);

    //! translates strings up to size 9 to size_t; needed for example for sending messages among processes -msgsnd, msgrcv cannot handle strings correctly on all OS
    size_t StringToSizeT( const std::string &STRING);

    //! translates size_t to string up to size 9
    std::string SizeTToString( const size_t &VALUE);

    //! searches for spaces and removes them from the string
    inline
    std::string RemoveSpacesFromString( const std::string &STRING)
    {
        std::string cleaned_string;
        for( size_t i = 0 ;i < STRING.size(); i++ )
        {
            if( STRING[i] != ' ') cleaned_string.push_back( STRING[i]);
        }
        return cleaned_string;
      }

    inline
    bool
    IsCapitolLetter( const char &CHAR)
    {
        size_t nr( CHAR);
        return nr >= 65 && nr <= 90;
    }

    inline
    bool
    IsSmallLetter( const char &CHAR)
    {
        size_t nr( CHAR);
        return nr >= 97 && nr <= 122;
    }

    inline
    bool
    IsLetter( const char &CHAR)
    {
        size_t nr( CHAR);
        return ( nr >= 65 && nr <= 90) || ( nr >= 97 && nr <= 122);
    }

    inline
    bool
    IsNumber( const char &CHAR)
    {
        size_t nr( CHAR);
        return nr >= 48 && nr <= 57;
    }


} // end namespace mystr

#endif
