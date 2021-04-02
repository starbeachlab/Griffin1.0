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


#include "../../include/string/string_functions.h"
#include <cmath>

namespace mystr
{

    //! translate a char array into a string, removing eventual non alphanumerical endings
    void BufferToFileName( const char BUFFER[], std::string &STR, const int &SIZE)
    {
	std::string str = BUFFER;
	
	STR = str.substr( 0, str.find( ".") + 4);

	std::cout << __FUNCTION__ << ": " << STR << std::endl;

/*
	str = str.substr( 0, str.size() - 1);
//	    const char * c = &BUFFER[0];
	for( std::string::iterator c = str.begin(); c != str.end(); ++c)
	    {
		if( isalnum( *c) || *c == '.' || *c == '-' || *c == '_' || *c == '+')
		{
		    STR += *c;
		}
		}
*/
    }



    //! trims spaces from beginning and end of a copied string
    std::string TrimString( const std::string &STRING)
    {
      //searches for last leading space
      const std::string::size_type pos1( STRING.find_first_not_of( " \n\t\r\0"));
      //searches for first tailing space
      const std::string::size_type pos2( STRING.find_last_not_of( " \n\t\r\0"));

      //returns substring from pos1 of length pos2 - pos1 + 1
      return STRING.substr( ( pos1 == std::string::npos ? 0 : pos1), ( pos2 == std::string::npos ? 0 : pos2 - pos1 + 1));
    }


   //! test whether string is numerical (float, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    bool IsNumerical( const std::string &STRING)
    {
      //remove leading and tailing spaces
      const std::string trimmed_string( TrimString( STRING));

      //check that string is not empty
      if( trimmed_string.empty())
      {
        return false;
      }
      bool contains_point( false);
      bool contains_eE( false);

      // test all characters
      for( size_t i( 0); i < trimmed_string.length(); ++i)
      {
        // character is between 0 and 9 -> ok
        if( trimmed_string[i] >= '0' && trimmed_string[i] <= '0' + 9)
        {
          continue;
        }
        //character is '.' and the only point-> ok
        else if( trimmed_string[i] == '.' && !contains_point)
        {
          contains_point = true;
          continue;
        }

        else if( ( trimmed_string[i] == 'e' || trimmed_string[i] == 'E') && !contains_eE)
        {
          contains_eE = true;
          continue;
        }

        //is charater '+' or '-' -> ok
        else if( ( ( trimmed_string[i] == '-') || ( trimmed_string[i] == '+')))
        {
          continue;
        }
        //if neither of above cases -> not a numerical value
        else
        {
          return false;
        }

      }

      // return true
      return true;
    }


    std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER)
    {
    	if( STRING.length() == 0)
    	{
    		return std::vector< std::string>();
    	}
    	std::vector< std::string> result;
    	size_t a(0),b(0);
    	do
    	{
    		a = STRING.find_first_not_of( SPLITTER,b);
    		b = STRING.find_first_of( SPLITTER, a);
    		if( b != std::string::npos)
    		{ result.push_back( STRING.substr( a, b-a));}
    		else
    		{ result.push_back( STRING.substr( a));}
    	} while( a != std::string::npos && b != std::string::npos && STRING.find_first_not_of( SPLITTER,b) != std::string::npos);

    	return result;
    }


    std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER_A, const std::string &SPLITTER_B)
    {
          std::vector< std::string>
              result;
          size_t
              first_of(0),
              first_not_of( 0),
              value_a,
              value_b;

          do
            {
//              DebugWrite( __FUNCTION__ << ": first of: " << first_of);
              value_a = STRING.find_first_not_of( SPLITTER_A, first_not_of);
              value_b = STRING.find_first_not_of( SPLITTER_B, first_not_of);

              if( value_a > first_not_of && value_b > first_not_of)
              {
                  first_of = std::min( value_a, value_b);
              }
              else if( value_a > first_not_of)
              {
                  first_of = value_a;
              }
              else if( value_b > first_not_of)
              {
                  first_of = value_b;
              }
//              DebugWrite( __FUNCTION__ << ": first of: " << first_of);
//              DebugWrite( __FUNCTION__ << ": first not of: " << first_not_of);


              first_not_of = std::min( STRING.find_first_of( SPLITTER_A, first_of), STRING.find_first_of( SPLITTER_B, first_of));
//              DebugWrite( __FUNCTION__ << ": first not of: " << first_not_of);

              if( first_not_of != std::string::npos)
              {
//                  DebugWrite( "<" << STRING.substr( first_of, first_not_of - first_of) << ">");
                  result.push_back( STRING.substr( first_of, first_not_of - first_of));
              }
              else if( first_of != std::string::npos)
              {
//                  DebugWrite( "<<" << STRING.substr( first_of) << ">>");
                  result.push_back( STRING.substr( first_of));
              }
            } while
                (
                        first_of != std::string::npos
                        && first_not_of != std::string::npos
                        && ( STRING.find_first_not_of( SPLITTER_A, first_not_of) != std::string::npos || STRING.find_first_not_of( SPLITTER_B, first_not_of) != std::string::npos)
                );
//          DebugWrite( __FUNCTION__ << ": done");

          return result;
    }



    std::string &AllToLowerCase( std::string &STRING)  // BE AWARE: this always mutates string, non const
    {
        std::transform( STRING.begin(), STRING.end(), STRING.begin(), ::tolower);
        return STRING;
    }

    std::string ToLower( const std::string &STRING)
    {
        std::string str( STRING);
        std::transform( str.begin(), str.end(), str.begin(), ::tolower);
        return str;
    }

    std::string &AllToUpperCase( std::string &STRING)
    {
        std::transform( STRING.begin(), STRING.end(), STRING.begin(), toupper);
        return STRING;
    }

    std::string ToUpper( const std::string &STRING)
    {
        std::string str( STRING);
        std::transform( str.begin(), str.end(), str.begin(), toupper);
        return str;
    }

    size_t StringToSizeT( const std::string &STRING) // no encryption
    {
        assert( STRING.size() < 10);
        int i(0);
        size_t result(0);
        for( std::string::const_reverse_iterator itr = STRING.rbegin(); itr != STRING.rend(); ++itr, ++i)
        {
            result += ( size_t( *itr) - 32) * size_t( pow( 100.0, i));
        }
        return result;
    }

    std::string SizeTToString( const size_t &VALUE)
    {
        size_t value( VALUE);
        std::string str( "");
        while( value > 0)
        {
            size_t rest = value % 100;
            char c ( rest + 32);
//            std::cout << c;
            str.insert( str.begin(), c);
            value = size_t( value / 100);
        }
//        std::cout << std::endl;
        return str;
    }



} // end namespace mystr
