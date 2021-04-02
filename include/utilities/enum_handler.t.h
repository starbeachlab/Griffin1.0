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


#ifndef ENUM_HANDLER_H
#define ENUM_HANDLER_H

#include <map>
#include <string>
#include <iostream>
#include "../macro/griffin_global.h"

namespace util
{
	template< typename t_ENUM>
	class EnumHandler
	{
	private:
		static std::map< t_ENUM, std::string>  sm_Id2String;
		static std::map< std::string, t_ENUM>  sm_String2Id;

	public:
		EnumHandler(){}

		~EnumHandler(){}

		void Insert( const t_ENUM &ID, const std::string &STR)
		{
		    if( sm_Id2String.find( ID) != sm_Id2String.end())
		    {
			std::cout << "\n====> ID: " << ID << " already exists in map (from " << __PRETTY_FUNCTION__ << ")" << std::endl;
		    }
		    if( sm_String2Id.find( STR) != sm_String2Id.end())
		    {
			std::cout << "\n====> STRING: " << STR << " already exists in map (from " << __PRETTY_FUNCTION__ << ")" << std::endl;
		    }

			sm_Id2String.insert( std::make_pair( ID, STR));
			sm_String2Id.insert( std::make_pair( STR, ID));
//			std::cout << "enum handler inserted: " << ID << "  " << STR << std::endl;
		}


		std::string String( const t_ENUM &ID) const
		{
			typename std::map< t_ENUM, std::string>::const_iterator itr = sm_Id2String.find( ID);
			if( itr == sm_Id2String.end())
			{
				std::cout << "\n====> ID: " << ID << " not found in EnumHandler";
				return std::string( "UNDEFINED");
			}
			return itr->second;
		}


		t_ENUM ID( const std::string &STR) const
		{
			typename std::map< std::string, t_ENUM>::const_iterator itr = sm_String2Id.find( STR);
			if( itr == sm_String2Id.end())
			{
				std::cout << "\n====> String: " << STR << " not found in EnumHandler";
				return t_ENUM( 88888888);
			}
	//		std::cout << std::endl <<  __FUNCTION__ << ": " << STR << " = " << itr->first << " => " << itr->second << std::endl;
			return itr->second;
		}


		std::ostream & WriteString2ID( std::ostream &STREAM = std::cout) const
		{
		    for( typename std::map< std::string, t_ENUM>::const_iterator itr = sm_String2Id.begin(); itr != sm_String2Id.end(); ++itr)
		    {
			STREAM << "<" << itr->first << ">\t<" << itr->second << ">" << std::endl;
		    }
		    return STREAM;
		}

	}; // end class EnumHandler


	// outside initialization of static members
	template< typename t_ENUM>
	std::map< t_ENUM, std::string> EnumHandler< t_ENUM>::sm_Id2String
		= std::map< t_ENUM, std::string>();
	template< typename t_ENUM>
	std::map< std::string, t_ENUM> EnumHandler< t_ENUM>::sm_String2Id
		= std::map< std::string, t_ENUM>();


	template< typename t_ENUM>
	void
	BuildEnumHandler()
	{
	    StandardWrite( "====> should not be called: " << __PRETTY_FUNCTION__);
	}

} // end namespace util

#endif

