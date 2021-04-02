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


#ifndef STREAM_FUNCTIONS_H_
#define STREAM_FUNCTIONS_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

inline
std::istream &JumpOverComments( std::istream &STREAM, const char &IDENTIFIER = '#')
{
    char c;
    std::string tmp;

    while( (c = STREAM.get()) == IDENTIFIER || c == '\n' || c == ' ')
    {
        if( c == IDENTIFIER)
        {
            std::getline( STREAM, tmp);
            DebugWriteNoFlush( "*delete line: #");
            DebugWrite( tmp);
        }
    }
    if( !STREAM.eof())
    {
        STREAM.putback( c);
    }
    return STREAM;
}

//inline
//void Check( std::ios &STREAM)
//{
//    if( !STREAM)
//    {
//        std::cerr << __PRETTY_FUNCTION__ << std::endl;
//        std::cerr << "stream not opened: " << std::endl;
//        exit( -1);
//    }
//}

inline
void Close( std::ifstream &STREAM)
{
    STREAM.close();
    STREAM.clear();
}

inline
void Close( std::ofstream &STREAM)
{
    STREAM.close();
    STREAM.clear();
}

inline
void Open( std::ifstream &STREAM, const std::string &FILE, std::ostream &LOG = std::cerr)
{
    STREAM.open( FILE.c_str());
    if( !STREAM)
    {
        LOG << __PRETTY_FUNCTION__ << std::endl;
        LOG << "stream not opened: <" << FILE << ">" << std::endl;
        exit( -1);
    }
}


inline
void Open( std::ofstream &STREAM, const std::string &FILE, std::ostream &LOG = std::cerr)
{
    STREAM.open( FILE.c_str());
    if( !STREAM)
    {
        LOG << __PRETTY_FUNCTION__ << std::endl;
        LOG << "stream not opened: <" << FILE << ">" << std::endl;
        exit( -1);
    }
}


inline
bool SoftOpen( std::ifstream &STREAM, const std::string &FILE, std::ostream &LOG = std::cerr)
{
    STREAM.open( FILE.c_str());
    if( !STREAM)
    {
        LOG <<  "===>" << __PRETTY_FUNCTION__ << std::endl;
        LOG <<  "===>" << " stream not opened: <" << FILE << ">" << std::endl;
        return false;
    }
    return true;
}


inline
bool SoftOpen( std::ofstream &STREAM, const std::string &FILE, std::ostream &LOG = std::cerr)
{
    STREAM.open( FILE.c_str());
    if( !STREAM)
    {
        LOG << "===>" << __PRETTY_FUNCTION__ << std::endl;
        LOG << "===>" << " stream not opened: <" << FILE << ">" << std::endl;
        return false;
    }
    return true;
}

template< typename T>
void WriteObject(  std::ofstream &STREAM, const T &OBJECT, const std::string &FILE)
{
    STREAM.open( FILE.c_str());
    if( !STREAM)
    {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        std::cerr << "stream not opened: <" << FILE << ">" << std::endl;
        exit( -1);
    }
    STREAM << OBJECT;
    STREAM.close();
    STREAM.clear();
}


inline
void WriteSystemVariable( std::ostream &STREAM, const std::string &NAME)
{
	char * descr = getenv( NAME.c_str());
	if (descr)
	{
		STREAM << NAME << " = " << descr << std::endl;
	}
	else
	{
		STREAM << "variable " << NAME << " is not defined" << std::endl;
	}
}


inline
std::string SystemVariable( const std::string &NAME)
{
	char * descr = getenv( NAME.c_str());
	if (descr)
	{
		return std::string( descr);
	}
	return "UNDEFINED";
}


#endif /* STREAM_FUNCTIONS_H_ */
