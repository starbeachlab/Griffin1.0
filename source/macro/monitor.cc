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


#include "../../include/macro/monitor.h"

#ifdef MONITORING

void Memory( const std::string &FILE, const int &STEP)
{
	static int call = 0;

	std::string
		pid_str,
		call_str,
		step_str;

	std::stringstream ss;

	ss << getpid();
	ss >> pid_str;

	ss.clear( std::stringstream::goodbit);
	ss << call;
	ss >> call_str;

	ss.clear( std::stringstream::goodbit);
	ss << STEP;
	ss >> step_str;

	if( call == 0)
	{
		system( ("echo memory for " + pid_str + " > " + FILE).c_str());
	}

	++call;


	// THIS SYNTAX MAY VARY ON DIFF SYSTEMS !!!
	system( ( "pmap -d " + pid_str + " | grep writ | awk '{ print \"" + call_str + " \"  \"" + step_str + " \" substr( $1, 1, length( $1) - 1) \"  \" substr( $3, 1, length( $3) - 1) ;}' >> " + FILE + "&").c_str());
}
#endif
