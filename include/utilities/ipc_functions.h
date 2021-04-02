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


#ifndef IPC_FUNCTIONS_H_
#define IPC_FUNCTIONS_H_

namespace util
{
	void CleanupQueues( const std::string &NODE)
	{
		std::cout << "\n======   " << __FUNCTION__ << " " << NODE << "   =====" << std::endl;
		std::cout << "this is a very basic implementation of a cleanup function.\n";
		std::cout << "a more fancy way would include either serious socket programming or another library like libssh2\n";
		std::cout << "password free access to your nodes should be enabled.\n";
		std::cout << "otherwise you will spend quite some time typing passwords.\n" << std::endl;

		std::ifstream in;
		std::string
			line,
			name = "ipcs_" + NODE + ".tmp",
			command = "ssh " + NODE + " \"ipcs\" >& " + name;
		bool
			correct_block = false;
		std::vector< std::string>
			cols;

		std::cout << command << std::endl;
		system( command.c_str());

		in.open( name.c_str());

		while( in)
		{
			std::getline( in, line);
			if( correct_block && line.size() > 1 && line.substr( 0, 3) != "key")
			{
				cols = mystr::SplitString( line);
				command = "ssh " + NODE + " \"ipcrm -q " + cols[1] + "\"";
				std::cout << "found: " << line << std::endl;
				std::cout << command << std::endl;
				system( command.c_str());
			}

			if( line.substr( 0, 1) == "-" && line.find( "Message Queues") != std::string::npos)
			{
				correct_block = true;
			}
		}

		in.close();
		command = "rm " + name;
		system( command.c_str());
	}

	void CleanupQueuesInFile( const std::string &FILE)
	{
		std::string file;
		std::ifstream in( FILE.c_str());
		if( !in)
		{ std::cerr << FILE << " not opened!" << std::endl; exit( -1);}

		while( in)
		{
			file.clear();
			in >> file;
			if( file.length() != 0)
			{
				CleanupQueues( file);
			}
		}
		in.close(); in.clear();
	}

} // end namespace

#endif /* IPC_FUNCTIONS_H_ */
