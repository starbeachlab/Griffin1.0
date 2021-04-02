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
//!	Tool to communicate with the daemon.
//!	Called by external programs.
//!	The interface between GRIFFIN and e.g. MD.
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include <sys/time.h>
#include "../include/communication/force_grid_messenger.h"
#include "../include/readwrite/stream_functions.h"
#include "../include/readwrite/quotes.h"

float g_PDBFactor = 0.0;


int main( const int ARGC, const char *ARGV[])
{
	if( ARGC < 2 || std::string( ARGV[1]) == "-help")
    {
		std::cout << std::endl << std::string( ARGV[0]) << std::endl << std::endl;

		std::cout << ThrowFullHeader();

        std::cout << "\n\n\nAbout this program:" << std::endl;
        std::cout << "This tool allows to link the Griffin daemon with external programs." << std::endl;
        std::cout << "It controls the basic functionality of the daemon, like setting up a new exlicit system (e.g. bilayer), calculating forces and termination\n\n\n" << std::endl;

        std::cout << "========================  HELP  ===========================\n\n\n";
        std::cout << "============  required flags  ============\n\n" << std::endl;
        std::cout << "============  optional flags: standard  ============\n\n" << std::endl;

        std::cout << "-help    TYPE" << std::endl;
        std::cout << "     TYPE can be 'standard' or 'debug' (default 'standard'), writes this help" << std::endl;
        std::cout << "     from 0 to 1 parameters\n" << std::endl;


        std::cout << "-setup  MOLECULES_INPUT_FILE  COORD_INPUT_FILE  FORCES_OUTPUT_FILE" << std::endl;
        std::cout << "     the set-up mode is used only once, to prepare the daemon for receiving jobs;\n     MOLECULES_INPUT_FILE is a Griffin-formatted file (.gfn)\n     with the coordinates and topology of the system used in the molecular dynamics simulation;\n     COORD_INPUT_FILE and FORCES_OUTPUT_FILE are filenames of the coordinate and forces files\n     that will be written at each step of the molecular dynamics simulation;\n     all filenames must be without paths, and should contain extensions of three characters\n     (e.g. namdf.tmp, forces.dat)\n";
        std::cout << "     3 parameters\n" << std::endl;

        std::cout << "-forces" << std::endl;
        std::cout << "     the forces mode is used repeatedly during a molecular dynamics simulation and should be called by the MD package at each step;\n     it tells the daemon to read the current coordinates in the working directory (COORD_INPUT_FILE),\n     to calculate forces for those coordinates using the force-grid loaded in memory,\n     and to write those forces into the force output file (FORCES_OUTPUT_FILE);\n     these files are deleted by default" << std::endl;
        std::cout << "     0 parameters\n" << std::endl;



        std::cout << "-terminate" << std::endl;
        std::cout << "     this mode terminates the running daemon" << std::endl;
        std::cout << "     0 parameters\n" << std::endl;

        std::cout << "-wait   VALUE" << std::endl;
        std::cout << "     in the wait mode, the messenger waits on ipc channel VALUE until receives a message, then terminates itself (default channel: 12342) " << std::endl;
        std::cout << "     1 parameter\n" << std::endl;
        std::cout << "-send   VALUE" << std::endl;
        std::cout << "     in the send mode, the messenger sends empty message on ipc channel VALUE, then terminates itself (default channel: 12342) " << std::endl;
        std::cout << "     NOTE: the VALUE given in '-wait' and '-send' have to match in order for them to communicate!" << std::endl;
        std::cout << "     1 parameter\n" << std::endl;
        if( ARGC == 3 && std::string( ARGV[2]) == "debug")
        {
            std::cout << "============  optional flags: debug  ============\n\n" << std::endl;
        }
        std::cout << "=====================  HELP done  ========================\n\n" << std::endl;
        std::cout << "=====================  No flags given!  ===================== \n\n"<< std::endl;
        return 0;
    }

    std::cout << ARGV[0] << "  " << ARGV[1] << std::endl;

    if( std::string( ARGV[1]) == "-forces")
    {
		comm::ForceGridMessenger messenger;

		if( messenger.InitializeConnection( 123) == -1)
		{ return 1;}

		time_t start, end;
		time( &start);

		messenger.SendCommand( comm::ForceGridJobInfo( comm::EXECUTE));

		std::cout << "execution message sent, waiting for deamon to finish calculation ...(user: " << getenv( "USER") << ")" << std::endl;
		comm::BriefMessage msg;
		int term_message_id = msgget( (key_t)234, 0666); // not creating queue here!
		if( term_message_id == -1)
		{
			std::cerr << "====> connection could not be created!" << std::endl;
		}
		else
		{
			std::cout << "connection for return message established on ipc channel 234, message queue id: " << term_message_id << std::endl;
		}

        if( msgrcv( term_message_id, ( void*)&msg    , sizeof( msg) - sizeof( long), 0, 0) == -1)
        {
            std::cerr << "====> message receiving failed ... but hey, why should you see this message then??" << std::endl;
        }

        time( &end);

        std::cout << "... job_gun done (" << difftime( end, start) << "s)" << std::endl;
        return 0;
    }
    else if( std::string( ARGV[1]) == "-terminate")
    {
		comm::ForceGridMessenger messenger;

		if( messenger.InitializeConnection( 123) == -1)
		{ return 1;}
    
        std::cout << ThrowFullHeader();
        std::cout << "send termination message to daemon ..." << std::endl;
        messenger.SendCommand( comm::ForceGridJobInfo( comm::TERMINATE));
        std::cout << "... job_gun done\n" << std::endl;
        return 0;
    }
    else if( std::string( ARGV[1]) == "-setup")
    {
		comm::ForceGridMessenger messenger;

		if( messenger.InitializeConnection( 123) == -1)
		{ return 1;}
    
        std::cout << ThrowFullHeader();
        time_t start, end;
        time( &start);

        // check arguments
        std::string
			s1 = ARGV[2],
			s2 = ARGV[3],
			s3 = ARGV[4];

        if( s1.find( ".") != s1.length() - 4 || s2.find( ".") != s2.length() - 4 || s3.find( ".") != s3.length() - 4)
        {
        	std::cout << "====> file name has to follow this pattern: STRING.abc (only a single dot, and three letters after)" << std::endl;
        	std::cout << "got: " << s1 << ", " << s2 << ", " << s3 << std::endl;
        }

        char
            input_file[ 120],
            coord_input_file[120],
            output_file[120];

        sprintf( input_file ,"%s\n", ARGV[2]);
        sprintf( coord_input_file ,"%s\n", ARGV[3]);
        sprintf( output_file ,"%s\n", ARGV[4]);

        std::cout << "set up grid for new explicit molecules" << std::endl;
        messenger.SendCommand( comm::ForceGridJobInfo( comm::NEW_MOL, input_file, coord_input_file, output_file));

        std::cout << "setup message sent, waiting for deamon to finish calculation ...(user: " << getenv( "USER") << ")" << std::endl;
        comm::BriefMessage msg;
        int term_message_id = msgget( (key_t)234, 0666); // not creating queue here!
        if( term_message_id == -1)
        {
            std::cerr << "====> connection to message queue could not be created!" << std::endl;
        }
		else
		{
			std::cout << "connection to queue for return message established on ipc channel 234, message queue id: " << term_message_id << std::endl;
		}

        if( msgrcv( term_message_id, ( void*)&msg    , sizeof( msg) - sizeof( long), 0, 0) == -1)
        {
            std::cerr << "====> message receiving failed ... but hey, why should you see this message then??" << std::endl;
        }
        time( &end);
        std::cout << "... job_gun done (" << difftime( end, start) << "s)\n" << std::endl;
        return 0;
    }
    else if( std::string( ARGV[1]) == "-wait")
    {
        time_t start, end;
        time( &start);
        int ipc_channel = 12342;
        if( ARGC == 3)
        {
        	ipc_channel = mystr::ConvertStringToNumericalValue< int>( std::string( ARGV[2]));
        }
        std::cout << "messenger waits for message on ipc channel " << ipc_channel << std::endl;


		comm::BriefMessage
			msg;
		int
			id = msgget( (key_t) ipc_channel, 0666 | IPC_CREAT);

		if( id == -1)
		{
			std::cout << "message queue could not be created! " << std::endl;
			return -1;
		}
		std::cout << "message queue: " << id << " created" << std::endl;

		msgrcv( id, (void*)&msg, sizeof( msg) - sizeof( long), 0, 0);

		if( msgctl( id, IPC_RMID, 0) == -1)
		{
			std::cout << "====> message queue " << id << " could not be deleted" << std::endl;
		}
		else
		{
			std::cout << "message queue ( " << id << "/" << ipc_channel << ") deleted " << std::endl;
		}
	
		time( &end);
		std::cout << "messenger received message, erased queue and terminates, waited: " << difftime( end, start) << "s" << std::endl;
		return 0;
	}
    else if( std::string( ARGV[1]) == "-send")
    {
		int
			value = mystr::ConvertStringToNumericalValue< int>( std::string( ARGV[2]));
		comm::BriefMessage
			msg;
		int
			id = msgget( (key_t) value, 0666); // not creating queue here!

		if( id == -1)
		{
			std::cerr << "====> connection could not be created!" << std::endl;
		}
		else
		{
			std::cout << "messenger sends message on ipc channel " << ARGV[2] << ", message queue id: " << id << std::endl;
		}
		msgsnd( (key_t) id, (void *) &msg, sizeof( msg) - sizeof( long), 0);
		return 0;
    }


//    else if( std::string( ARGV[1]) == "-rescale_sforces")
//    {
//    	std::cout << ThrowFullHeader();
//    	char
//			empty[0],
//			value[12];
//
//    	sprintf( value, "%s\n", ARGV[2]);
//
//    	std::cout << "rescale sforces and energies inside implicit volume" << std::endl;
//    	messenger.SendCommand( comm::ForceGridJobInfo( comm::RESCALE, value, empty, empty, empty));
//    	std::cout << "... job gun done" << std::endl;
//    	return 0;
//    }

    std::cout << "nothing to do, no valid option given\n" << std::endl;
    return 1;
};
