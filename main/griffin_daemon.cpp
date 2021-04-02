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
//!	The daemon 'carrying' the implicit potential force field.
//!	Is designed to collaborate with MD package.
//! Communication between GRIFFIN and MD is provided by the griffin_messenger.exe.
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#include <ctime>

// order of inclusion is crucial !


#include "../include/macro/griffin_definitions.h"
#include "../include/macro/macro_functions_read_write.h"
#include "../include/macro/griffin_includes.h"

float g_PDBFactor = 0.0;



int main( int ARGC, char *ARGV[])
{
    std::cout << std::endl << std::string( ARGV[0]) << std::endl << std::endl;

    std::cout << ThrowFullHeader() << std::endl;


    /////////////////
    // DECLARATIONS
    /////////////////
    std::ofstream
        write;
    std::ifstream
        read;
    CommandLineManager
        cmd;
    boost::shared_ptr< mol::MoleculeForceGrid>
        molecule_grid;

    std::time_t
        start,
        begin,
        now;


    std::time( &start);

    //////////////////
    // READ OPTIONS
    //////////////////
    command_line_factory::SetUpManagerForDaemon( cmd, ARGC, ( const char **) ARGV);

    util::BuildEnumHandler< util::FunctorEnum>();
    util::BuildEnumHandler< util::EnergyEnum>();

	/////////////////////////////
	/////   FORK PROCESS    /////
	/////////////////////////////
	//  fork to enable to send the process to background

    comm::BriefMessage msg;

    int
		this_process = 0;



    // split process
#ifdef DAEMON
    int
		ipc_channel = 0;

    ipc_channel = this_process + 1000;

    std::cout << "ipc channel for forked main processes: " << ipc_channel << std::endl;


    pid_t pid;
	
    std::cout << "\nfork process" << std::endl;

    if ((pid = fork()) < 0)
    {
		perror("fork failed");
		return 1;
    }
    
    // child process
    if (pid == 0)
    {
#endif
	
//////    DISCONNECT PROCESS FROM SHELL (no output into log) not necessary    //////
////    chdir("/");
//    // disconnect process from shell
//    setsid();  // set session id
//	  // close standard in- and output
//    close(STDIN_FILENO);
//    close(STDOUT_FILENO);
//    close(STDERR_FILENO);
//    // open ouput to 'trash'
//    open("/dev/null", O_RDWR);
//    dup(STDIN_FILENO);

		///////////////////////////////////////////////////////////////
		/////            CONSTRUCT GRID                     ///////////
		///////////////////////////////////////////////////////////////

#ifdef MPI_PARALLEL

	int 
	    nr_processes;
	
	MPI_Init( &ARGC, &ARGV);
	MPI_Comm_rank( MPI_COMM_WORLD, &this_process);
	MPI_Comm_size( MPI_COMM_WORLD, &nr_processes);
	std::cout << "process " << this_process << " out of " << nr_processes << std::endl;
	
	std::string
	    this_process_str = mystr::NumericalValueToString( this_process);
#endif

	// READ GRID FROM FILE //
	StandardWrite( "\nread force grid from file ...");


		molecule_grid = boost::shared_ptr< mol::MoleculeForceGrid>( new mol::MoleculeForceGrid());

		Open( read, cmd.GetArgumentStringForFlag( "force_grid"));
		read >> molecule_grid;
		Close( read);

//		Open( write, "test_grid.txt");
//		write << molecule_grid;
//		Close( write);

		StandardWrite( "grid reading done\n");

		if( cmd.IsFlagSet( "sforce_scale"))
		{
			float sforce_length = cmd.GetArgumentForFlag< float>( "sforce_scale");
			StandardWrite(  "rescale sforces: " + mystr::NumericalValueToString( sforce_length));
			mol::factory::RescaleConstantForces( molecule_grid, sforce_length);
/*
			Open( write, "rescaled_grid_" + mystr::NumericalValueToString( sforce_length) + ".txt");
			molecule_grid->Write( write);
			Close( write);
			std::cout << "done" << std::endl;
			return 0;
*/

		}

		if( cmd.IsFlagSet( "sforce_reset_min_angle"))
		{
			float value = cmd.GetArgumentForFlag< float>( "sforce_reset_min_angle");
			StandardWrite(  "sforce redirection minimum angle: " + mystr::NumericalValueToString( value));
			molecule_grid->GetSetVectorAngle().ResetAngle( math::DegreesToRadians( value));
		}


		std::string format = cmd.GetArgumentStringForFlag( "unit_system");
		StandardWrite( "coordinate input file expected in " + format + " format");
		molecule_grid->SetForceFieldType( format);


		if( cmd.IsFlagSet( "rescale"))
		{
			// eg rescale coulomb 0.98
			StandardWrite(  "rescale: ");
			std::vector< std::string> arguments = cmd.GetAllOptions( "rescale");
			assert( arguments.size() % 2 == 0);
			for( int i = 0; i < arguments.size() * 0.5; ++i)
			{
				StandardWrite( arguments[ 2*i] << "  " << arguments[ 2*i+1]);
				mol::factory::RescaleForces( molecule_grid, arguments[ 2 * i], mystr::ConvertStringToNumericalValue< float>(arguments[ 2 * i + 1]));
			}
		}

//	if( cmd.IsFlagSet( "smooth_force_transition"))
//	{
//	    StandardWrite( "smooth force transition layer");
//	    std::vector< float> values =  cmd.GetArgumentsForFlag< float>( "smooth_force_transition");
//	    store::TransientInteractionGridPoint().SetStartValue( values[0]);
//	    store::TransientInteractionGridPoint().SetTotalSteps( math::Round< int>( values[1]));
//	    mol::factory::ReplaceInteractionWithTransientGridPoints( molecule_grid);
//	}

		std::vector< std::string> lipid_names = cmd.GetAllOptions( "lipid_types");
		std::cout << "lipid names: ";
		std::copy( lipid_names.begin(), lipid_names.end(), std::ostream_iterator< std::string>( std::cout, " "));
		std::cout << std::endl;
		molecule_grid->SetLipidNames( lipid_names);
		store::TypeMappedGridPoint< mol::Atom, math::Vector3N>().SetMolTypes( molecule_grid->GetSurfGrid().GetMoleculeTypes());
		store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N>().SetMolTypes( molecule_grid->GetSurfGrid().GetMoleculeTypes());

		if( cmd.IsFlagSet( "sforce_exclude_hydrogens"))
		{
			StandardWrite( "no sforces for hydrogen atoms");
			molecule_grid->SetIgnoreHydrogens();
		}
		StandardWrite(  "read grid done");

		StandardWrite(  "insert interpolation layer");
		molecule_grid->GetSetAtomForceGrid()->GetSetPositionGrid() = boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
		    ( new store::InterpolGrid( molecule_grid->GetAtomForceGrid()->GetPositionGrid()));



#ifdef MPI_PARALLEL
		if( cmd.IsFlagSet( "logfile"))
		{
		    std::string log = cmd.GetArgumentStringForFlag( "logfile");
		    log = log.substr( 0, log.find( "."));
		    log += mystr::NumericalValueToString< int>( this_process) + ".log";
		    std::cout << "logfile: " << log << std::endl;
		    Open( write, log);
		}
		else
		{
		    std::cout << "logfile: griffin_daemon_" << this_process_str + ".log" << std::endl;
		    Open( write, "griffin_daemon_" + this_process_str + ".log");
		}
//		system( ( "ipcs > ipcs_before_" + this_process_str + ".log").c_str());
//		system( ( "echo $HOST > host_" + this_process_str + ".txt").c_str());

#else
		if( cmd.IsFlagSet( "logfile"))
		{
			// open logfile
			Open( write, cmd.GetArgumentStringForFlag( "logfile"));
			system( ( "ipcs >> " + cmd.GetArgumentStringForFlag( "logfile")).c_str());
		}
		else
		{
			Open( write, "griffin_daemon.log");
		}
		system( "ipcs > ipcs_before.log");

#endif


		// pass grid to messenger grid
		boost::shared_ptr< mol::MessengerMoleculeForceGrid>
			messenger_grid( new mol::MessengerMoleculeForceGrid( molecule_grid));

		StandardWrite( "molecule iterator");

		boost::shared_ptr< mol::MoleculeIterator>
			molecule_iterator( new mol::MoleculeIterator());

		if( cmd.IsFlagSet( "neighbor_list_update"))
		{
			int frequency = cmd.GetArgumentForFlag< int>( "neighbor_list_update");
			StandardWrite( ( "periodic molecule iterator: recalculate zero forces every " + mystr::NumericalValueToString( frequency) + " steps").c_str());
			molecule_iterator = boost::shared_ptr< mol::MoleculeIterator>( new mol::PeriodicMoleculeIterator( frequency));
		}

#ifdef MPI_PARALLEL

		StandardWrite( "parallel mpi molecule iterator, daemon is copied, watch memory!");
		molecule_iterator = boost::shared_ptr< mol::MoleculeIterator>( new mol::MPIMoleculeIterator( *molecule_iterator)); 

		if( this_process == 0)
		{
		    write << "head ";
#endif
		    write << "clear connection" << std::endl;
		    messenger_grid->ClearConnection( 123);
		    // set up connection
		    write << "initialize connection" << std::endl;
		    if( messenger_grid->InitializeConnection( 123) == -1)
		    { exit( -1);}
#ifdef MPI_PARALLEL
		}
#endif

		write << ThrowFullHeader() << std::endl;
		write << "\nthis is your friendly daemon most humbly waiting for your commands ...\n\n\n" << std::endl;

		StandardWrite(  "\n...grid is daemonized, ready to receive jobs\n");
		StandardWrite(  "shoot job at daemon with griffin_messenger.exe \n");
		    
#ifdef DAEMON
		std::cout << "child process: sending message to terminate parent process" << std::endl;
		    
		// release parent process from fork
		int parent_term_id = msgget( (key_t)ipc_channel, 0666|IPC_CREAT);
		if( -1 ==  msgsnd(
			parent_term_id,
			(void *)&msg,
			sizeof( msg) - sizeof( long),
			0))
		{
		    std::cerr << "====> parent process termination msgsnd failed!" << std::endl;
		    return -1;
		}
#endif

		// execution loop
		messenger_grid->ExecuteMessages( molecule_iterator, write);

		// termination
		std::time( &now);
		write << "daemon terminated" << std::endl;
		write << "time: " << ctime( &now) << std::endl;

#ifdef MPI_PARALLEL
		system( ( "ipcs > ipcs_after_" + this_process_str + ".log").c_str());
		if( this_process == 0)
		{
		    // disconnect
		    messenger_grid->CutConnection();
		}
		MPI_Finalize();
		StandardWrite( "\nMPI connection closed");
#else
		system( "ipcs > ipcs_after.log");
#endif

		Close( write);

#ifdef DAEMON
    } // end child
    else
    {
    	std::cout << "parent process: waiting for message from child (on ipc channel" << ipc_channel << ") ..." << std::endl;

		int term_message_id = msgget( (key_t)ipc_channel, 0666|IPC_CREAT);

		if( msgrcv( term_message_id, ( void*)&msg, sizeof( msg) - sizeof( long), 0, 0) == -1)
		{
			std::cerr << "====> parent process: message receiving failed ... but hey, why should you see this message then?? (did you remove the message queue?)" << std::endl;
		}
		std::cout << "parent process: received message from child, terminates itself (free shell) and the message queue\n" << std::endl;

		if( msgctl( term_message_id, IPC_RMID, 0) == -1)
		{
		    std::cout << "=====> termination of parent closing message queue failed" << std::endl;
		}
		std::cout << "parent process: Daemon has PID: " << pid << std::endl << std::endl;
    }
#endif

    std::time( &now);
    std::cout << "\nGriffin: total time: " << std::difftime( now, start) << "s\n\n" << std::endl;
    std::time( &begin);

    return 0;
};  // end main
