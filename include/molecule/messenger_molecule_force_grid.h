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


#ifndef MESSENGER_MOLECULE_GRID_H_
#define MESSENGER_MOLECULE_GRID_H_

//#include <ctime>

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"

#include "molecule_force_grid.h"
#include "molecule_force_grid_factory.h"
#include "../communication/messenger.h"
#include "../communication/force_grid_job_info.h"
#include "../utilities/time.h"

#include "../macro/griffin_definitions.h"

namespace mol
{

    class MessengerMoleculeForceGrid
    : public comm::Messenger
    {
    protected:
        boost::shared_ptr< MoleculeForceGrid>  m_Grid;
        std::string                            m_InputFile;
        std::string                            m_OutputFile;
    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        MessengerMoleculeForceGrid()
        : m_Grid(),
        m_InputFile(),
        m_OutputFile()
        {}

        //! construct from MoleculeForceGrid
        MessengerMoleculeForceGrid( const boost::shared_ptr< MoleculeForceGrid> &GRID)
        : m_Grid( GRID),
        m_InputFile(),
        m_OutputFile()
        {}

        //! copy constructor
        MessengerMoleculeForceGrid( const MessengerMoleculeForceGrid &ORIGINAL)
        : m_Grid( ORIGINAL.m_Grid),
        m_InputFile( ORIGINAL.m_InputFile),
        m_OutputFile( ORIGINAL.m_OutputFile)
        {
        	std::cout << __FUNCTION__ << " copy constructor" << std::endl;
        }

        //! virtual destructor
        virtual ~MessengerMoleculeForceGrid(){}

        //! virtual copy constructor
        virtual MessengerMoleculeForceGrid *Clone() const{ return new MessengerMoleculeForceGrid( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////


        int ExecuteMessages
        (
                boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
                std::ostream &LOG
//#ifdef MPI_PARALLEL
//		, std::string SYS_CALL
//#endif
        )
        {
            bool is_finished( false);

            struct timeval  start, end;

//	 	    time_t  begin, final;

            comm::ForceGridJobInfo message;
	    
            comm::BriefMessage brief_message;
//            comm::BriefMessage brief_message = { 0, ""};

	    int term_message_id = 0;

#ifdef MPI_PARALLEL
            // only head process is messaging via IPC
	    // head process is the LAST one!
            int
		message_type,
		nr_processes,
		this_process,
		tag = 1;

	    MPI_Status status;

	    //	    int empty[0];

	    MPI_Comm_size( MPI_COMM_WORLD, &nr_processes);
            MPI_Comm_rank( MPI_COMM_WORLD, &this_process);

	    const int
		head_id = 0,
		first_kid = 1,
		last_kid = nr_processes - 1;

            LOG << "mpi setup, process " << this_process << " out of " << nr_processes << std::endl;
            if( this_process == head_id)
            {
            	LOG << "head: cleans up old confirmation queues" << std::endl;
#else
				LOG << "clean up old confirmation message queues" << std::endl;
#endif

				int cc = 0;

				// if confirmation message queue is already open, close it and reopen, not to catch old messages
				while( term_message_id != -1 && cc++ < 10)
				{
					term_message_id =  msgget( (key_t)234, 0666);
					if( term_message_id != -1)
					{
						LOG << "====> messenger: old confirmation message queue found in ipc channel 234! id: " << term_message_id << std::endl;
						if( msgctl( term_message_id , IPC_RMID, 0) == -1)
						{
								LOG << "====> " << __FUNCTION__ << " message queue could not be closed: " << term_message_id  << std::endl;
						}
						else
						{
						LOG << "confirmation message queue closed" << std::endl;
						}
					}
				}

				// open || reopen confirmation message queue
				term_message_id = msgget( (key_t)234, 0666|IPC_CREAT);
		//                    int term_message_id = msgget( IPC_PRIVATE,0600|IPC_CREAT|IPC_EXCL);
				if( term_message_id == -1)
				{
					LOG  << "====> connection could not be created! (" << __FUNCTION__ << ")" << std::endl;
				}
				LOG << __FUNCTION__ << ": connection for execution confirmation established, id: " << term_message_id << std::endl;


#ifdef MPI_PARALLEL

				// SEND THE CONFIRMATION TO WAITING MD PROCESS

				int
					id = msgget( (key_t) 12342, 0666);  // sender does not create queue, queue is opened and closed by waiting process

				if( id == -1)
				{
					std::cout << "message queue could not be opened! " << std::endl;
					return -1;
				}
				std::cout << "message queue: " << id << " opened" << std::endl;

				msgsnd( (key_t) id, (void *) &brief_message, sizeof( brief_message) - sizeof( long), 0);
		
		
//		if( SYS_CALL.size() != 0)
//		{
//		    LOG << "head process calls external script: <" << SYS_CALL << ">" << std::endl;
//		    if( SYS_CALL.find("&") != std::string::npos)
//		    {
//			LOG << "====> system call should not contain '&'" << std::endl;
//		    }
//		    SYS_CALL += "&";
//		    system( SYS_CALL.c_str());
//		}
            }
#endif


            while( !is_finished)
            {

#ifdef MPI_PARALLEL
            	if( this_process == head_id)
            	{
//            		LOG << "head proc waits for message from job-gun" << std::endl;
#endif

				// IPC receive of job
				// if parallel and daddy
				// or if not parallel
				if( msgrcv( m_MessageID, ( void*)&message, sizeof( message) - sizeof( long), 0, 0) == -1)
				{
				LOG << "====> msgrcv failed! (" << __FUNCTION__ << ")" << std::endl;
				return -1;
				}



#ifdef MPI_PARALLEL
					gettimeofday( &start, NULL);
//					LOG << "head: received message, send order to kids" << std::endl;

					message_type = (int) message.m_Type;

					// if parallel and daddy send message to kids to act
					for( int i = first_kid; i <= last_kid; ++i)
					{
//							LOG << "head: send to kid " << i << " type " << message_type << std::endl;
							MPI_Send( &message_type, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					}
            	}
            	// if parallel and child: wait for message from daddy
            	else
            	{
    				gettimeofday( &start, NULL);
//            		LOG << "kid " << this_process << " waits" << std::endl;
            		MPI_Recv( &message_type, 1, MPI_INT, head_id, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//            		LOG << "kid " << this_process << " got the message "  << " type " << message_type << std::endl;
            		message.m_Type = comm::e_JobType( message_type);
            	}
#else
				gettimeofday( &start, NULL);
#endif

//		time( &begin);
//		LOG << ctime( &begin );

                switch( message.m_Type)
                {
                case comm::EXECUTE:
                {
                	//                    LOG << "execute" << std::endl;
					m_Grid->Execute
					(
									m_InputFile,
									m_OutputFile,
									MOL_ITR,
									LOG
					);

#ifdef MPI_PARALLEL
#ifdef FULL_PARALLEL
					if( this_process == head_id)
					{
#else // FULL_PARALLEL
					if( this_process != head_id)
					{
						// waiting for confirmation
						LOG << "send confirmation to head process" << std::endl;
						MPI_Send( 0, 0, MPI_INT, head_id, 99999, MPI_COMM_WORLD);
					}
					else
					{
						int empty[0];
						// send termination confirmation to kids
						for( int i = first_kid; i <= last_kid; ++i)
						{
						LOG << "waiting for confirmation of child: " << i << std::endl;
						MPI_Recv( empty, 0, MPI_INT, i, 99999, MPI_COMM_WORLD, &status);
						}
						LOG << "send confirmation to messenger" << std::endl;
#endif // FULL_PARALLEL
#endif // MPI_PARALLEL
					if( msgsnd( term_message_id, ( void *)&brief_message, 0, 0) == -1)
					{
					LOG << "====> msgsnd failed! (" << __FUNCTION__ << ")" << std::endl;
					}

//                    time( &end);
//                    LOG << "execution time: " << difftime( end, now)<< std::endl;

#ifdef MPI_PARALLEL
				}
#endif
				}
					break;
                case comm::TERMINATE:

                	LOG << "terminate grid daemon on request ..." << std::endl;
                	is_finished = true;

                	break;
                case comm::NEW_MOL:
                {
                	std::string
						mol_file;

#ifdef MPI_PARALLEL
					if( this_process == head_id)
					{
#endif
					mol_file = message.m_MolsFile;
					m_InputFile = message.m_CoordInputFile;
					m_OutputFile = message.m_ForceOutputFile;

					mol_file = mol_file.substr( 0, mol_file.find( "\n"));
					m_InputFile = m_InputFile.substr( 0, m_InputFile.find( "\n"));
					m_OutputFile = m_OutputFile.substr( 0, m_OutputFile.find( "\n"));

#ifdef MPI_PARALLEL
						for( int i = first_kid; i <= last_kid; ++i)
						{
						MPI_Send( ( void *) mol_file.c_str(), mol_file.size(), MPI_CHAR, i, 100, MPI_COMM_WORLD);
						}
						for( int i = first_kid; i <= last_kid; ++i)
						{
						MPI_Send( ( void *) m_InputFile.c_str(), m_InputFile.size(), MPI_CHAR, i, 101, MPI_COMM_WORLD);
						}
						for( int i = first_kid; i <= last_kid; ++i)
						{
						MPI_Send( ( void *) m_OutputFile.c_str(), m_OutputFile.size(), MPI_CHAR, i, 102, MPI_COMM_WORLD);
						}
						LOG << "daddy sent setup information" << std::endl;
					}
					else
					{
						char buffer1[120], buffer2[120],buffer3[120];
		//			    MPI_Recv( ( void *) mol_file.c_str(), 120, MPI_CHAR, 0, 100, MPI_COMM_WORLD, &status); // this is unstable
						MPI_Recv( buffer1, 120, MPI_CHAR, head_id, 100, MPI_COMM_WORLD, &status);
						mystr::BufferToFileName( buffer1, mol_file, 120);
		//			    sprintf( buffer, "");
						MPI_Recv( buffer2, 120, MPI_CHAR, head_id, 101, MPI_COMM_WORLD, &status);
						mystr::BufferToFileName( buffer2, m_InputFile, 120);
		//			    sprintf( buffer, "");
						MPI_Recv( buffer3, 120, MPI_CHAR, head_id, 102, MPI_COMM_WORLD, &status);
						mystr::BufferToFileName( buffer3,  m_OutputFile, 120);
						LOG << "kid: " << this_process << " received setup information " << std::endl;

					}

#ifndef FULL_PARALLEL
					m_OutputFile = m_OutputFile.substr( 0, m_OutputFile.find_last_of( ".") + 1) + mystr::NumericalValueToString( this_process) + ".tmp";
#endif   //  FULL_PARALLEL

#endif   //  MPI_PARALLEL
			
				LOG << __FUNCTION__ << ": setup new molecule:  <" <<  mol_file << ">  coor-in:  <" << m_InputFile << ">  forces-out:  <" <<  m_OutputFile << ">" << std::endl;
			
				std::ifstream
					read;
			
				read.open( mol_file.c_str());

				if( read)
				{
					LOG << "read mols" << std::endl;
					MOL_ITR->SetMols( *mol::file::ReadMoleculesInGriffinFormat( read, LOG));

//					LOG << "******check content of molitr: " << MOL_ITR << "\n*********" << std::endl;

#ifdef CPPIO
					std::ofstream write;
#ifdef MPI_PARALLEL			    
					if( this_process == head_id && !SoftOpen( write, m_OutputFile))
					{
					LOG << "====> out file could not be opened" << std::endl;
					}
#else
					if( !SoftOpen( write, m_OutputFile))
					{
					LOG << "====> out file could not be opened" << std::endl;
					}
#endif

#else // CPPIO
					FILE *write = NULL; // = fopen( m_OutputFile.c_str(), "w");
#ifdef MPI_PARALLEL			    
					if( this_process == head_id)
					{
					    write = fopen( m_OutputFile.c_str(), "w");
					}
#else  // MPI
					write = fopen( m_OutputFile.c_str(), "w");
#endif // MPI
#endif // CPPIO


					LOG << "calculate and write forces, energy and virial for initial conformation" << std::endl;
					m_Grid->CalculateAndWriteForcesEnergyVirial( write, MOL_ITR, LOG);



#ifdef CPPIO
#ifdef MPI_PARALLEL			    
					if( this_process == head_process)
					{
					    Close( write);
					}
#else
					Close(write);
#endif 
#else // CPPIO
#ifdef MPI_PARALLEL			    
					if( this_process == head_id)
					{
					    fclose( write);
					}
#else  // MPI
					fclose( write);
#endif // MPI
#endif // CPPIO

				}
				else
				{
					LOG << "=====>  explicit system file not opened" << std::endl;
				}
				read.close();
				read.clear();
			
			
			
#ifdef MPI_PARALLEL
//#ifdef FULL_PARALLEL
	//			if( this_process == 0)  // todo: probably unnecesary!!
	//			{
	//			    // daddy wait for kids to be done
	//			    for( int i = 1; i < nr_processes; ++i)
	//			    {
	//				LOG << "daddy: waits for kid " << i << " to verify that it is done (probably waisting time with  this!!!!)" << std::endl;
	//				MPI_Recv( empty, 0, MPI_INT, i, 999, MPI_COMM_WORLD, &status);
	//				LOG << "daddy: kid " << i << " said its done" << std::endl;
	//			    }
//
//#else // FULL_PARALLEL
				if( this_process == head_id)
				{
/*  CHECK WHETHER NEEDED:
				for( int i = first_kid; i < head_id; ++i)
			    {
					LOG << "head process waits for process " << i << " to verify that it is done" << std::endl;
					MPI_Recv( empty, 0, MPI_INT, i, 999, MPI_COMM_WORLD, &status);
					LOG << "head: process " << i << " is done" << std::endl;
			    }
*/

//#endif // FULL_PARALLEL
#endif // MPI_PARALLEL
			    
			    if( msgsnd( term_message_id, ( void *)&brief_message, 0, 0) == -1)
			    {
				LOG << "====> msgsnd failed! (" << __FUNCTION__ << ")" << std::endl;
			    }

#ifdef MPI_PARALLEL
//#ifdef FULL_PARALLEL
//			}
//			else
//			{
//			// all kids send message to daddy that they're done with their job
//			    for( int i = 1 ; i < nr_processes; ++i)
//			    {
//				LOG << "kid: " << i << " is done" << std::endl;
//				MPI_Send( 0, 0, MPI_INT, 0, 999, MPI_COMM_WORLD);
//			    }
//			}
//#else // FULL_PARALLEL
            }
/*  CHECK WHETHER NEEDED:
			else
			{
			    for( int i = first_kid; i < head_id; ++i)
			    {
				LOG << "process " << i << " is done and sends mpi confirmation" << std::endl;
				MPI_Send( 0, 0, MPI_INT, head_id, 999, MPI_COMM_WORLD);
			    }
			}
*/
//#endif // FULL_PARALLEL
#endif // MPI_PARALLEL
			
            }
                break;
            default: LOG << "no useful option" << std::endl;
			
            } // end switch
		
		
//		time( &final);
//		LOG << "oldschool: " << difftime( final, begin) << std::endl;
		
            gettimeofday( &end, NULL);
            LOG << "execution time: " << TimeDiff( NULL, &end, &start) << "s\n"<< std::endl;
        } // end 'endless' while
        LOG << "grid daemon terminated" << std::endl;
#ifdef MPI_PARALLEL
	    if( this_process == head_id)
	    {
#endif
		if( msgctl( term_message_id , IPC_RMID, 0) == -1)
		{
		    LOG << "====> " << __FUNCTION__ << " message queue could not be closed: " << term_message_id  << std::endl;
		}
		else
		{
		    LOG << "confirmation message queue closed" << std::endl;
		}
#ifdef MPI_PARALLEL
	    }
#endif
            return 0;
        } // end ExecuteMessages()


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << std::endl;
            STREAM << m_Grid;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class MessengerMoleculeForceGrid
} // end namespace mol




#endif /* MESSENGER_MOLECULE_GRID_H_ */
