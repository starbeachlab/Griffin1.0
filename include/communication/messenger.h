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


#ifndef Messenger_H_
#define Messenger_H_

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>


#include <csignal>

#include <unistd.h>


#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"



namespace comm
{

    class Messenger
    : public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        int m_MessageID;

    public:

        virtual ~Messenger(){}

	void ClearConnection( const int &VALUE)
	    {
		int 
		    message_id = 0,
		    cc = 0;
		
		while( message_id != -1 && cc++ < 10)
		{
		message_id = msgget( (key_t) VALUE, 0666);
		if( message_id != -1)
		{
		    std::cout << "====> " << __FUNCTION__ << ": old message queue found, id: " << message_id << " in ipc channel: " << VALUE << std::endl;
		    std::cout << "====>  did you start this program multiple times on the same machine or did you have other programs running that might have used a message queue with the same key? if so, that one is in trouble!" << std::endl;

		    if( msgctl( message_id, IPC_RMID, 0) == -1)
		    {
			std::cout << "====> message queue " << message_id << " could not be deleted" << std::endl;
		    }
		    else
		    {
			std::cout << "===> message queue ( " << message_id << "/" << VALUE << ") deleted " << std::endl;
		    }
		}
		else
		{
		    std::cout << "nothing to clear in ipc channel: " << VALUE << std::endl;
		}
		}
	    }


//        inline int InitializeConnection()
//        {
//            m_MessageID = msgget((key_t)123, 0666 | IPC_CREAT);
//
//            if( m_MessageID == -1)
//            {
//                std::cerr << "====> connection could not be created!" << std::endl;
//                return -1;
//            }
//            std::cout << "connection is established. messageID: " << m_MessageID << std::endl;
//            return 0;
//        }

        inline int InitializeConnection( const int &VALUE)
        {
            m_MessageID = msgget((key_t)VALUE, 0666 | IPC_CREAT);

            if( m_MessageID == -1)
            {
                std::cerr << "====> connection could not be created!" << std::endl;
                return -1;
            }
            std::cout << "connection is established. message queue id: " << m_MessageID << " ipc channel: " << VALUE << std::endl;
            return 0;
        }

        inline int CutConnection()
        {
            if (msgctl(m_MessageID, IPC_RMID, 0) == -1)
            {
                std::cerr << "====> msgctl() failed!" << std::endl;
                return -1;
            }
            std::cout << "connection is cut. messageID: " << m_MessageID << std::endl;
            return 0;
        }

        const int &
        GetMessageID() const
        {
        	return m_MessageID;
        }
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
            STREAM << "messageID: " << m_MessageID << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Messenger
} // end namespace comm




#endif /* Messenger_H_ */
