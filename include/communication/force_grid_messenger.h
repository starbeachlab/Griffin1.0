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


#ifndef MESSENGER_H_
#define MESSENGER_H_


#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "messenger.h"
#include "force_grid_job_info.h"

namespace comm
{

    class ForceGridMessenger
    : public Messenger
    {

    public:
        //! virtual destructor
        virtual ~ForceGridMessenger(){}

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////
        int SendCommand( const ForceGridJobInfo INFO) // test with and without reference here
        {
            int check = msgsnd
            (
                    m_MessageID,
                    (void *)&INFO,
                    sizeof( INFO) - sizeof( long),
                    0
            );
            if( check == -1)
            {
                std::cerr << "====> msgsnd failed! (" << __FUNCTION__ << ")" << std::endl;
                return -1;
            }
            return 0;
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
            STREAM << m_MessageID << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class ForceGridMessenger
} // end namespace comm




#endif /* MESSENGER_H_ */
