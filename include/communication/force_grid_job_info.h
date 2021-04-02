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


#ifndef GRID_FORCE_JOB_INFO_H_
#define GRID_FORCE_JOB_INFO_H_

//#include "../readwrite/stream_operator.h"
//#include "../string/io_string_functions.h"

namespace comm
{
    enum e_JobType { NEW_MOL, EXECUTE, TERMINATE, RESCALE};

    struct BriefMessage
    {
        long    m_ID;
//        int		m_X;

        BriefMessage()
        : m_ID( 1)//, // NEVER USE A ZERO FOR THE LONG, DUDE !!!!!
 //         m_X( 0)
          {}
    };

    struct ForceGridJobInfo
    {
        long          m_ID;
        e_JobType     m_Type;
        char          m_MolsFile[120];
        char          m_CoordInputFile[120];
        char          m_ForceOutputFile[120];

        ForceGridJobInfo
        (
                const e_JobType &TYPE = EXECUTE,
                const char MOLS_INPUT_FILE[] = "",
                const char COORD_FILE[] = "",
                const char FORCE_FILE[] = ""
        )
        : m_ID( 1),
        m_Type( TYPE)
        {
            sprintf( m_MolsFile ,"%s\n", MOLS_INPUT_FILE);
            sprintf( m_CoordInputFile ,"%s\n", COORD_FILE);
            sprintf( m_ForceOutputFile ,"%s\n", FORCE_FILE);
        }

    }; // end struct ForceGridJobInfo
} // end namespace comm





#endif /* GRID_FORCE_JOB_INFO_H_ */
