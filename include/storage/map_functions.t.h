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


#ifndef MAP_FUNCTIONS_T_H_
#define MAP_FUNCTIONS_T_H_


namespace store
{
    template< typename t_KEY, typename t_VALUE>
    inline
    bool IsValueElementofMap( const t_VALUE &VALUE, const boost::shared_ptr< Map< t_KEY, t_VALUE> > &MAP)
    {
        for
        (
                typename std::map< t_KEY, t_VALUE>::const_iterator itr( MAP->begin());
                itr != MAP->end();
                ++itr
        )
        {
//            DebugWrite( itr->first << " " << itr->second);
            if( itr->second == VALUE)
            {
                return true;
            }
        }
        return false;
    }



    inline
    float
    ReadFromMap
    (
            const boost::shared_ptr< store::Map< std::string, float> > &MAP,
            const std::string &RESIDUE,
            const std::string &ATOM
    )
    {
        std::map< std::string, float>::const_iterator itr = MAP->find( RESIDUE + ":" + ATOM);
        if( itr == MAP->end())
        {
            itr = MAP->find( "NTER:" + ATOM);
            if( itr == MAP->end())
            {
                itr = MAP->find( "CTER:" + ATOM);
                if( itr == MAP->end())
                {
                    std::cout << "===> " << __FUNCTION__ << ": "<< RESIDUE << ":" << ATOM << " NOT found in map!" << std::endl;
                    return std::numeric_limits< float>::max();
                }
            }
        }
        return itr->second;
    }


    inline
    std::string
    ReadFromMap
    (
            const boost::shared_ptr< store::Map< std::string, std::string> > &MAP,
            const std::string &RESIDUE,
            const std::string &ATOM
    )
    {
        std::map< std::string, std::string>::const_iterator itr = MAP->find( RESIDUE + ":" + ATOM);
        if( itr == MAP->end())
        {
            itr = MAP->find( "NTER:" + ATOM);
            if( itr == MAP->end())
            {
                itr = MAP->find( "CTER:" + ATOM);
                if( itr == MAP->end())
                {
                    std::cout << "===> " << __FUNCTION__ << ": "<< RESIDUE << ":" << ATOM << " NOT found in map!" << std::endl;
                    return "";
                }
            }
        }
        return itr->second;
    }

} // end namespace store




#endif /* MAP_FUNCTIONS_T_H_ */
