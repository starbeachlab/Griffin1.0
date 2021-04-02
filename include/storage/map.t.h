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


#ifndef MAP_H_
#define MAP_H_

#include <map>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include "../external/boost_functions.h"

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"


namespace store
{
    template< typename t_KEY, typename t_VALUE>
    class Map
    : public std::map< t_KEY, t_VALUE>,
    public StreamOperator
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Map()
        : std::map< t_KEY, t_VALUE>()
        {}

        //! construct from default values
        Map( const std::vector< t_KEY> &DEFAULT)
        : std::map< t_KEY, t_VALUE>()
        {
            for( typename std::vector< t_KEY>::const_iterator itr = DEFAULT.begin(); itr != DEFAULT.end(); ++itr)
            { std::map< t_KEY, t_VALUE>::operator []( *itr);}
        }
        //! copy constructor
        Map( const Map &ORIGINAL)
        : std::map< t_KEY, t_VALUE>( ORIGINAL)
        {}


        //! virtual destructor
        virtual ~Map(){}

        //! virtual copy constructor
        virtual Map *Clone() const{ return new Map( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! safe but strict access by key, program dies if key is not contained!
        virtual const t_VALUE &operator () ( const t_KEY &KEY) const
        {
            typename std::map< t_KEY, t_VALUE>::const_iterator itr( std::map< t_KEY, t_VALUE>::find( KEY));
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            std::cout << "====> key: <" << KEY << "> not found in map!" << std::endl;
//            if( std::map< t_KEY, t_VALUE>::size() < 50)
            {
            	std::cout << "\nmap content:" << std::endl;

				for( itr = std::map< t_KEY, t_VALUE>::begin(); itr != std::map< t_KEY, t_VALUE>::end(); ++itr)
				{
					std::cout << itr->first << ": " << itr->second << std::endl;
				}
	            std::cout << "====> key: <" << KEY << "> not found in map! END OF MAP" << std::endl;
            }
            exit( 1);
            return itr->second;
        }

        //! safe but strict access by key, program dies if key is not contained!
        virtual t_VALUE &operator () ( const t_KEY &KEY)
        {
            typename std::map< t_KEY, t_VALUE>::iterator itr( std::map< t_KEY, t_VALUE>::find( KEY));
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            std::cout << "Key: <" << KEY << "> not found in map!" << std::endl;
            exit( 1);
            return itr->second;
        }

        virtual const t_VALUE &operator []( const t_KEY &KEY) const
        { return operator()( KEY);}

        virtual t_VALUE &operator []( const t_KEY &KEY)
        { return operator()( KEY);}


//        virtual void operator += ( const std::pair< t_KEY, t_VALUE> &DATA)
//        {
//            typename std::map< t_KEY, t_VALUE>::iterator itr( std::map< t_KEY, t_VALUE>::find( DATA.first));
//            if( itr == std::map< t_KEY, t_VALUE>::end())
//            {
//                std::cout << DATA.first << " is not a defined key!" << std::endl;
//                exit( -1);
////                std::map< t_KEY, t_VALUE>::insert( DATA);
//            }
//            itr->second += DATA.second;
//        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const std::vector< t_KEY> GetKeys() const
        {
            std::vector< t_KEY> keys( std::map< t_KEY, t_VALUE>::size());
            typename std::vector< t_KEY>::iterator key_itr = keys.begin();

            for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::begin(); itr != std::map< t_KEY, t_VALUE>::end(); ++itr, ++key_itr)
            {
                *key_itr = itr->first;
            }
            return keys;
        }

        virtual bool IsValidKey( const t_KEY &KEY) const
        {
            DebugWrite( __FUNCTION__);
//            typename std::map< t_KEY, t_VALUE>::iterator itr( std::map< t_KEY, t_VALUE>::find( KEY);
            if( std::map< t_KEY, t_VALUE>::find( KEY) != std::map< t_KEY, t_VALUE>::end())
            {
                return true;
            }
            return false;
        }


        virtual void InsertNewKeyAndValue( const t_KEY &KEY, const t_VALUE &VALUE)
        {
            if( std::map< t_KEY, t_VALUE>::find( KEY) != std::map< t_KEY, t_VALUE>::end())
            {
                std::cout << __FUNCTION__ << " " << KEY << " already exists in map (value: " << VALUE << ")" << std::endl;
                exit( -1);
            }
            std::map< t_KEY, t_VALUE>::insert( std::make_pair( KEY, VALUE));
        }



        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
        	if(  !STREAM)
        	{
        		std::cout << "===> stream not open in " << __PRETTY_FUNCTION__ << std::endl;
        		exit( -1);
        	}
            t_KEY key;
            t_VALUE value;
            std::string str;
            size_t nr;
            STREAM >> str;
            if( str != GetClassName())
            {
            	std::cout << "===> class name was <" << str << "> expected: " << GetClassName() << std::endl;
            	exit( -1);
            }

            STREAM >> nr;
            for( size_t i = 0; i < nr; ++i)
            {
                STREAM >> key;
                STREAM >> value;
                std::map< t_KEY, t_VALUE>::insert( std::make_pair( key, value));
            }
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << std::endl;
            STREAM << std::map< t_KEY, t_VALUE>::size() << std::endl;
            for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::begin(); itr != std::map< t_KEY, t_VALUE>::end(); ++itr)
            {
                STREAM << itr->first << "  " << itr->second << std::endl;
            }
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
//            std::cout << __FUNCTION__ << std::endl;
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class Map
} // end namespace store




#endif /* SHARED_POINTER_MAP_H_ */
