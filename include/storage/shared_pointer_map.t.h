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


#ifndef SHARED_POINTER_MAP_H_
#define SHARED_POINTER_MAP_H_

#include <map>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"


namespace store
{
    template< typename t_KEY, typename t_VALUE>
    class ShPtrMap
    : public std::map< t_KEY, boost::shared_ptr< t_VALUE> >,
    public StreamOperator
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ShPtrMap()
        : std::map< t_KEY, boost::shared_ptr< t_VALUE> >()
        {}

        //! construct from default values
        ShPtrMap( const std::vector< t_KEY> &DEFAULT)
        : std::map< t_KEY, boost::shared_ptr< t_VALUE> >()
        {
            for( typename std::vector< t_KEY>::const_iterator itr = DEFAULT.begin(); itr != DEFAULT.end(); ++itr)
            { std::map< t_KEY, boost::shared_ptr< t_VALUE> >::insert( std::make_pair( *itr, boost::shared_ptr< t_VALUE>()));}
            std::cout << GetClassName() << " build" << std::endl;
        }
        //! copy constructor
        ShPtrMap( const ShPtrMap &ORIGINAL)
        : std::map< t_KEY, boost::shared_ptr< t_VALUE> >( ORIGINAL)
        {}


        //! virtual destructor
        virtual ~ShPtrMap(){}

        //! virtual copy constructor
        virtual ShPtrMap *Clone() const{ return new ShPtrMap( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! safe but strict access by key, program dies if key is not contained!
        virtual const boost::shared_ptr< t_VALUE> &operator () ( const t_KEY &KEY) const
        {
            typename std::map< t_KEY, boost::shared_ptr< t_VALUE> >::const_iterator itr( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::find( KEY));
            if( itr != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end())
            { return itr->second;}
            std::cout << "Key: " << KEY << " not found in map!" << std::endl;
            exit( 1);
            return itr->second;
        }

        //! safe but strict access by key, program dies if key is not contained!
        virtual boost::shared_ptr< t_VALUE>  &operator () ( const t_KEY &KEY)
        {
            typename std::map< t_KEY, boost::shared_ptr< t_VALUE> >::iterator itr( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::find( KEY));
            if( itr != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end())
            { return itr->second;}
            std::cout << "Key: " << KEY << " not found in map!" << std::endl;
            exit( 1);
            return itr->second;
        }

        virtual const boost::shared_ptr< t_VALUE>  &operator []( const t_KEY &KEY) const
        { return operator()( KEY);}

        virtual boost::shared_ptr< t_VALUE>  &operator []( const t_KEY &KEY)
        { return operator()( KEY);}

//        virtual void operator += ( const std::pair< t_KEY, boost::shared_ptr< t_VALUE> > &DATA)
//        {
//            typename std::map< t_KEY, boost::shared_ptr< t_VALUE> >::iterator itr( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::find( DATA.first));
//            if( itr == std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end())
//            {
//                std::cout << DATA.first << " is not a defined key!" << std::endl;
//                exit( -1);
////                std::map< t_KEY, boost::shared_ptr< t_VALUE> >::insert( DATA);
//            }
//            else
//            {
//                *itr->second += *( DATA.second);
//            }
//        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const std::vector< t_KEY> GetKeys() const
        {
            std::vector< t_KEY> keys( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::size());
            typename std::vector< t_KEY>::iterator key_itr = keys.begin();

            for( typename std::map< t_KEY, boost::shared_ptr< t_VALUE> >::const_iterator itr = std::map< t_KEY, boost::shared_ptr< t_VALUE> >::begin(); itr != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end(); ++itr, ++key_itr)
            {
                *key_itr = itr->first;
            }
            return keys;
        }

        virtual void InsertNewKeyAndValue( const t_KEY &KEY, const boost::shared_ptr< t_VALUE> &VALUE)
        {
        	DebugWrite( __PRETTY_FUNCTION__ << " " << KEY);
            if( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::find( KEY) != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end())
            {
                std::cout << KEY << " already exists in map" << std::endl;
                exit( -1);
            }
            std::map< t_KEY, boost::shared_ptr< t_VALUE> >::insert( std::make_pair( KEY, VALUE));
        }


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
//        	std::string str;
//        	int size;
//        	t_KEY key;
//        	boost::shared_ptr< t_VALUE> value = new t_VALUE();
//
//        	STREAM >> str;
//        	if( str != GetClassName())
//        	{
//        		std::cout << "===> invalid type id (" << str << ") in " << __PRETTY_FUNCTION__ << std::endl;
//        	}
//        	STREAM >> size;
//        	for( int i = 0; i < size; ++i)
//        	{
//        		STREAM >> key >> *value;
//                if( std::map< t_KEY, boost::shared_ptr< t_VALUE> >::find( KEY) != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end())
//                {
//                    std::cout << "===> " << KEY << " already exists in " << __PRETTY_FUNCTION__ << std::endl;
//                    exit( -1);
//                }
//                std::map< t_KEY, boost::shared_ptr< t_VALUE> >::insert( std::make_pair( KEY, VALUE));
//        	}
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << " " << std::map< t_KEY, boost::shared_ptr< t_VALUE> >::size() << std::endl;
            for( typename std::map< t_KEY, boost::shared_ptr< t_VALUE> >::const_iterator itr = std::map< t_KEY, boost::shared_ptr< t_VALUE> >::begin(); itr != std::map< t_KEY, boost::shared_ptr< t_VALUE> >::end(); ++itr)
            {
                STREAM << itr->first << "  " << *itr->second << std::endl;
            }
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class ShPtrMap
} // end namespace store




#endif /* SHARED_POINTER_MAP_H_ */
