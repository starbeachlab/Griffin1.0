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


#ifndef DEFAULT_UNIQUE_MAP_H
#define DEFAULT_UNIQUE_MAP_H

#include <map>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include "../external/boost_functions.h"

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"


namespace store
{
    template< typename t_KEY, typename t_VALUE>
    class DefaultUniqueMap
    : public std::map< t_KEY, t_VALUE>,
    public StreamOperator
    {
    private:
	typename std::map< t_KEY, t_VALUE>::iterator         m_Default;
//	t_VALUE  *        m_Default;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

	DefaultUniqueMap()
	: std::map< t_KEY, t_VALUE>(),
	  m_Default()
	{}

        //! construct from default values
        DefaultUniqueMap( const t_KEY &KEY, const t_VALUE &VALUE)
        : std::map< t_KEY, t_VALUE>(),
	  m_Default()
        {
	    std::map< t_KEY, t_VALUE>::insert( std::make_pair( KEY, VALUE));
	    m_Default = std::map< t_KEY, t_VALUE>::find( KEY);
//	    m_Default = &VALUE;
        }

        //! copy constructor
        DefaultUniqueMap( const DefaultUniqueMap &ORIGINAL)
        : std::map< t_KEY, t_VALUE>( ORIGINAL)
        {}


        //! virtual destructor
        virtual ~DefaultUniqueMap(){}

        //! virtual copy constructor
        virtual DefaultUniqueMap *Clone() const{ return new DefaultUniqueMap( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! safe but strict access by key, program dies if key is not contained!
        virtual const t_VALUE &operator () ( const t_KEY &KEY) const
        {
            typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::find( KEY);
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            return m_Default->second;
        }

        //! safe but strict access by key, program dies if key is not contained!
        virtual t_VALUE &operator () ( const t_KEY &KEY)
        {
            typename std::map< t_KEY, t_VALUE>::iterator itr = std::map< t_KEY, t_VALUE>::find( KEY);
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            return m_Default->second;
        }

        virtual const t_VALUE &operator []( const t_KEY &KEY) const
        {
            typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::find( KEY);
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            return m_Default->second;
        }


        virtual t_VALUE &operator []( const t_KEY &KEY)
        {
            typename std::map< t_KEY, t_VALUE>::iterator itr = std::map< t_KEY, t_VALUE>::find( KEY);
            if( itr != std::map< t_KEY, t_VALUE>::end())
            { return itr->second;}
            return m_Default->second;
        }



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


	// both key and value are unique!
        virtual void InsertNewKeyAndValue( const t_KEY &KEY, const t_VALUE &VALUE)
        {
            for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::begin(); itr != std::map< t_KEY, t_VALUE>::end(); ++itr)
            {
                if( KEY == itr->first || VALUE == itr->second)
		{
		    std::cout << "====> KEY or VALUE already exists in DefaultUniqueMap: (" << KEY << "," << VALUE << ")" << std::endl;
		    return;
		}
            }
            std::map< t_KEY, t_VALUE>::insert( std::make_pair( KEY, VALUE));
        }

	void SetDefault( const t_KEY &KEY)
	{
            typename std::map< t_KEY, t_VALUE>::iterator itr = std::map< t_KEY, t_VALUE>::find( KEY);
            if( itr != std::map< t_KEY, t_VALUE>::end())
	    {
		m_Default = itr;
	    }
	    std::cout << "====> key was not found in map, therefore could not be set as default: <" << KEY << ">\nthe map:" <<std::endl;
            for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::begin(); itr != std::map< t_KEY, t_VALUE>::end(); ++itr)
	    {
		std::cout << itr->first << " " << itr->second << std::endl;
	    }
	}

	typename std::map< t_KEY, t_VALUE>::iterator &FindValue( const t_VALUE &VALUE)
	{
	    typename std::map< t_KEY, t_VALUE>::const_iterator itr = std::map< t_KEY, t_VALUE>::begin();
	    for( ; itr != std::map< t_KEY, t_VALUE>::end(); ++itr)
            {
                if( VALUE == itr->second)
		{
		    return itr;
		}
            }
	    return itr;
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
	    std::cout << "size: " << nr << std::endl;
            for( size_t i = 0; i < nr; ++i)
            {
                STREAM >> key;
                STREAM >> value;
                std::map< t_KEY, t_VALUE>::insert( std::make_pair( key, value));
            }
	    STREAM >> str >> key >>  value;

	    std::cout << str << "  " << key << " "<<value << std::endl;

	    m_Default = std::map< t_KEY, t_VALUE>::find( key);

	    if( m_Default == std::map< t_KEY, t_VALUE>::end())
	    {
		std::cout << "====> default values could not be assigned to default unique map when reading from file" << std::endl;
		Write( std::cout);
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
	    STREAM << "default: " << m_Default->first << " " << m_Default->second << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
//            std::cout << __FUNCTION__ << std::endl;
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class DefaultUniqueMap
} // end namespace store




#endif /* DEFAULT_UNIQUE_MAP_H */
