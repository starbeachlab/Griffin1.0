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


#ifndef GEOMETRIC_OBJECT_FUNCTIONS_H
#define GEOMETRIC_OBJECT_FUNCTIONS_H

#include "surface_factory.h"


namespace geom
{
    void DirectedGeometricObjects( const int ARGC, const char *ARGV[])
    {
    	std::cout << "write vmd commands for user defined objects, both for the overall shape and the point surface" << std::endl;
//	        assert( ARGC == 4);
        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > >
            surf_object_map( new store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >());

        surf_object_map->InsertNewKeyAndValue( "all", store::ShPtrVec< geom::DirectedPointSurfaceObject>());

        StandardWrite( "read user defined objects");
        std::ifstream read;
        read.open( ARGV[2]);
        if( !read)
        {
        	std::cout << "<" << ARGV[2] << "> not opened" << std::endl;
        }
        geom::factory::ReadIntoSurfObjectMap( read, surf_object_map);
        read.close();
        read.clear();
        StandardWrite( "write vmd commands for " << surf_object_map->size() << " keys");
        std::ofstream write;
        write.open( ARGV[3]);
        for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator type_itr = surf_object_map->begin(); type_itr != surf_object_map->end(); ++type_itr)
        {
            StandardWrite(  "write vmd commands for: " << type_itr->first);

            for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator itr = type_itr->second.begin(); itr != type_itr->second.end(); ++itr)
            {
                ( *itr)->WriteVmdCommands( write);
                write << std::endl;
                write << std::endl;
            }
        }
        write.clear();
        write.close();
    }

    void GeometricObjects( const int ARGC, const char *ARGV[])
    {
        std::cout << "write vmd commands for user defined objects, both for the overall shape and the point surface" << std::endl;
//        assert( ARGC == 4);
        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
            surf_object_map( new store::Map< std::string, store::ShPtrVec< geom::Object> >());

        surf_object_map->InsertNewKeyAndValue( "all", store::ShPtrVec< geom::Object>());

        StandardWrite( "read user defined objects");
        std::ifstream read;
        read.open( ARGV[2]);
        if( !read)
        {
        	std::cout << ARGV[2] << std::endl;
        }
        surf_object_map = geom::factory::ReadIntoSurfObjectMap( read);
        read.close();
        read.clear();

        StandardWrite( "write vmd commands for " << surf_object_map->size() << " keys");
        std::ofstream write;
        write.open( ARGV[3]);
	if( !write)
	  {
	    std::cout << "output file not opened! \n"<< std::endl;
	    exit(-1);
	  }
        write << "mol new" << std::endl;
        for( std::map< std::string, store::ShPtrVec< geom::Object> >::const_iterator type_itr = surf_object_map->begin(); type_itr != surf_object_map->end(); ++type_itr)
        {
            StandardWrite(  "write vmd commands for: " << type_itr->first);

            for( std::vector< boost::shared_ptr< geom::Object> >::const_iterator itr = type_itr->second.begin(); itr != type_itr->second.end(); ++itr)
            {
                ( *itr)->WriteVmdCommands( write);
                write << std::endl;
            }
        }
        write.clear();
        write.close();

    }

} // end namespace geom

#endif //  GEOMETRIC_OBJECT_FUNCTIONS_H
