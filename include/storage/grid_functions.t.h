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


#ifndef GRID_FUNCTIONS_T_H_
#define GRID_FUNCTIONS_T_H_

#ifdef SQLITE
#include <sqlite3.h>
#endif
	
#include "grid1D.t.h"
#include "grid3D.t.h"
#include "../molecule/molecule_force_grid.h"
#include "../molecule/molecule_force_grid_factory.h"
#include "../utilities/enum_handler_energy.t.h"
#include "../utilities/enum_handler_functor.t.h"

namespace store
{

template< typename t_TYPE>
void
Convert
(
		const boost::shared_ptr< GridVector3D< t_TYPE> > &GRID3D,
		boost::shared_ptr< Grid1D< t_TYPE> > &GRID1D
)
{
	StandardWrite( __FUNCTION__ << " grid3D to grid1D");
	std::cout <<  __FUNCTION__ << " grid3D to grid1D" << std::endl;

	std::cout << "pass position grid: 3D => 1D" << std::endl;
	GRID1D->SetPositionGrid( GRID3D->GetPositionGrid());

	std::cout << "allocate" << std::endl;
	GRID1D->Allocate();

    store::Vector3N<size_t>
		bins = GRID3D->GetNrBins( );

    std::cout << "bins: ";
    std::copy( bins.begin(), bins.end(), std::ostream_iterator< size_t>( std::cout, "  "));
    std::cout << std::endl;

    int
		count = bins[0];

    std::cout << "layers remaining to convert:  ";
    for( size_t i = 0; i < bins[0]; ++i)
    {
    	std::cout << count-- << " ";
    	std::cout.flush();

    	for( size_t j = 0; j < bins[1]; ++j)
    		for( size_t k = 0; k < bins[2]; ++k)
    		{
    			store::Vector3N< int> indices( i, j, k);
    			GRID1D->SetGridPoint( indices, GRID3D->GetGridPoint( indices));
    		}
    }
    std::cout << std::endl;

} // end Convert()


template< typename t_TYPE>
void
Convert3DgridInto1Dgrid( const int ARGC, const char *ARGV[])
{
    std::cout << __FUNCTION__ << std::endl;
    	assert( ARGC == 4);

    	std::string
			str,
			in = ARGV[ 2],
			out = ARGV[ 3];

    	std::ifstream
			read( in.c_str());
    	std::ofstream
			write( out.c_str());

       	boost::shared_ptr< store::GridVector3D< t_TYPE> >
        	 grid3D( new store::GridVector3D< t_TYPE>());

       	boost::shared_ptr< store::Grid1D< t_TYPE> >
        	 grid1D( new store::Grid1D< t_TYPE>());

       	std::vector< std::string>
			buffer;

       	store::DefaultUniqueMap< std::string, int>
			mol_types;

	std::cout << "build enum handler " << std::endl;
	util::BuildEnumHandler< util::FunctorEnum>();
//	util::EnumHandler< util::FunctorEnum>().WriteString2ID( std::cout);
	util::BuildEnumHandler< util::EnergyEnum>();
//	util::EnumHandler< util::EnergyEnum>().WriteString2ID( std::cout);

	std::cout << "read grid3D " << std::endl;

	// molgrid
       	read >> str;
       	write << str << std::endl;
	// atomgrid
       	read >> str;
       	write << str << std::endl;
	// grid
       	read >> str;
//     	write << str << std::endl;


       	grid3D->Read( read);

       	bool
			check = false;

       	int
			c = 0;

    	while( read)
    	{
    	    std::getline( read, str);

    	    if( check && str.find( "surfpoints:") != str.npos)
    	    {
    	    	check = false;
    	    }

    	    if( check)
    	    {
    	    	mol_types.InsertNewKeyAndValue( str, c);
    	    	++c;
    	    }

    	    if( !check && str.find( "moltypes:") != str.npos)
    	    {
    	    	check = true;
    	    }

    	    buffer.push_back( str);
//    	    write << str << std::endl;
    	}

	mol_types.SetDefault( "all");

    	grid1D->SetMolTypes( mol_types);

//	std::cout << "grid read" << std::endl;
//
//	std::ofstream test( "test.txt");
//  	grid3D->Write( test);
//	test.close();
//	test.clear();

	std::cout << "convert 3D => 1D" << std::endl;

       	store::Convert( grid3D, grid1D);

	std::cout << "write grid1D" << std::endl;

    	grid1D->Write( write);

	std::cout << "write molgrid" << std::endl;

	std::copy( buffer.begin(), buffer.end(), std::ostream_iterator< std::string>( write, "\n"));
	write << std::endl;

    	write.close();
    	write.clear();
    	read.close();
    	read.clear();
}


template< typename t_TYPE>
void
Convert
(
		const boost::shared_ptr< Grid1D< t_TYPE> > &GRID1D,
		boost::shared_ptr< GridVector3D< t_TYPE> > &GRID3D
)
{
    StandardWrite( __FUNCTION__ << " grid1D to grid3D");

    std::cout << "pass position grid: 1D => 3D" << std::endl;
    GRID3D->SetPositionGrid( GRID1D->GetPositionGrid());
    std::cout << "allocate" << std::endl;
    GRID3D->InitializeDataStructure();

    store::Vector3N<size_t>
	bins = GRID3D->GetNrBins( );

    std::cout << "bins: ";
    std::copy( bins.begin(), bins.end(), std::ostream_iterator< size_t>( std::cout, "  "));
    std::cout << std::endl;

    int
	count = bins[0];

    std::cout << "layers remaining to convert:  ";
    for( size_t i = 0; i < bins[0]; ++i)
    {
    	std::cout << count-- << " ";
//	std::cout << std::endl;
    	std::cout.flush();

    	for( size_t j = 0; j < bins[1]; ++j)
	    for( size_t k = 0; k < bins[2]; ++k)
	    {
//		std::cout << i << " " << j << " " << k << std::endl;
    			store::Vector3N< int> indices( i, j, k);
			GRID3D->SetGridPoint( indices, GRID1D->GetGridPoint( indices));
	    }
    }
    std::cout << std::endl;

} // end Convert()




template< typename t_TYPE>
void
Convert1DgridInto3Dgrid( const int ARGC, const char *ARGV[])
{
    std::cout << __FUNCTION__ << std::endl;
    	assert( ARGC == 4);

    	std::string
			str,
			in = ARGV[ 2],
			out = ARGV[ 3];

    	std::ifstream
			read( in.c_str());
    	std::ofstream
			write( out.c_str());

       	boost::shared_ptr< store::GridVector3D< t_TYPE> >
        	 grid3D( new store::GridVector3D< t_TYPE>());

       	boost::shared_ptr< store::Grid1D< t_TYPE> >
        	 grid1D( new store::Grid1D< t_TYPE>());

	std::cout << "build enum handler " << std::endl;
	util::BuildEnumHandler< util::FunctorEnum>();
//	util::EnumHandler< util::FunctorEnum>().WriteString2ID( std::cout);
	util::BuildEnumHandler< util::EnergyEnum>();
//	util::EnumHandler< util::EnergyEnum>().WriteString2ID( std::cout);

	std::cout << "read grid1D " << std::endl;

	// molgrid
       	read >> str;
       	write << str << std::endl;
	// atomgrid
       	read >> str;
       	write << str << std::endl;
	// grid
       	read >> str;
//     	write << str << std::endl;

       	grid1D->Read( read);

//	std::ofstream test( "test.exe");
//  	grid1D->Write( test);
//	test.close();
//	test.clear();

	std::cout << "convert 1D => 3D" << std::endl;

       	store::Convert( grid1D, grid3D);

	std::cout << "write grid3D" << std::endl;

    	grid3D->Write( write);

	std::cout << "write molgrid" << std::endl;

	while( read)
	{
	    std::getline( read, str);
	    write << str << std::endl;
	}

    	write.close();
    	write.clear();
    	read.close();
    	read.clear();
}



template< typename t_TYPE>
void
ConvertNamdGridIntoGromacsGrid( const int ARGC, const char *ARGV[])
{
    std::cout << std::endl <<  __FUNCTION__ << std::endl;

    std::cout << "\n==============================================================================" << std::endl;
    std::cout << "\nPLEASE NOTE: this transforms a grid calculated within NAMD by a simple scaling." << std::endl;
    std::cout << "This may not be realistic because within GROMACS  parameters are different.\n" << std::endl; 
    std::cout << "==============================================================================\n\n" << std::endl;

    assert( ARGC == 4);

    std::string
	str,
	in = ARGV[ 2],
	out = ARGV[ 3];
    
    std::ifstream
	read (in.c_str());

    std::ofstream
	write;

    math::Vector3N
	vec;

    float
	scale =  4.1868 * 0.1; // (kcal/mol,A) -> (kjoule/mol,nm)
    
    boost::shared_ptr <mol::MoleculeForceGrid>
	molecule_grid (new  mol::MoleculeForceGrid ());
 


    util::BuildEnumHandler< util::FunctorEnum>();
    util::BuildEnumHandler< util::EnergyEnum>();
   
    read >> molecule_grid;

    read.close();
    read.clear();

    std::cout << "changing limits" << std::endl;

    // set limits
    vec = molecule_grid->GetAtomForceGrid()->GetPositionGrid()->GetMin();
    molecule_grid->GetSetAtomForceGrid()->GetSetPositionGrid()->SetMin (0.1 *  vec);

    vec = molecule_grid->GetAtomForceGrid()->GetPositionGrid()->GetMax();
    molecule_grid->GetSetAtomForceGrid()->GetSetPositionGrid()->SetMax (0.1 *  vec);

    vec = molecule_grid->GetAtomForceGrid()->GetPositionGrid()->GetDelta();
    molecule_grid->GetSetAtomForceGrid()->GetSetPositionGrid()->SetDelta (0.1 * vec);

//    sforce_scale = molecule_grid->GetSForceScale();
//    molecule_grid->SetSForceScale( sforce_scale * scale);

    vec = molecule_grid->GetSurfGrid().GetMinimum();
    molecule_grid->GetSetSurfGrid().SetMinimum (0.1 * vec);
   
    vec = molecule_grid->GetSurfGrid().GetDelta();
    molecule_grid->GetSetSurfGrid().SetDelta (0.1 * vec);
   

    // rescale grid points
    StandardWrite ("rescale coulomb");
    mol::factory::RescaleForces( molecule_grid, "coulomb", scale);

// TODO: check whether all contributions need to rescaled with same factor (i.e. check vdw parameter)
    StandardWrite ("rescale attractive_vdw");
    mol::factory::RescaleForces (molecule_grid, "vdw-attractive", scale);

    StandardWrite ("rescale repulsive_vdw");
    mol::factory::RescaleForces (molecule_grid, "vdw-repulsive", scale);

    StandardWrite( "rescale sforces");
    mol::factory::RescaleConstantForces (molecule_grid, molecule_grid->GetSForceScale() * scale);

    write.open (out.c_str());

    StandardWrite( "write new grid");

    write << molecule_grid;

    write.close();
    write.clear();
}





#ifdef SQLITE

void BuildSqliteDB( const std::string &NAME, const std::vector< std::string> &MOL_TYPES, const std::vector< std::string> &ENERGY_TYPES)
{

	// Conventions:
	// ID: entry id, no information about other tables
	// Locator: id to find entry in other table

    sqlite3 *db;
    char *db_err;

    sqlite3_open( NAME.c_str(), &db);

    std::string
		str,
		command = "create table 'tGridPoints' (ID integer primary key, FunctionID integer, Locator integer);";

    std::cout << command << std::endl;
    if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
    {
    	std::cout << "sql bug: " << db_err << std::endl;
    }

    command = "create table 'tForceEnergyNSP' ( ID integer primary key, Fx float, Fy float, Fz float, Energy float, NSP integer);";
    std::cout << command << std::endl;
    if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
    {
    	std::cout << "sql bug: " << db_err << std::endl;
    }


    command = "create table 'tEnergyTypes' ( Type text primary key, Locator integer);";
//    command = "create table 'EnergyMap' ( ID integer primary key, Type text, Locator integer);";
//    command = "create table 'EnergyMap' ( Locator integer primary key, Type text);";
    std::cout << command << std::endl;
    sqlite3_exec( db, command.c_str(), NULL, 0, &db_err);

    for( unsigned int i = 0; i < ENERGY_TYPES.size(); ++i)
    {
    	command = "insert into tEnergyTypes values('" + ENERGY_TYPES[i] + "'," + mystr::NumericalValueToString< int>( i) + ")";
        std::cout << command << std::endl;
        if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
        {
        	std::cout << "sql bug: " << db_err << std::endl;
        }
    }

    command = "create table 'tInteractions' ( ID integer primary key, ";
    for( unsigned int i = 0; i < ENERGY_TYPES.size(); ++i)
    {
    	command += "Locator" + mystr::NumericalValueToString( i) + " integer";
		if( i + 1 < ENERGY_TYPES.size())
		{
			command += ", ";
		}
		else
		{
			command += ");";
		}
    }
    std::cout << command << std::endl;
    if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
    {
    	std::cout << "sql bug: " << db_err << std::endl;
    }

    command = "create table 'tMolTypes' ( Type text primary key, Locator integer);";
    std::cout << command << std::endl;
    sqlite3_exec( db, command.c_str(), NULL, 0, &db_err);

    for( unsigned int i = 0; i < MOL_TYPES.size(); ++i)
    {
    	command = "insert into tMolTypes values('" + MOL_TYPES[i] + "'," + mystr::NumericalValueToString< int>( i) + ")";
        std::cout << command << std::endl;
        if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
        {
        	std::cout << "sql bug: " << db_err << std::endl;
        }
    }

    command = "create table 'tRecursive' ( ID integer primary key, MolType text, ";
    for( unsigned int i = 0; i < MOL_TYPES.size(); ++i)
    {
    	str = mystr::NumericalValueToString( i);
    	command += "Functor" + str + " integer, Locator" + str + " integer";
		if( i + 1 < MOL_TYPES.size())
		{
			command += ", ";
		}
		else
		{
			command += ");";
		}
    }
    std::cout << command << std::endl;
    if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
    {
    	std::cout << "sql bug: " << db_err << std::endl;
    }


    command = "create table 'tTypeMapped' ( ID integer primary key,  MolType text, ";
    for( unsigned int i = 0; i < MOL_TYPES.size(); ++i)
    {
    	str = mystr::NumericalValueToString( i);
    	command += "Functor" + str + " integer, Locator" + str + " integer";
		if( i + 1 < MOL_TYPES.size())
		{
			command += ", ";
		}
		else
		{
			command += ");";
		}
    }
    std::cout << command << std::endl;
    if( sqlite3_exec( db, command.c_str(), NULL, 0, &db_err))
    {
    	std::cout << "sql bug: " << db_err << std::endl;
    }



    sqlite3_close( db);
}


void BuildSqliteDB( const int ARGC, const char *ARGV[])
{
    std::cout << __FUNCTION__ << std::endl;

    int
		nr_mol_types = mystr::ConvertStringToNumericalValue< int>( ARGV[3]),
		start_energy_types = 4 + nr_mol_types,
		nr_energy_types = mystr::ConvertStringToNumericalValue< int>( ARGV[start_energy_types]);

    std::vector< std::string>
		mol_types( nr_mol_types),
		energy_types( nr_energy_types);

    std::vector< std::string>::iterator
		mol_itr = mol_types.begin(),
		e_itr = energy_types.begin();

    for( int i = 0; i < nr_mol_types; ++i, ++mol_itr)
    {
    	*mol_itr = (std::string) ARGV[ 4 + i];
    }
    for( int i = 0; i < nr_energy_types; ++i, ++e_itr)
    {
    	*e_itr = (std::string) ARGV[ start_energy_types + i + 1];
    }

    std::cout << "mol-types: ";
	std::copy( mol_types.begin(), mol_types.end(), std::ostream_iterator< std::string>( std::cout, " "));
	std::cout << std::endl;

    std::cout << "energy-types: ";
	std::copy( energy_types.begin(), energy_types.end(), std::ostream_iterator< std::string>( std::cout, " "));
	std::cout << std::endl;

    BuildSqliteDB( ARGV[2], mol_types, energy_types);
}



template< typename t_TYPE>
void ConvertGrid1DIntoGridSqlite( const int ARGC, const char *ARGV[])
{
	std::cout << __FUNCTION__ << std::endl;
	assert( ARGC == 4);

	std::string
		str,
		in = ARGV[ 2],
		out = ARGV[ 3];

	std::ifstream
		read( in.c_str());

    boost::shared_ptr< store::GridVector3D< t_TYPE> >
		 grid_sqlite( new store::GridSqlite< t_TYPE>());

    grid_sqlite->ReadGrid1D( read);
}

#endif




} // end namespace store


#endif /* GRID_FUNCTIONS_T_H_ */
