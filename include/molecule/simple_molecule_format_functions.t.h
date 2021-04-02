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
//!	Translation functions between different formats.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef SIMPLE_MOLECULE_FORMAT_FUNCTIONS_T_H
#define SIMPLE_MOLECULE_FORMAT_FUNCTIONS_T_H

#include "simple_molecule_file_handler.t.h"

namespace mol
{

    void NamdIntoGriffin( const int ARGC, const char *ARGV[])
    {
    	std::ifstream
			read;
    	std::ofstream
			write;
		std::string
			file = ARGV[5];
    	std::cout << "creating a griffin formated molcule file from NAMD pdb and xpsf file" << std::endl;

    	if( ARGC != 6)
    	{
    		std::cerr << "\nyou did not provide the correct number of input strings." << std::endl;
    		std::cerr << "call executable without options to get help.\nbye\n\n\n" << std::endl;
    		exit( -1);
    	}

		boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
		  all_molecules( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());

		boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
			vdw_epsilon_and_radius_map( new store::Map< std::string, std::pair< float, float> >());

		std::cout << "read xpsf from " << ARGV[3] << std::endl;

		read.open( ARGV[3]);

//		assert( ending.substr( ending.size() - 5, 5) == ".xpsf");

		if( !read){ std::cout << "xpsf file not opened" << std::endl; return;}

		mol::file::ReadMoleculesFromPSF( read, *all_molecules);

		read.close();
		read.clear();


		std::cout << "reading pdb from " << ARGV[2] << std::endl;
//		ending = ARGV[2];
//		assert( ending.substr( ending.size() - 4, 4) == ".pdb");

		read.open( ARGV[2]);
		if( !read){ std::cout << "pdb file not opened" << std::endl; return;}

		mol::file::ReadNewCoordinatesFromPdb( read, *all_molecules);

		read.close();
		read.clear();

//		ending = ARGV[4];
//			assert( ending.substr( ending.size() - 4, 4) == ".prm");

		std::cout << "read vdw values into map from " << ARGV[4] << std::endl;

		read.open( ARGV[4]);
		if( !read){ std::cout << "charmm parameter file not opened" << std::endl; return;}

		mol::file::ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( read, vdw_epsilon_and_radius_map);

		read.close();
		read.clear();

		std::cout << "read vdw parameters into molecule" << std::endl;

		mol::file::ReadVdwEpsilonAndRadiusIntoMolecules( vdw_epsilon_and_radius_map, all_molecules);


		if( file.substr( file.length() - 4, 4) != ".gfn")
		{
			file += ".gfn";
		}
		write.open( file.c_str());
	
		std::cout << "write griffin formatted molecule to " << file << std::endl;
	
		if( !write){ std::cout << "output file not opened!" << std::endl; return;}
		mol::file::WriteInGriffinFormat( write, all_molecules);
		write.close();
		write.clear();
    }



    void CharmmIntoGriffin( const int ARGC, const char *ARGV[])
    {
    std::ifstream read;
    std::ofstream write;
/*
	store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >
	    all_molecules;

    	StandardWrite(  "read vdw radii and epsilon from charmm parameter file ");

		Open( read, ARGV[2]);
		mol::file::ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( read, vdw_epsilon_and_radius_map);
//            StandardWrite(  partial_charge_map;
		Close( read);


		//            StandardWrite(  "read charm force field parameter and topology files");
		StandardWrite(  "read atom type translating map from topology file");

		Open( read, ARGV[3]);
		mol::file::ReadAtomTypeTranslator( read, atom_type_map);
		Close( read);

		DebugWrite( atom_type_map);

		StandardWrite(  "read partial charges from charmm topology file");

		Open( read, ARGV[3]);
		mol::file::ReadPartialChargeMapFromCharmmTopologyFile( read, partial_charge_map);
//            StandardWrite(  partial_charge_map);
		Close( read);

		StandardWrite(  "read mass values from charmm topology file");
		Open( read, ARGV[3]);
		mol::file::ReadMassFromCharmm( read, mass_map);
		Close( read);

		Open( read, ARGV[4]);
		all_moleculess = file::ReadMoleculesFromPdb( read);
		Close( read);


		// TODO: getting rid of reading partial charges externally?
		StandardWrite(  "read partial charges into molecule ...");
		mol::file::ReadPartialChargesIntoMolecule( partial_charge_map, implicit_molecule);
		DebugWrite( "write partial_charge_map");
		DebugWrite( partial_charge_map);

                for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::iterator mol_itr = all_molecules->begin(); mol_itr != all_molecules->end(); ++mol_itr)
                {
                    mol::file::ReadMassIntoMolecule( mass_map, atom_type_map, *mol_itr);
                }


	write.open( ARGV[4]);
	if( !write){ std::cout << "output file not opened!" << std::endl; return 0;}
	mol::file::WriteInGriffinFormat( write, all_molecules);
	write.close();
	write.clear();
*/
    }


    void GmxIntoGriffin( const int ARGC, const char *ARGV[])
    {    

	std::cout << "Not yet asked FAQs:\n" << std::endl;
	std::cout << "Nothing works: check the order you provide the files!\n(was that helpful?:-)\n" << std::endl;


	assert( ARGC == 5 || ARGC == 6);

	std::ifstream 
	    read;
	std::ofstream 
	    write;
	std::string
	    file;
	boost::shared_ptr< store::Map< std::string, float> > 
	    mass_map( new store::Map< std::string, float>());
	boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
	    vdw_epsilon_and_radius_map( new store::Map< std::string, std::pair< float, float> >());
	boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >
	    implicit_molecule;
	bool
	    is_charmm_tip3p_set = false;
	
	if( ARGC > 5 && std::string( ARGV[5]) == "tip3p")
	{
	    std::cout << "use hydrogens of charmm-tip3p defintions" << std::endl;
	    is_charmm_tip3p_set = true;
	}

	Open( read, ARGV[2]);
	mol::file::ReadAtomNameAndTypeAndResidueNameAndPartialChargeFromGmxDump( read, implicit_molecule);
	mol::file::ReadNewCoordinatesFromGmxDump( read, implicit_molecule);
	Close( read);
	
	Open( read, ARGV[3]);
	mol::file::ReadMassAndVdwRadiiAndEpsilonMapsFromGromacsCharmm( read, mass_map, vdw_epsilon_and_radius_map, is_charmm_tip3p_set);
	Close( read);	

	DebugWrite( "write mass_map");
	DebugWrite( mass_map);
	DebugWrite( "write vdw map");
	DebugWrite( vdw_epsilon_and_radius_map);

	mol::file::ReadMassIntoMolecule( mass_map, implicit_molecule);

	mol::file::ReadVdwEpsilonAndRadiusIntoMolecule( vdw_epsilon_and_radius_map, implicit_molecule);

	Open( write, ARGV[4]);
	mol::file::WriteInGriffinFormat( write, implicit_molecule);
	Close( write);	
    }



    void XyzIntoPdb( const int ARGC, const char *ARGV[])
    {  
    	assert( ARGC == 5);
    	std::string
	    pdb_file = ARGV[2],
	    xyz_file = ARGV[3],
	    out = ARGV[4],
	    pdb_line,
	    xyz_line;

    	std::vector< std::string>
			value_strings;

    	std::vector< double>
			v( 3);

    	if( pdb_file.substr( pdb_file.length() - 4, 4) != ".pdb")
    	{
    		std::cerr << "second argument has to be .pdb file! (" << pdb_file.substr( pdb_file.length() - 5, 4) << ")" << std::endl;
    		exit( -1);
    	}

    	std::ifstream
			read_pdb( pdb_file.c_str()),
			read_xyz( xyz_file.c_str());

    	std::ofstream
			write( out.c_str());

    	if( ! read_pdb || ! read_xyz)
    	{
			std::cerr << "input files not correctly opened!" << std::endl;
			exit( -1);
    	}

    	int count = 0;

    	while( read_pdb && read_xyz)
    	{
			pdb_line.clear();
			std::getline( read_pdb, pdb_line);

			if( pdb_line.substr( 0, 4) != "ATOM")
			{
				continue;
			}

			xyz_line.clear();
			std::getline( read_xyz, xyz_line);

			if( pdb_line.length() > 0 && xyz_line.length() > 0)
			{
				value_strings = mystr::SplitString( xyz_line);
				pdb_line.replace( 30, 8, value_strings[2].substr( 0, 8));
				pdb_line.replace( 38, 8, value_strings[3].substr( 0, 8));
				pdb_line.replace( 46, 8, value_strings[4].substr( 0, 8));
				++count;
				write << pdb_line << std::endl;
			}
    	}
    	std::cout << count << " lines" << std::endl;
	    
    	read_pdb.close();
    	read_xyz.close();
    	write.close();
    	read_pdb.clear();
    	read_xyz.clear();
    	write.clear();
    }






    void GmxIntoXyz( const int ARGC, const char *ARGV[])
    {  

//        std::cout << "convert gmx dump file to xzy format" << std::endl;
        std::ifstream read( ARGV[ 2]);
        if( !read)
        {
            std::cout << ARGV[2] << " could not be opened!\n" << std::endl;
            return;
        }

        std::string
            line;
        bool
            found( false);

        while( !found && !read.eof())
        {
            line.clear();
            std::getline( read, line);
            if( line.substr( 0, 3) == "x (")
            {
                found = true;
//                std::cout << "block found" << std::endl;
            }
        }

        found = false;

        while( !read.eof() && !found)
        {
            line.clear();
            std::getline( read, line);
            if( line.substr( 0, 3) == "v (")
            {
                found = true;
            }
            else
            {
                size_t
                    first( line.find( "{") + 1),
                    second( line.find( ",", first)),
                    third( line.find( ",", second + 2)),
                    last( line.find( "}", third));
                first = line.find_first_not_of( " ", first);
//                std::cout << first << " " << second << "  " << third << " " << last << std::endl;
                std::cout << line.substr( first, second-first) << "   " << line.substr( second + 1, third - second - 1) << "   " << line.substr( third + 1, last - third - 1) << std::endl;
            }
        }

    }



} // end namespace

#endif // SIMPLE_MOLECULE_FORMAT_FUNCTIONS_T_H
