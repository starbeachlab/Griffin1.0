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
//!	A collection of tools for GRIFFIN.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#include <string>
#include <iostream>

#include "../include/geometry/geometric_object_functions.h"
#include "../include/molecule/simple_molecule_format_functions.t.h"
#include "../include/molecule/density_map_functions.h"
#include "../include/storage/grid_functions.t.h"
#include "../include/utilities/ipc_functions.h"
#include "../include/readwrite/quotes.h"

float g_PDBFactor = 1.0;

int main( const int ARGC, const char *ARGV[])
{
    std::cout << "\n\n" << __FILE__ << "\n\n" << std::string( ARGV[0]) << std::endl;

    std::cout << ThrowFullHeader() << std::endl;

    std::cout << "\ncommand line:" << std::endl;

    for( int i = 0; i < ARGC; ++i)
    {
    	std::cout << ARGV[i] << "  ";
    }
    std::cout << "\n\n";

    if( ARGC == 1 || std::string( ARGV[1]) == "-help")
    {
    	std::cout << "About this program:" << std::endl;
    	std::cout << "This is a collection of tools for Griffin.\n\n" << std::endl;

    	std::cout << "========================  HELP  ===========================\n\n\n";
    	std::cout << "============  required flags  ============\n\n" << std::endl;
    	std::cout << "============  optional flags: standard  ============\n\n" << std::endl;

    	std::cout << "-help    TYPE" << std::endl;
    	std::cout << "     TYPE can be 'standard' or 'debug' (default 'standard'), writes this help" << std::endl;
    	std::cout << "     from 0 to 1 parameters\n" << std::endl;


         std::cout << "-namd2griffin     PDB   XPSF   CHARMM_PARAMETER_FILE   OUT" << std::endl;
        std::cout << "     combines PDB, X-PLOR PSF, and CHARMM force field parameter file (.prm) into GRIFFIN input file (.gfn)\n" << std::endl;
        std::cout << "     4 parameters\n" << std::endl;

        std::cout << "-gmx2griffin   GMXDUMP  CHARMM  OUT  (tip3p)" << std::endl;
		std::cout << "     combines Gromacs files into Griffin input file (.gfn) and requires the following files:" << std::endl;
		std::cout << "     GMXDUMP = Gromacs file from running 'gmxdump' on .tpr file," << std::endl;
		std::cout << "     CHARMM = CHARMM force field parameter file in Gromacs format (.itp)" << std::endl;
		std::cout << "     i.e. ffcharmm27nb.itp for older versions of Gromacs" << std::endl;
		std::cout << "     or share/gromacs/top/charmm27.ff/ffnb.itp or ffnonbonded.itp for Gromacs versions 4+" << std::endl;
		std::cout << "     add 'tip3p if you want to use TIP3P waters from CHARMM" << std::endl;
        std::cout << "     From 3 to 4 parameters\n" << std::endl;




        std::cout << "-2D_density    GFN  X_MIN  X_MAX  Y_MIN  Y_MAX  Z_MIN  Z_MAX  VOXELSIZE TYPE(e.g.'POP','TIP','all')  MODE('atom_count','mass','gaussian')  OUT " << std::endl;
		std::cout << "     to analyze the density of a (lipid) system projected onto the X,Y plane, using a grid-based method" << std::endl;
		std::cout << "     requires an input Griffin file (.gfn), the dimensions of the region of interest," << std::endl;
		std::cout << "     the grid-point spacing (VOXELSIZE)" << std::endl;
		std::cout << "     residue name (TYPE) and the type of calculation to be made (MODE) of which there are three to choose:" << std::endl;
		std::cout << "     'atom_count' average number of atoms per voxel (3D bin)" << std::endl;
		std::cout << "     'mass' average mass per voxel" << std::endl;
		std::cout << "     'gaussian' atoms in voxel are counted with gaussian function (smoother)" << std::endl;
        std::cout << "     11 parameters\n" << std::endl;

        std::cout << "-cleanup_queues     FILE   " << std::endl;
        std::cout << "     remove all ipc message queues that can remain after jobs crash " << std::endl;
        std::cout << "     mostly if simulations do not start properly old ipc message queues are floating around on the node where the head-process is running." << std::endl;
        std::cout << "     1 parameter\n" << std::endl;

        std::cout << "-geometric_objects  SURF_FILE VMD_CMD_FILE" << std::endl;
        std::cout << "     translates a Griffin SURF_FILE into a list of vmd commands that allow to visualize the geometric objects in VMD" << std::endl;
        std::cout << "     call: 'vmd -e VMD_CMD_FILE' or from the console in VMD: 'play VMD_CMD_FILE'" << std::endl;
        std::cout << "     2 parameters\n" << std::endl;

        std::cout << "-3Dgrid2grid1D      GRID3D    OUT1D" << std::endl;
        std::cout << "     for development purposes" << std::endl;
        std::cout << "     2 parameters\n" << std::endl;

        std::cout << "-1Dgrid2grid3D      GRID1D    OUT3D" << std::endl;
        std::cout << "     for development purposes" << std::endl;
        std::cout << "     2 parameters\n" << std::endl;

        std::cout << "-gfn_format" << std::endl;
        std::cout << "     write description of griffin format" << std::endl;
        std::cout << "     0 parameters\n" << std::endl;


#ifdef SQLITE
        std::cout << "\n\n===>  "  << ARGV[ 0] << "  -build_sqlite_db     NAME   #MOL-TYPES  MOLTYPE_1 .... MOLTYPE_n   #ENERGY-TYPES  ET_1 ... ET_n \n\n" << std::endl;

        std::cout << "\n\n===>  "  << ARGV[ 0] << "  -1Dgrid2sqliteDB     NAME\n\n" << std::endl;
#endif

        if( ARGC == 3 && std::string( ARGV[2]) == "debug")
        {
            std::cout << "============  optional flags: debug  ============\n\n" << std::endl;
        }
        std::cout << "=====================  HELP done  ========================\n\n" << std::endl;
        std::cout << "=====================  No flags given!  ===================== \n\n"<< std::endl;
        return 0;
    }

    std::string mode( ARGV[ 1]);
    std::ifstream read;
    std::ofstream write;

    if( ARGC == 2 && ( mode == "hidden" || mode == "-hidden"))
    {
    	std::cout << ARGV[ 0] << "  -spherical2xyz       R PHI PSI\n" << std::endl;
		std::cout << ARGV[ 0] << "  -xyz2spherical       X Y Z\n" << std::endl;
		std::cout << ARGV[ 0] << "  -grad2rad            GRAD\n" << std::endl;
		std::cout << ARGV[ 0] << "  -rad2grad            RAD\n" << std::endl;
		std::cout << ARGV[ 0] << "  -gmx2xyz             GMX_FILE\n" << std::endl;
		std::cout << ARGV[ 0] << "  -dist_check          FILE CONDITION(SMALLER/EQUAL/LARGER) THRESHOLD\n" << std::endl;
		std::cout << ARGV [0] << "  -namd2gromacs        INPUT_GRID  OUTPUT_GRID   |  unit transformation: angstroem -> nm, kcal/mol -> kjoul/mol \n" << std::endl;
		std::cout << ARGV[ 0] << "  -norm                X1 ... XN \n" << std::endl;
		std::cout << ARGV[ 0] << "  -directed_geometric_objects     SURF_FILE  VMD_CMD_FILE\n" << std::endl;
		std::cout << ARGV[ 0] << "  -clean               NODE [NODE]\n" << std::endl;
		std::cout << "\n\n\n" << std::endl;
    }



    if( mode == "-namd2griffin")
    {
    	mol::NamdIntoGriffin( ARGC, ARGV);
    }
    else if( mode == "-gmx2griffin")
    {
    	mol::GmxIntoGriffin( ARGC, ARGV);
    }
    else if( mode == "-XyzIntoPdb")
    {
    	mol::XyzIntoPdb( ARGC, ARGV);
    }
    else if( mode == "-gmxdump2xyz")
    {
    	mol::GmxIntoXyz(  ARGC, ARGV);
    }
    else if( mode == "-2D_density")
    {
    	mol::DensityMap( ARGC, ARGV);
    }
    else if( mode == "-directed_geometric_objects")
    {
    	geom::DirectedGeometricObjects( ARGC, ARGV);
    }
    else if( mode == "-geometric_objects")
    {
    	geom::GeometricObjects( ARGC, ARGV);
    }
    else if( mode == "-3Dgrid2grid1D")
    {
    	store::Convert3DgridInto1Dgrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( ARGC, ARGV);
    }
    else if( mode == "-1Dgrid2grid3D")
    {
    	store::Convert1DgridInto3Dgrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( ARGC, ARGV);
    }
    else if( mode == "-namd2gromacs")
    {
    	store::ConvertNamdGridIntoGromacsGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( ARGC, ARGV);
    }
    else if( mode == "-clean")
    {
    	std::string node;
    	for( int i = 2; i < ARGC; ++i)
    	{
    		node = ARGV[i];
    		util::CleanupQueues( node);
    	}
    }
    else if( mode == "-cleanup_queues")
    {
    	util::CleanupQueuesInFile( ARGV[2]);
    }
#ifdef SQLITE
    else if( mode == "-1Dgrid2sqliteDB")
    {
    	store::ConvertGrid1DIntoGridSqlite< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( ARGC, ARGV);
    }
    else if( mode == "-build_sqlite_db")
    {
    	store::BuildSqliteDB( ARGC, ARGV);
    }
#endif
    else if( mode == "-spherical2xyz")
    {
        StandardWrite( "transform spherical coordinates in grad to cartesian\n\n");
        float r, phi, psi;
        r = mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 2]));
        phi = math::Pi / 180.0 * mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 3]));
        psi = math::Pi / 180.0 * mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 4]));
        std::cout << "x: " << r * sin( psi) * cos( phi) << std::endl;
        std::cout << "y: " << r * sin( psi) * sin( phi) << std::endl;
        std::cout << "z: " << r * cos( psi) << std::endl << std::endl << std::endl;
    }
    else if( mode == "-xyz2spherical")
    {
        StandardWrite( "transform cartesian coordinates to spherical in grad \n\n");
        float x, y, z;
        x = mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 2]));
        y = mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 3]));
        z = mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 4]));
        std::cout << "r: " << sqrt( x*x + y*y + z*z) << std::endl;
        float phi( acos( x / ( sqrt( x*x + y*y))));
        if( y < 0)
        {
            phi = math::TwoPi - phi;
        }
        std::cout << "phi: " << phi * math::RadToGrad << std::endl;
        std::cout << "psi: " << ( math::HalfPi - atan( z / ( sqrt( x*x + y*y)))) * math::RadToGrad << std::endl;
    }
    else if( mode == "-grad2rad")
    {
        StandardWrite( "transform grad to radians");
        float value( mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 2])));
        std::cout << "rad: " << math::GradToRad * value << std::endl << std::endl;
    }
    else if( mode == "-rad2grad")
    {
        StandardWrite( "transform radians to grad");
        float value( mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ 2])));
        std::cout << "grad: " << math::RadToGrad * value << std::endl << std::endl;
    }
    else if( mode == "-norm")
    {
        std::vector< float> v( ARGC - 2);
        float sum = 0.0;
        for( unsigned int i = 0; i < v.size(); ++i)
        {
        	sum += math::Square( mystr::ConvertStringToNumericalValue< float>( std::string( ARGV[ i+2])));
        }
        std::cout << "\nnorm: " << sqrt( sum) << std::endl << std::endl;
    }
    else if( mode == "-dist_check")
    {
		std::string
			file = ARGV[2],
			condition = ARGV[3];
		float
			dist,
			threshold = mystr::ConvertStringToNumericalValue< float>( ARGV[4]);
		std::ifstream
			read( file.c_str());

			boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
				all_molecules = mol::file::ReadMoleculesInGriffinFormat( read);

		int
			size = all_molecules->size();

		for( store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >::const_iterator mol_a = all_molecules->begin(); mol_a != all_molecules->end(); ++mol_a, --size)
		{
			std::cout << size << " mols remaining" << std::endl;
			for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_a = ( *mol_a)->GetAtoms().begin(); atom_a != ( *mol_a)->GetAtoms().end(); ++atom_a)
			for( store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >::const_iterator mol_b = all_molecules->begin(); mol_b != all_molecules->end(); ++mol_b)
				for(std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_b = ( *mol_b)->GetAtoms().begin(); atom_b != ( *mol_b)->GetAtoms().end(); ++atom_b)
				if( atom_a != atom_b)
				{
					dist = math::Distance( ( *atom_a)->GetPosition(), ( *atom_b)->GetPosition());
	//			    std::cout << dist << "  ";
					if( condition == "SMALLER" && dist < threshold)
					{
					std::cout << dist << "is smaller than the threshold: " << threshold << std::endl;
					}
					else if( condition == "LARGER" && dist > threshold)
					{
					std::cout << dist << "is larger than the threshold: " << threshold << std::endl;
					}
					else if( condition == "EQUAL" && dist == threshold)
					{
					std::cout << dist << "is equal the threshold: " << threshold << std::endl;
					}
				}
		}
    }
    else if( mode == "-gfn_format")
    {
        std::cout << "GRIFFIN format: each of the following lines describes one column in the .gfn files, which are separated by one or more spaces:" << std::endl;
        std::cout << "     atom-id        [integer]" << std::endl;
        std::cout << "     atom-name      [string]" << std::endl;
        std::cout << "     residue-id     [integer]" << std::endl;
        std::cout << "     residue-name   [string]" << std::endl;
        std::cout << "     mol-id         [integer]" << std::endl;
        std::cout << "     pos-x          [angstroem/nm]" << std::endl;
        std::cout << "     pos-y          [angstroem/nm]" << std::endl;
        std::cout << "     pos-z          [angstroem/nm]" << std::endl;
        std::cout << "     mass           [atomic mass]" << std::endl;
        std::cout << "     charge         [electron charge]" << std::endl;
        std::cout << "     vdw-radius     [angstroem/nm]" << std::endl;
        std::cout << "     vdw-epsilon    [none]" << std::endl;
        std::cout << "     temperature    [K]\n" << std::endl;
        std::cout << "default units of coordinates is Angstroem, can be changed into nm (mind setting flags in the other executables accordingly!)\n" << std::endl;

    }
    else
    {
        std::cout << "no allowed option given" << std::endl;
    }

    std::cout << "toolbox done\n\n" << std::endl;
    return 0;
}

