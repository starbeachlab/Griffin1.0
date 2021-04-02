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
//!	Carves a lipid bilayer to fit a protein.
//!	It will leave some lipid tails sticking into the protein volume.
//! This returns the initial conformation for the optimization using the daemon.
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#include <ctime>

#include "../include/macro/griffin_definitions.h"
#include "../include/macro/macro_functions_read_write.h"
#include "../include/command/membrane_carver_command_line_factory.h"
#include "../include/storage/map.t.h"
#include "../include/storage/limits3D.h"
#include "../include/geometry/surf_grid.h"
#include "../include/geometry/surface_factory.h"
#include "../include/molecule/simple_molecule_file_handler.t.h"
#include "../include/molecule/molecule_factory.h"
#include "../include/geometry/volume_functions.h"
#include "../include/readwrite/quotes.h"


float g_PDBFactor = 1.0;


int main( const int ARGC, const char *ARGV[])
{
    std::cout << "\n\n" << std::string( ARGV[0]) << std::endl;

    std::cout << ThrowFullHeader() << std::endl;



    /////////////////
    // DECLARATIONS
    /////////////////
    std::ofstream
        write;
    std::ifstream
        read;
    CommandLineManager
        cmd;

    boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
        explicit_molecules;
    boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >
        implicit_molecule;

    std::time_t
        start,
        begin,
        now;

    math::Vector3N
        min,
        max,
        delta;
    store::Vector3N< size_t>
        nr_bins;
    std::string
        file;

    float
      surf_voxel_size = 0.2,
      probe_radius = 1.4,
      offset = 0.0,
//      offset = 1.4,
      factor = 1.0,
	surf_grid_margin = 3.0;

    std::string
        str;
    std::vector< std::string>
        lipid_names,
        sol_names,
        nr_lipids( 2, ""),
        tmp_vec;
    std::vector< std::pair< std::string, int> >
        keep_mols,
        remove_mols;
    bool
      smooth = true;

    store::Limits3D
      surf_grid_limits;


    std::time( &start);
    std::time( &begin);

    //////////////////
    // READ OPTIONS
    //////////////////

    command_line_factory::SetUpManagerForMembraneCarver( cmd, ARGC, ARGV);

    std::cout << std::endl;

    if( cmd.IsFlagSet( "kjoule_nm_units"))
    {
    	std::cout << "\nunits set to kjoule and nm !" << std::endl;
     	std::cout << "====> recall to provide membrane, solute, voxel-size, delta, limits and probe radii in nm !!!" << std::endl;
     	std::cout << "====> default values are adjusted to nm !!! \n" << std::endl;
     	std::cout << "====> set this flag also during force-grid calculation \n" << std::endl;
     	g_PDBFactor = 10;            // for writing pdb files in angstroem
     	surf_voxel_size *= 0.1;      // adjust default to nm
     	probe_radius *= 0.1;         // adjust default to nm
	surf_grid_margin *= 0.1;
    }

    if( cmd.IsFlagSet( "surf_voxel_size"))
    {
        surf_voxel_size = cmd.GetArgumentForFlag< float>( "surf_voxel_size");
        StandardWrite( "surf voxel size set to " << surf_voxel_size);
    }
    if( cmd.IsFlagSet( "smooth_surf"))
    {
        probe_radius = cmd.GetArgumentForFlag< float>( "smooth_surf");
        StandardWrite( "probe radius: " << probe_radius);
    }
    if( cmd.IsFlagSet( "vdw_surf"))
    {
	  smooth = false;
	  probe_radius = cmd.GetArgumentForFlag< float>( "vdw_surf");
	  StandardWrite( "water radius: " << probe_radius);
    }

    if( cmd.IsFlagSet( "magnify_radii"))
    {
        std::vector< float> vec( cmd.GetArgumentsForFlag< float>( "magnify_radii"));
        offset = vec[ 0];
        if( vec.size() > 1)
        {
            factor = vec[ 1];
        }
    }
    StandardWrite( "probe radius: " << probe_radius);
    StandardWrite( "vdw radii offset: " << offset);
    StandardWrite( "vdw radii factor: " << factor);


    std::cout << "\nlipid names: ";
	std::vector< std::string> arguments = cmd.GetAllOptions( "lipid_types");
	for( size_t i = 0; i < arguments.size(); ++i)
	{
	    std::cout << arguments[i] << "  ";
		lipid_names.push_back( arguments[ i]);
	}
	std::cout << std::endl;

	std::cout << "\nsolvent names: ";
	arguments = cmd.GetAllOptions( "solvent_names");
	for( size_t i = 0; i < arguments.size(); ++i)
	{
	    std::cout << arguments[i] << "  ";
		sol_names.push_back( arguments[ i]);
	}
	std::cout << std::endl;

    if( cmd.IsFlagSet( "surf_grid_limits"))
    {
		float
		  xmin, xmax,
		  ymin, ymax,
		  zmin, zmax;
		read >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
		std::cout << "surf grid limits: xmin: " << xmin << " xmax: " << xmax << " ymin: " << ymin << " ymax: " << ymax << " zmin: " << zmin << " zmax: " << zmax << std::endl; 
        surf_grid_limits = store::Limits3D( xmin, xmax, ymin, ymax, zmin, zmax);
    }



    if( cmd.IsFlagSet( "nr_lipids_to_be_removed"))
    {
    	nr_lipids = cmd.GetAllOptions( "nr_lipids_to_be_removed");
	std::cout << "\nnumber lipids to be removed (per leaflet): ";
	std::copy( nr_lipids.begin(), nr_lipids.end(), std::ostream_iterator< std::string>( std::cout, "  "));
	std::cout << std::endl;
    	assert( nr_lipids.size() == 2);
    }

    if( cmd.IsFlagSet( "keep_molecules"))
    {
    	StandardWrite( "keep molecules: ");
    	tmp_vec = cmd.GetAllOptions( "keep_molecules");
    	assert( tmp_vec.size() % 2 == 0);
    	for( int i = 0; i < 0.5 * tmp_vec.size(); ++i)
    	{
    		StandardWrite( tmp_vec[ 2*i] << "  " << tmp_vec[ 2*i+1]);
    		keep_mols.push_back( std::make_pair( tmp_vec[ 2*i], mystr::ConvertStringToNumericalValue< int>( tmp_vec[ 2*i+1])));
    	}
    }

    if( cmd.IsFlagSet( "remove_molecules"))
    {
    	StandardWrite( "remove molecules:");
    	tmp_vec = cmd.GetAllOptions( "remove_molecules");
    	assert( tmp_vec.size() % 2 == 0);
    	for( int i = 0; i < 0.5 * tmp_vec.size(); ++i)
    	{
    		StandardWrite( tmp_vec[ 2*i] << "  " << tmp_vec[ 2*i+1]);
    		remove_mols.push_back( std::make_pair( tmp_vec[ 2*i], mystr::ConvertStringToNumericalValue< int>( tmp_vec[ 2*i+1])));
    	}
    }

    ////////////////////////////////////////////////////
    /////            MAIN                    ///////////
    ////////////////////////////////////////////////////

    boost::shared_ptr< mol::Membrane>
        membrane( new mol::Membrane( lipid_names, sol_names));

    StandardWrite( "\n\ncheck and build membrane");

    membrane->AnalyseAndBuild( cmd);

    std::time( &now);
    std::cout << "membrane reading: " << std::difftime( now, begin) << "s" << std::endl;
    std::time( &begin);


//    StandardWrite( "membrane limits: " << membrane->Limits());

//    membrane->Write( write, cmd);


//    if( cmd.IsFlagSet( "write_membrane"))
//    {
//        StandardWrite( "write membrane ...");
//        WriteObject( write, membrane, cmd.GetArgumentStringForFlag( "write_membrane"));
//    }


    if( cmd.IsFlagSet( "check_dimensions_only"))
    {
        std::cout << "limit check only, no carving, bye" << std::endl << std::endl;
        return 0;
    }


    if( !cmd.IsFlagSet( "read_surf_from_file"))
    {
    	StandardWrite(  "build implicit molecule ...");

    	read.open( cmd.GetArgumentStringForFlag( "solute").c_str());

    	implicit_molecule = mol::file::ReadMoleculeInGriffinFormat( read);

		read.close();
		read.clear();

        DebugWrite( "implicit molecule: " << implicit_molecule);

        std::time( &now);
        std::cout << "molecule built: " << std::difftime( now, begin) << "s" << std::endl;
        std::time( &begin);
    }

    StandardWrite( "build surface grid ...");

    boost::shared_ptr< geom::SurfGrid>
      surf_grid( new geom::SurfGrid( surf_voxel_size, surf_grid_limits));

    int size = implicit_molecule->GetAtoms().size();

    float x_coo[ size];
    float y_coo[ size];
    float z_coo[ size];
    float vdw_radii[ size];
    std::vector< float> center(3);

    surf_grid->Prepare( implicit_molecule, x_coo, y_coo, z_coo, vdw_radii, center, factor, offset, surf_grid_margin);

    std::cout << "\nlimits of the surf grid based on the implicit molecule:" << std::endl;
    std::cout << "min: ";
    std::copy( surf_grid->GetMinimum().begin(), surf_grid->GetMinimum().end(), std::ostream_iterator< float>( std::cout, "  "));
    std::cout << std::endl;

    std::cout << "max: ";
    std::vector< float> maxy = surf_grid->CalcMax();
    std::copy( maxy.begin(), maxy.end(), std::ostream_iterator< float>( std::cout, "  "));
    std::cout << std::endl;

    std::cout << "nr grid-points: ";
    std::vector< int> nrelem = surf_grid->GetNrElements();
    std::copy( nrelem.begin(), nrelem.end(), std::ostream_iterator< float>( std::cout, "  "));
    std::cout << std::endl;
    std::cout << std::endl;


    boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
        user_defined_surface_objects( new store::Map< std::string, store::ShPtrVec< geom::Object> >());

    // adjust limits of surf grid for user defined surface objects
    if( cmd.IsFlagSet( "add_surf_objects"))
    {
        StandardWrite( "adjust limits of surf grid for surf objects ...");

        Open( read, cmd.GetArgumentStringForFlag( "add_surf_objects"));
        user_defined_surface_objects = geom::factory::ReadIntoSurfObjectMap( read);
        Close( read);

        DebugWrite( "objects: " << user_defined_surface_objects);

        surf_grid->AdjustLimitsForSurfaceObjects( user_defined_surface_objects);

        std::cout << "\nlimits after adding surf objects: " << std::endl;
        std::cout << "min: ";
        std::copy( surf_grid->GetMinimum().begin(), surf_grid->GetMinimum().end(), std::ostream_iterator< float>( std::cout, "  "));
        std::cout << std::endl;
        std::cout << "max: ";
        std::vector< float> maxb = surf_grid->CalcMax();
        std::copy( maxb.begin(), maxb.end(), std::ostream_iterator< float>( std::cout, "  "));
        std::cout << std::endl;
        std::cout << "grid-points: ";
        std::vector< int> nrelemb = surf_grid->GetNrElements();
        std::copy( nrelemb.begin(), nrelemb.end(), std::ostream_iterator< float>( std::cout, "  "));
        std::cout << std::endl;
    }

    center = surf_grid->CalcCenter();

    surf_grid->Build( smooth, probe_radius, size, x_coo, y_coo, z_coo, vdw_radii, center);


    // add user defined surface objects
    if( cmd.IsFlagSet( "add_surf_objects"))
    {
        StandardWrite( "add surface objects ...");

        surf_grid->AddSurfaceObjects( user_defined_surface_objects);
    }

    StandardWrite( "find surf points within grid ...");

    surf_grid->IndicesOfSurfaceGridPoints();

    if( cmd.IsFlagSet( "write_surf"))
    {
    	StandardWrite( "write surf grid");
    	Open( write, cmd.GetArgumentStringForFlag( "write_surf"));
    	surf_grid->Write( write);
    	Close( write);
    }


    if( cmd.IsFlagSet( "write_surf_as_pdb"))
    {
       StandardWrite( "write surf as pdb");
       Open( write, cmd.GetArgumentStringForFlag( "write_surf_as_pdb"));
       surf_grid->WriteSurfAsPdb( write);
       Close( write);
    }

    std::time( &now);
    std::cout << "surface_building: " << std::difftime( now, begin) << " s" << std::endl;
    std::time( &begin);

    StandardWrite( "\ncarve membrane");

    if( cmd.IsFlagSet( "write_erased_pdb"))
    {
    	Open( write, cmd.GetArgumentStringForFlag( "write_erased_pdb"));
    }
    explicit_molecules = geom::RemoveMoleculesInImplicitMembraneVolume( write, membrane, surf_grid, nr_lipids[0], nr_lipids[1], keep_mols, remove_mols);
    Close( write);

    std::time( &now);
    std::cout << "membrane carving: " << std::difftime( now, begin) << " s" << std::endl;
    std::time( &begin);

//    if( cmd.IsFlagSet( "write_explicit_molecules"))
//    {
//        StandardWrite( "write explicit molecules ...");
//        write.open( ( cmd.GetArgumentStringForFlag( "write_explicit_molecules")).c_str());
//        Check( write);
//        write << explicit_molecules;
//        write.close();
//        write.clear();
//        StandardWrite( "... explicit molecules written");
//    }
    if( cmd.IsFlagSet( "write_carved_pdb"))
    {
        StandardWrite( "write explicit molecules in pdb format ...");
        Open( write, cmd.GetArgumentStringForFlag( "write_carved_pdb"));
        mol::file::SortByTypeAndID( explicit_molecules);
        mol::file::WriteToPdb< mol::Atom>( write, explicit_molecules);
        Close( write);
        StandardWrite( "... explicit molecules written");
    }

    if( cmd.IsFlagSet( "write_carved"))
    {
        StandardWrite( "write explicit molecules in griffin format ...");
        Open( write, cmd.GetArgumentStringForFlag( "write_carved"));
        mol::file::SortByTypeAndID( explicit_molecules);
        mol::file::WriteInGriffinFormat( write, explicit_molecules);
        Close( write);
        StandardWrite( "... explicit molecules written");
    }
    StandardWrite( "... membrane carving done");

    std::time( &now);
    std::cout << "carving: total time: " << std::difftime( now, start) << " s" << std::endl;
    std::time( &begin);

    return 0;
};  // end main
