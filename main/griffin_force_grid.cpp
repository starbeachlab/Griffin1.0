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
//!	Calculates the implicit potential force field of a protein.
//! Force field will be used from daemon.
//!									 
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#include <ctime>

// order of inclusion is crucial !


#include "../include/macro/griffin_definitions.h"
#include "../include/macro/macro_functions_read_write.h"
#include "../include/macro/griffin_includes.h"


float g_PDBFactor = 1.0;

int main( const int ARGC, const char *ARGV[])
{
    std::cout << std::string( ARGV[0]) << std::endl;

    std::cout << ThrowFullHeader() << std::endl;
//    WriteSystemVariable( std::cout, "USER");
//    WriteSystemVariable( std::cout, "HOSTNAME");
//    WriteSystemVariable( std::cout, "PWD");


    /////////////////
    // DECLARATIONS
    /////////////////
    std::ofstream
        write;
    std::ifstream
        read;
    CommandLineManager
        cmd;
    boost::shared_ptr< mol::MoleculeForceGrid>
        molecule_grid;

    size_t
        nr_processes( 1),
        process_id( 0);

    std::time_t
        start,
        begin,
        now;

    std::time( &start);

    //////////////////
    // READ OPTIONS
    //////////////////
    command_line_factory::SetUpManagerForForceGrid( cmd, ARGC, ARGV);

    util::BuildEnumHandler< util::EnergyEnum>();
    util::BuildEnumHandler< util::FunctorEnum>();

    util::EnumHandler< util::FunctorEnum>().WriteString2ID();

    ///////////////////////////////////////////////////////////////
    /////            CONSTRUCT GRID                     ///////////
    ///////////////////////////////////////////////////////////////

	if( cmd.IsFlagSet( "combine_subgrids"))
	{
		StandardWrite(  "combine subgrid files from parallel calculation ...");
		mol::factory::BuildGridFromParallelFiles( molecule_grid, cmd);
	}
	else
	{
		StandardWrite(  "build grid from scratch");

		//  DECLARATIONS AND DEFAULTS  //
		boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >
			implicit_molecule;
		boost::shared_ptr< phys::ForceContainer>
			forces( new phys::ForceContainer());
		math::Vector3N
			min,
			max,
			delta( std::vector< float>( 3, 0.5));
		store::Vector3N< size_t>
			nr_bins;
		boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
			position_grid;
		float
			surface_factor = 1.0,
			surface_offset = 1.4,
			probe_radius = 1.4,
		        surf_grid_margin = 3.0;
		std::string
//			force_field_type( "charmm"),
			file;
		store::Limits3D
		  surf_grid_limits;
//		std::vector< std::string>
//			lipid_names;

		bool
		  smooth = true;
		float
			max_cutoff,
			min_force_magnitude = 1e-11,
			sforce_length = 1.0,
			sforce_reset_min_angle = math::DegreesToRadians( 110.0);

		//  READ OPTIONS  //
		if( cmd.IsFlagSet( "kjoule_nm_units"))
		{
			StandardWrite( "switching to kjoule/nm units (e.g. for gromacs)");
			StandardWrite( "====> recall to provide solute, voxel-size, limits and probe radii in nm !!!");
			StandardWrite( "====> default values are adjusted to nm !!! \n");
			phys::ElectrostaticForce().SetKjouleNmUnits();
			probe_radius *= 0.1;
			delta *= 0.1;
			sforce_length *= 41.868;
			g_PDBFactor = 10.0;  // to write surf as pdb
			surf_grid_margin *= 0.1;
		}
//		if( cmd.IsFlagSet( "lipid_types")) // todo: FINISH THIS
//		{
//			std::vector< std::string>
//				lipid_names = cmd.GetAllOptions( "lipid_types");
//		}
//		else
//		{
//			lipid_names.push_back( "POPE");
//			lipid_names.push_back( "POPC");
//		}
		if( cmd.IsFlagSet( "parallel"))
		{
			std::vector< size_t> values( cmd.GetArgumentsForFlag< size_t>( "parallel"));
			nr_processes = values[0];
			process_id = values[1];
		}
		if( cmd.IsFlagSet( "sforce_scale"))
		{
			sforce_length = cmd.GetArgumentForFlag< float>( "sforce_scale");
		}
		if( cmd.IsFlagSet( "sforce_reset_min_angle"))
		{
			sforce_reset_min_angle = math::DegreesToRadians( cmd.GetArgumentForFlag< float>( "sforce_reset_min_angle"));
		}
		if( cmd.IsFlagSet( "minimum_potential_force_magnitude"))
		{
			min_force_magnitude = cmd.GetArgumentForFlag< float>( "min_force_magnitude");
		}
		if( cmd.IsFlagSet( "voxel_size"))
		{
			std::vector< float> binning( cmd.GetArgumentsForFlag< float>( "voxel_size"));
			if( binning.size() == 1)
			{ delta.SetAll( binning[0]);}
			else if( binning.size() == 3)
			{ delta = binning;}
			else
			{
				std::cout << "pass either 1 or 3 values to \'-voxel_size\' flag" << std::endl;
				exit( -1);
			}
		}
		if( cmd.IsFlagSet( "surf_grid_limits"))
		{
			Open( read, cmd.GetArgumentStringForFlag( "surf_grid_limits"));
			read >> surf_grid_limits;
			Close( read);
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
			StandardWrite( "probe radius: " << probe_radius);
		}
		if( cmd.IsFlagSet( "magnify_radii"))
		{
			std::vector< float> vec( cmd.GetArgumentsForFlag< float>( "magnify_radii"));
			surface_offset = vec[ 0];
			if( vec.size() > 1)
			{
				surface_factor = vec[ 1];
			}
		}
		else
		{
			surface_offset = probe_radius;
		}
		StandardWrite( "vdw radii offset: " << surface_offset);
		StandardWrite( "vdw radii factor: " << surface_factor);
		StandardWrite( "probe radius: " << probe_radius);
		StandardWrite( "delta (voxel size): " << delta);
		StandardWrite( "sforce length: " << sforce_length);
		StandardWrite( "g_PDBFactor: " << g_PDBFactor);
		StandardWrite( "surf_grid_margin: " << surf_grid_margin);


		std::time( &now);
		std::cout << "\nsetup: " << std::difftime( now, start) << "s\n" << std::endl;
		std::time( &begin);

		Open( read, cmd.GetArgumentStringForFlag( "solute"));

		implicit_molecule = mol::file::ReadMoleculeInGriffinFormat( read);

		Close( read);

		std::time( &now);
		std::cout << "\nmolecule reading: " << std::difftime( now, begin) << "s" << std::endl << std::endl;
		std::time( &begin);


		//////     BUILD SURFACE OBJECT   //////
		StandardWrite( "build surface grid ...");

		boost::shared_ptr< geom::SurfGrid>
		  surf_grid( new geom::SurfGrid( /*implicit_molecule,*/ delta( 0), surf_grid_limits));


		if( cmd.IsFlagSet( "read_surf"))
		{
			Open( read, cmd.GetArgumentStringForFlag( "read_surf"));
			read >> surf_grid;
			Close( read);
		}
		else
		{

			int size = implicit_molecule->GetAtoms().size();

			float x_coo[ size];
			float y_coo[ size];
			float z_coo[ size];
			float vdw_radii[ size];
			std::vector< float> center(3);

			surf_grid->Prepare( implicit_molecule, x_coo, y_coo, z_coo, vdw_radii, center, surface_factor, surface_offset, surf_grid_margin);

			std::cout << "\nlimits of the surf grid based on the implicit molecule:" << std::endl;
			std::cout << "min: ";
			std::copy( surf_grid->GetMinimum().begin(), surf_grid->GetMinimum().end(), std::ostream_iterator< float>( std::cout, "  "));
			std::cout << std::endl;

			std::cout << "max: ";
			std::vector< float> maxy = surf_grid->CalcMax();
			std::copy( maxy.begin(), maxy.end(), std::ostream_iterator< float>( std::cout, "  "));
			std::cout << std::endl;

			std::cout << "grid-points: ";
			std::vector< int> nrelem = surf_grid->GetNrElements();
			std::copy( nrelem.begin(), nrelem.end(), std::ostream_iterator< float>( std::cout, "  "));
			std::cout << std::endl;
			std::cout << std::endl;


			boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
				user_defined_surface_objects( new store::Map< std::string, store::ShPtrVec< geom::Object> >());

			// adjust limits of surf grid for user defined surface objects
			if( cmd.IsFlagSet( "add_surf_objects"))
			{
				StandardWrite( "adjust limits of surf grid for surf objects");

				if( cmd.IsFlagSet( "kjoule_nm_units"))
				{
					StandardWrite( "==> make sure that the geometric objects are defined in nm <==");
				}

				Open( read, cmd.GetArgumentStringForFlag( "add_surf_objects"));
				user_defined_surface_objects = geom::factory::ReadIntoSurfObjectMap( read);
				Close( read);

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
		} // end else // surf grid calc

/*		if( cmd.IsFlagSet( "write_surf"))
		{
			StandardWrite( "write surf grid");
			Open( write, cmd.GetArgumentStringForFlag( "write_surf"));
//			Open( write, "surf.txt");
			write << surf_grid;
			Close( write);
		}
*/

		if( cmd.IsFlagSet( "write_surf_as_pdb"))
		{
			StandardWrite( "write surf as pdb");
//			Open( write, cmd.GetArgumentStringForFlag( "write_surf_as_pdb"));
			Open( write, "surf.pdb");
			surf_grid->WriteSurfAsPdb( write);
			Close( write);  exit(-1);
		}

		std::time( &now);
		std::cout << "\nsurface_building: " << std::difftime( now, begin) << "s\n" << std::endl;
		std::time( &begin);

		StandardWrite( "build forces ...");
		////      BUILD FORCE OBJECT  ////
		if( cmd.IsFlagSet( "force_file"))
		{
			StandardWrite( "read forces from file");
			Open( read, cmd.GetArgumentStringForFlag( "force_file"));
			phys::factory::ReadForces( read, forces, implicit_molecule);
			max_cutoff = forces->GetMaximumCutoff();
			Close( read);
		}
		else
		{
			float
				coulomb_cutoff = 18.0,
				attractive_vdw_cutoff = 12.0,
				repulsive_vdw_cutoff = 8.0;

			if( cmd.IsFlagSet( "kjoule_nm_units"))
			{
				coulomb_cutoff *= 0.1;
				attractive_vdw_cutoff *= 0.1;
				repulsive_vdw_cutoff *= 0.1;
			}

			StandardWrite( "mount default forces"); // \ncutoffs:\ncoulomb: " << coulomb_cutoff << "\nattractive vdw: " << attractive_vdw_cutoff << "\nrepulsive: " << repulsive_vdw_cutoff);
			

			forces->InsertNewKeyAndValue( "vdw-attractive", boost::shared_ptr< phys::AttractiveVanDerWaalsForce>( new phys::AttractiveVanDerWaalsForce( implicit_molecule, attractive_vdw_cutoff)));
			forces->InsertNewKeyAndValue( "vdw-repulsive", boost::shared_ptr< phys::RepulsiveVanDerWaalsForce>( new phys::RepulsiveVanDerWaalsForce( implicit_molecule, repulsive_vdw_cutoff)));
			forces->InsertNewKeyAndValue( "coulomb", boost::shared_ptr< phys::ElectrostaticForce>( new phys::ElectrostaticForce( implicit_molecule, coulomb_cutoff)));
			max_cutoff = std::max( coulomb_cutoff, std::max( repulsive_vdw_cutoff, attractive_vdw_cutoff));
		}

		if( cmd.IsFlagSet( "force_grid_limits"))
		{
			StandardWrite( "read grid limits");
			Open( read, cmd.GetArgumentStringForFlag( "force_grid_limits"));
			read >> min;
			read >> max;
			StandardWrite( "min: " << min << std::endl << "max: " << max);
			Close( read);
		}
		else
		{
			StandardWrite(  "calculate grid limits");
			mol::factory::CalculateGridLimits( implicit_molecule, min, max, delta, max_cutoff);
		}
		StandardWrite( "adjust limits to match definitions of surf grid");
		mol::factory::AdjustLimits( min, max, delta, math::Vector3N( surf_grid->GetMinimum()), math::Vector3N( surf_grid->CalcMax()), math::Vector3N( surf_grid->GetDelta()));

#ifdef STANDARD
		std::cout << "\nmin:  ";
		std::copy( min.begin(), min.end(), std::ostream_iterator< float>( std::cout, " "));
		std::cout << std::endl;
		std::cout << "max: ";
		std::copy( max.begin(), max.end(), std::ostream_iterator< float>( std::cout, " "));
		std::cout << std::endl;
		std::cout << "delta: ";
		std::copy( delta.begin(), delta.end(), std::ostream_iterator< float>( std::cout, " "));
		std::cout << std::endl << std::endl;
#endif

		if( nr_processes > 1)
		{
			StandardWrite( "parallize grid calculation");
			util::ParallelizeProcess( min, max, delta, nr_processes, process_id);
			StandardWrite( "new min and max: " << min << " " << max);
		}

		nr_bins = mol::factory::CalculateNumberOfBins( min, max, delta);

/*		if( cmd.IsFlagSet( "grid1D"))
		{
			store::Map< std::string, int>
				energy_types;
			store::DefaultUniqueMap< std::string, int>
			    mol_types;

			int
				i = 0;

			for( std::map< std::string, boost::shared_ptr< phys::Force> >::const_iterator itr = forces->begin(); itr != forces->end(); ++itr, ++i)
			{
				energy_types.insert( std::make_pair( itr->first, i));
			}

			i = 0;
			for( std::vector< std::string>::const_iterator itr = surf_grid->GetMoleculeTypes().begin(); itr != surf_grid->GetMoleculeTypes().end(); ++itr, ++i)
			{
				mol_types.insert( std::make_pair( *itr, i));
			}

			StandardWrite(  "initialize position grid");
			position_grid = boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
			( new store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( min, max, delta, nr_bins, mol_types, energy_types));
		}
		else   */
		{
			StandardWrite(  "initialize position grid");
			position_grid = boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
			( new store::GridVector3D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( min, max, delta, nr_bins));
		}

		StandardWrite(  "build atom force grid ...");
		boost::shared_ptr< mol::AtomForceGrid>
		atom_grid
		(
				new mol::AtomForceGrid
				(
						position_grid,
						surf_grid,
						forces,
						min_force_magnitude,
						sforce_length
				)
		);

		if( cmd.IsFlagSet( "fill_surf_layer"))
		{
			atom_grid->FillSurfaceLayer
			(
					surf_grid,
					forces,
					min_force_magnitude
			);
		}
		else
		{
			atom_grid->EmptySurfaceLayer
			(
					surf_grid
			);
		}

		StandardWrite( "pass atom grid to molecule grid ...");
		molecule_grid = boost::shared_ptr< mol::MoleculeForceGrid>
		(
				new mol::MoleculeForceGrid
				(
						atom_grid,
						*surf_grid,
						sforce_reset_min_angle,
						sforce_length
//						lipid_names
				)
		);

		store::TypeMappedGridPoint< mol::Atom, math::Vector3N>().SetMolTypes( surf_grid->GetMoleculeTypes());
		store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N>().SetMolTypes( surf_grid->GetMoleculeTypes());

		StandardWrite( "... build molecule grid done");
		std::time( &now);
		std::cout << "\nbuilding grid: " << std::difftime( now, begin) << "s\n" << std::endl;
		std::time( &begin);

	} // end else // if not "collect parallel"

	////////////////////////
	// WRITE GRID INTO FILE
	////////////////////////

	if( cmd.IsFlagSet( "write_force_grid"))
	{
		std::string grid_file( cmd.GetArgumentStringForFlag( "write_force_grid"));
		if( grid_file.substr( grid_file.length() - 4) == ".txt")
		{ grid_file = grid_file.substr( 0, grid_file.length() - 4);}
		if( cmd.IsFlagSet( "parallel"))
		{
			grid_file += mystr::NumericalValueToString( nr_processes) + "_" + mystr::NumericalValueToString( process_id);
		}
		grid_file += ".txt";

		StandardWrite(  "write grid to: " << grid_file);
		Open( write, grid_file);
		write.setf( std::ios::fixed, std::ios::floatfield);
		write.precision( 6);
		write.unsetf( std::ios::showpoint);
		write << molecule_grid;
		Close( write);

		// terminate for parallel processes
		if( cmd.IsFlagSet( "parallel"))
		{
			std::cout << "sub-grid " << process_id << " finished" << std::endl;
			return 0;
		}
	}

    std::time( &now);
    std::cout << "\nGriffin: total time: " << std::difftime( now, start) << "s\n\n" << std::endl;
    std::time( &begin);

    return 0;
};  // end main
