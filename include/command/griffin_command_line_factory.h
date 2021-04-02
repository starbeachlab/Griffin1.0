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
//!	Definitions of flags for force grid calculation and daemon.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef COMMAND_LINE_FACTORY_H_
#define COMMAND_LINE_FACTORY_H_

#include "command_line_manager.h"



namespace command_line_factory
{
    inline
    void SetUpManagerForForceGrid( CommandLineManager &MANAGER, const size_t &ARGC, const char *ARGV[])
    {
        /////////////////////////
        //  DEFINE ALLOWED FLAGS
        /////////////////////////

        MANAGER.DefineFlag( "write_force_grid", "FILE", "output force-field grid file", 1, e_Required);
        MANAGER.DefineFlag( "solute", "FILE", "Griffin formatted input (.gfn) of protein molecule", 1, e_Required);


        MANAGER.DefineFlag( "parallel", "N  SUBGRID", "the total number of subgrids to be calculated is N**3;\n     the current job will compute one subgrid, identified with SUBGRID (integer between 0 and N**3-1)", 2);
        MANAGER.DefineFlag( "combine_subgrids", "FILE_PREFIX  N", "combines the subgrids that have been calculated in parallel into one complete grid;\n     FILE_PREFIX is the subgrid filename (before “_SUBGRID.txt”)", 2);
        MANAGER.DefineFlag( "voxel_size", "VALUE (default 0.5)", "bin sizes of the force and surf grid", 1);
        MANAGER.DefineFlag( "kjoule_nm_units", "set this flag when input coordinates are in nanometers and energies are to be calculated in kJoules\n     (e.g. in Gromacs) rather than the default kcal and Angstroem units");
        MANAGER.DefineFlag( "add_surf_objects", "FILE", "file specifying additional geometric objects for altering the protein surface", 1);
        MANAGER.DefineFlag( "magnify_radii", "OFFSET  [FACTOR]", "offset (in Angstroems) and optional scaling factor for increasing VdW radius for calculating protein surface\n     (defaults: 1.4, 1.0)", 1, 2);
        MANAGER.DefineFlag( "force_file", "FILE", "file specifying which physical interactions to include, and their cutoff values (in Angstroems)\n     by default these are: coulomb 18, vdw-attractive 12, vdw-repulsive 8");
        MANAGER.DefineFlag( "force_grid_limits", "FILE", "allows to define a smaller force grid than bilayer box;\n     FILE contains force grid dimensions; useful for reducing memory", 1);
        MANAGER.DefineFlag( "write_surf_as_pdb", "PDB", "writes surface points as atoms into PDB file", 1);

        MANAGER.DefineFlag( "smooth_surf", "PROBE_RADIUS", "calculates smooth protein surface (used by default) with variable probe radius (default: 1.4 Angstroems)", 0, 1, e_Advanced);
        MANAGER.DefineFlag( "vdw_surf", "PROBE_RADIUS", "calculates vdW molecular rather than smooth surface, with variable probe radius (default: 1.4 Angstroems)", 0, 1, e_Advanced);
        MANAGER.DefineFlag( "read_surf", "FILE", "reads Griffin-formatted, pre-calculated protein surface; by default the surface is calculated from the input PDB", 1, e_Advanced);
        MANAGER.DefineFlag( "fill_surf_layer", "allow layer of grid points at protein surface to contain physical interaction forces;\n     default is to have no physical forces at this layer", 0, e_Advanced);
        MANAGER.DefineFlag( "minimum_potential_force_magnitude", "VALUE", "grid forces that have a magnitude lower than the specified value will be skipped (default: 1e-11)", 1, e_Advanced);


//        MANAGER.DefineFlag( "surf_grid_limits", "FILE", "containing xmin, xmax, ..., zmax of grid describing the surface/volume of the implicit molecule", 1, e_Advanced);
//        MANAGER.DefineFlag( "smooth_force_transition", "START_MAGNITUDE  NR_STEPS limits interaction forces to maximum magnitude and slowly increases them", 2, e_Advanced);

        ///////////////////////////////////////////////////////////
        //  READ OPTIONS AND CHECK FOR CONSISTENCY WITH DEFINTIONS
        ///////////////////////////////////////////////////////////

        MANAGER.SetGeneralIntroduction( "This program computes a force field created by a protein structure, and maps this force field on a 3D grid. This force field derives from the electrostatic and van der Waals interactions between the protein atoms and a probe particle, outside of the protein volume; it also includes so-called surface forces inside of the protein volume. The grid also includes energies associated with those forces. To speed up the calculation of this force grid, subsections of the grid can be calculated on different computers, in parallel. This program can also combine a set of such sub-grids into a complete grid.\n\n");

        std::cout << MANAGER.GetGeneralIntroduction() << std::endl;

        MANAGER.IfNoFlagIsSetWriteHelpAndExit( ARGC, std::cout);

        if( !MANAGER.ReadAndCheckFlags( ARGC, ARGV))
        {
            MANAGER.WriteHelp( std::cout);
            std::cout << ">>>>  Errors in reading flags  <<<<" << std::endl;
            exit( -1);
        }

        //////////////////////
        // CONSISTENCY CHECKS
        //////////////////////
        if( MANAGER.IsFlagSet( "parallel") && MANAGER.IsFlagSet( "combine_subgrids"))
        {
        	std::cout << "Use either 'parallel' or 'combine_subgrids'. Both together make no sense.\n" << std::endl;
        	exit(-1);
        }
    }



    inline
    void SetUpManagerForDaemon( CommandLineManager &MANAGER, const size_t &ARGC, const char *ARGV[])
    {

		/////////////////////////
		//  DEFINE ALLOWED FLAGS
		/////////////////////////

		MANAGER.DefineFlag( "force_grid", "FILE", "reads Griffin-formatted force-grid from FILE ", 1, e_Required);
        MANAGER.DefineFlag( "unit_system", "TYPE", "coordinate and energy unit system used in the input .gfn file\n     TYPE can be 'namd' or 'gromacs'", 1, e_Required);
		MANAGER.DefineFlag( "lipid_types", "RESNAME  [RESNAME]", "list of lipid types in the system (needed for redirection)", 1, 9, e_Required);


		MANAGER.DefineFlag( "logfile", "FILE","standard output when daemon is running (default: 'griffin_daemon.log') ", 1);
		MANAGER.DefineFlag( "sforce_scale", "VALUE", "rescale surface forces by VALUE (default: 1)", 1);
		MANAGER.DefineFlag( "rescale", "TERM  VALUE  [TERM  VALUE]", "rescale one or more interaction terms by VALUE;\n     TERM can be 'vdw-attractive', 'vdw-repulsive' or 'coulumb'\n     e.g.: -rescale vdw-attractive 1.5 coulomb 10", 2, 6);
		MANAGER.DefineFlag( "sforce_exclude_hydrogens", "surface forces will not be calculated for buried hydrogens");
		MANAGER.DefineFlag( "neighbor_list_update", "VALUE", "a list of atoms with non-negligible forces is updated every VALUE steps;\n     subsequent force calculations are restricted to this list only.", 1);
		MANAGER.DefineFlag( "sforce_reset_min_angle", "VALUE", "the surface force acting on a lipid atom inside the protein volume is redirected\n     towards the center of the corresponding molecule if the angle (in degrees)\n     between the surface-force vector and the vector to the molecular center exceeds the threshold VALUE\n     (default: 110 degrees) ", 1);

		//		MANAGER.DefineFlag( "interpolate", "interpolate forces from neighboring grid points");
		//     	 MANAGER.DefineFlag( "smooth_force_transition", "START_MAGNITUDE  NR_STEPS limits interaction forces to maximum magnitude and slowly increases them", 2, e_Advanced);
		//  	 MANAGER.DefineFlag( "", "", 1, 1);

		///////////////////////////////////////////////////////////
		//  READ OPTIONS AND CHECK FOR CONSISTENCY WITH DEFINTIONS
		///////////////////////////////////////////////////////////

		MANAGER.SetGeneralIntroduction( "This program creates a daemon to mediate access to the GRIFFIN force-grid in the course of a molecular dynamics simulation. Use griffin_messenger.exe to communicate with the daemon.\n\n");

        std::cout << MANAGER.GetGeneralIntroduction() << std::endl;

		MANAGER.IfNoFlagIsSetWriteHelpAndExit( ARGC, std::cout);

		if( !MANAGER.ReadAndCheckFlags( ARGC, ARGV))
		{
			MANAGER.WriteHelp( std::cout);
			std::cout << ">>>>  Errors in reading flags  <<<<" << std::endl;
			exit( -1);
		}

    }


} // end namespace command_line_factory




#endif /* COMMAND_LINE_FACTORY_H_ */
