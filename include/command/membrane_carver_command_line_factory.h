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
//!	Flag definitions for the carver.exe
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
    void SetUpManagerForMembraneCarver( CommandLineManager &MANAGER, const size_t &ARGC, const char *ARGV[])
    {
        /////////////////////////
        //  DEFINE ALLOWED FLAGS
        /////////////////////////
        MANAGER.DefineFlag( "membrane", "FILE", "Griffin formatted input (.gfn) of the membrane/solvent system to be carved", 1, e_Required);
        MANAGER.DefineFlag( "solute", "FILE", "Griffin formatted input (.gfn) of molecule to be embedded in membrane, e.g. protein", 1, e_Required);
        MANAGER.DefineFlag( "lipid_types", "RESNAME [RESNAME]", "lipid types in membrane, e.g.: POPE POPC", 1, 20, e_Required);
        MANAGER.DefineFlag( "solvent_names", "RESNAME [RESNAME]", "non-lipid molecule types, e.g.: TIP3, NACL", 1, 20, e_Required);
        MANAGER.DefineFlag( "write_carved", "FILE", "Griffin-formatted output (.gfn) with carved membrane/solvent", 1, e_Required);


        MANAGER.DefineFlag( "dimensions", "XMIN  XMAX  YMIN  YMAX  ZMIN  ZMAX  Z_BILAYER_CENTER  Z_BILAYER_THICKNESS", "dimensions of the periodic box containing the membrane, membrane center and approximate thickness", 8);
        MANAGER.DefineFlag( "add_surf_objects", "FILE", "file specifying additional geometric objects for altering the protein surface", 1);
        MANAGER.DefineFlag( "kjoule_nm_units", "set this flag when input coordinates are in nanometers and not Angstroem (e.g. in Gromacs) AND want to write a PDB file", 0);
        MANAGER.DefineFlag( "write_carved_pdb", "PDB", "PDB-formatted output file with carved membrane/solvent molecules", 1);
        MANAGER.DefineFlag( "write_erased_pdb", "PDB", "PDB-formatted output file with deleted membrane/solvent molecules", 1);
        MANAGER.DefineFlag( "smooth_surf", "PROBE_RADIUS", "calculates smooth protein surface: used by default with default probe radius of 1.4 Angstroems", 1);
        MANAGER.DefineFlag( "surf_voxel_size", "VALUE", "resolution of calculated protein surface (default: 0.2 Angstroem)", 1);
        MANAGER.DefineFlag( "write_surf_as_pdb", "PDB", "writes surface points as atoms into PDB file", 1);
        MANAGER.DefineFlag( "keep_molecules", "RESNAME  RESNUM  [RESNAME  RESNUM]", "specify individual molecules that shall NOT be removed during carving, e.g.: POPE 232 POPE 123 TIP3 876", 2, 200);
        MANAGER.DefineFlag( "remove_molecules", "RESNAME  RESNUM  [RESNAME  RESNUM]", "specify individual molecules that shall be removed during carving, e.g.: POPE 232 POPE 123 TIP3 876", 2, 200);


        MANAGER.DefineFlag( "nr_lipids_to_be_removed", "VALUE_LOWER_LAYER  VALUE_UPPER_LAYER", "modify number of lipids that shall be removed in each layer.\n     if '-n' is given, the number calculated automatically is reduced by n,\n     if '+n' is given, the number is increased by n\n     if 'n' is given, n lipids will be removed", 2, e_Advanced);
        MANAGER.DefineFlag( "write_histograms", "FILE", "creates all-atom distribution profiles in X, Y, and Z", 1, e_Advanced);
        MANAGER.DefineFlag( "check_dimensions_only", "exit after all system dimensions are calculated, no carving", e_Advanced);
        MANAGER.DefineFlag( "surf_grid_limits", "XMIN  XMAX  YMIN  YMAX  ZMIN  ZMAX", "grid describing the surface/volume of the protein (plus objects); overrides values calculated automatically", 6, e_Advanced);
        MANAGER.DefineFlag( "auto_dimensions",  "MODE", "automated calculation of the system dimensions, which overrides values specified in '-dimensions' flag. MODE can be\n     'max' (maximal values), 'average' (mean values), 'standard-deviation' (mean plus standard deviation)", 1, e_Advanced);
        MANAGER.DefineFlag( "magnify_radii", "OFFSET  FACTOR", "offset (in Angstroems) and scaling factor for atomic radii (default: 0.0, 1.0)", 2, e_Advanced);
        MANAGER.DefineFlag( "use_generic_radii", "", "alternate radii based on  element types for comparing surfaces with e.g. VMD", e_Advanced);
        MANAGER.DefineFlag( "vdw_surf", "PROBE_RADIUS", "calculates vdw surfaces rather than smooth (default 1.4 Angstroems)", 0, 1, e_Advanced);
        MANAGER.DefineFlag( "write_surf", "FILE", "writes grid-based surface of the protein to be embedded in a Griffin-specific format", 1, e_Advanced);
        MANAGER.DefineFlag( "read_surf_from_file", "FILE", "reads Griffin-specific formatted pre-calculated grid-based surface of the protein to be embedded", 1, e_Advanced);


        MANAGER.SetGeneralIntroduction( "This algorithm carves out lipid and water molecules from a membrane system to make room for a protein. The number and selection of molecules to be removed is based on the lipid density in the upper and lower leaflets, the volume of the protein in these layers, and the degree to which the molecules overlap with the protein.\n\n");

        std::cout << MANAGER.GetGeneralIntroduction() << std::endl;


        MANAGER.IfNoFlagIsSetWriteHelpAndExit( ARGC, std::cout);

        if( !MANAGER.ReadAndCheckFlags( ARGC, ARGV))
        {
            MANAGER.WriteFlagsWithOptions( std::cout);
            std::cout << std::endl << std::endl;
            MANAGER.WriteHelp( std::cout);
            std::cout << ">>>>  Errors in reading flags  <<<<" << std::endl;
            exit( -1);
        }

        //////////////////////
        // CONSISTENCY CHECKS
        //////////////////////
        if( MANAGER.IsFlagSet( "dimensions") && MANAGER.IsFlagSet( "auto_dimensions"))
        {
        	std::cout << "ERROR: it does not make sense to have both '-dimensions' and '-auto_dimensions' set!" << std::endl;
        	exit( -1);
        }
        if( !MANAGER.IsFlagSet( "dimensions") && !MANAGER.IsFlagSet( "auto_dimensions"))
        {
        	std::cout << "ERROR: either '-dimensions' or '-auto_dimensions' has to be set!" << std::endl;
        	exit( -1);
        }

    } // end Setup...

} // end namespace command_line_factory




#endif /* COMMAND_LINE_FACTORY_H_ */
