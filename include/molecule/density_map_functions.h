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


#ifndef DENSITY_MAP_FUNCTIONS_H
#define DENSITY_MAP_FUNCTIONS_H

#include "../math/distributionXD.t.h"
#include "density_map.h"

namespace mol
{

    void DensityMap(  const int ARGC, const char *ARGV[])
    {
    	float
			x_min = mystr::ConvertStringToNumericalValue< float>( ARGV[ 3]),
			x_max = mystr::ConvertStringToNumericalValue< float>( ARGV[ 4]),
			y_min = mystr::ConvertStringToNumericalValue< float>( ARGV[ 5]),
			y_max = mystr::ConvertStringToNumericalValue< float>( ARGV[ 6]),
			z_min = mystr::ConvertStringToNumericalValue< float>( ARGV[ 7]),
			z_max = mystr::ConvertStringToNumericalValue< float>( ARGV[ 8]),
			voxelsize = mystr::ConvertStringToNumericalValue< float>( ARGV[ 9]);

    	std::string
			in = ARGV[ 2],
			type = ARGV[ 10],
			mode = ARGV[ 11],
			out = ARGV[ 12];

    	std::ifstream read( in.c_str());
    	std::ofstream write( out.c_str());

    	if( !read)
    	{
    		std::cerr << "input file not opened! " << std::endl;
    		return;
    	}
    	if( !write)
    	{
    		std::cerr << "output file not opened!" << std::endl;
    		return;
    	}

        boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
            all_molecules( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());

        all_molecules = mol::file::ReadMoleculesInGriffinFormat( read);

    	read.close();
    	read.clear();

    	math::DistributionXD< double, 3> density3D;

    	if( mode == "atom_count")
    	{
    		density3D = mol::AtomCountDensityMap( all_molecules, x_min, x_max, y_min, y_max, z_min, z_max, voxelsize, type);
    	}
    	else if( mode == "mass")
    	{
    		density3D = mol::MassDensityMap( all_molecules, x_min, x_max, y_min, y_max, z_min, z_max, voxelsize, type);
    	}
    	else if( mode == "gaussian")
    	{
    		density3D = mol::GaussianMassDensityMap( all_molecules, x_min, x_max, y_min, y_max, z_min, z_max, voxelsize, type);
    	}
    	else
    	{
    		std::cout << "no defined mode (atom_count, mass, gaussian)" << std::endl;
    		return;
    	}

    	write.close();
    	write.clear();

//    	write.open( "density3D.txt");
//    	write << density3D;
//    	write.close();
//    	write.clear();

    	write.open( out.c_str());

    	math::DistributionXD< double, 2> density2D( density3D, 2);
	std::cout << "write distribution: " << out << std::endl;
    	write << density2D;

    	write.close();
    	write.clear();

    	out = out.substr( 0, out.find("."));
    	out += ".gnu";
    	write.open( out.c_str());
	std::cout << "write gnuplot style: " << out << std::endl;


    	density2D.WriteGnuplot( write);

    	write.close();
    	write.clear();

    	std::cout << "2D_density done\n" << std::endl;
    	std::cout << "NOTE: you will find an additional file 'OUT.gnu' containing the density map in the format needed to plot in gnuplot (commands: 'set pm3d (map)' and 'splot \"OUT.gnu\" with pm3d)\n" << std::endl;
	std::cout << "NOTE: altering the color scale: set palette defined (0 \"black\", 4 \"blue\", 40 \"red\")" << std::endl;

    }
    

} // end namespace mol

#endif //  DENSITY_MAP_FUNCTIONS_H
