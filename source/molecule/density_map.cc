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


#include "../../include/molecule/density_map.h"


namespace mol
{
	math::DistributionXD< double, 3>
	AtomCountDensityMap
	(
			const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &ALL_MOLECULES,
			const float &X_MIN,
			const float &X_MAX,
			const float &Y_MIN,
			const float &Y_MAX,
			const float &Z_MIN,
			const float &Z_MAX,
			const float &VOXELSIZE,
			const std::string &TYPE
	)
	{
		std::cout << __FUNCTION__ << " for: " << TYPE << std::endl;
		std::vector< float>
			min( 3, 0.0),
			max( 3, 0.0);

		min[ 0] = X_MIN;
		min[ 1] = Y_MIN;
		min[ 2] = Z_MIN;

		max[ 0] = X_MAX;
		max[ 1] = Y_MAX;
		max[ 2] = Z_MAX;

		math::DistributionXD< double, 3> density( VOXELSIZE, min, max);

//		std::cout << "mols: " << ALL_MOLECULES->size() << std::endl;

		for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = ALL_MOLECULES->begin(); mol_itr != ALL_MOLECULES->end(); ++mol_itr)
		{
//			std::cout << "atoms: " << ( *mol_itr)->GetAtoms().size() << std::endl;
			for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = ( *mol_itr)->GetAtoms().begin(); atom_itr != ( *mol_itr)->GetAtoms().end(); ++atom_itr)
			{
				if( TYPE == "all" || ( *atom_itr)->GetResidueType().find( TYPE) != std::string::npos)
				{
					density.Insert( ( *atom_itr)->GetPosition());
				}
			}
		}
		return density;
	}

	math::DistributionXD< double, 3>
	MassDensityMap
	(
			const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &ALL_MOLECULES,
			const float &X_MIN,
			const float &X_MAX,
			const float &Y_MIN,
			const float &Y_MAX,
			const float &Z_MIN,
			const float &Z_MAX,
			const float &VOXELSIZE,
			const std::string &TYPE
	)
	{
		std::cout << __FUNCTION__ << std::endl;
		std::vector< float>
			min( 3, 0.0),
			max( 3, 0.0);

		min[ 0] = X_MIN;
		min[ 1] = Y_MIN;
		min[ 2] = Z_MIN;

		max[ 0] = X_MAX;
		max[ 1] = Y_MAX;
		max[ 2] = Z_MAX;

		math::DistributionXD< double, 3> density( VOXELSIZE, min, max);

		for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = ALL_MOLECULES->begin(); mol_itr != ALL_MOLECULES->end(); ++mol_itr)
			for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = ( *mol_itr)->GetAtoms().begin(); atom_itr != ( *mol_itr)->GetAtoms().end(); ++atom_itr)
			{
				if( TYPE == "all" || ( *atom_itr)->GetResidueType().find( TYPE) != std::string::npos)
				{
					density.Insert( ( *atom_itr)->GetPosition(), ( *atom_itr)->GetMass());
				}
			}
		return density;
	}

	math::DistributionXD< double, 3>
	GaussianMassDensityMap
	(
			const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &ALL_MOLECULES,
			const float &X_MIN,
			const float &X_MAX,
			const float &Y_MIN,
			const float &Y_MAX,
			const float &Z_MIN,
			const float &Z_MAX,
			const float &VOXELSIZE,
			const std::string &TYPE
	)
	{
		std::cout << __FUNCTION__ << std::endl;
		std::vector< float>
			min( 3, 0.0),
			max( 3, 0.0);

		min[ 0] = X_MIN;
		min[ 1] = Y_MIN;
		min[ 2] = Z_MIN;

		max[ 0] = X_MAX;
		max[ 1] = Y_MAX;
		max[ 2] = Z_MAX;

		math::DistributionXD< double, 3> density( VOXELSIZE, min, max);

		float d, dist;

		math::Vector3N
			min_pos,
			max_pos,
			pos;

		for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = ALL_MOLECULES->begin(); mol_itr != ALL_MOLECULES->end(); ++mol_itr)
			for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = ( *mol_itr)->GetAtoms().begin(); atom_itr != ( *mol_itr)->GetAtoms().end(); ++atom_itr)
			{
				if( TYPE == "all" || ( *atom_itr)->GetResidueType().find( TYPE) != std::string::npos)
				{
					density.Insert( ( *atom_itr)->GetPosition(), ( *atom_itr)->GetMass());
					d = 2.0 * ( *atom_itr)->GetVanDerWaalsRadius();
					d -= math::Modulo( d, VOXELSIZE);
					min_pos = ( *atom_itr)->GetPosition() - d * math::Vector3N( 1.0, 1.0, 1.0);
					max_pos = ( *atom_itr)->GetPosition() + d * math::Vector3N( 1.0, 1.0, 1.0);
					for( float x = min_pos[0]; x <= max_pos[0]; x += VOXELSIZE)
						for( float y = min_pos[1]; y <= max_pos[1]; y += VOXELSIZE)
							for( float z = min_pos[2]; z <= max_pos[2]; z += VOXELSIZE)
							{
								pos = math::Vector3N( x, y, z);
								dist = math::Distance( pos, ( *atom_itr)->GetPosition());
								density.Insert( pos, math::GaussDistribution( dist, ( *atom_itr)->GetVanDerWaalsRadius()));
							}
				}
			}
		return density;
	}






}
