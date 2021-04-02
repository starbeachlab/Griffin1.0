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


#ifndef DENSITY_MAP_H_
#define DENSITY_MAP_H_

#include "../math/distributionXD.t.h"
#include "../storage/shared_pointer_vector.t.h"
#include "../molecule/simple_molecule.t.h"
#include "../molecule/atom.h"

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
	);



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
	);



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
	);





//	math::TupelXD< float, 3> InterpolateTupleXD( const math::TupelXD< float, 3> &DENSITY)
//	{
//		math::TupelXD< float, 3> interpolated;
//		return interpolated;
//	}
//
//	math::TupelXD< float, 3> GaussianSmoothedTupleXD( const math::TupelXD< float, 3> &DENSITY);
//
//	template< int N>
//	math::TupleXD< float, N-1> CondensingTupleXD( const math::TupelXD< float, N> &DENSITY, const int &XYZ_LAYER);

} // end namespace mol


#endif /* DENSITY_MAP_H_ */
