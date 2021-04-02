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


#ifndef VOLUME_FUNCTIONS_H_
#define VOLUME_FUNCTIONS_H_

#include <iterator>

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../storage/vector3N.t.h"
#include "../storage/limits3D.h"
#include "../molecule/membrane.h"
#include "sphere.h"
#include "cube.h"
#include "directed_point_surface_object.h"
#include "../external/std_functions.h"
#include "../molecule/simple_molecule_file_handler.t.h"
#include "../geometry/surf_grid.h"

namespace geom
{

    enum GeometricRelation { NO_OVERLAP, PARTLY_OVERLAPPED, A_BURIED_IN_B, B_BURIED_IN_A, COMPLETE_OVERLAP, UNKNOWN_RELATION};


    float VolumeOfCapsOfOverlappingSpheres( const Sphere &A, const Sphere &B);


    bool IsSphereInsideCube( const Cube &CUBE, const Sphere &SPHERE);

    GeometricRelation IdentifyGeometricRelation( const Cube &CUBE, const Sphere &SPHERE);


    store::Limits3D CalcLimits( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE);

    store::Limits3D CalcLimits( const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLECULES);



    float
    Volume
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const store::Limits3D &LIMITS,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS,
            float &VOLUME
    );


    float
    Volume
    (
            const boost::shared_ptr< mol::Membrane> MEMBRANE,
            const boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURFACES,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS,
            float &VOLUME
    );


    float
    Volume
    (
            const boost::shared_ptr< geom::DirectedPointSurface> &SURFACE,
            const store::Limits3D &LIMITS,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS
    );


    float
    Volume
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const float &VOXELSIZE = 1.0,
            const size_t &NR_LEVELS = 5
    );


    boost::shared_ptr< geom::DirectedPointSurface>
    SurfacePointsWithinCube
    (
            const Cube &CUBE,
            const boost::shared_ptr< std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> > > &DISTANCE_MAP
    );


    bool
    IsMoleculeWithinLimits
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const store::Limits3D &LIMITS
    );


    size_t
    NumberMoleculesInLimits
    (
            const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
            const store::Limits3D &LIMITS
    );


    float
    SimpleDensity
    (
            const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
            const store::Limits3D &LIMITS
    );



    float
    SumOverSquaredDepthOfPenetretion
    (
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > MOLECULE
    );


    float
    SumOverSquaredDepthOfPenetretion
    (
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            boost::shared_ptr< geom::DirectedPointSurface> SURF
    );


    float
    SquaredDepthOfHeavyAtomPenetration
    (
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            geom::DirectedPointSurface SURFACE
    );

    boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
    RemoveMoleculesInImplicitMembraneVolume
    (
            std::ostream &STREAM,
            boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > SURF_MAP,
            const float &VOXELSIZE = 1.0,
            const size_t &NR_LEVELS = 4
    );

    std::vector< store::Map< std::string, float> >
    Volume
    (
            const boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID
    );


    boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
    RemoveMoleculesInImplicitMembraneVolume
    (
            std::ostream &STREAM,
            boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID,
            const std::string &NR_LIPIDS_LOWER_LAYER,
            const std::string &NR_LIPIDS_UPPER_LAYER,
            const std::vector< std::pair< std::string, int> > &KEEP_MOLS,
            const std::vector< std::pair< std::string, int> > &REMOVE_MOLS
    );

} // end namespace geom




#endif /* VOLUME_FUNCTIONS_H_ */
