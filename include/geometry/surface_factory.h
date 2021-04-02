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


#ifndef GEOM_SURFACE_FACTORY_H
#define GEOM_SURFACE_FACTORY_H

#include "sphere.h"
#include "point_surface.h"
#include "directed_point_surface.h"
#include "point_surface_object.h"
#include "directed_point_surface_object.h"
#include "../storage/shared_pointer_vector.t.h"
#include "../storage/shared_pointer_map.t.h"
#include "../storage/map.t.h"
#include "../molecule/simple_molecule.t.h"
#include "../molecule/surf_atom.h"
#include "point_surface_cube.h"
#include "point_surface_cylinder.h"
#include "torus_slice.h"
#include "torus_functions.h"
#include "../command/command_line_manager.h"
#include "../readwrite/stream_functions.h"

#include <boost/shared_ptr.hpp>
//#include "disc.h"
//#include "sheet.h"
//#include "cylinder.h"
//#include "cube.h"


namespace geom
{

    namespace factory
    {
        bool PDBToSurface( std::istream &STREAM, SurfaceIF &SURF, const float &RESOLUTION);

        //   boost::shared_ptr< SurfaceIF> NonOverlappingSurface( const std::vector< boost::shared_ptr< PointSurfaceObject> > & SURFS);

        bool DoObjectsOverlap( const boost::shared_ptr< Object> &FIRST, const boost::shared_ptr< Object> &SECOND, const float &PROBE_RADIUS = 0.0);

        boost::shared_ptr< PointSurface> BuildNonOverlappingPointSurface( const store::ShPtrVec< PointSurfaceObject> &SURFS);

        std::vector< std::vector< size_t> > ListOverlappingObjects( const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS, const float &PROBE_RADIUS = 0.0);

        boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> >
        RemoveOverlappingSurfaces( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &MAP);

        store::ShPtrVec< DirectedPointSurfaceObject> &RemoveOverlappingSurface( store::ShPtrVec< DirectedPointSurfaceObject> &SURFS);

        store::ShPtrVec< PointSurfaceObject> &SmoothSurface( store::ShPtrVec< PointSurfaceObject> &SURFS, const float &PROBE_RADIUS);

        boost::shared_ptr< store::Map< std::string, geom::PointSurface> >
        SmoothSurfaces( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::PointSurfaceObject> > > &MAP, const float &PROBE_RADIUS);


        template< typename t_ATOM>
        mol::SimpleMolecule< t_ATOM> &RemoveOverlappingSurface( mol::SimpleMolecule< t_ATOM> &MOL)
        {
            store::ShPtrVec< t_ATOM>
                atoms( MOL.GetAtoms());
            store::ShPtrVec< DirectedPointSurfaceObject>
                surfs( atoms.size());
            std::vector< boost::shared_ptr< DirectedPointSurfaceObject> >::iterator
                pstr = surfs.begin();
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator satr = atoms.begin(); satr != atoms.end(); ++satr, ++pstr)
            {
                *pstr = ( *satr)->GetSurf();
            }
            RemoveOverlappingSurface( surfs);
            return MOL;
        }

        PointSurface PointSurfaceFromPointSurfaceObjects( const store::ShPtrVec< geom::PointSurfaceObject> &OBJECTS);

        DirectedPointSurface DirectedPointSurfaceFromPointSurfaceObjects( const store::ShPtrVec< geom::DirectedPointSurfaceObject> &OBJECTS);

        DirectedPointSurface PointSurfaceFromSimpleMolecule( const mol::SimpleMolecule< mol::SurfAtom> &MOL);


        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        SurfaceObjectsFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &RESOLUTION = 10.0,
                const float &FACTOR = 1.0,
                const float &OFFSET = 0.3
        );

        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        BuildSmoothSurfacesFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &DELTA,
                const float &FACTOR,
                const float &OFFSET,
                const float &PROBE_RADIUS
        );

        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        IntegrateSmoothSurfacesFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &DELTA,
                const float &FACTOR,
                const float &OFFSET,
                const float &PROBE_RADIUS
        );

        void AddToSurfaceMapIfNoOverlapWithOtherAtoms
        (
                const float &PROBE_RADIUS,
                const math::Vector3N &POSITION,
                const math::Vector3N &NORMAL,
                const float &SURF_SIZE,
                boost::shared_ptr< geom::DirectedPointSurfaceObject> &SURF,
                const std::multimap< size_t, boost::shared_ptr< TorusSlice> > &TORUS
        );

        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > >
        ReadIntoSurfObjectMap( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &SURF);

        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
        ReadIntoSurfObjectMap
        (
                std::istream &STREAM
        );

        bool IsPointInside( const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT, const math::Vector3N &POS);

        boost::shared_ptr< std::pair< boost::shared_ptr< DirectedPointSurfaceObject>, boost::shared_ptr< DirectedSurfacePoint> > >
        ClosestObjectAndPoint( const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS, const math::Vector3N &POSITION);

        boost::shared_ptr< std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> > >
        ObjectDistanceMap
        (
                const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS,
                const math::Vector3N &POSITION
        );

        boost::shared_ptr< geom::DirectedSurfacePoint>
        ClosestSurfacePoint
        (
                const boost::shared_ptr< std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> > > MAP,
                const math::Vector3N &POS
        );


        std::vector< std::vector< size_t> >
        FuseClusters
        (
            std::vector< std::vector< size_t> > &CLUSTERS,
            const std::vector< size_t> &TO_FUSE
        );


        bool
        IsElementOfCluster
        (
              const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT,
              const std::vector< size_t> &CLUSTER_IDS,
              const geom::DirectedPointSurface &SURF,
              const float &CLUSTER_THRESHOLD
        );


        std::vector< size_t>
        ElementOfExistingClusters
        (
              const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT,
              const std::vector< std::vector< size_t> >    &CLUSTER_IDS,
              const geom::DirectedPointSurface &SURF,
              const float &CLUSTER_THRESHOLD
        );


        void
        FilterSurfArtefacts
        (
                boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURF_MAP,
                const size_t &CLUSTER_MIN_SIZE,
                const float &CLUSTER_THRESHOLD
        );


        void
        BuildSurface
        (
                const CommandLineManager &COMMAND,
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE,
                boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURF_MAP,
                boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &SURF_OBJECT_MAP
        );

    } // end namespace factory

} // end namespace geom

#endif
