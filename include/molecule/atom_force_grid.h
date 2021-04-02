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


#ifndef ATOM_FORCE_GRID_H
#define ATOM_FORCE_GRID_H

#include "../macro/monitor.h"
#include "../storage/grid1D.t.h"
#include "../storage/grid3D.t.h"

#ifdef SQLITE
#include "../storage/grid_sqlite.t.h"
#endif

#include "../molecule/surf_atom.h"
#include "../molecule/simple_molecule.t.h"
#include "../phys/force_container.h"
#include "../storage/void_gridpoint.t.h"
#include "../storage/surf_gridpoint.t.h"
#include "../storage/constant_gridpoint.t.h"
#include "../storage/recursive_typemapped_gridpoint.t.h"
#include "../storage/typemapped_gridpoint.t.h"
//#include "../storage/transient_interaction_gridpoint.h"
//#include "../storage/map.t.h"
#include "../storage/position_grid.t.h"
#include "../utilities/parallelization.h"
#include "../geometry/surface_factory.h"
#include "../phys/force_factory.h"
#include "../storage/triplet.h"
#include "../geometry/surf_grid.h"

namespace mol
{

    class AtomForceGrid
    : public math::Function< mol::Atom, math::Vector3N >
    {
    protected:
        boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >    m_Grid;
    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        AtomForceGrid()
        : math::Function< mol::Atom, math::Vector3N >(),
        m_Grid()
        {}

        //! construct from data
        AtomForceGrid( const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID)
        : math::Function< mol::Atom, math::Vector3N >(),
        m_Grid( GRID)
        {}

        //! copy constructor
        AtomForceGrid( const AtomForceGrid &ORIGINAL)
        : math::Function< mol::Atom, math::Vector3N >( ORIGINAL),
        m_Grid( ORIGINAL.m_Grid)
        {
        	std::cout << __FUNCTION__ << " copy constructor!" << std::endl;
        }

        AtomForceGrid
        (
                const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID,
                const boost::shared_ptr< geom::SurfGrid> &SURF_GRID,
                const boost::shared_ptr< phys::ForceContainer> &FORCES,
                const float &MIN_FORCE_MAGNITUDE,
                const float &SFORCE_MAX_LENGTH = 1.0
        )
        : math::Function< mol::Atom, math::Vector3N >(),
          m_Grid(    GRID)
          { CalculateGrid( m_Grid, SURF_GRID, FORCES, MIN_FORCE_MAGNITUDE, SFORCE_MAX_LENGTH);}


        //! virtual destructor
        virtual ~AtomForceGrid(){}

        //! virtual copy constructor
        virtual AtomForceGrid *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////


        boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >&
        GetSetPositionGrid()
        {
        	return m_Grid;
        }

        const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >&
        GetPositionGrid() const
        {
        	return m_Grid;
        }


        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        GetGridPoint( const math::Vector3N &POSITION);

        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        GetGridPoint( const store::Vector3N< int> &INDICES);

        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE);

#ifdef FORCE_INDICES
        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        GetGridPoint
        (
	    const math::Vector3N &POSITION,
	    const std::string &MOL_TYPE, 
#ifdef CPPIO
	std::ostream &STREAM
#else
	FILE *STREAM
#endif
	);
#endif

        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        GetGridPoint( const store::Vector3N< int> &INDICES, const std::string &MOL_TYPE);

        virtual
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > &
        GetSetGridPoint( const store::Vector3N< int> &INDICES);

//        virtual
//        std::pair< std::string, int>
//        GetGridPointTypeAndNSP( const math::Vector3N &POSITION, const std::string &MOL_TYPE);

        virtual
        util::FunctorEnum
        GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE);

        virtual math::Vector3N & operator()( const Atom &ATOM) const;

//        virtual boost::shared_ptr< ForceAtom> operator[]( const Atom &ATOM) const
//        {
//            return boost::shared_ptr< ForceAtom>( new ForceAtom( ATOM, math::Vector3N( new math::Vector3N( operator()( ATOM)))));
//        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        //! LOGIC MAP for CalculateGrid
        //!
        //! if: single mol type
        //!	  iterate grid
        //!     get SURF_GRID value for this point
        //!     if value < 1  (inside)
        //!       set GP to ConstantGP
        //!     elsif dist > max: set GP to VoidGP
        //!     else GP to InteractionGP
        //! elsif multiple mol types
        //!   (pass mol types to recursive- and type mapped GP)
        //!   iterate grid
        //!     get SURF_GRID value for this point
        //!     main IF: value == 0 (inside for all)
        //!       check whether all mol types point to same NSP
		//!         if yes: insert ConstantGP
        //!         if no: insert TypeMappedGP
        //!     IF value > 0 (outside for all)
        //!       is dist within max cutoff:
        //!         if yes: insert InteractionGP
        //!         if no:  insert VoidGP
        //!     IF value < 0  (inside for some mol types)
        //!       build RecursiveGP:
        //!         insert InteractionGP for "all"
        //!         (translate value => vec of indices of mol_types being inside)
        //!         iterate mol types
        //!         if mol type represented by value (indices):
        //!            yes: insert ConstantGP for mol type
        //!            no:  insert InteractionGP for mol type (copy GP from "all" ==> should one recalculate NSP?? (not really needed))
        //!
        virtual const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
        &CalculateGrid
        (
                boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID,
                const boost::shared_ptr< geom::SurfGrid> &SURF_GRID,
                const boost::shared_ptr< phys::ForceContainer> &FORCES,
                const float &MIN_FORCE_MAGNITUDE,
                const float &SFORCE_MAX_LENGTH
        );


       //! LOGIC MAP for FillSurfaceLayer:
       //! if: single mol type (mol indep)
       //!   iterate: surf points (SP)
       //!     if: SP within grid limits
       //!       if: grid point (GP) type == ConstantGP
       //!          set GP to Interaction
       //! elsif: mol type dependent
       //!   iterate: mol types
       //!     iterate: surf points
       //!       if SP within grid limits
       //!         if GP type == ConstantGP or TypeMapped
       //!           set GP to Interaction
       //!         elsif GP type == RecursiveTypeMappedGP
       //!		     extract GP for this mol type and if GP type == ConstantGP set to Interaction
       virtual void
        FillSurfaceLayer
        (
                const boost::shared_ptr< geom::SurfGrid> &SURF_GRID,
                const boost::shared_ptr< phys::ForceContainer> &FORCES,
                const float &MIN_FORCE_MAGNITUDE
        );


        //! same LOGIC MAP as in FillSurfLayer
        virtual void
        EmptySurfaceLayer
        (
                const boost::shared_ptr< geom::SurfGrid> &SURF_GRID
        );




//        virtual boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
//        CalculateGridPoint
//        (
//                const math::Vector3N &POSITION,
//                const store::ShPtrVec< geom::DirectedPointSurfaceObject> &SURF_VEC,
//                const std::vector< std::vector< size_t> > &NEIGHBOR_LIST,
//                const boost::shared_ptr< phys::ForceContainer> &FORCES,
//                const float &CUTOFF,
//                const float &SFORCE_MAX_LENGTH = 1.0,
//                const float PROBE_RADIUS = 1.4
//        );


        virtual math::Vector3N
        SmoothedSurfaceRepulsion
        (
                const boost::shared_ptr< geom::DirectedPointSurfaceObject> &FIRST,
                const boost::shared_ptr< geom::DirectedPointSurfaceObject> &SECOND,
                const math::Vector3N &POS,
                const float &PROBE_RADIUS,
                const float &LENGTH
        );

//        virtual void SetGridLimits( mol::SimpleMolecule< mol::Atom> &MOL, const math::Vector3N &DELTA, const float &CUTOFF);


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::string GetClassName() const;

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::ostream &WriteAsPdb( std::ostream &STREAM) const;

//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const;




    }; // end class AtomForceGrid
} // end namespace mol




#endif /* ATOM_FORCE_GRID_H */
