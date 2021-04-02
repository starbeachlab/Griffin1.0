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


#ifndef MOLECULE_FORCE_GRID_H_
#define MOLECULE_FORCE_GRID_H_

#include <cstdio>


#include "../storage/position_grid.t.h"
#include "../readwrite/stream_operator.h"
#include "../readwrite/stream_functions.h"
#include "../string/io_string_functions.h"
#include "atom_force_grid.h"
#include "simple_molecule.t.h"
#include "simple_molecule_file_handler.t.h"
#include "molecule_factory.h"
#include "../math/vector_angle.h"
#include "../math/matrix3x3N.h"
#include "../math/matrix_vector_functions.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"
#include "../macro/monitor.h"
#include "../utilities/time.h"
#include "periodic_molecule_iterator.h"

#ifdef MPI_PARALLEL
	#include "../molecule/mpi_molecule_iterator.h"
#endif


namespace mol
{

    class MoleculeForceGrid
    : public math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >
    {
        ///////////////
        //    DATA     //
        ///////////////

    protected:
        boost::shared_ptr< AtomForceGrid>                                           m_Grid;
//        boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >    m_VdwEpsilonAndRadiusMap;
//        boost::shared_ptr< store::Map< std::string, float> >                        m_PartialChargeMap;
//        boost::shared_ptr< store::Map< std::string, float> >                        m_MassMap;
        std::string                                                                 m_ForceFieldType;
        geom::SurfGrid                                                              m_SurfGrid;
        math::VectorAngle                                                           m_AngleChecker;
        float                                                                       m_SforceLength;
        std::vector< std::string>													m_LipidNames;
        bool																		m_IgnoreHydrogens;
#ifdef TIMERS
        SumTimer																	m_ForceTimer;
        SumTimer																	m_GetGPTimer;
        SumTimer																	m_SimpleGetGPTimer;
        SumTimer																	m_RedirectTimer;
        SumTimer                                                                    m_CalcForceVirialEnergyTimer;
        SumTimer                                                                    m_CheckTimer;
        SumTimer																	m_SimpleTimer;
#endif


    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        MoleculeForceGrid()
        : math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >(),
        m_Grid(),
//        m_VdwEpsilonAndRadiusMap(),
//        m_PartialChargeMap(),
//        m_MassMap(),
        m_ForceFieldType(),
        m_SurfGrid(),
        m_AngleChecker(),
        m_SforceLength(),
        m_LipidNames(),
	m_IgnoreHydrogens( false)
#ifdef TIMERS
	    ,m_ForceTimer("MoleculeForceGrid:force"),
	    m_GetGPTimer("MoleculeForceGrid:Get GridPoint and type and geometric center"),
	    m_SimpleGetGPTimer("MoleculeForceGrid:Simple gridpoint and type and geometric center"),
	    m_RedirectTimer("MoleculeForceGrid:redirect"),
	    m_CalcForceVirialEnergyTimer( "MoleculeForceGrid:CalcForceVirialEnergyTimer per mol sum"),
	    m_CheckTimer( "MoleculeForceGrid:check and write force energy virial"),
	    m_SimpleTimer( "MoleculeForceGrid:simple")
#endif
        {}


//	//! construct passing positiongrid
//        MoleculeForceGrid( const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID)
//        : math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >(),
//        m_Grid( new AtomForceGrid( GRID)),
////        m_VdwEpsilonAndRadiusMap(),
////        m_PartialChargeMap(),
////        m_MassMap(),
//        m_ForceFieldType(),
//        m_SurfGrid(),
//        m_AngleChecker(),
//        m_SforceLength(),
//        m_LipidNames(),
//	m_IgnoreHydrogens( false)
//#ifdef TIMERS
//	    ,m_ForceTimer("force"),
//	    m_GetGPTimer("getGP"),
//	    m_RedirectTimer("redirect")
//#endif
//        {}




        //! construct from data
        MoleculeForceGrid
        (
                const boost::shared_ptr< AtomForceGrid> &GRID,
//                const boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >  &VDW_EPSILON_AND_RADIUS_MAP,
//                const boost::shared_ptr< store::Map< std::string, float> > &PARTIAL_CHARGE_MAP,
//                const boost::shared_ptr< store::Map< std::string, float> > &MASS_MAP,
//                const std::string &FORCE_FIELD,
                const geom::SurfGrid &SURF_GRID,
                const float &ANGLE_THRESHOLD,
                const float &SFORCE_LENGTH,
                const std::vector< std::string> &LIPID_NAMES = std::vector< std::string>()
        )
        : math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >(),
        m_Grid( GRID),
//        m_VdwEpsilonAndRadiusMap( VDW_EPSILON_AND_RADIUS_MAP),
//        m_PartialChargeMap( PARTIAL_CHARGE_MAP),
//        m_MassMap( MASS_MAP),
        m_ForceFieldType( ""),
        m_SurfGrid( SURF_GRID),
        m_AngleChecker( ANGLE_THRESHOLD),
        m_SforceLength( SFORCE_LENGTH),
        m_LipidNames( LIPID_NAMES),
        m_IgnoreHydrogens( false)
#ifdef TIMERS
        , m_ForceTimer( "MoleculeForceGrid: ForceCalculation"),
        m_GetGPTimer( "MoleculeForceGrid: Get GridPoint and type and geometric center"),
        m_SimpleGetGPTimer("MoleculeForceGrid:Simple gridpoint and type and geometric center"),
        m_RedirectTimer( "MoleculeForceGrid: Redirecting"),
        m_CalcForceVirialEnergyTimer( "MoleculeForceGrid: CalcForceVirialEnergyTimer per mol sum"),
	    m_CheckTimer( "MoleculeForceGrid: check and write force energy virial"),
	    m_SimpleTimer( "MoleculeForceGrid:simple")
#endif
        {
            m_SurfGrid.ClearData();
        }

        //! copy constructor
        MoleculeForceGrid( const MoleculeForceGrid &ORIGINAL)
        : math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >( ORIGINAL),
        m_Grid( ORIGINAL.m_Grid),
//        m_VdwEpsilonAndRadiusMap( ORIGINAL.m_VdwEpsilonAndRadiusMap),
//        m_PartialChargeMap( ORIGINAL.m_PartialChargeMap),
//        m_MassMap( ORIGINAL.m_MassMap),
        m_ForceFieldType( ORIGINAL.m_ForceFieldType),
        m_SurfGrid( ORIGINAL.m_SurfGrid),
        m_AngleChecker( ORIGINAL.m_AngleChecker),
        m_SforceLength( ORIGINAL.m_SforceLength),
	    m_LipidNames( ORIGINAL.m_LipidNames),
	    m_IgnoreHydrogens( ORIGINAL.m_IgnoreHydrogens)
#ifdef TIMERS
	    , m_ForceTimer( ORIGINAL.m_ForceTimer),
	    m_GetGPTimer( ORIGINAL.m_GetGPTimer),
	    m_SimpleGetGPTimer( ORIGINAL.m_SimpleGetGPTimer),
	    m_RedirectTimer( ORIGINAL.m_RedirectTimer)
		, m_CalcForceVirialEnergyTimer( ORIGINAL.m_CalcForceVirialEnergyTimer),
		m_CheckTimer( ORIGINAL.m_CheckTimer),
		m_SimpleTimer( ORIGINAL.m_SimpleTimer)
#endif
        {
        	StandardWrite( __PRETTY_FUNCTION__ << " copy constructor");
        }

        //! virtual destructor
        virtual ~MoleculeForceGrid(){}

        //! virtual copy constructor
        virtual MoleculeForceGrid *Clone() const;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! returns a vector of forces for molecule
        virtual boost::shared_ptr< std::vector< math::Vector3N> >  & operator()( const boost::shared_ptr< SimpleMolecule< Atom> > &MOL) const;

        // THIS SHOULD NOT BE HERE, THIS SHOULD NOT BE ABOUT
        //! returns a vector of forces for shared pointer vector of molecules
        virtual boost::shared_ptr< std::vector< math::Vector3N> > operator()( const boost::shared_ptr< store::ShPtrVec< SimpleMolecule< Atom> > > &MOLS) const;


        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////


//        virtual int Execute
//        (
//                const std::string &INPUT_FILE,
//                const std::string &OUTPUT_FILE,
//                boost::shared_ptr< SimpleMolecule< Atom> > &MOL,
//                std::ostream &WRITE
//        );

        virtual int Execute
        (
                const std::string &INPUT_FILE,
                const std::string &OUTPUT_FILE,
//                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
                boost::shared_ptr< mol::MoleculeIterator> &MOLS,
                std::ostream &LOG
        );

//        virtual int Execute
//        (
//                const std::string &INPUT_FILE,
//                const std::string &INPUT_FILE_FORMAT,
//                const std::string &OUTPUT_FILE,
//                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
//                std::ostream &WRITE
//        );


//        std::ostream &CalculateForcesAndWriteToStream( std::ostream &STREAM, const boost::shared_ptr< SimpleMolecule< Atom> > &MOL);
//
//        std::ostream &CalculateAndWriteForces( std::ostream &STREAM, const boost::shared_ptr< SimpleMolecule< Atom> > &MOL, int &COUNT, std::ostream &LOG = std::cout);
//
//        std::ostream &SimpleCalcAndWriteForces( std::ostream &STREAM, const boost::shared_ptr< SimpleMolecule< Atom> > &MOL, int &COUNT, std::ostream &LOG = std::cout);


        void CalculateAndWriteForcesEnergyVirial
        (
#ifdef CPPIO
    		std::ostream &WRITE,
#else
    		FILE *WRITE,
#endif
              boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
      		  std::ostream &LOG
        );


        void CalculateAndWriteForcesEnergyVirial
        (
#ifdef CPPIO
    		std::ostream &WRITE,
#else
    		FILE *WRITE,
#endif
              boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
      		  math::Matrix3x3N &VIRIAL,
      		  float &ENERGY,
      		  float &MAX_DIST,
      		  float &MAX_REDIRECTED_DIST,
      		  int &SFORCE_COUNT,
      		  int &REDIRECTED_COUNT,
      		  std::ostream &LOG
#ifdef ANALYZE_BURIAL 
      		  , math::Histogram< float> & REDIRECTED_DISTANCES,
      		  math::Histogram< float> & SURF_DISTANCES
#endif
        );

        void SimpleCalcAndWriteForcesEnergyVirial
        (
#ifdef CPPIO
    		std::ostream &WRITE,
#else
    		FILE *WRITE,
#endif
              boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
      		  math::Matrix3x3N &VIRIAL,
      		  float &ENERGY,
      		  int &SFORCE_COUNT,
      		  float &MAX_DIST,
      		  std::ostream &LOG = std::cout
        );


        virtual const std::string &GetForceFieldType() const
        {
            return m_ForceFieldType;
        }

        virtual void SetForceFieldType( const std::string &TYPE)
        {
        	m_ForceFieldType = TYPE;
        }

//        virtual const boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
//        &GetVdwEpsilonAndRadiusMap() const
//        {
//            return m_VdwEpsilonAndRadiusMap;
//        }
//        virtual const boost::shared_ptr< store::Map< std::string, float> >
//        &GetPartialChargeMap() const
//        {
//            return m_PartialChargeMap;
//        }
//        virtual const boost::shared_ptr< store::Map< std::string, float> >
//        &GetMassMap() const
//        {
//            return m_MassMap;
//        }

        boost::shared_ptr< AtomForceGrid> &
        GetSetAtomForceGrid()
        {
        	return m_Grid;
        }

        const boost::shared_ptr< AtomForceGrid> &
        GetAtomForceGrid() const
        {
        	return m_Grid;
        }

//        const std::string &GetForceField() const
//        {
//        	return m_ForceFieldType;
//        }


        const geom::SurfGrid &GetSurfGrid() const
        {
        	return m_SurfGrid;
        }

        geom::SurfGrid &GetSetSurfGrid()
        {
        	return m_SurfGrid;
        }

        const float &GetSForceScale() const
	{
	    return m_SforceLength;
	}


        void SetSForceScale( const float &SCALE)
	{
	    m_SforceLength = SCALE;
	}


         const math::VectorAngle &GetVectorAngle() const
 		{
         	return m_AngleChecker;
 		}

        math::VectorAngle &GetSetVectorAngle()
	    {
		return m_AngleChecker;
	    }

        void SetVectorAngle( const math::VectorAngle &VA)
 		{
         	m_AngleChecker = VA;
 		}

        const std::vector< std::string> &GetLipidNames() const
		{
        	return m_LipidNames;
		}

        void SetLipidNames( const std::vector< std::string> &LIPID_NAMES)
        {
        	m_LipidNames = LIPID_NAMES;
        }

        void SetIgnoreHydrogens( const bool &BOOL = true)
        {
        	m_IgnoreHydrogens = BOOL;
        }

#ifdef TIMERS
        void ResetTimers()
        {
        	m_ForceTimer.Reset();
        	m_GetGPTimer.Reset();
		m_SimpleGetGPTimer.Reset();
        	m_RedirectTimer.Reset();
        	m_CalcForceVirialEnergyTimer.Reset();
        	m_CheckTimer.Reset();
        	m_SimpleTimer.Reset();
        }
#endif


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        std::ostream &WriteAsPdb( std::ostream &STREAM) const;

#ifdef TIMERS
        std::ostream &WriteTimerStatus( std::ostream &STREAM) const
        {
        	m_ForceTimer.WriteStatus( STREAM);
        	m_GetGPTimer.WriteStatus( STREAM);
		m_SimpleGetGPTimer.WriteStatus( STREAM);
        	m_RedirectTimer.WriteStatus( STREAM);
        	m_CalcForceVirialEnergyTimer.WriteStatus( STREAM);
        	m_CheckTimer.WriteStatus( STREAM);
        	m_SimpleTimer.WriteStatus( STREAM);
		return STREAM;
        }
#endif


//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;



    }; // end class MoleculeForceGrid


    inline
    std::istream &operator >> ( std::istream &STREAM, MoleculeForceGrid& GRID)
    {
        return GRID.Read( STREAM);
    }


    inline
    std::ostream &operator << ( std::ostream &STREAM, const MoleculeForceGrid &GRID)
    {
        return GRID.Write( STREAM);
    }

} // end namespace mol




#endif /* MOLECULE_FORCE_GRID_H_ */
