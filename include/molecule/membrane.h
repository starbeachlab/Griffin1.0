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


#ifndef MEMBRANE_H_
#define MEMBRANE_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../storage/limits3D.h"
#include "../storage/map.t.h"
#include "../storage/shared_pointer_vector.t.h"
#include "simple_molecule.t.h"
#include "simple_molecule_file_handler.t.h"
#include "atom.h"
#include "atom_file_handler.t.h"
#include "../storage/multiplets.h"
#include "../command/command_line_manager.h"
#include "../math/histogram.t.h"
#include "../math/distribution.t.h"

namespace mol
{

    class Membrane
    : public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        std::vector< std::string>
            m_LipidKeys;
        std::vector< std::string>
            m_SolKeys;
        store::Limits3D
            m_TotalLimits;
        float
            m_LipidMin;
        float
            m_Center;
        float
            m_LipidMax;
        std::vector< store::Map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > > >
            m_MoleculeLayer;
        int
            m_PreviousID;
        std::string
            m_PreviousType;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Membrane( const std::vector< std::string> &LIPID_KEYS, const std::vector< std::string> &SOL_KEYS)
        :
            m_LipidKeys( LIPID_KEYS),
            m_SolKeys( SOL_KEYS),
            m_TotalLimits(),
            m_LipidMin(),
            m_Center(),
            m_LipidMax(),
            m_MoleculeLayer( 4),
            m_PreviousID( std::numeric_limits< int>::max()),
            m_PreviousType( "")
            {}

        //! copy constructor
        Membrane( const Membrane &ORIGINAL)
        :
            m_LipidKeys( ORIGINAL.m_LipidKeys),
            m_SolKeys( ORIGINAL.m_SolKeys),
            m_TotalLimits( ORIGINAL.m_TotalLimits),
            m_LipidMin( ORIGINAL.m_LipidMin),
            m_Center( ORIGINAL.m_Center),
            m_LipidMax( ORIGINAL.m_LipidMax),
            m_MoleculeLayer( ORIGINAL.m_MoleculeLayer),
            m_PreviousID( ORIGINAL.m_PreviousID),
            m_PreviousType( ORIGINAL.m_PreviousType)
            {}

        //! virtual destructor
        virtual ~Membrane(){}

        //! virtual copy constructor
        virtual Membrane *Clone() const{ return new Membrane( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const store::Limits3D &GetTotalLimits() const;

        virtual std::pair< float, float> GetXLimits() const;

        virtual std::pair< float, float> GetYLimits() const;

        virtual std::vector< float> GetZLimits() const;

        virtual void SetXLimits( const float &MIN, const float &MAX);

        virtual void SetYLimits( const float &MIN, const float &MAX);

        virtual void SetZLimits( const float &LOWER_SOLUTION, const float &LOWER_LIPID, const float &UPPER_LIPID, const float &UPPER_SOLUTION);

        virtual void SetZLimits( const float &LOWER_SOLUTION, const float &LOWER_LIPID, const float &CENTER, const float &UPPER_LIPID, const float &UPPER_SOLUTION);

        virtual std::vector< std::string> GetMoleculeTypesInLayer( const size_t &LAYER) const;

        virtual std::vector< std::string> GetMoleculeTypes() const;

        virtual
        const store::Map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > > &
        GetMoleculeLayer( const size_t &LAYER) const;

        virtual
        store::ShPtrVec< SimpleMolecule< Atom> >
        GetMolecules() const;

        virtual
        const std::vector< std::string> &
        GetLipidKeys() const;

        virtual
        const std::vector< std::string> &
        GetSolKeys() const;

        virtual
        store::Map< std::string, float>
        Densities( const size_t &LAYER) const;

        virtual
        store::Limits3D
        FourLayerBoxLimits( const size_t &LAYER) const;

        virtual
        bool
        IsValidKey( const std::string &KEY) const;

        virtual
        size_t
        LayerID( const float &Z_COORDINATE) const;

        virtual
        void
        AdjustLimits( const std::string &MOLECULE_TYPE, const math::Vector3N &POSITION/*, const float &RADIUS*/);

        virtual
        void
        AdjustLimits( const store::Limits3D &LIMITS);

        virtual
        std::vector< float>
        ZLayerLimits() const;

        virtual
        bool
        LineChecker
        (
                const boost::shared_ptr< Atom> &ATOM,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOLECULE
        );
        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM){ return STREAM;}

        virtual float BilayerThickness() const;

        virtual void AnalyseAndBuild(const CommandLineManager &COMMAND);

        virtual std::istream &Read
        (
                std::istream &STREAM,
                const boost::shared_ptr< store::Map< std::string, float> > &MASS_MAP,
                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR
        );

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::ostream &WriteLimits( std::ostream &STREAM) const;

        virtual std::ostream &Write( std::ostream &STREAM, const CommandLineManager &COMMAND) const;

        virtual std::ostream &WriteLayerOverview( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class Membrane
} // end namespace mol




#endif /* MEMBRANE_H_ */
