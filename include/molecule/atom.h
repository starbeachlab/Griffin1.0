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


#ifndef SIMPLE_ATOM_H_
#define SIMPLE_ATOM_H_


#include "../math/vector3N.h"
#include "../readwrite/stream_operator.h"
#include "moleculeIF.h"
#include "../string/io_string_functions.h"


namespace mol
{

    class Atom // TODO: clean mess with atom name and type!!!
    : public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        math::Vector3N                  m_Position;
        std::string                     m_AtomType;      // e.g.: CT1,CT3
        std::string                     m_AtomName;     // e.g.: CA,CB
        std::string                     m_ResidueType;
        MoleculeIF< Atom>*              m_Mol;
        float                           m_PartialCharge;
        float                           m_Epsilon;
        float                           m_VdwRadius;
        float                           m_Mass;
        int                         	m_ResidueID;
        int                          	m_AtomID;
        float							m_AttractiveVdwFactor;
        float							m_RepulsiveVdwFactor;
	int                             m_MolID;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        Atom()
        : m_Position(),
        m_AtomType( "UNDEFINED"),
        m_AtomName( "UNDEFINED"),
        m_ResidueType( "UNDEFINED"),
        m_Mol(),
        m_PartialCharge(),
        m_Epsilon(),
        m_VdwRadius(),
        m_Mass(),
        m_ResidueID(),
        m_AtomID(),
        m_AttractiveVdwFactor( std::numeric_limits< float>::max()),
	    m_RepulsiveVdwFactor( std::numeric_limits< float>::max()),
	    m_MolID()
        {}

        //! construct from data
        Atom( const std::string &TYPE, const math::Vector3N &POS = math::Vector3N())
        : m_Position( POS),
        m_AtomType( TYPE),
        m_AtomName( "UNDEFINED"),
        m_ResidueType( "UNDEFINED"),
        m_Mol(),
        m_PartialCharge(),
        m_Epsilon(),
        m_VdwRadius(),
        m_Mass(),
        m_ResidueID(),
        m_AtomID(),
        m_AttractiveVdwFactor( std::numeric_limits< float>::max()),
	    m_RepulsiveVdwFactor( std::numeric_limits< float>::max()),
	    m_MolID()
        {} 


        //! construct from data
        Atom( const std::string &ATOM_TYPE, const std::string & RESIDUE_TYPE, const math::Vector3N &POS = math::Vector3N())
        : m_Position( POS),
        m_AtomType( ATOM_TYPE),
        m_AtomName( "UNDEFINED"),
        m_ResidueType( RESIDUE_TYPE),
        m_Mol(),
        m_PartialCharge(),
        m_Epsilon(),
        m_VdwRadius(),
        m_Mass(),
        m_ResidueID(),
        m_AtomID(),
        m_AttractiveVdwFactor( std::numeric_limits< float>::max()),
	    m_RepulsiveVdwFactor( std::numeric_limits< float>::max()),
	    m_MolID()
        {} 


        //! construct from data
        Atom
        (
                const std::string &ATOM_TYPE,
                const std::string &ATOM_NAME,
                const std::string & RESIDUE_TYPE,
                const math::Vector3N &POS = math::Vector3N(),
                const float &PARTIAL_CHARGE = 0.0,
                const int &RESIDUE_ID = 0,
                const int &ATOM_ID = 0,
		const int &MOL_ID = 0
        )
        : m_Position( POS),
        m_AtomType( ATOM_TYPE),
        m_AtomName( ATOM_NAME),
        m_ResidueType( RESIDUE_TYPE),
        m_Mol(),
        m_PartialCharge( PARTIAL_CHARGE),
        m_Epsilon(),
        m_VdwRadius(),
        m_Mass( /*ElementDataMap()( ATOM_NAME)->GetMass()*/),
        m_ResidueID( RESIDUE_ID),
        m_AtomID( ATOM_ID),
        m_AttractiveVdwFactor( std::numeric_limits< float>::max()),
	    m_RepulsiveVdwFactor( std::numeric_limits< float>::max()),
	    m_MolID( MOL_ID)
	      {}



        //! construct from data
        Atom
        (
                const std::string &ATOM_TYPE,
                const std::string &ATOM_NAME,
                const std::string & RESIDUE_TYPE,
                const math::Vector3N &POS,
                const int &RESIDUE_ID,
                const int &ATOM_ID,
                const float &MASS,
                const float &PARTIAL_CHARGE,
                const float &EPSILON,
                const float &VDW_RADIUS,
		const int &MOL_ID = 0
        )
        : m_Position( POS),
        m_AtomType( ATOM_TYPE),
        m_AtomName( ATOM_NAME),
        m_ResidueType( RESIDUE_TYPE),
        m_Mol(),
        m_PartialCharge( PARTIAL_CHARGE),
        m_Epsilon( EPSILON),
        m_VdwRadius( VDW_RADIUS),
        m_Mass( MASS),
        m_ResidueID( RESIDUE_ID),
        m_AtomID( ATOM_ID),
        m_AttractiveVdwFactor(),
	    m_RepulsiveVdwFactor(),
	    m_MolID( MOL_ID)
        {
                	CalcVdwFactors();
        }



        //! copy constructor
        Atom( const Atom &ORIGINAL)
        : m_Position( ORIGINAL.m_Position),
        m_AtomType( ORIGINAL.m_AtomType),
        m_AtomName( ORIGINAL.m_AtomName),
        m_ResidueType( ORIGINAL.m_ResidueType),
        m_Mol( ORIGINAL.m_Mol),
        m_PartialCharge( ORIGINAL.m_PartialCharge),
        m_Epsilon( ORIGINAL.m_Epsilon),
        m_VdwRadius( ORIGINAL.m_VdwRadius),
        m_Mass( ORIGINAL.m_Mass),
        m_ResidueID( ORIGINAL.m_ResidueID),
        m_AtomID( ORIGINAL.m_AtomID),
        m_AttractiveVdwFactor( ORIGINAL.m_AttractiveVdwFactor),
	    m_RepulsiveVdwFactor( ORIGINAL.m_RepulsiveVdwFactor),
	    m_MolID( ORIGINAL.m_MolID)
        {}

        //! virtual destructor
        virtual ~Atom(){}

        //! virtual copy constructor
        virtual Atom *Clone() const{ return new Atom( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////


        virtual const math::Vector3N &GetPosition() const
        { return m_Position;}

        virtual void SetPosition( const math::Vector3N &POSITION)
        { m_Position = POSITION;}

		virtual void SetPosition( const float &X, const float &Y, const float &Z)
		{
			std::vector< float>::iterator
				itr = m_Position.begin();
			*itr = X;
			++itr;
			*itr = Y;
			++itr;
			*itr = Z;
		}

        virtual const std::string &GetType() const
        { return m_AtomType;}

        virtual const std::string &GetAtomName() const
        { return m_AtomName;}

        virtual const std::string &GetResidueType() const
        { return m_ResidueType;}

        virtual void SetResidueType( const std::string &RESIDUE_TYPE)
        { m_ResidueType = RESIDUE_TYPE;}

        virtual const int &GetResidueID() const
        { return m_ResidueID;}

        virtual void SetResidueID( const int &ID)
        { m_ResidueID = ID;}

        virtual const int &GetMolID() const
        { return m_MolID;}

        virtual void SetMolID( const int &ID)
        { m_MolID = ID;}

        virtual void SetOwningMolecule( MoleculeIF< Atom> *MOL)
        { m_Mol = MOL;}

        virtual MoleculeIF< Atom> *GetOwningMolecule() const
        { return m_Mol;}

        virtual const float &GetPartialCharge() const
        { return m_PartialCharge;}

        virtual void SetPartialCharge( const float &CHARGE)
        { m_PartialCharge = CHARGE;}

        virtual const float &GetEpsilon() const
        { return m_Epsilon;}

        virtual void SetEpsilon( const float &VALUE)
        { m_Epsilon = VALUE;}

        virtual const float &GetVanDerWaalsRadius() const
        { return m_VdwRadius;}

        virtual void SetVanDerWaalsRadius( const float &RADIUS)
        { m_VdwRadius = RADIUS;}

        virtual const float &GetMass() const
        { return m_Mass;}

        virtual void SetMass( const float &MASS)
        {
        	m_Mass = MASS;
        }

        const float &GetAttractiveVdwFactor() const
      	{
        	return m_AttractiveVdwFactor;
      	}

        const float &GetRepulsiveVdwFactor() const
      	{
        	return m_RepulsiveVdwFactor;
      	}

        void CalcVdwFactors()
        {
        	float
				radius = pow( m_VdwRadius, 3);
        	m_AttractiveVdwFactor = sqrt( fabs( m_Epsilon)) * radius;
        	m_RepulsiveVdwFactor = m_AttractiveVdwFactor * radius;
        	DebugWrite( "vdw-factor: " << m_AttractiveVdwFactor << " epsilon: " << m_Epsilon << " radius: " << m_VdwRadius);
        }

        virtual const int &GetAtomID() const
        { return m_AtomID;}

        virtual void SetAtomID( const int &ID)
        { m_AtomID = ID;}

        /////////////////////////
        //      Read/Write     //
        /////////////////////////


        virtual std::istream &Read( std::istream &STREAM)
        {
        	std::string
				tmp;
        	STREAM >> tmp;
        	if( tmp != GetClassName())
        	{
        		std::cout << __PRETTY_FUNCTION__ << "====> wrong class name";
        	}
            STREAM >> m_AtomID;
            STREAM >> m_AtomType;
            STREAM >> m_AtomName;
            STREAM >> m_ResidueID;
            STREAM >> m_ResidueType;
            STREAM >> m_Position[0];
            STREAM >> m_Position[1];
            STREAM >> m_Position[2];
            STREAM >> m_PartialCharge;
            STREAM >> m_Epsilon;
            STREAM >> m_VdwRadius;
            STREAM >> m_Mass;
            STREAM >> m_AttractiveVdwFactor;
            STREAM >> m_RepulsiveVdwFactor;
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
//            STREAM << GetClassName() << " atom-name atom-type residue-name id (x y z) partial-charge epsilon vdw-radius / element-data " << std::endl;
            STREAM << GetClassName() << "  ";
            STREAM << m_AtomID << " ";
            STREAM << m_AtomType << "  ";
            STREAM << m_AtomName << "  ";
            STREAM << m_ResidueID << " ";
            STREAM << m_ResidueType << "  ";
            STREAM << m_Position[ 0] << "  ";
            STREAM << m_Position[ 1] << "  ";
            STREAM << m_Position[ 2] << "  ";
            STREAM << m_PartialCharge << "  ";
            STREAM << m_Epsilon << "  ";
            STREAM << m_VdwRadius << "  ";
            STREAM << m_Mass << "  ";
            STREAM << m_AttractiveVdwFactor << " ";
            STREAM << m_RepulsiveVdwFactor << " ";
            STREAM << std::endl;
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
        {
            STREAM << "graphics 0 color red" <<  std::endl;
            STREAM << "graphics 0 sphere {" << GetPosition()(0) << " " << GetPosition()(1) << " " << GetPosition()(2) << "} radius 0.1 resolution 21" << std::endl;
            return STREAM;
        }

    }; // end class Atom


	mol::Atom CreateNeutralAtom();


} // end namespace mol




#endif /* SIMPLE_ATOM_H_ */
