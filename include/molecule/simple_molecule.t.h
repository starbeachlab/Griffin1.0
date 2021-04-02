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


#ifndef SIMPLE_MOLECULE_H_
#define SIMPLE_MOLECULE_H_

#include <algorithm>

#include "../storage/shared_pointer_vector.t.h"
#include "../utilities/id.h"
#include "../string/io_string_functions.h"
#include "../readwrite/stream_operator.h"
#include "moleculeIF.h"
#include "../math/vector3N.h"
#include "../storage/vector3N.t.h"

namespace mol
{
    template< typename t_ATOM>
    class SimpleMolecule
    : public MoleculeIF< t_ATOM>, public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        store::ShPtrVec< t_ATOM> m_Atoms;
        std::string              m_Type;
        size_t                   m_ID;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        SimpleMolecule( const std::string &TYPE = "", const size_t &ID =  util::ID().GetID())
        : m_Atoms(),
        m_Type( TYPE),
        m_ID( ID)
        {}

        //! copy constructor
        SimpleMolecule( const SimpleMolecule &ORIGINAL)
        : m_Atoms( ORIGINAL.m_Atoms),
        m_Type( ORIGINAL.m_Type),
        m_ID( ORIGINAL.m_ID)
        {
        	std::cout << __PRETTY_FUNCTION__ << std::endl;
        }

        //! construct from data
        SimpleMolecule( const store::ShPtrVec< t_ATOM> &ATOMS, const std::string &TYPE = "", const size_t &ID =  util::ID().GetID())
        : m_Atoms( ATOMS),
        m_Type( TYPE),
        m_ID( ID)
        { LinkAtomsToThisMolecule();}

        //! virtual destructor
        virtual ~SimpleMolecule(){}

        //! virtual copy constructor
        //!(needed for copying pointers to derived classes)
        virtual SimpleMolecule *Clone() const{ return new SimpleMolecule( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const store::ShPtrVec< t_ATOM> &GetAtoms() const{ return m_Atoms;}

        virtual store::ShPtrVec< t_ATOM> &Atoms(){ return m_Atoms;}

        //! returns all atoms matching TYPE
        store::ShPtrVec< t_ATOM> GetAtoms( const std::string &TYPE) const
        {
            store::ShPtrVec< t_ATOM> atoms;
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                if( ( *itr)->GetType() == TYPE)
                { atoms.push_back( *itr);}
            }
            return atoms;
        }

        virtual const size_t &GetID() const{ return m_ID;}

        virtual void SetID( const size_t &ID)
        {
        	m_ID = ID;
        }

        virtual const std::string &GetType() const{ return m_Type;}

        virtual void SetType( const std::string &TYPE){ m_Type = TYPE;}


//        math::Vector3N GetCms() const;
//
//        float MaxRadius() const;
//
//        float TotalCharge() const;
//
//        math::Vector3N TotalDipole() const;

        virtual void AddAtom( const boost::shared_ptr< t_ATOM> &ATOM)
        {
//            DebugWrite( __PRETTY_FUNCTION__);
            m_Atoms.push_back( ATOM);
            ATOM->SetOwningMolecule( this);
        }

        virtual void LinkAtomsToThisMolecule()
        {
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                ( *itr)->SetOwningMolecule( this);
            }
        }

        virtual math::Vector3N CMS() const
        {
            math::Vector3N mass_times_location;
            float mass_sum( 0.0);
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                DebugWrite( ( *itr)->GetType() << " mass: " << ( *itr)->GetMass());
                mass_sum += ( *itr)->GetMass();
                mass_times_location += ( *itr)->GetMass() * ( *itr)->GetPosition();
            }
            if( mass_sum > 0.0)
            {
                return mass_times_location / mass_sum;
            }
            return math::Vector3N( std::numeric_limits< float>::max());
        }


        virtual math::Vector3N GeometricCenter() const
        {
            math::Vector3N center;
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                center += ( *itr)->GetPosition();
            }
            return center / float( std::max( size_t( 1), m_Atoms.size()));
        }


        virtual store::Vector3N< std::pair< float, float> > Limits() const
        {
            store::Vector3N< std::pair< float, float> > limits;
            std::vector< std::pair< float, float> >::iterator limit_itr;

            typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator atom_itr = m_Atoms.begin();

            math::Vector3N pos = ( *atom_itr)->GetPosition();

            limits( 0) = std::make_pair( pos( 0), pos( 0));
            limits( 1) = std::make_pair( pos( 1), pos( 1));
            limits( 2) = std::make_pair( pos( 2), pos( 2));

            std::vector< float>::const_iterator pos_itr;

            for( ; atom_itr != m_Atoms.end(); ++atom_itr)
            {
                pos = ( *atom_itr)->GetPosition();
                pos_itr = pos.begin();
                limit_itr = limits.begin();
                for( ; pos_itr != pos.end(); ++pos_itr, ++limit_itr)
                {
                    if( *pos_itr < limit_itr->first)
                    {
                        limit_itr->first = *pos_itr;
                    }
                    else if( *pos_itr > limit_itr->second)
                    {
                        limit_itr->second = *pos_itr;
                    }
                }
            }

            return limits;
        }

        virtual float TotalMass() const
        {
            float mass_sum( 0.0);
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                mass_sum += ( *itr)->GetMass();
            }
            return mass_sum;
        }
        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            STREAM >> m_Atoms;
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << " ";
            STREAM << "type: <" << m_Type << ">" << std::endl;
            STREAM << "id: " << m_ID << std::endl;
            STREAM << "nr_atoms: " << m_Atoms.size() << std::endl;
            STREAM << "#atomlines: atom-name atom-type residue-name id (x y z) partial-charge epsilon vdw-radius mass atom-id residue-id(molecule-id)" << std::endl;

            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                ( *itr)->Write( STREAM);
            }
            return STREAM;
        }


        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }


        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
        {
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = m_Atoms.begin(); itr != m_Atoms.end(); ++itr)
            {
                ( *itr)->WriteVmdCommands( STREAM);
            }
            return STREAM;
        }

    }; // end class SimpleMolecule
} // end namespace mol




#endif /* SIMPLE_MOLECULE_H_ */
