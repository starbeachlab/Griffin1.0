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


#ifndef POTENTIAL_FORCE_H_
#define POTENTIAL_FORCE_H_


#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"
#include "../math/vector3N.h"
#include "../molecule/simple_molecule.t.h"
#include "../molecule/atom.h"
#include "../external/boost_functions.h"

namespace phys
{

    class PotentialForce
    : public StreamOperator
    {
    protected:
        math::Vector3N   m_Force;
        float			 m_Energy;


    public:

        PotentialForce()
        : m_Force( 0.0, 0.0, 0.0),
        m_Energy( 0)
        {}

        PotentialForce( const math::Vector3N &FORCE, const float &ENERGY)
        : m_Force( FORCE),
        m_Energy( ENERGY)
        {}

        PotentialForce( const PotentialForce &FORCE)
        : m_Force( FORCE.m_Force),
		m_Energy( FORCE.m_Energy)
        {}

        //! virtual copy constructor
        virtual PotentialForce *Clone() const;

        math::Vector3N &
        GetSetVector()
        {
        	return m_Force;
        }

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        //! returns the potential force without requiring any actual information from atom; this function needs to be overwritten in derived classes for interactions requiring specific information of the atom
        virtual math::Vector3N operator()( const mol::Atom &ATOM) const;

//        virtual math::Vector3N &operator()( const mol::Atom &ATOM);

        virtual void operator = ( const math::Vector3N &FORCE);

//        virtual void operator += ( const math::Vector3N &FORCE, const float &ENERGY);  /// TODO: CHECK THIS!!!

        virtual PotentialForce operator += ( const PotentialForce &FORCE)   /// TODO: CHECK THIS!!!
        {
            m_Force += FORCE.m_Force;
            m_Energy += FORCE.m_Energy;
            return *this;
        }

        virtual float Energy( const mol::Atom &ATOM) const;

        float &GetSetEnergy()
        {
        	return m_Energy;
        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        void AddToPotentialForce( const math::Vector3N &FORCE, const float &ENERGY);

        void AddToPotentialForce( const math::Vector3N &FORCE);

        void AddToPotentialForce( const float &ENERGY);

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        //! writes forces without class identifier; to be called from derived classes
        virtual std::ostream &WriteData( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class
} // end namespace $




#endif /* PotentialForce_H_ */
