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


#ifndef FORCE_H_
#define FORCE_H_


#include "../readwrite/stream_operator.h"
#include "potential_force.h"


namespace phys
{

    // forward declaration
    class PotentialForce;


    class Force
    : public StreamOperator
    {
    protected:
        float    m_Cutoff;
    public:
        Force( const float CUTOFF = g_FloatMax)
        : m_Cutoff( CUTOFF)
        {}

        Force( const Force &FORCE)
        : m_Cutoff( FORCE.m_Cutoff)
        {}

        //! virtual copy constructor
        virtual Force *Clone() const = 0;

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////
        virtual math::Vector3N operator()( const mol::Atom &ATOM_A, const mol::Atom &ATOM_B) const = 0;

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////
        virtual boost::shared_ptr< phys::PotentialForce> PotentialForce( const math::Vector3N & POS) const
        {
            std::cout << GetClassName() << ": should not call this function" << std::endl;
            return boost::shared_ptr< phys::PotentialForce>();
        }

//        virtual math::Vector3N PotentialForce( const math::Vector3N & POS, const mol::Atom &ATOM) const
//        { return math::Vector3N();}

        virtual boost::shared_ptr< phys::PotentialForce> GetAssociatedPotentialForceObject() const
        { return boost::shared_ptr< phys::PotentialForce>( new phys::PotentialForce());}

        virtual const float &GetCutoff() const
        { return m_Cutoff;}

        virtual void SetCutoff( const float &CUTOFF)
        { m_Cutoff = CUTOFF;}
        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM) = 0;

        virtual std::ostream &Write( std::ostream &STREAM)const = 0;

        virtual std::string GetClassName() const = 0;

    }; // end class Force
} // end namespace phys




#endif /* FORCE_H_ */
