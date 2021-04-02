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


#ifndef ELECTROSTATIC_ENERGY_H_
#define ELECTROSTATIC_ENERGY_H_


#include "../readwrite/class_name.h"


namespace phys
{

    class ElectrostaticEnergy
    : public math::Function< store::ShPtrVec< 2, mol::Atom>, float>
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        data m_Y;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ElectrostaticEnergy()
        : math::Function< store::ShPtrVec< 2, mol::Atom>, float>(),
        m_Y()
        {}

        //! construct from data
        ElectrostaticEnergy( const data &DATA)
        : math::Function< store::ShPtrVec< 2, mol::Atom>, float>(),
        m_Y( DATA)
        {}

        //! copy constructor
        ElectrostaticEnergy( const ElectrostaticEnergy &ORIGINAL)
        : math::Function< store::ShPtrVec< 2, mol::Atom>, float>( ORIGINAL),
        m_Y( ORIGINAL.m_Y)
        {}

        //! virtual destructor
        virtual ~ElectrostaticEnergy(){}

        //! virtual copy constructor
        virtual ElectrostaticEnergy *Clone() const{ return new ElectrostaticEnergy( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

    }; // end class ElectrostaticEnergy
} // end namespace phys




#endif /* ELECTROSTATIC_ENERGY_H_ */
