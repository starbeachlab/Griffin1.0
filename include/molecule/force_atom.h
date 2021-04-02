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


#ifndef FORCE_ATOM_H_
#define FORCE_ATOM_H_

#include "../readwrite/stream_operator.h"
#include "../string/io_string_functions.h"



namespace mol
{

    class ForceAtom
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        boost::shared_ptr< Atom>            m_Atom;
        math::Vector3N  m_Force;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ForceAtom()
        : m_Atom(),
        m_Force()
        {}

        //! construct from data
        ForceAtom( const boost::shared_ptr< Atom> &ATOM, const math::Vector3N &FORCE = math::Vector3N())
        : m_Atom( ATOM),
        m_Force( FORCE)
        {}

        //! copy constructor
        ForceAtom( const ForceAtom &ORIGINAL)
        : m_Atom( ORIGINAL.m_Atom),
        m_Force( ORIGINAL.m_Force)
        {}

        //! virtual destructor
        virtual ~ForceAtom(){}

        //! virtual copy constructor
        virtual ForceAtom *Clone() const{ return new ForceAtom( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual const math::Vector3N &GetForce() const{ return m_Force;}

        virtual const boost::shared_ptr< Atom> &GetAtom() const{ return m_Atom;}


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

    }; // end class ForceAtom
} // end namespace mol




#endif /* FORCE_ATOM_H_ */
