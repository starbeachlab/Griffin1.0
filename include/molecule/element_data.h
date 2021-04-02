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
//! Default element type information. 									 
//!									 
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef ELEMENT_DATA_H_
#define ELEMENT_DATA_H_


#include <limits>
#include <string>

#include <boost/shared_ptr.hpp>

#include "../readwrite/stream_operator.h"



namespace mol
{

    class ElementData
    : public StreamOperator
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
        std::string    m_Name;
        std::string    m_Symbol;
        char           m_Period;
        size_t         m_AtomNr;
        float         m_Mass;
        float         m_CovalentRadius;
        float         m_VdwRadius;


    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ElementData()
        : m_Name( "UNDEFINED"),
        m_Symbol( "UNDEFINED"),
        m_Period( 'X'),
        m_AtomNr( std::numeric_limits< size_t>::max()),
        m_Mass( std::numeric_limits< float>::max()),
        m_CovalentRadius( std::numeric_limits< float>::max()),
        m_VdwRadius( std::numeric_limits< float>::max())//,
        {}

        //! construct from data
        ElementData
        (
                const std::string &NAME,
                const std::string &SYMBOL =  "UNDEFINED",
                const char &PERIOD = 'X',
                const size_t &ATOM_NR = std::numeric_limits< size_t>::max(),
                const float &MASS =  std::numeric_limits< float>::max(),
                const float &COVALENT_RADIUS =  std::numeric_limits< float>::max(),
                const float &VDW_RADIUS =  std::numeric_limits< float>::max()
        )
        : m_Name(NAME),
        m_Symbol(SYMBOL),
        m_Period( PERIOD),
        m_AtomNr(ATOM_NR),
        m_Mass(MASS),
        m_CovalentRadius(COVALENT_RADIUS),
        m_VdwRadius(VDW_RADIUS)
        {}

        //! copy constructor
        ElementData( const ElementData &ORIGINAL)
        : m_Name( ORIGINAL.m_Name),
        m_Symbol( ORIGINAL.m_Symbol),
        m_Period( ORIGINAL.m_Period),
        m_AtomNr( ORIGINAL.m_AtomNr),
        m_Mass( ORIGINAL.m_Mass),
        m_CovalentRadius( ORIGINAL.m_CovalentRadius),
        m_VdwRadius( ORIGINAL.m_VdwRadius)
        {}

        //! virtual destructor
        virtual ~ElementData(){}

        //! virtual copy constructor
        virtual ElementData *Clone() const{ return new ElementData( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        const std::string &GetName() const{ return m_Name;}

        const std::string &GetSymbol() const{ return m_Symbol;}

        const size_t &GetAtomNr() const{ return m_AtomNr;}

        const float &GetVanDerWaalsRadius() const { return m_VdwRadius;}

        const float &GetMass() const{ return m_Mass;}


        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            STREAM >>  m_Name;
            STREAM >>  m_Symbol;
            STREAM >>  m_Period;
            STREAM >>  m_AtomNr;
            STREAM >>  m_Mass;
            STREAM >>  m_CovalentRadius;
            STREAM >>  m_VdwRadius;
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM <<  m_Name << "  ";
            STREAM <<  m_Symbol <<  "  ";
            STREAM <<  m_Period << "  ";
            STREAM <<  m_AtomNr << "  ";
            STREAM <<  m_Mass << "  ";
            STREAM <<  m_CovalentRadius << "  ";
            STREAM <<  m_VdwRadius << "  ";
            return STREAM;
        }

    }; // end class ElementData
} // end namespace mol




#endif /* ELEMENT_DATA_H_ */
