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


#ifndef CONSTANT_GRIDPOINT_H
#define CONSTANT_GRIDPOINT_H

#include "../math/function.t.h"
#include "../string/io_string_functions.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"

namespace store
{

    template< typename t_INPUT, typename t_RETURN>
    class ConstantGridPoint
    : public math::Function< t_INPUT, t_RETURN>
    {
    protected:
        t_RETURN  m_Return;
        int       m_NSP;       //!< nearest surface point id
        float     m_Energy;

    public:
        ConstantGridPoint()
        : m_Return(),
          m_NSP( std::numeric_limits< int>::max()),
          m_Energy( 0.0)
        {}

        ConstantGridPoint( const t_RETURN &RETURN, const int &NSP = std::numeric_limits< int>::max(), const float &ENERGY = 0.0)
        : m_Return( RETURN),
          m_NSP( NSP),
          m_Energy( ENERGY)
        {}

        ConstantGridPoint( const ConstantGridPoint & CGP)
        : m_Return( CGP.m_Return),
          m_NSP( CGP.m_NSP),
          m_Energy( CGP.m_Energy)
        {}

        virtual ~ConstantGridPoint(){}

        virtual ConstantGridPoint *Clone() const
        { return new ConstantGridPoint( *this);}

        virtual const int &GetNSP() const
        {
            return m_NSP;
        }

        t_RETURN &ReturnValue()
        {
        	return m_Return;
        }

        const t_RETURN &GetReturnValue() const
        {
        	return m_Return;
        }

        virtual float Energy( const mol::Atom &ATOM) const
        {
        	return m_Energy;
        }

        virtual const float & GetEnergy() const
        {
        	return m_Energy;
        }

        float &Energy()
        {
        	return m_Energy;
        }

//        virtual void SetReturnData( const t_RETURN &DATA)
//        {
//        	m_Energy *= DATA;
//
//        }


//        virtual
//        ConstantGridPoint *
//        MolTypeFilter( const std::string &MOL_TYPE)
//        {
//	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
//	    return new ConstantGridPoint( *this);
//        }

        virtual t_RETURN & operator()( const t_INPUT &DATA) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return math::Function< t_INPUT, t_RETURN>::s_Tmp = m_Return;
        }


        virtual std::ostream & Write( std::ostream &STREAM) const
        {
	    STREAM << util::EnumHandler< util::FunctorEnum>().String( util::e_ConstantGridPoint) << "  ";
            STREAM.flush();
            STREAM << m_Return;
            STREAM << "nsp: " << m_NSP << std::endl;
            STREAM << "energy: " << m_Energy << std::endl;
            return STREAM;
        }

//        virtual std::ostream &WriteVmdCommands( std::ostream &STREAM) const
//        {
//            return STREAM;
//        }


        virtual std::istream &Read( std::istream &STREAM)
        {
        	DebugWrite(  __PRETTY_FUNCTION__ );
            std::string str;
            if( !STREAM)
            {
				std::cout << "===> stream closed in " << __PRETTY_FUNCTION__ << std::endl;
				exit( -1);
            }
//  	          STREAM >> str;
//            assert( str == GetClassName());
            m_Return = factory::BuildNeutralObject< t_RETURN>();
            STREAM >> m_Return >> str >> m_NSP;
            if( str != "nsp:")
            {
            	std::cout << "===> should be nsp: " << str << " in " << __PRETTY_FUNCTION__ << std::endl;
            	std::cout << "nsp: " << m_NSP << std::endl;
            	exit( -1);
            }

            STREAM >> str >> m_Energy;
            if( str != "energy:")
            {
            	std::cout << "===> string <" << str << "> should be \'energy:\'" << std::endl;
            	exit( -1);
            }

            return STREAM;
        }

        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_ConstantGridPoint;
        }


    }; // end class ConstantGridPoint
} // end namespace store
#endif
