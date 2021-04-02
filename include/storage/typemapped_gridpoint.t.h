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


#ifndef TYPEMAPPED_GRIDPOINT_H_
#define TYPEMAPPED_GRIDPOINT_H_


#include "../math/function.t.h"
#include "../storage/map.t.h"
#include "../string/io_string_functions.h"
#include "constant_gridpoint.t.h"

namespace store
{
    template< typename t_INPUT, typename t_RETURN>
    class TypeMappedGridPoint
    : public math::Function< t_INPUT, t_RETURN>
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////

    	store::Map< std::string, ConstantGridPoint< t_INPUT, t_RETURN> >  m_TypeMap;  // instead of three maps below ???

    	static std::vector< std::string>   m_MolTypes;

//    	store::Map< std::string, t_RETURN> m_TypeMap;
//        store::Map< std::string, int>      m_NSP;       //!< nearest surface point ids
//        store::Map< std::string, float>    m_Energy;       //!< nearest surface point ids


    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
         TypeMappedGridPoint()
         : math::Function< t_INPUT, t_RETURN>(),
         m_TypeMap()
//         m_NSP(),
//         m_Energy()
         {}

//         //! construct from store::Map< std::string, t_RETURN>
//		TypeMappedGridPoint
//		(
//				const store::Map< std::string, t_RETURN> &DATA,
//				const store::Map< std::string, int> &NSP,
//				const store::Map< std::string, float> &ENERGY
//		)
//		: math::Function< t_INPUT, t_RETURN>(),
//		m_TypeMap( DATA),
//		m_NSP( NSP),
//		m_Energy( ENERGY)
//		{}

            //! copy constructor
        TypeMappedGridPoint( const TypeMappedGridPoint &ORIGINAL)
        : math::Function< t_INPUT, t_RETURN>( ORIGINAL),
        m_TypeMap( ORIGINAL.m_TypeMap)
//        m_NSP( ORIGINAL.m_NSP),
//        m_Energy( ORIGINAL.m_Energy)
        {}

        //! virtual destructor
        virtual ~TypeMappedGridPoint(){}

        //! virtual copy constructor
        virtual TypeMappedGridPoint *Clone() const{ return new TypeMappedGridPoint( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        virtual t_RETURN & operator()( const t_INPUT &DATA) const
        {
	    DebugWrite( __PRETTY_FUNCTION__);
            return math::Function< t_INPUT, t_RETURN>::s_Tmp = m_TypeMap( MolNameToMolType( ( *DATA.GetOwningMolecule()).GetType())).GetReturnValue();
//            return m_TypeMap( DATA->GetTypeString())( DATA);           /// BAD: atom will have to provide molecule information... || have typestring: mol:atom and split at ':' !!!
        }

//        virtual
//        math::Function< t_INPUT, t_RETURN> *
//        MolTypeFilter( const std::string &MOL_TYPE)
//        {
//	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
//	    return new ConstantGridPoint< t_INPUT, t_RETURN>( m_TypeMap( MolNameToMolType( MOL_TYPE)));
//        }

        const std::vector< std::string> &GetMolTypes() const
		{
        	return m_MolTypes;
		}

        void SetMolTypes( const std::vector< std::string> &TYPES)
        {
        	m_MolTypes = TYPES;
        }

        virtual float Energy( const t_INPUT &DATA) const
        {
        	return m_TypeMap( MolNameToMolType( ( *DATA.GetOwningMolecule()).GetType())).GetEnergy();
        }

        virtual std::string MolNameToMolType( const std::string &NAME) const
        {
	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << NAME);
        	if( std::find( m_MolTypes.begin(), m_MolTypes.end(), NAME) == m_MolTypes.end())
        	{
		    DebugWrite( "returns 'all'");
        		return "all";
        	}
        	return NAME;
        }

        virtual 
	store::Map< std::string, ConstantGridPoint< t_INPUT, t_RETURN> > &
	TypeMap()
        {
            return m_TypeMap;
        }

	ConstantGridPoint< t_INPUT, t_RETURN> &
	GetFromTypeMap( const std::string &TYPE)
	{
	    return m_TypeMap( MolNameToMolType( TYPE));
	}
		
//        virtual const int &GetNSP() const
//        {
//        	std::cout << "===> WHY THIS? " << __PRETTY_FUNCTION__ << std::endl;
//            return m_NSP.begin()->second;
//        }
        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

//        virtual void Insert( const t_INPUT &INPUT, const t_RETURN &OUT)
//        {
//            if( !m_TypeMap.IsValidKey( ( *INPUT->GetOwningMolecule()).GetType()))
//            {
//                m_TypeMap.InsertNewKeyAndValue( ( *INPUT->GetOwningMolecule()).GetType(), OUT);
//            }
//        }

        virtual void InsertNewKeyAndValues( const std::string &KEY, const t_RETURN &OUT, const int &NSP, const float &ENERGY)
        {
            if( !m_TypeMap.IsValidKey( KEY))
            {
                m_TypeMap.InsertNewKeyAndValue( KEY, ConstantGridPoint< t_INPUT, t_RETURN>( OUT, NSP, ENERGY));
//                m_NSP.InsertNewKeyAndValue( KEY, NSP);
//                m_Energy.InsertNewKeyAndValue( KEY, ENERGY);
            }
        }

//        virtual void InsertTypeAndNSP( const std::string &TYPE, const int &NSP)
//        {
//        	m_NSP.InsertNewKeyAndValue( TYPE, NSP);
//        }

//        void SetMolTypes( const std::vector< std::string> &MOLTYPES)
//        {
//        	for( std::vector< std::string>::const_iterator itr = MOLTYPES.begin(); itr != MOLTYPES.end(); ++itr)
//        		if( *itr != "all")
//        		{
//        			m_MolTypes.push_back( *itr);
//        		}
//        }

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
	    DebugWrite(  __PRETTY_FUNCTION__ );
            std::string str, key;
            if( !STREAM)
            {
				std::cout << "===> stream closed in " << __PRETTY_FUNCTION__ << std::endl;
				exit( -1);
            }

            float energy;
            int nsp, size;
            STREAM >> size;

            for( int i = 0; i < size; ++i)
            {
                t_RETURN object = factory::BuildNeutralObject< t_RETURN>();
            	STREAM >> str >> key;
            	STREAM >> object;
            	STREAM >> str >> energy;
            	STREAM >> str >> nsp;

            	m_TypeMap.InsertNewKeyAndValue( key, ConstantGridPoint< t_INPUT, t_RETURN>( object, nsp, energy));
//            	m_Energy.InsertNewKeyAndValue( key, energy);
//            	m_NSP.InsertNewKeyAndValue( key, nsp);
             }


            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << util::EnumHandler< util::FunctorEnum>().String( util::e_TypeMappedGridPoint) << " " << m_TypeMap.size() << std::endl;

//            typename store::Map< std::string, t_RETURN>::const_iterator ret_itr =  m_TypeMap.begin();
//            store::Map< std::string, int>::const_iterator nsp_itr = m_NSP.begin();
//            store::Map< std::string, float>::const_iterator en_itr = m_Energy.begin();
//
//			while( ret_itr != m_TypeMap.end())
//			{
//				STREAM << "type: " << ret_itr->first << "  " << ret_itr->second;
//				STREAM << "energy: " << en_itr->second << std::endl;
//				STREAM << "nsp: " << nsp_itr->second << std::endl;
//				++ret_itr;
//				++en_itr;
//				++nsp_itr;
//			}

	    for( typename Map< std::string, ConstantGridPoint< t_INPUT, t_RETURN> >::const_iterator itr = m_TypeMap.begin(); itr != m_TypeMap.end(); ++itr)
	    {
		STREAM << "type: " << itr->first << " " << itr->second.GetReturnValue();
		STREAM << "energy: " << itr->second.GetEnergy() << std::endl;
		STREAM << "nsp: " << itr->second.GetNSP() << std::endl;
	    }

            return STREAM;
        }

 
        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_TypeMappedGridPoint;
        }


    }; // end class TypeMappedGridPoint

    template< typename t_INPUT, typename t_RETURN>
    std::vector< std::string> TypeMappedGridPoint< t_INPUT, t_RETURN>::m_MolTypes = std::vector< std::string>();

} // end namespace store

#endif /* TYPEMAPPED_GRIDPOINT_H_ */
