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


#ifndef RECURSIVE_TYPEMAPPED_GRIDPOINT_H_
#define RECURSIVE_TYPEMAPPED_GRIDPOINT_H_


#include "../storage/map.t.h"
#include "../math/function.t.h"

#include "interaction_gridpoint.h"
//#include "transient_interaction_gridpoint.h"

namespace  store
{

    template< typename t_INPUT, typename t_RETURN>
    class RecursiveTypeMappedGridPoint
    : public math::Function< t_INPUT, t_RETURN>
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////
    	static std::vector< std::string>  		   						     					   s_MolTypes;
        store::Map< std::string, boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > >         m_TypeMap;

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        RecursiveTypeMappedGridPoint( const int &NSP = std::numeric_limits< int>::max())
        : math::Function< t_INPUT, t_RETURN>(),
        m_TypeMap()
        {}


        //! construct from map< std::string, t_RETURN>
        RecursiveTypeMappedGridPoint( const store::Map< std::string, boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > > &DATA, const int &NSP = std::numeric_limits< int>::max())
        : math::Function< t_INPUT, t_RETURN>(),
        m_TypeMap( DATA)
        {}


        //! copy constructor
        RecursiveTypeMappedGridPoint( const RecursiveTypeMappedGridPoint &ORIGINAL)
        : math::Function< t_INPUT, t_RETURN>( ORIGINAL),
        m_TypeMap( ORIGINAL.m_TypeMap)
        {}

        //! virtual destructor
        virtual ~RecursiveTypeMappedGridPoint(){}

        //! virtual copy constructor
        virtual RecursiveTypeMappedGridPoint *Clone() const{ return new RecursiveTypeMappedGridPoint( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////


//        virtual
//        math::Function< t_INPUT, t_RETURN> *
//        MolTypeFilter( const std::string &MOL_TYPE)
//        {
//	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << MOL_TYPE);
//	    return  m_TypeMap( MolNameToMolType( MOL_TYPE))->Clone();
//        }


        virtual t_RETURN & operator()( const std::string &TYPE, const t_INPUT &DATA) const
        {
            return ( math::Function< t_INPUT, t_RETURN>::s_Tmp = m_TypeMap( MolNameToMolType( TYPE))->operator()( DATA));
        }

        virtual t_RETURN & operator()( const t_INPUT &DATA) const
        {
	    DebugWrite( __PRETTY_FUNCTION__);
            return ( math::Function< t_INPUT, t_RETURN>::s_Tmp = m_TypeMap( MolNameToMolType( DATA.GetOwningMolecule()->GetType()))->operator()( DATA));  // not general!!!
        }

        virtual float Energy(  const t_INPUT &DATA) const
        {
        	return m_TypeMap( MolNameToMolType( DATA.GetOwningMolecule()->GetType()))->Energy( DATA);
        }

        virtual std::string MolNameToMolType( const std::string &NAME) const
        {
	    DebugWrite( __PRETTY_FUNCTION__ << " type: " << NAME);
        	if( std::find( s_MolTypes.begin(), s_MolTypes.end(), NAME) == s_MolTypes.end())
        	{
		    DebugWrite( "returns 'all'");
		    return "all";
        	}
        	return NAME;
        }

        virtual const int &GetNSP() const
        {
        	std::cout << "===> should not be called (individual grid point functions should be called) " << __PRETTY_FUNCTION__ << std::endl;
            return g_IntMax;
        }
        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual void InsertNewKeyAndValue( const std::string &TYPE, const boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > &FUNCTION)
        {
        	DebugWrite( __FUNCTION__ << " " << TYPE);
            m_TypeMap.InsertNewKeyAndValue( TYPE, FUNCTION);
        }

        store::Map< std::string, boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > > &
        TypeMap()
        {
        	return m_TypeMap;
        }

		boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > GetFromTypeMap( const std::string &TYPE)
		{
			return m_TypeMap( MolNameToMolType( TYPE));
		}
		
	    
		const std::vector< std::string> &GetMolTypes() const
		{
			return s_MolTypes;
		}

        void SetMolTypes( const std::vector< std::string> &MOLTYPES)
        { // todo: check whether elements exist already?
        	for( std::vector< std::string>::const_iterator itr = MOLTYPES.begin(); itr != MOLTYPES.end(); ++itr)
        		if( *itr != "all")
				{
        			s_MolTypes.push_back( *itr);
				}
        }

        boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >
        GetGridPoint( const std::string &MOLNAME)
		{
        	return m_TypeMap( MOLNAME);
		}

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
           	DebugWrite( __PRETTY_FUNCTION__);
        	std::string
				str,
				moltype;
        	int
				size;
        	math::Function< t_INPUT,t_RETURN>
				func;
        	boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >
				ptr;
        	boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >
				void_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N >());
        	boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >
				surf_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N >());

        	util::FunctorEnum
				gp_type;

//		util::EnumHandler< util::FunctorEnum>().WriteString2ID();

        	STREAM >> str;
//        	if( str != GetClassName())
//        	{
//        		std::cout << "===> invalid type id (" << str << ") in " << __PRETTY_FUNCTION__ << std::endl;
//        	}
        	STREAM >> size;
        	for( int i = 0; i < size; ++i)
        	{
        		STREAM >> moltype >> str;
        		DebugWrite( "now: " << moltype << " " << str << " in " << __FUNCTION__);
                if( m_TypeMap.find( moltype) != m_TypeMap.end())
                {
                    std::cout << "===> moltype " << moltype << " already exists in " << __PRETTY_FUNCTION__ << std::endl;
                    exit( -1);
                }

                gp_type = util::EnumHandler< util::FunctorEnum>().ID( str);

//		std::cout << "now: <" << moltype << "> <" << str << ">  <" << gp_type << ">" << std::endl;

                if( gp_type == util::e_VoidGridPoint)
                {
                    ptr = boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >( void_grid_point);
                    m_TypeMap.InsertNewKeyAndValue( moltype, void_grid_point);
                }
                else if( gp_type == util::e_SurfGridPoint)
                {
                    ptr = boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >( surf_grid_point);
                    m_TypeMap.InsertNewKeyAndValue( moltype, void_grid_point);
                }
                else if( gp_type == util::e_ConstantGridPoint)
                {
                	store::ConstantGridPoint< mol::Atom, math::Vector3N > grid_point;
                	STREAM >> grid_point;
                    ptr = boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >( new store::ConstantGridPoint< mol::Atom, math::Vector3N >( grid_point));
                    m_TypeMap.InsertNewKeyAndValue( moltype, ptr);
                }
                else if( gp_type == util::e_InteractionGridPoint)
                {
                    store::InteractionGridPoint grid_point;
                    STREAM >> grid_point;
                    ptr = boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >( new store::InteractionGridPoint( grid_point));
                    m_TypeMap.InsertNewKeyAndValue( moltype, ptr);
                }
//                else if( str == store::TransientInteractionGridPoint().GetClassName())
//                {
//                    store::TransientInteractionGridPoint grid_point;
//                    STREAM >> grid_point;
//                    ptr = boost::shared_ptr< math::Function< t_INPUT, t_RETURN> >( new store::TransientInteractionGridPoint( grid_point));
//                    m_TypeMap.InsertNewKeyAndValue( moltype, ptr);
//                }
                else
                {
                	std::cout << "\n====> undefined gridpoint type (" << str << "/" << gp_type << ") in " << __PRETTY_FUNCTION__ << std::endl;
                }

#ifndef NO_INTERACTION_FORCES
#endif

//                Memory();
       	}
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
	    STREAM << util::EnumHandler< util::FunctorEnum>().String( util::e_RecursiveTypeMappedGridPoint) << " ";
            STREAM << m_TypeMap << std::endl;
            return STREAM;
        }


        virtual util::FunctorEnum GetClassID() const
        {
        	return util::e_RecursiveTypeMappedGridPoint;
        }


    }; // end class

    template< typename t_INPUT, typename t_RETURN>
    std::vector< std::string> RecursiveTypeMappedGridPoint< t_INPUT, t_RETURN>::s_MolTypes = std::vector< std::string>();

} // end namespace $




#endif /* RECURSIVE_TYPEMAPPED_GRIDPOINT_H_ */
