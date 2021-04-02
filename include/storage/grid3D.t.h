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


#ifndef GRID_VECTOR_H
#define GRID_VECTOR_H

#include "position_grid.t.h"
#include "atom_grid_factory.t.h"
#include "typemapped_gridpoint.t.h"
#include "recursive_typemapped_gridpoint.t.h"
#include "surf_gridpoint.t.h"

namespace store
{
    template< typename t_RETURN>
    class GridVector3D
    : public PositionGrid< t_RETURN>
    {
    protected:
        std::vector< std::vector< std::vector< t_RETURN> > >   m_Data;

    public:

        GridVector3D()
        : PositionGrid< t_RETURN>(),
        m_Data()
        { /*std::cerr << "calling this constructor does not make much sense" << std::endl; exit( -1);*/}


        GridVector3D(  const math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA, const store::Vector3N< size_t> &NR_BINS)
        : PositionGrid< t_RETURN>( MIN, MAX, DELTA, NR_BINS),
        m_Data()
        {
            InitializeDataStructure();
        }

        GridVector3D( const GridVector3D & GV)
        : PositionGrid< t_RETURN>( GV),
        m_Data( GV.m_Data)
        { std::cerr << "why do you want to copy the entire grid?" << std::endl; exit( -1);}


        virtual GridVector3D *Clone() const{ return new GridVector3D( *this);}

        const PositionGrid< t_RETURN> &
        GetPositionGrid() const
        {
        	return *this;
        }

        PositionGrid< t_RETURN> &
        GetSetPositionGrid()
        {
        	return *this;
        }

        virtual void SetPositionGrid( const PositionGrid< t_RETURN> &POSGRID)
	{
//	    ( PositionGrid< t_RETURN>) *this = POSGRID;
//	    GetSetPositionGrid() = POSGRID;
	    PositionGrid< t_RETURN>::operator = (POSGRID);
        }

        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION));
        }


        virtual t_RETURN GetGridPoint( const store::Vector3N< int> &IDS) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")");

            if( IDS[0] < 0 || IDS[0] >= (int) m_Data.size() || IDS[1] < 0 || IDS[1] >= (int) m_Data[0].size() || IDS[2] < 0 || IDS[2] >= (int) m_Data[0][0].size())
            {
            	DebugWrite( "return neutral object");
                return t_RETURN( factory::BuildNeutralObject< t_RETURN>());
            }
			DebugWrite( "return value");
            return m_Data[ IDS[0]][ IDS[1]][ IDS[2]];
        }


#ifdef FORCE_INDICES
        virtual t_RETURN GetGridPoint
	    ( 
		const math::Vector3N &POSITION, 
		const std::string &MOL_TYPE, 
#ifdef CPPIO
		std::ostream &STREAM
#else
		FILE *STREAM
#endif
		) const
        {
            StandardWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION), MOL_TYPE, STREAM);
        }
#endif


        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION), MOL_TYPE);
        }




        virtual t_RETURN 
	GetGridPoint
	( 
	    const store::Vector3N< int> &IDS, 
	    const std::string &MOL_TYPE
	) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")" << " moltype: " << MOL_TYPE);

            if( IDS[0] < 0 || IDS[0] >= (int) m_Data.size() || IDS[1] < 0 || IDS[1] >= (int) m_Data[0].size() || IDS[2] < 0 || IDS[2] >= (int) m_Data[0][0].size())
            {
            	DebugWrite( "return neutral object");
                return t_RETURN( factory::BuildNeutralObject< t_RETURN>());
            }
	
            t_RETURN result = m_Data[ IDS[0]][ IDS[1]][ IDS[2]];

			if( result->GetClassID() == util::e_TypeMappedGridPoint)
			{
				TypeMappedGridPoint<  mol::Atom, math::Vector3N> * cptr = ( TypeMappedGridPoint<  mol::Atom, math::Vector3N> *) result.get();
				DebugWrite( "return value from typemapped: " << cptr->GetFromTypeMap( MOL_TYPE));
				return t_RETURN( new ConstantGridPoint<  mol::Atom, math::Vector3N>( cptr->GetFromTypeMap( MOL_TYPE)));
			}
			else if(  result->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
			{
				RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) result.get();
				DebugWrite( "return value from recursive: " << cptr->GetFromTypeMap( MOL_TYPE));
				return cptr->GetFromTypeMap( MOL_TYPE);
			}

			DebugWrite( "return value: " << result);
			return result;
        }




#ifdef FORCE_INDICES
        virtual t_RETURN 
	GetGridPoint
	( 
	    const store::Vector3N< int> &IDS, 
	    const std::string &MOL_TYPE,
#ifdef CPPIO
		std::ostream &STREAM
#else
		FILE *STREAM
#endif
	) const
        {
            StandardWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")" << " moltype: " << MOL_TYPE);

#ifdef CPPIO
	    STREAM << IDS[0] << "  " << IDS[1] << "  " << IDS[2] << "    ";
#else
	    fprintf( STREAM, "%d  %d  %d  ", IDS[0], IDS[1], IDS[2]);
#endif

            if( IDS[0] < 0 || IDS[0] >= (int) m_Data.size() || IDS[1] < 0 || IDS[1] >= (int) m_Data[0].size() || IDS[2] < 0 || IDS[2] >= (int) m_Data[0][0].size())
            {
            	DebugWrite( "return neutral object");
                return t_RETURN( factory::BuildNeutralObject< t_RETURN>());
            }
	
            t_RETURN result = m_Data[ IDS[0]][ IDS[1]][ IDS[2]];

			if( result->GetClassID() == util::e_TypeMappedGridPoint)
			{
				TypeMappedGridPoint<  mol::Atom, math::Vector3N> * cptr = ( TypeMappedGridPoint<  mol::Atom, math::Vector3N> *) result.get();
				DebugWrite( "return value from typemapped: " << cptr->GetFromTypeMap( MOL_TYPE));
				return t_RETURN( new ConstantGridPoint<  mol::Atom, math::Vector3N>( cptr->GetFromTypeMap( MOL_TYPE)));
			}
			else if(  result->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
			{
				RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) result.get();
				DebugWrite( "return value from recursive: " << cptr->GetFromTypeMap( MOL_TYPE));
				return cptr->GetFromTypeMap( MOL_TYPE);
			}

			DebugWrite( "return value: " << result);
			return result;
        }
#endif



        virtual
        util::FunctorEnum
        GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
        {
	        DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
        	std::pair< std::string, int> result;
        	store::Vector3N< int> ids = PositionGrid< t_RETURN>::IDsFromPosition( POSITION);
			if( ids[0] < 0 || ids[0] >= (int) m_Data.size() || ids[1] < 0 || ids[1] >= (int) m_Data[0].size() || ids[2] < 0 || ids[2] >= (int) m_Data[0][0].size())
			{
				DebugWrite( __FUNCTION__<< ": return void information. moltype: " << MOL_TYPE << " pos: " << POSITION);
				return util::e_VoidGridPoint;
			}
        	t_RETURN value( m_Data[ ids[0]][ ids[1]][ ids[2]]);
        	util::FunctorEnum
				gp_type = value->GetClassID();
			if( gp_type == util::e_TypeMappedGridPoint)
			{
				return util::e_ConstantGridPoint;
			}
			else if( gp_type == util::e_RecursiveTypeMappedGridPoint)
			{
				RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) value.get();
				DebugWrite( "nsp from recursive: " << cptr->GetFromTypeMap( MOL_TYPE)->GetClassID());
				return cptr->GetFromTypeMap( MOL_TYPE)->GetClassID();
			}
        	return gp_type;
        }


        virtual
        std::pair< util::FunctorEnum, int>
        GetGridPointTypeAndNSP( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
        {
	        DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
        	std::pair< util::FunctorEnum, int> result;
        	store::Vector3N< int> ids = PositionGrid< t_RETURN>::IDsFromPosition( POSITION);
			if( ids[0] < 0 || ids[0] >= (int) m_Data.size() || ids[1] < 0 || ids[1] >= (int) m_Data[0].size() || ids[2] < 0 || ids[2] >= (int) m_Data[0][0].size())
			{
				DebugWrite( __FUNCTION__<< ": return void information. moltype: " << MOL_TYPE << " pos: " << POSITION);
				return std::make_pair( util::e_VoidGridPoint, g_IntMax);
			}
        	t_RETURN value( m_Data[ ids[0]][ ids[1]][ ids[2]]);

			if( value->GetClassID() == util::e_TypeMappedGridPoint)
			{
				TypeMappedGridPoint<  mol::Atom, math::Vector3N> * cptr = ( TypeMappedGridPoint<  mol::Atom, math::Vector3N> *) value.get();
				DebugWrite( "nsp from typemapped: " << cptr->GetFromTypeMap( MOL_TYPE).GetNSP());
				result.first = util::e_ConstantGridPoint;
				result.second =  cptr->GetFromTypeMap( MOL_TYPE).GetNSP();
			}
			else if( value->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
			{
				RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) value.get();
				DebugWrite( "nsp from recursive: " << cptr->GetFromTypeMap( MOL_TYPE)->GetNSP());
				result.first = cptr->GetFromTypeMap( MOL_TYPE)->GetClassID();
				result.second =  cptr->GetFromTypeMap( MOL_TYPE)->GetNSP();
			}
			else
			{
				result.first = value->GetClassID();
				result.second = value->GetNSP();
			}
        	return result;
        }


        virtual void SetGridPoint( const math::Vector3N &POSITION, const t_RETURN &DATA)
        {
	  //	    std::cout << __FUNCTION__ << std::endl;
            store::Vector3N< int> indices( PositionGrid< t_RETURN>::IDsFromPosition( POSITION));  /// todo: return safety if out of limits!

//	    std::cout << m_Data.size() << "  ";
//	    std::cout << m_Data[indices[0]].size() << "  ";
//	    std::cout << m_Data[indices[0]][indices[1]].size() << std::endl;
            m_Data[ indices[0]][ indices[1]][ indices[2]] =  DATA;
	    //	    std::cout << __FUNCTION__ << "  " << DATA << std::endl;
        }

        virtual void SetGridPoint( const store::Vector3N< int> &INDICES, const t_RETURN &DATA)
        {
	  //	    std::cout << __FUNCTION__ << std::endl;
//	    std::cout << DATA << std::endl;
            m_Data[ INDICES[0]][ INDICES[1]][ INDICES[2]] =  DATA;
        }

        virtual void SetGridPoint( const int &I, const int &J, const int &K, const t_RETURN &DATA)
        {
	  //	    std::cout << __FUNCTION__ << std::endl;
            m_Data[ I][ J][ K] = DATA;
        }

        virtual void SetGridPoint(  t_RETURN *PTR_TO_GRIDPOINT, const t_RETURN &DATA_FOR_GRIDPOINT)
        { // use with iterator: &*itr
            *PTR_TO_GRIDPOINT = DATA_FOR_GRIDPOINT;  // not exactly safe
        }

        virtual t_RETURN &GetSetGridPoint(  const store::Vector3N< int> &INDICES)
        {
        	return m_Data[ INDICES[0]][ INDICES[1]][ INDICES[2]];
        }


        virtual const math::Vector3N &GetMax() const{ return PositionGrid< t_RETURN>::m_Max;}

        virtual const math::Vector3N &GetMin() const{ return PositionGrid< t_RETURN>::m_Min;}

        virtual const math::Vector3N &GetDelta() const{ return PositionGrid< t_RETURN>::m_Delta;}

        virtual void SetMax( const math::Vector3N &MAX)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            PositionGrid< t_RETURN>::m_Max = MAX;
        }

        virtual void SetMin( const math::Vector3N &MIN)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            PositionGrid< t_RETURN>::m_Min = MIN;
        }

        virtual void SetDelta( const math::Vector3N &DELTA)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            PositionGrid< t_RETURN>::m_Delta = DELTA;
        }


        virtual void InitializeDataStructure()
        {
            store::Vector3N< size_t> dims( PositionGrid< t_RETURN>::GetNrBins());
            m_Data = std::vector< std::vector< std::vector< t_RETURN> > >( dims( 0), std::vector< std::vector< t_RETURN> >( dims( 1), std::vector< t_RETURN>( dims( 2))));
        }


        std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << std::endl;
            size_t i = 0;
            PositionGrid< t_RETURN>::Write( STREAM);
            for( typename std::vector< std::vector< std::vector< t_RETURN> > >::const_iterator x_itr = m_Data.begin(); x_itr != m_Data.end(); ++x_itr, ++i)
            {
                size_t j = 0;
                for( typename std::vector< std::vector< t_RETURN> >::const_iterator y_itr = x_itr->begin(); y_itr != x_itr->end(); ++y_itr, ++j)
                {
                    size_t k = 0;
                    for( typename std::vector< t_RETURN>::const_iterator z_itr = y_itr->begin(); z_itr != y_itr->end(); ++z_itr, ++k)
                    {
#ifdef DEBUG_XXL
                        STREAM << i << "  " << j << "  " << k << "  ";
                        STREAM.flush();
                        math::Vector3N pos( PositionGrid< t_RETURN>::PositionFromIDs( store::Vector3N< size_t>( i, j, k)));
                        STREAM << pos( 0) << "  " << pos( 1) << "  " << pos( 2) << "  ";
                        STREAM.flush();
#endif
                        STREAM << *z_itr;
                    }
                }
            }
            return STREAM;
        }





        std::istream &Read( std::istream &STREAM) // todo: use iterators!!!
        {
        	DebugWrite( __PRETTY_FUNCTION__);
        	std::string str;
        	PositionGrid< t_RETURN>::Read( STREAM);



    	math::Vector3N
			min = GetMin(),
			max = GetMax(),
			delta = GetDelta();

	// predefinition of void grid point for speed reasons // all void gridpoint actually link to same object
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
			void_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N >());
	// same as void grid point but conserving surface identity
        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
			surf_grid_point( new store::SurfGridPoint< mol::Atom, math::Vector3N >());

    	// Memory();

        InitializeDataStructure();

        float
            x_start( min( 0) + 0.5 * delta( 0)),
            x( x_start),
            y_start( min( 1) + 0.5 * delta( 1)),
            y( y_start),
            z_start( min( 2) + 0.5 * delta( 2)),
            z( z_start);

        int nr_layers = math::Round< int>( fabs( x_start - max( 0)) / delta( 0));

    	std::cout << "read gridpoints" << std::endl;

        int
			i = 0,
			j, k;


    	// Memory();

			std::cout << "layers remaining in reading: ";

            for( ; x < max( 0); x += delta( 0), ++i)
            {
            	std::cout << nr_layers - i << " ";
		std::cout.flush();
                for( y = y_start, j = 0; y < max( 1); y += delta( 1), ++j)
                    for( z = z_start, k = 0; z < max( 2); z += delta( 2), ++k)
                    {
		      //		      Memory( "monitor.txt", z +  1000 * y + 1e6 * x);
                        math::Vector3N grid_position( x, y, z);

                        STREAM >> str;

                        util::FunctorEnum
							type = util::EnumHandler< util::FunctorEnum>().ID( str);

			//			std::cout << x << " : " << y << " : " << z << std::endl;
                        if( type == util::e_VoidGridPoint)
                        {

			  //			  std::cout << "void grid point in " << __PRETTY_FUNCTION__ << std::endl;
                            SetGridPoint
                            (
                                    grid_position,
                                    void_grid_point
                            );
                            // Memory( "memory.txt", 1);
                        }
                        else if( type == util::e_SurfGridPoint)
                        {

			  //			  std::cout << "surf grid point in " << __PRETTY_FUNCTION__ << std::endl;
                            SetGridPoint
                            (
                                    grid_position,
                                    surf_grid_point
                            );
                            // Memory( "memory.txt", 1);
                        }
                        else if( type == util::e_ConstantGridPoint)
                        {
                        	store::ConstantGridPoint< mol::Atom, math::Vector3N > grid_point;
				//			  std::cout << "const grid point in " << __PRETTY_FUNCTION__ << std::endl;
                        	STREAM >> grid_point;
                            SetGridPoint
                            (
                            		grid_position,
                            		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::ConstantGridPoint< mol::Atom, math::Vector3N >( grid_point))
                            );
                            // Memory( "memory.txt", 2);
                        }
                        else if( type == util::e_InteractionGridPoint)
                        {
			  //			  std::cout << "interaction grid point in " << __PRETTY_FUNCTION__ << std::endl;
#ifndef NO_INTERACTION_FORCES
                            store::InteractionGridPoint grid_point;
                            STREAM >> grid_point;
                            SetGridPoint
                            (
                            		grid_position,
                            		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::InteractionGridPoint( grid_point))
                            );
                            // Memory( "memory.txt", 3);
#else
                            SetGridPoint
                            (
                                    grid_position,
                                    void_grid_point
                            );
                            // Memory( "memory.txt", 1);
#endif
                        }
                        else if( type == util::e_TypeMappedGridPoint)
                        {
                        	store::TypeMappedGridPoint< mol::Atom, math::Vector3N > grid_point;
                            STREAM >> grid_point;
                            SetGridPoint
                            (
                            		grid_position,
                            		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::TypeMappedGridPoint< mol::Atom, math::Vector3N >( grid_point))
                            );
                            // Memory( "memory.txt", 4);
                        }
                        else if( type == util::e_RecursiveTypeMappedGridPoint)
                        {
                            store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > grid_point;
                            STREAM >> grid_point;
                            SetGridPoint
                            (
                            		grid_position,
                            		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >( grid_point))
                            );
                            // Memory( "memory.txt", 5);
                        }
                        else
                        {
                        	std::cout << "===> undefined grid point type in " << __PRETTY_FUNCTION__ << ": <" << str << "> type: <" << type << ">" << std::endl;
                        }
			if( !STREAM)
			{
			    std::cout << "====> stream closed after reading gridpoint" << std::endl;
			    exit(-1);
			}
			//	    	Memory();
			//			std::cout << std::endl;
                    }
            }
            std::cout << std::endl;


        	return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }







    }; // end class GridVector3D


} // end namespace store

#endif
