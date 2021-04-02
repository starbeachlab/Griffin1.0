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


#ifndef INTERPOLGRID_T_H
#define INTERPOLGRID_T_H

#include "position_grid.t.h"
#include "../external/boost_functions.h"
#include "../math/function_functions.t.h"
#include <cmath>
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"
#include "atom_grid_factory.t.h"
#include "void_gridpoint.t.h"
#include "../molecule/atom.h"

#include "sum_gridpoint.t.h"


namespace store
{
	typedef boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >  td_FUNC;

	class InterpolGrid
	: public PositionGrid< td_FUNC>
	{
	protected:
		boost::shared_ptr< PositionGrid< td_FUNC> >     m_Grid;  // TODO: try without member !!!!!!
		float                                           m_Energy;

	public:

		InterpolGrid()
		: PositionGrid< td_FUNC>(),
		m_Grid(),
		m_Energy( std::numeric_limits< float>::max())
		{
	    	DebugWrite( __PRETTY_FUNCTION__);
		}


		InterpolGrid( const InterpolGrid &GRID)
		: PositionGrid< td_FUNC>(),
		m_Grid( GRID.m_Grid),
		m_Energy( std::numeric_limits< float>::max())
		{
			StandardWrite( __PRETTY_FUNCTION__ << " - copy constructor");
			Memory();
		}


		InterpolGrid( const boost::shared_ptr< PositionGrid< td_FUNC> >&GRID)
		: PositionGrid< td_FUNC>(),
		m_Grid( GRID),
		m_Energy( std::numeric_limits< float>::max())
		{
			DebugWrite( __PRETTY_FUNCTION__ << " grid: " << GRID->GetClassName());
			Memory();
		}


		virtual ~InterpolGrid(){}


		virtual InterpolGrid *Clone() const
		{
			return new InterpolGrid( *this);
		}


		virtual void SetGridPoint( const math::Vector3N &POSITION, const td_FUNC &DATA)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetGridPoint( POSITION, DATA);
		}

		virtual void SetGridPoint( const store::Vector3N< int> &INDICES, const td_FUNC &DATA)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetGridPoint( INDICES, DATA);
		}

		virtual void SetGridPoint( const int &I, const int &J, const int &K, const td_FUNC &DATA)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetGridPoint( I, J, K, DATA);
		}

        virtual
        util::FunctorEnum
        GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
        {
        	return m_Grid->GetGridPointType( POSITION, MOL_TYPE);
        }

		virtual const math::Vector3N &GetMax() const{ return m_Grid->GetMax();}

		virtual const math::Vector3N &GetMin() const{ return m_Grid->GetMin();}

		virtual const math::Vector3N &GetDelta() const{ return m_Grid->GetDelta();}

		virtual void SetMax( const math::Vector3N &MAX)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetMax( MAX);
		}

		virtual void SetMin( const math::Vector3N &MIN)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetMin( MIN);
		}

		virtual void SetDelta( const math::Vector3N &DELTA)
		{
			DebugWrite( __PRETTY_FUNCTION__);
			m_Grid->SetDelta( DELTA);
		}

		virtual std::ostream &Write( std::ostream &STREAM) const
		   {
//			   STREAM << GetClassName()<< std::endl;
			   m_Grid->Write( STREAM);
			   return STREAM;
		   }


		virtual std::istream &Read( std::istream &STREAM)
		{
			return STREAM;
		}


		virtual std::string GetClassName() const
		{
			return mystr::GetClassName( __PRETTY_FUNCTION__);
		}


		virtual
		td_FUNC
		GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE) const
		{
			DebugWrite( "* " << __PRETTY_FUNCTION__ << " pos: " << POSITION[0] << "  " << POSITION[1] << "   " << POSITION[2] << " moltype: " << MOL_TYPE);

			// if position is outside of grid limits return neutral object
			for( size_t i = 0; i < 3; ++i)  // TODO use iterators instead of indices
			{
				if
				(
						POSITION[ i] < m_Grid->GetMin()[ i] + 0.5 * m_Grid->GetDelta()[ i]
						|| POSITION[ i] > m_Grid->GetMax()[ i] - 0.5 * m_Grid->GetDelta()[ i]
				)
				{
					DebugWrite( "==> interpolation grid applied to points outside grid points (outer frame)" << " pos: " << POSITION[0] << "  " << POSITION[1] << "   " << POSITION[2] << " moltype: " << MOL_TYPE);
					return td_FUNC( new VoidGridPoint< mol::Atom, math::Vector3N >());
				}
			}

			td_FUNC
				center_gridpoint = m_Grid->GetGridPoint( POSITION, MOL_TYPE);

			util::FunctorEnum
				gridpoint_type = center_gridpoint->GetClassID();

			int
				nsp = center_gridpoint->GetNSP();

			// determine indices of zero point (front lower left corner)
			store::Vector3N< size_t>
				zero_point_indices = m_Grid->NeighborWithLowestIndices( POSITION);


			// translate coordinates into unit cell
			math::Vector3N
				relative_position = POSITION - m_Grid->PositionFromIDs( zero_point_indices);

			relative_position /= m_Grid->GetDelta();

			DebugWrite( "position: " << POSITION << std::endl << "relative position:  " << relative_position);

			boost::shared_ptr< SumGridPoint< mol::Atom, math::Vector3N> >   // TODO: TEMPLATE!
				sum_gridpoint( new SumGridPoint< mol::Atom, math::Vector3N>( gridpoint_type, nsp));  // SumGP has type and nsp of center
#ifdef DEBUG
			float sum = 0;
#endif
			float
				factor;

			//		DebugWrite( m_Grid->GetClassName() << " min: " << m_Grid->GetMin() << " max: " << m_Grid->GetMax() << " delta: " << m_Grid->GetDelta());


			// weighted sum of grid points
			for( size_t i = 0; i < 2; ++i)
				for( size_t j = 0; j < 2; ++j)
					for( size_t k = 0; k < 2; ++k)
					{
						factor = std::fabs( ( 1 - i - relative_position( 0)) * ( 1 - j - relative_position( 1)) * ( 1 - k - relative_position( 2)));
						// todo: use moltype dependent get gp function instead?:
						td_FUNC gp = m_Grid->GetGridPoint( store::Vector3N< int>( zero_point_indices( 0) + i, zero_point_indices( 1) + j, zero_point_indices( 2) + k));

						if( gp->GetClassID() == util::e_TypeMappedGridPoint)
						{
						    TypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( TypeMappedGridPoint< mol::Atom, math::Vector3N> *) gp.get();

						    sum_gridpoint->Insert
						    (
						    		factor,
						    		td_FUNC( new ConstantGridPoint< mol::Atom, math::Vector3N>( cptr->GetFromTypeMap( MOL_TYPE)))
						    );
						}
						else if( gp->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
						{
						    
						    RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) gp.get();

						    sum_gridpoint->Insert
						    (
						    		factor,
						    		cptr->GetFromTypeMap( MOL_TYPE) // check for type to omit void??
						    );
						}
						else if( gp->GetClassID() != util::e_VoidGridPoint && gp->GetClassID() != util::e_SurfGridPoint)
						{
						    sum_gridpoint->Insert
						    (
						    		factor,
						    		gp
						    );
						}
#ifdef DEBUG
						else{ std::cout << "=> void gridpoint is not inserted to sumgridpoint" << std::endl;}

						sum += factor;
						std::cout << "factor: " << factor;
						std::cout << "    sum: " << sum << std::endl;
					}
			std::cout << "passed sum_gridpoint: " << sum_gridpoint << std::endl;
#else
		            }
#endif
			return sum_gridpoint;
		}
















#ifdef FORCE_INDICES
		virtual
		td_FUNC
		GetGridPoint
		( 
		    const math::Vector3N &POSITION, 
		    const std::string &MOL_TYPE
#ifdef CPPIO
		    ,std::ostream &STREAM
#else
		    , FILE *STREAM
#endif
		    ) const
		{
			DebugWrite( "* " << __PRETTY_FUNCTION__ << " pos: " << POSITION[0] << "  " << POSITION[1] << "   " << POSITION[2] << " moltype: " << MOL_TYPE);

			// if position is outside of grid limits return neutral object
			for( size_t i = 0; i < 3; ++i)  // TODO use iterators instead of indices
			{
				if
				(
						POSITION[ i] < m_Grid->GetMin()[ i] + 0.5 * m_Grid->GetDelta()[ i]
						|| POSITION[ i] > m_Grid->GetMax()[ i] - 0.5 * m_Grid->GetDelta()[ i]
				)
				{
					DebugWrite( "==> interpolation grid applied to points outside grid points (outer frame)" << " pos: " << POSITION[0] << "  " << POSITION[1] << "   " << POSITION[2] << " moltype: " << MOL_TYPE);
					return td_FUNC( new VoidGridPoint< mol::Atom, math::Vector3N >());
				}
			}

			td_FUNC
				center_gridpoint = m_Grid->GetGridPoint( POSITION, MOL_TYPE, STREAM);

			util::FunctorEnum
				gridpoint_type = center_gridpoint->GetClassID();

			int
				nsp = center_gridpoint->GetNSP();

			// determine indices of zero point (front lower left corner)
			store::Vector3N< size_t>
				zero_point_indices = m_Grid->NeighborWithLowestIndices( POSITION);


			// translate coordinates into unit cell
			math::Vector3N
				relative_position = POSITION - m_Grid->PositionFromIDs( zero_point_indices);

			relative_position /= m_Grid->GetDelta();

			DebugWrite( "position: " << POSITION << std::endl << "relative position:  " << relative_position);

			boost::shared_ptr< SumGridPoint< mol::Atom, math::Vector3N> >   // TODO: TEMPLATE!
				sum_gridpoint( new SumGridPoint< mol::Atom, math::Vector3N>( gridpoint_type, nsp));  // SumGP has type and nsp of center
#ifdef DEBUG
			float sum = 0;
#endif
			float
				factor;

			//		DebugWrite( m_Grid->GetClassName() << " min: " << m_Grid->GetMin() << " max: " << m_Grid->GetMax() << " delta: " << m_Grid->GetDelta());


			// weighted sum of grid points
			for( size_t i = 0; i < 2; ++i)
				for( size_t j = 0; j < 2; ++j)
					for( size_t k = 0; k < 2; ++k)
					{
						factor = std::fabs( ( 1 - i - relative_position( 0)) * ( 1 - j - relative_position( 1)) * ( 1 - k - relative_position( 2)));
						// todo: use moltype dependent get gp function instead?:
						td_FUNC gp = m_Grid->GetGridPoint( store::Vector3N< int>( zero_point_indices( 0) + i, zero_point_indices( 1) + j, zero_point_indices( 2) + k));

						if( gp->GetClassID() == util::e_TypeMappedGridPoint)
						{
						    TypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( TypeMappedGridPoint< mol::Atom, math::Vector3N> *) gp.get();

						    sum_gridpoint->Insert
						    (
						    		factor,
						    		td_FUNC( new ConstantGridPoint< mol::Atom, math::Vector3N>( cptr->GetFromTypeMap( MOL_TYPE)))
						    );
						}
						else if( gp->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
						{
						    
						    RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> * cptr = ( RecursiveTypeMappedGridPoint<  mol::Atom, math::Vector3N> *) gp.get();

						    sum_gridpoint->Insert
						    (
						    		factor,
						    		cptr->GetFromTypeMap( MOL_TYPE) // check for type to omit void??
						    );
						}
						else if( gp->GetClassID() != util::e_VoidGridPoint && gp->GetClassID() != util::e_SurfGridPoint)
						{
						    sum_gridpoint->Insert
						    (
						    		factor,
						    		gp
						    );
						}
#ifdef DEBUG
						else{ std::cout << "=> void gridpoint is not inserted to sumgridpoint" << std::endl;}

						sum += factor;
						std::cout << "factor: " << factor;
						std::cout << "    sum: " << sum << std::endl;
					}
			std::cout << "passed sum_gridpoint: " << sum_gridpoint << std::endl;
#else
		            }
#endif
			return sum_gridpoint;
		}
#endif






		virtual td_FUNC GetGridPoint( const store::Vector3N< int> &INDICES) const
		{
			std::cout << "=====> " << __PRETTY_FUNCTION__ << " not implemented yet" << std::endl;
			exit( -1); // TODO implement
			return td_FUNC();
		}

		virtual td_FUNC GetGridPoint( const store::Vector3N< int> &INDICES, const std::string &MOL_TYPE) const
		{
			std::cout << "=====> " << __PRETTY_FUNCTION__ << " not implemented yet" << std::endl;
			exit( -1); // TODO implement
			return td_FUNC();
		}

		virtual td_FUNC GetGridPoint( const math::Vector3N &POS) const
		{
			std::cout << "===> " << __PRETTY_FUNCTION__ << " not implemented yet" << std::endl;
			return m_Grid->GetGridPoint( POS);
		}



	}; // end class InterpolGrid



} // end namespace store

#endif /* INTERPOLGRID_T_H */



















