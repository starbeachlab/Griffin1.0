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


#ifndef DISTRIBUTIONXD_T_H_
#define DISTRIBUTIONXD_T_H_

#include "tupleXD.t.h"
#include "iteration_vector.t.h"
#include "../external/std_functions.h"

namespace math
{
	template< typename t_DATA, int DIM>
	class DistributionXD
	: public TupleXD< t_DATA, DIM>
	{
	private:
		t_DATA    m_Outside;


	public:

		DistributionXD()
		: TupleXD< t_DATA, DIM>(),
		m_Outside()
		{}

		DistributionXD(  const float &DELTA, const std::vector< float> &MINIMUM, const std::vector< float> &MAXIMUM, const t_DATA &DEFAULT = 0)
		: TupleXD< t_DATA, DIM>( DELTA, MINIMUM, MAXIMUM, DEFAULT),
		m_Outside( t_DATA( 0))
		{}

		//! sum over one dimension
		DistributionXD( const DistributionXD< t_DATA, DIM+1> &HIGHER_DIM_DISTR, const unsigned int &CONTRACT_AXIS)
		: TupleXD< t_DATA, DIM>(),
		m_Outside( HIGHER_DIM_DISTR.GetOutsideValue())
		{
			std::cout << __PRETTY_FUNCTION__ << std::endl;

			std::vector< int>
				sizes = HIGHER_DIM_DISTR.GetNrElements(),
				indices,
				all_indices( DIM, 0);

			TupleXD< t_DATA, DIM>::m_Default = HIGHER_DIM_DISTR.GetDefault();

			int id;
			t_DATA sum;

			int total_size = 1;
			int cc = 0;
			for( unsigned int i = 0; i < sizes.size(); ++i)
				if( i != CONTRACT_AXIS)
				{
					std::cout << "insert layer " << i << std::endl;
					TupleXD< t_DATA, DIM>::m_Min[ cc] =  HIGHER_DIM_DISTR.GetMinimum()[i];
					TupleXD< t_DATA, DIM>::m_Delta[ cc] = HIGHER_DIM_DISTR.GetDelta()[i];
					TupleXD< t_DATA, DIM>::m_NrElements[ cc] = HIGHER_DIM_DISTR.GetNrElements()[i];
					total_size *= HIGHER_DIM_DISTR.GetNrElements()[ i];
					++cc;
				}

			TupleXD< t_DATA, DIM>::m_Data = std::vector< t_DATA>( total_size, t_DATA( 0));

			std::cout << "nr elements: ";
			std::copy( HIGHER_DIM_DISTR.GetNrElements().begin(), HIGHER_DIM_DISTR.GetNrElements().end(), std::ostream_iterator<int>( std::cout, " "));
			std::cout << std::endl;

			std::cout << "nr elements: ";
			std::copy( TupleXD< t_DATA, DIM>::m_NrElements.begin(), TupleXD< t_DATA, DIM>::m_NrElements.end(), std::ostream_iterator<int>( std::cout, " "));
			std::cout << std::endl;

			IterationVector< int> iterator( std::vector< int>( DIM, 0), TupleXD< t_DATA, DIM>::m_NrElements, std::vector< int>( DIM, 1));

			typename std::vector< t_DATA>::iterator data_itr;

			std::cout << "iterate: " << iterator.GetValues().size() << std::endl;

			do
			{
				sum = t_DATA( 0);
				indices = iterator.GetValues();
				id = TupleXD< t_DATA, DIM>::CalcID( indices);
//				std::cout << "iteration: " << indices[0] << " " << indices[1]  << "   " << id << std::endl;

				all_indices = indices;
				all_indices.insert( all_indices.begin() + CONTRACT_AXIS, t_DATA( 0));

//				std::copy( all_indices.begin(), all_indices.end(), std::ostream_iterator<int>( std::cout, " "));
//				std::cout << std::endl;

				data_itr = TupleXD< t_DATA, DIM>::m_Data.begin() + id;
				*data_itr = t_DATA( 0);

				for( int k = 0; k < sizes[ CONTRACT_AXIS]; ++k)
				{
					all_indices[ CONTRACT_AXIS] = k;
					*data_itr += HIGHER_DIM_DISTR( all_indices);
//					std::cout << "now: ";
//					std::copy( all_indices.begin(), all_indices.end(), std::ostream_iterator< int>( std::cout, " "));
//					std::cout << *data_itr << std::endl;
				}
//				std::cout << "sum: ";
//				std::copy( indices.begin(), indices.end(), std::ostream_iterator< int>( std::cout, " "));
//				std::cout << "   (" << id << "): " << *data_itr << std::endl;
			} while( iterator.Iterate());
		}


		const t_DATA &GetOutsideValue() const
		{
			return m_Outside;
		}

		void Insert( const std::vector< float> &POS, const t_DATA &VALUE = t_DATA( 1))
		{
//			std::cout << __FUNCTION__ << " ";
			assert( POS.size() == DIM);  // TODO: math/store::VectorXD< t_DATA, DIM> : public std::vector< t_DATA>
			if( !TupleXD< t_DATA, DIM>::IsInsideGrid( POS))
			{
//				std::cout << "outside" << std::endl;
				m_Outside += VALUE;
			}
			else
			{
//				std::cout <<  " inside " <<  TupleXD< t_DATA, DIM>::CalcID( TupleXD< t_DATA, DIM>::CalcIndices( POS)) << std::endl;
				TupleXD< t_DATA, DIM>::m_Data[ TupleXD< t_DATA, DIM>::CalcID( TupleXD< t_DATA, DIM>::CalcIndices( POS))] += VALUE;
			}
		}

        virtual
        std::ostream &
        Write( std::ostream& STREAM) const
        {
        	STREAM << "DistributionXD\noutside: " << m_Outside << std::endl;
        	return TupleXD< t_DATA, DIM>::Write( STREAM);
        }


        virtual
        std::ostream &
        WriteGnuplot( std::ostream& STREAM) const
        {
         	return TupleXD< t_DATA, DIM>::WriteGnuplot( STREAM);
        }

        virtual
        std::ostream &
        WriteAsArray( std::ostream& STREAM) const
        {
        	return TupleXD< t_DATA, DIM>::WriteAsArray( STREAM);
        }


	}; // end class DistributionXD

} // end namspace math


#endif /* DISTRIBUTIONXD_T_H_ */
