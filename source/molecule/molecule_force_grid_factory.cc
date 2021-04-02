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


#include "../../include/molecule/molecule_force_grid_factory.h"



namespace mol
{
    namespace factory
    {

		void
		RescaleConstantForces
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const float &VALUE
		)
		{

			float
				old_scale = GRID->GetSForceScale(),
				rescale = VALUE / old_scale;

			GRID->SetSForceScale( VALUE);

			StandardWrite( __FUNCTION__ << " old size: " << old_scale);
			StandardWrite( __FUNCTION__ << " new size: " << GRID->GetSForceScale());

			math::Vector3N
				min( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMin()),
				max( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMax()),
				delta( GRID->GetAtomForceGrid()->GetPositionGrid()->GetDelta()),
				tmp( ( max - min) / delta);

			store::Vector3N< int>
				max_index,
				indices;

			util::FunctorEnum
				grid_point_type;

			max_index( 0) = int( ( max( 0) - min( 0)) / delta( 0));
			max_index( 1) = int( ( max( 1) - min( 1)) / delta( 1));
			max_index( 2) = int( ( max( 2) - min( 2)) / delta( 2));

			boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				grid_point;

			bool
				is_grid1D = false;
			if( GRID->GetAtomForceGrid()->GetPositionGrid()->GetClassName() == store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >().GetClassName())
			{
			    std::cout << __FUNCTION__ << " grid1D" << std::endl;
				is_grid1D = true;
			}
 
			for( int i = 0; i < max_index( 0); ++i)
				for( int j = 0; j < max_index( 1); ++j)
					for( int k = 0; k < max_index( 2); ++k)
					{
						indices = store::Vector3N< int>( i, j, k);
						grid_point = GRID->GetAtomForceGrid()->GetGridPoint( indices);
						grid_point_type = grid_point->GetClassID();
						if( grid_point_type == util::e_ConstantGridPoint)
						{
							store::ConstantGridPoint< mol::Atom, math::Vector3N > *
								cptr = ( store::ConstantGridPoint< mol::Atom, math::Vector3N > *) grid_point.get();
							cptr->ReturnValue() *= rescale;
							cptr->Energy() *= rescale;
						}
						else if( grid_point_type == util::e_RecursiveTypeMappedGridPoint)
						{
						    store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > * 
							rptr = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *) grid_point.get();

						    for( store::Map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> > >::iterator itr = rptr->TypeMap().begin(); itr != rptr->TypeMap().end(); ++itr)
						    {
							if( itr->second->GetClassID() == util::e_ConstantGridPoint)
							{
							    store::ConstantGridPoint< mol::Atom, math::Vector3N > *
								cptr = ( store::ConstantGridPoint< mol::Atom, math::Vector3N > *) itr->second.get();
							    cptr->ReturnValue() *= rescale;
							    cptr->Energy() *= rescale;
							}
						    }
						}
						else if( grid_point_type == util::e_TypeMappedGridPoint)
						{
						    store::TypeMappedGridPoint< mol::Atom, math::Vector3N > * 
							tptr = ( store::TypeMappedGridPoint< mol::Atom, math::Vector3N > *) grid_point.get();

						    for( store::Map< std::string, store::ConstantGridPoint< mol::Atom, math::Vector3N> >::iterator itr = tptr->TypeMap().begin(); itr != tptr->TypeMap().end(); ++itr)
						    {
							    itr->second.ReturnValue() *= rescale;
							    itr->second.Energy() *= rescale;							
						    }
						}

						if( is_grid1D)
						{ 
//						    std::cout << "it is grid1D " << i << " " << j << " " << k << std::endl;
							store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > * 
							    grid_ptr = (store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> > > *) GRID->GetSetAtomForceGrid()->GetSetPositionGrid().get();
//							std::cout << "grid  type: " << grid_ptr->GetClassName() << std::endl;
//							std::cout << "grid point type: " << util::EnumHandler< util::FunctorEnum>().String( grid_point->GetClassID()) << std::endl;
							grid_ptr->UpdateGridPoint( indices, grid_point);
						}

					}
		}



		void
		RescaleForces
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const std::string &FORCE_TYPE,
//				const std::string &SCALING_STYLE,
				const float &VALUE
		)
		{
		    StandardWrite( __FUNCTION__ << " " << FORCE_TYPE << " value: " << VALUE);

			math::Vector3N
				min( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMin()),
				max( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMax()),
				delta( GRID->GetAtomForceGrid()->GetPositionGrid()->GetDelta()),
				tmp( ( max - min) / delta),
				before;

			store::Vector3N< int>
				max_index,
				indices;


			max_index( 0) = int( ( max( 0) - min( 0)) / delta( 0));
			max_index( 1) = int( ( max( 1) - min( 1)) / delta( 1));
			max_index( 2) = int( ( max( 2) - min( 2)) / delta( 2));

			boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				grid_point;

			boost::shared_ptr< store::VoidGridPoint< mol::Atom, math::Vector3N> >
			  void_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N>());

			util::FunctorEnum
				grid_point_type;

			for( int i = 0; i < max_index( 0); ++i)
				for( int j = 0; j < max_index( 1); ++j)
					for( int k = 0; k < max_index( 2); ++k)
					{
						indices = store::Vector3N< int>( i, j, k);
						grid_point = GRID->GetSetAtomForceGrid()->GetSetGridPoint( indices);
						grid_point_type = grid_point->GetClassID();
						if( grid_point_type == util::e_InteractionGridPoint)
						{
//						    if( VALUE < 1e-5)
//						    {
//						  grid_point = void_grid_point;
//						      GRID->GetSetAtomForceGrid()->GetSetPositionGrid()->SetGridPoint( indices, void_grid_point);
//						    }
//						    else
//						    {
							store::InteractionGridPoint *interaction_grid_point = ( store::InteractionGridPoint *) grid_point.get();
							//	phys::PotentialForceContainer potential_forces = interaction_grid_point->GetSetPotentialForcesContainer();
							for( std::map< std::string, boost::shared_ptr< phys::PotentialForce> >::iterator itr = interaction_grid_point->GetSetPotentialForcesContainer().begin(); itr != interaction_grid_point->GetSetPotentialForcesContainer().end(); ++itr)
							{
								if( itr->first == FORCE_TYPE)
								{
									DebugWrite( "found matching interaction gridpoint");
//									std::cout << "before: " << itr->second->GetSetVector()[0] << "  " << itr->second->GetSetVector()[1] << "  " << itr->second->GetSetVector()[2] << "  " << FORCE_TYPE << "  " << VALUE << std::endl;
									itr->second->GetSetVector() *= VALUE;
									itr->second->GetSetEnergy() *= VALUE;
//									std::cout << "after:  " << itr->second->GetSetVector()[0] << "  " << itr->second->GetSetVector()[1] << "  " << itr->second->GetSetVector()[2] << std::endl;
								}
							}
//						    }
						}
						else if( grid_point_type == util::e_RecursiveTypeMappedGridPoint)
						{
						    store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > * rptr = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *) grid_point.get();
						    for( store::Map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> > >::iterator rec_itr = rptr->TypeMap().begin(); rec_itr != rptr->TypeMap().end(); ++rec_itr)
						    {
							if( rec_itr->second->GetClassID() == util::e_InteractionGridPoint)
							{
//							    if( VALUE < 1e-5)
//							    {
//								rec_itr->second = void_grid_point;
//								StandardWrite( "=> check forces (" << __FILE__ << ":" << __LINE__ << ")");
//							    }
//							    else
//							    {
								store::InteractionGridPoint *interaction_grid_point = ( store::InteractionGridPoint *) rec_itr->second.get();
								for( std::map< std::string, boost::shared_ptr< phys::PotentialForce> >::iterator pot_itr = interaction_grid_point->GetSetPotentialForcesContainer().begin(); pot_itr != interaction_grid_point->GetSetPotentialForcesContainer().end(); ++pot_itr)
								{
								    if( pot_itr->first == FORCE_TYPE)
								    {
									pot_itr->second->GetSetVector() *= VALUE;
									pot_itr->second->GetSetEnergy() *= VALUE;
								    }
								}
//							    }
							}
						    }
						}
						else if( grid_point_type == util::e_TypeMappedGridPoint)
						{
							DebugWrite( "==> no interactions in type mapped gridpoints");
						}
					}
		}



		void
		BuildGridFromParallelFiles
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID,
				const CommandLineManager &CMD
		)
		{
		  std::vector< std::string> prefix_devisions = CMD.GetAllOptions( "combine_subgrids");

		  if( mystr::IsNumerical( prefix_devisions[0]) || !mystr::IsNumerical( prefix_devisions[1]))
		  {
		    std::cout << "ERROR in '-combine_subgrids  PREXIX  N' values were not passed correctly" << std::cout;
		    exit( -1);
		  }

		  int nr_devisions = mystr::ConvertStringToNumericalValue< int>( prefix_devisions[1]);
		  int nr_subgrids = math::Round< int>( pow( ( float) nr_devisions,3));
		  std::string
		    prefix = prefix_devisions[0] + mystr::NumericalValueToString( nr_devisions) + "_",
		    name;
		  		  
			std::ifstream read;
			std::vector< boost::shared_ptr< mol::MoleculeForceGrid> > sub_grids( nr_subgrids);
			std::vector< boost::shared_ptr< mol::MoleculeForceGrid> >::iterator grid_itr = sub_grids.begin();
			std::cout << __FUNCTION__ << " collects now: " << std::endl;
			int count = nr_subgrids;
			for( int i = 0; i < nr_subgrids; ++i, ++grid_itr)
			{
			  name = prefix + mystr::NumericalValueToString( i) + ".txt";
			  std::cout << name << " (remaining: " << --count << ")" << std::endl;
			  *grid_itr = boost::shared_ptr< mol::MoleculeForceGrid>( new mol::MoleculeForceGrid());
			  Open( read, name);
			  read >> *grid_itr;
			  Close( read);
			}
			GRID = MergeGrids( sub_grids, CMD);
		} // end BuildGridFromParallelFiles


		boost::shared_ptr< mol::MoleculeForceGrid>
		MergeGrids
		(
				const std::vector< boost::shared_ptr< mol::MoleculeForceGrid> > &SUB_GRIDS,
				const CommandLineManager &CMD
		)
		{

//	        boost::shared_ptr< mol::MoleculeForceGrid> grid( new mol::MoleculeForceGrid());

			std::vector< boost::shared_ptr< mol::MoleculeForceGrid> >::const_iterator grid_itr = SUB_GRIDS.begin();

			math::Vector3N
				total_min = ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMin(),
				total_max = ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMax(),
				delta = ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetDelta();

//			boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
//				vdw_epsilon_and_radius_map = ( *grid_itr)->GetVdwEpsilonAndRadiusMap();
//			boost::shared_ptr< store::Map< std::string, float> >
//				partial_charge_map = ( *grid_itr)->GetPartialChargeMap(),
//				mass_map = ( *grid_itr)->GetMassMap();
			std::string
				force_field_type = ( *grid_itr)->GetForceFieldType();
			geom::SurfGrid
				surf_grid = ( *grid_itr)->GetSurfGrid();
			math::VectorAngle
				angle_checker = ( *grid_itr)->GetVectorAngle();
			float
				sforce_scale = ( *grid_itr)->GetSForceScale();
			int 
			    count = SUB_GRIDS.size();

			std::cout << __FUNCTION__ << " consistency checks " << std::endl;

			// determine total limits of grid
			// assert that other members match
			for( ; grid_itr != SUB_GRIDS.end(); ++grid_itr)
			{
				total_min = math::Min( total_min, ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMin());
				total_max = math::Max( total_max, ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMax());
				if( delta != ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetDelta())
				{
					std::cout << "===> delta values don't match: " << delta << " vs: " << ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetDelta() << std::endl;
					exit( -1);
				}
//				if( *vdw_epsilon_and_radius_map != *( *grid_itr)->GetVdwEpsilonAndRadiusMap())
//				{
//					std::cout << "===> vdw epsilon and radius maps do not match!" << std::endl;
////					std::cout << vdw_epsilon_and_radius_map;
////					std::cout << ( *grid_itr)->GetVdwEpsilonAndRadiusMap();
////					exit( -1);
//				}
//				if( *partial_charge_map != *( *grid_itr)->GetPartialChargeMap())
//				{
//					std::cout << "===> partial charge maps don't match" << std::endl;
////					std::cout << partial_charge_map;
////					std::cout << ( *grid_itr)->GetPartialChargeMap();
////					exit( -1);
//				}
//				if( *mass_map != *( *grid_itr)->GetMassMap())
//				{
//					std::cout << "===> mass maps don't match" << std::endl;
////					std::cout << mass_map;
////					std::cout << ( *grid_itr)->GetMassMap();
////					exit( -1);
//				}
				if( angle_checker.GetCosThreshold() != ( *grid_itr)->GetVectorAngle().GetCosThreshold())
				{
					std::cout << "===> angle checker don't match" << std::endl;
					std::cout << angle_checker;
					std::cout << ( *grid_itr)->GetVectorAngle();
					exit( -1);
				}
				if( sforce_scale != ( *grid_itr)->GetSForceScale())
				{
					std::cout << "===> sforce scales don't match" << std::endl;
					std::cout << sforce_scale;
					std::cout << ( *grid_itr)->GetSForceScale();
					exit( -1);
				}
			}
//			grid->GetAtomForceGrid()->GetPositionGrid()->SetMax( total_max);
//			grid->GetAtomForceGrid()->GetPositionGrid()->SetMin( total_min);
//			grid->GetAtomForceGrid()->GetPositionGrid()->SetDelta( delta);

			store::Vector3N< size_t>
				nr_bins;

			nr_bins(0) = math::Round<size_t>( ( total_max( 0) - total_min( 0)) / delta( 0));
			nr_bins(1) = math::Round<size_t>( ( total_max( 1) - total_min( 1)) / delta( 1));
			nr_bins(2) = math::Round<size_t>( ( total_max( 2) - total_min( 2)) / delta( 2));
			
			StandardWrite( "total min: " << total_min);
			StandardWrite( "total max: " << total_max);
			StandardWrite( "delta: " << delta);
			StandardWrite( "nr bins: " << nr_bins[0] << " " << nr_bins[1] << " " << nr_bins[2]);

	        StandardWrite( "initialize position grid");
	        boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
				position_grid( new store::GridVector3D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( total_min, total_max, delta, nr_bins));



 //	        //  add interpolation layer
 //	        if( CMD.IsFlagSet( "interpolate"))
 //	        {
 //	            StandardWrite( "insert interpolation layer");
 //	            position_grid = boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
 //	            ( new store::InterpolGrid( position_grid));
 ////	            ( new store::InterpolGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >( position_grid));
 //	        }

 			// hook-up gridpoints
 			store::Vector3N< size_t>
 				max_index,
 				indices;

 			math::Vector3N
 				local_min,
 				local_max,
 				diff_min,
 				diff_max;

 			std::cout << "remaining sub grids: ";
			std::cout.flush();

 			for( grid_itr = SUB_GRIDS.begin(); grid_itr != SUB_GRIDS.end(); ++grid_itr)
 			{
			    std::cout<< count-- << " ";


 				local_min = ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMin();
				local_max = ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetMax();
				diff_min = local_min - total_min;
				diff_max = local_max - total_min;

//				std::cout << "min:   " << local_min << std::endl;
//				std::cout << "max:   " << local_max << std::endl;
//
//				std::cout << "total_min:   " << total_min << std::endl;
//				std::cout << "total_max:   " << total_max << std::endl;
//
//				std::cout << "diff_min:   " << diff_min << std::endl;
//				std::cout << "diff_max:   " << diff_max << std::endl;
//
//				std::cout << "delta:  " << delta << std::endl;
//
//				std::cout << "min_ids: "  << math::Round< size_t>( diff_min( 0) / delta( 0)) << " "  << math::Round< size_t>( diff_min( 1) / delta( 1)) << " "  << math::Round< size_t>( diff_min( 2) / delta( 2)) << std::endl;
//				std::cout << "max_ids: "  << math::Round< size_t>( diff_max( 0) / delta( 0)) << " "  << math::Round< size_t>( diff_max( 1) / delta( 1)) << " "  << math::Round< size_t>( diff_max( 2) / delta( 2)) << std::endl;
//
//				std::cout << "sub grid bins: " << ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetNrBins()[0] << "  " << ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetNrBins()[1] << "  " << ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->GetNrBins()[2] << "  " << std::endl;


				for( int i = (int) math::Round< size_t>( diff_min( 0) / delta( 0)), ii = 0; i < (int) math::Round< size_t>( diff_max( 0) / delta( 0)); ++i, ++ii)
					for( int j = (int) math::Round< size_t>( diff_min( 1) / delta( 1)), jj = 0; j < (int) math::Round< size_t>( diff_max( 1) / delta( 1)); ++j, ++jj)
						for( int k = (int) math::Round< size_t>( diff_min( 2) / delta( 2)), kk = 0; k < (int) math::Round< size_t>( diff_max( 2) / delta( 2)); ++k, ++kk)
						{
							indices = store::Vector3N< size_t>( i, j, k);
							position_grid->SetGridPoint( i, j, k, ( *grid_itr)->GetAtomForceGrid()->GetGridPoint( store::Vector3N< int>( ii, jj, kk)));
//							std::cout << "====>  " << i << "  " << j << "  " << k << "   " << ii << "  " << jj << "  " << kk << "   " << ( *grid_itr)->GetAtomForceGrid()->GetGridPoint( store::Vector3N< int>( ii, jj, kk)) << std::endl;
//							std::cout << "====>  " << i << "  " << j << "  " << k << "   " << ii << "  " << jj << "  " << kk << "   " << position_grid->GetGridPoint( store::Vector3N< int>( i, j, k)) << std::endl;
//							std::cout << "pos in total grid:  " << position_grid->PositionFromIDs( store::Vector3N< size_t>( i, j, k)) << std::endl;
//							std::cout << "pos in sub grid:    " << ( *grid_itr)->GetAtomForceGrid()->GetPositionGrid()->PositionFromIDs( store::Vector3N< size_t>( ii, jj, kk)) << std::endl;
						}
			}

			std::cout << std::endl;

			boost::shared_ptr< mol::MoleculeForceGrid>
			mol_grid
			(
					new mol::MoleculeForceGrid
					(
							boost::shared_ptr< AtomForceGrid>( new AtomForceGrid( position_grid)),
//							vdw_epsilon_and_radius_map,
//							partial_charge_map,
//							mass_map,
//							force_field_type,
							surf_grid,
							0.0,
							sforce_scale
					)
			);
			mol_grid->SetForceFieldType( force_field_type);
			mol_grid->SetVectorAngle( angle_checker);
			return mol_grid;
		}

		/*
		void
		ReplaceInteractionWithTransientGridPoints  // TODO: also replace within type mapped GPs
		(
				boost::shared_ptr< mol::MoleculeForceGrid> &GRID
		)
		{
			StandardWrite( __FUNCTION__ );

			math::Vector3N
				min( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMin()),
				max( GRID->GetAtomForceGrid()->GetPositionGrid()->GetMax()),
				delta( GRID->GetAtomForceGrid()->GetPositionGrid()->GetDelta()),
				tmp( ( max - min) / delta),
				before;

			store::Vector3N< int>
				max_index,
				indices;

			max_index( 0) = int( ( max( 0) - min( 0)) / delta( 0));
			max_index( 1) = int( ( max( 1) - min( 1)) / delta( 1));
			max_index( 2) = int( ( max( 2) - min( 2)) / delta( 2));

			boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				grid_point;

			std::string 
			    grid_point_type;

			for( int i = 0; i < max_index( 0); ++i)
				for( int j = 0; j < max_index( 1); ++j)
					for( int k = 0; k < max_index( 2); ++k)
					{
						indices = store::Vector3N< int>( i, j, k);
						grid_point = GRID->GetSetAtomForceGrid()->GetSetGridPoint( indices);
						grid_point_type = grid_point->GetClassID();
						if( grid_point->GetClassID() == store::InteractionGridPoint().GetClassID()) // todo: expand here to the type mapped GPs
						{
							grid_point = boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
							(
									new store::TransientInteractionGridPoint( *((store::InteractionGridPoint *)grid_point.get()))
							);
							assert( grid_point);
							GRID->GetSetAtomForceGrid()->GetSetPositionGrid()->SetGridPoint( i, j, k, grid_point);
//							DebugWrite( "now: " << GRID->GetSetAtomForceGrid()->GetSetGridPoint( indices)->GetClassID());
						}
						else if( grid_point_type == store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >().GetClassID())
						{
							std::cout << "====> replace with transient GP not written yet if moltype dependent surface objects are added" << std::endl;
							exit( -1);
						}
						else if( grid_point_type == store::TypeMappedGridPoint< mol::Atom, math::Vector3N >().GetClassID())
						{
							std::cout << "====> replace with transient GP not written yet if moltype dependent surface objects are added" << std::endl;
							exit( -1);
						}
					}
		}
*/

    } // end namespace factory
} // end namespace mol
