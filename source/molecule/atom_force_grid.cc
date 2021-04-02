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


#include "../../include/molecule/atom_force_grid.h"


namespace mol
{

    //! virtual copy constructor
    AtomForceGrid *AtomForceGrid::Clone() const{ return new AtomForceGrid( *this);}

    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    math::Vector3N & AtomForceGrid::operator()( const Atom &ATOM) const
    {
        DebugWrite( __PRETTY_FUNCTION__ << " (" << ATOM.GetPosition()( 0) << " " << ATOM.GetPosition()( 1) << " " << ATOM.GetPosition()( 2) << ")");
        return math::Function< mol::Atom, math::Vector3N>::s_Tmp = m_Grid->GetGridPoint( ATOM.GetPosition(), ATOM.GetOwningMolecule()->GetType())->operator()( ATOM);
    }

//        boost::shared_ptr< ForceAtom> operator[]( const Atom &ATOM) const
//        {
//            return boost::shared_ptr< ForceAtom>( new ForceAtom( ATOM, math::Vector3N( new math::Vector3N( operator()( ATOM)))));
//        }

    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////

    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
    {
    	DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
        return m_Grid->GetGridPoint( POSITION, MOL_TYPE);
    }


#ifdef FORCE_INDICES
    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::GetGridPoint
    ( 
	const math::Vector3N &POSITION, 
	const std::string &MOL_TYPE, 
#ifdef CPPIO
	std::ostream &STREAM
#else
	FILE *STREAM
#endif
    )
    {
    	DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
        return m_Grid->GetGridPoint( POSITION, MOL_TYPE, STREAM);
    }

#endif

    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::GetGridPoint( const store::Vector3N< int> &INDICES, const std::string &MOL_TYPE)
    {
    	DebugWrite( __PRETTY_FUNCTION__ << " indices: " << INDICES << " moltype: " << MOL_TYPE);
        return m_Grid->GetGridPoint( INDICES, MOL_TYPE);
    }

    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::GetGridPoint( const math::Vector3N &POSITION)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
        return m_Grid->GetGridPoint( POSITION);
    }

    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::GetGridPoint( const store::Vector3N< int> &INDICES)
    {
    	DebugWrite( __PRETTY_FUNCTION__ << " indices: " << INDICES);
        return m_Grid->GetGridPoint( INDICES);
    }

    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > &
    AtomForceGrid::GetSetGridPoint( const store::Vector3N< int> &INDICES)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
        return m_Grid->GetSetGridPoint( INDICES);
    }

//    std::pair< std::string, int>
//    AtomForceGrid::GetGridPointTypeAndNSP( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
//    {
//    	DebugWrite( __PRETTY_FUNCTION__);
//        return m_Grid->GetGridPointTypeAndNSP( POSITION, MOL_TYPE);
//    }

    util::FunctorEnum
    AtomForceGrid::GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
        return m_Grid->GetGridPointType( POSITION, MOL_TYPE);
    }

    const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
    &AtomForceGrid::CalculateGrid
    (
            boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID,
            const boost::shared_ptr< geom::SurfGrid> &SURF_GRID,
            const boost::shared_ptr< phys::ForceContainer> &FORCES,
            const float &MIN_FORCE_MAGNITUDE,
            const float &SFORCE_LENGTH = 1.0
    )
    {
//        SetGridLimits( *MOL, DELTA, CUTOFF);  // maybe split and don't redo in case of parallelization???

        math::Vector3N
            min( GRID->GetMin()),
            max( GRID->GetMax()),
            delta( GRID->GetDelta());

        DebugWrite( __FUNCTION__ << " min: " << GRID->GetMin() << " max: " << GRID->GetMax());

        boost::shared_ptr< store::VoidGridPoint< mol::Atom, math::Vector3N > >
            void_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N >());

        // iterate through grid
        float
            x_start( min( 0) + 0.5 * delta( 0)), x( x_start),
            y_start( min( 1) + 0.5 * delta( 1)), y( y_start),
            z_start( min( 2) + 0.5 * delta( 2)), z( z_start),
            max_cutoff = FORCES->GetMaximumCutoff();

        std::pair< float, int>
            dist_and_id;

        math::Vector3N
            grid_position;

        std::vector< int>
            surf_ids;

        short
            value;

        store::Map< std::string, store::Triplet< std::string, std::pair< float, int>, math::Vector3N> >
            map;

//        std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator
//            surf_itr( SURF_MAP->begin());

//        store::Map< float, boost::shared_ptr< geom::PointSurfaceObject> > object_distance_map( geom::factory::ObjectDistanceMap( surf_itr));

        int count = GRID->GetNrBins( )( 0);

        if( SURF_GRID->GetMoleculeTypes().size() == 1)
        {
            StandardWrite( "uniform grid for all molecule types");
            std::cout << "layers remaining to be calculated: ";
            for( ; x < max( 0); x += delta( 0))  // todo: use indices and avoid calculation of position to indices
            {
                std::cout << count-- << " ";
                std::cout . flush();
                for( y = y_start; y < max( 1); y += delta( 1))
                    for( z = z_start; z < max( 2); z += delta( 2))
                    {
                        grid_position = math::Vector3N( x, y, z);

                        dist_and_id = SURF_GRID->BasicClosestSurfPoint( grid_position, "all");

                        math::Vector3N closest_connection( SURF_GRID->CalcPositionFromID( dist_and_id.second) - grid_position);


                        math::Vector3N idpos = SURF_GRID->CalcPositionFromID( dist_and_id.second);

                        DebugWrite( x << " " << y << " " << z << " dist: " << dist_and_id.first << " id: " << dist_and_id.second << " pos: " << idpos( 0) << " " << idpos( 1) << " " << idpos( 2));

                        DebugWrite( closest_connection( 0) << " " << closest_connection( 1) << " " << closest_connection( 2));

                        DebugWrite( "value: " << ( *SURF_GRID)( grid_position));

                        if( dist_and_id.first < 0.001)
                        {
                        	dist_and_id.first = 0.0;
                        }


                        if( ( *SURF_GRID)( grid_position) < 1)
                        {
                        	DebugWrite( "sforce: " << closest_connection.Length());
                        	GRID->SetGridPoint
                            (
                                    grid_position,
                                    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                    (
                                            new store::ConstantGridPoint< mol::Atom, math::Vector3N >
                                            (
                                                    closest_connection.SetToLength( SFORCE_LENGTH),
                                                    dist_and_id.second,
                                                    dist_and_id.first * SFORCE_LENGTH
                                            )
                                    )
                            );
                        }
                        else if( closest_connection.Length() > max_cutoff)
                        {
                            DebugWrite( "out of limits");
                            GRID->SetGridPoint
                            (
                                    grid_position,
                                    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( void_grid_point)
                            );
                        }
                        else
                        {
                        	DebugWrite( "interaction range");
                            GRID->SetGridPoint
                            (
                                    grid_position,
                                    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                    (
                                            new store::InteractionGridPoint
                                            (
                                                    phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, dist_and_id.first),
                                                    dist_and_id.second
                                            )
                                    )
                            );
                        }
                    }
            }
            std::cout << std::endl;
            DebugWrite( "");
        }
        else
        {
        	StandardWrite( "grid depends on molecule type ");

            std::vector< std::string>
                mol_types = SURF_GRID->GetMoleculeTypes();

            std::cout << "mol-types: ";
            std::copy( mol_types.begin(), mol_types.end(), std::ostream_iterator< std::string>( std::cout, " "));
            std::cout << std::endl;

            store::TypeMappedGridPoint< mol::Atom, math::Vector3N >().SetMolTypes( mol_types);
            store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >().SetMolTypes( mol_types);

            std::string
                type;
            math::Vector3N
                closest_connection;

            for( ; x < max( 0); x += delta( 0))
            {
                std::cout << count-- <<  " ";
                std::cout . flush();
                for( y = y_start; y < max( 1); y += delta( 1))
                    for( z = z_start; z < max( 2); z += delta( 2))
                    {
                        map.clear();

                        grid_position = math::Vector3N( x, y, z);

                        value = ( *SURF_GRID)( grid_position);

                        dist_and_id = SURF_GRID->BasicClosestSurfPoint( grid_position, "all");
                        closest_connection = SURF_GRID->CalcPositionFromID( dist_and_id.second) - grid_position;


#ifdef DEBUG_XXL
                    	DebugWrite( "\n\npos: " << x << " " << y << " " << z);
                        math::Vector3N idpos = SURF_GRID->CalcPositionFromID( dist_and_id.second);
                        DebugWrite( x << " " << y << " " << z << " dist: " << dist_and_id.first << " id: " << dist_and_id.second << " pos: " << idpos( 0) << " " << idpos( 1) << " " << idpos( 2));
                        DebugWrite( closest_connection( 0) << " " << closest_connection( 1) << " " << closest_connection( 2));
                        DebugWrite( "value: " << ( *SURF_GRID)( grid_position));
#endif


                        if( dist_and_id.first < 0.001)
                        {
                        	dist_and_id.first = 0.0;
                        }

                        // checking which molecule type experiences which situation for the determination of the gridpoint type
                        if( value == 0)
                        {
                        	DebugWrite( "sforces for all mol types");
                            map.InsertNewKeyAndValue( "all", store::Triplet< std::string, std::pair< float, int>, math::Vector3N>( "inside", dist_and_id, closest_connection));
                            // collect connections for other mol-types
                            for( std::vector< std::string>::const_iterator mol_itr = mol_types.begin(); mol_itr != mol_types.end(); ++mol_itr)
								if( *mol_itr !=  "all")
								{
			                        dist_and_id = SURF_GRID->BasicClosestSurfPoint( grid_position, *mol_itr);
			                        closest_connection = SURF_GRID->CalcPositionFromID( dist_and_id.second) - grid_position;
			                        map.InsertNewKeyAndValue( *mol_itr, store::Triplet< std::string, std::pair< float, int>, math::Vector3N>( "inside", dist_and_id, closest_connection));
								}
                            // check whether all point to same surf point
                            bool same( true);
                            std::map< std::string, store::Triplet< std::string, std::pair< float, int>, math::Vector3N> >::const_iterator itr = map.begin();
                            int id = itr->second.Second().second;

                            while( ++itr != map.end())
                            {
                                if( id !=  itr->second.Second().second)
                                {
                                    DebugWrite( "different connectors");
                                    same = false;
                                    break;
                                }
                            }

                            // if all sforces point to same id use a ConstantGP
                            if( same)
                            {
                                DebugWrite( "same sforce for all mol types");
                                GRID->SetGridPoint
                                (
                                        grid_position,
                                        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                        (
                                                new store::ConstantGridPoint< mol::Atom, math::Vector3N >
                                                (
                                                        map.begin()->second.Third().SetToLength( SFORCE_LENGTH),
                                                        dist_and_id.second,
                                                        dist_and_id.first * SFORCE_LENGTH
                                                )
                                        )
                                );
                            }
                            // else use TypeMappedGP
                            else
                            {
                                DebugWrite( "insert type mapped grid point: same type of (constant-sforce) gridpoint but different values");
                                store::TypeMappedGridPoint< mol::Atom, math::Vector3N > type_mapped_grid_point;
                                for( std::map< std::string, store::Triplet< std::string, std::pair< float, int>, math::Vector3N> >::const_iterator m_itr = map.begin(); m_itr != map.end(); ++m_itr)
                                {
                                    type_mapped_grid_point.InsertNewKeyAndValues
                                    (
                                    		m_itr->first,
                                            math::Vector3N( m_itr->second.Third()).SetToLength( SFORCE_LENGTH),
                                            m_itr->second.Second().second,
                                            m_itr->second.Second().first * SFORCE_LENGTH
                                    );
                                }
                                GRID->SetGridPoint
                                (
                                        grid_position,
                                        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                        (
                                                new store::TypeMappedGridPoint< mol::Atom, math::Vector3N >
                                                (
                                                        type_mapped_grid_point
                                                )
                                        )
                                );
                            }
                        }
                        else if( value > 0)
                        {
                        	// the only unique situation => directly insert grid points and exit loop
                            if( closest_connection.Length() <= max_cutoff)
                            {
//                            	DebugWrite( "interaction range for all");
//                                map.InsertNewKeyAndValue( "all", store::Triplet< std::string, std::pair< float, int>, math::Vector3N>( "force", dist_and_id, closest_connection));
                                DebugWrite( "interaction grid point inserted, same for all mol types");
                                GRID->SetGridPoint
                                (
                                        grid_position,
                                        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                        (
                                                new store::InteractionGridPoint
                                                (
                                                        phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, dist_and_id.first),
                                                        dist_and_id.second
                                                )
                                        )
                                );
                            }
                            else
                            {
//                            	DebugWrite( "void for all");
//                                map.InsertNewKeyAndValue( "all", store::Triplet< std::string, std::pair< float, int>, math::Vector3N>( "void", dist_and_id, closest_connection));
                                DebugWrite( "void grid point inserted, same for all mol types");
                                 GRID->SetGridPoint
                                 (
                                         grid_position,
                                         void_grid_point
                                 );
                            }
//                            continue;
                        }
                        else // if( value < 0)
                        {
                        	store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > recursive_grid_point;
                            recursive_grid_point.InsertNewKeyAndValue
                            (
                                    "all",
                                    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                                    ( new store::InteractionGridPoint( phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, dist_and_id.first), dist_and_id.second))
                            );


                            surf_ids = SURF_GRID->KeyToIndices( value);

                            std::vector< std::string> sforce_types( surf_ids.size());
                            std::vector< std::string>::iterator s_itr = sforce_types.begin();

                            for( std::vector< int>::const_iterator itr = surf_ids.begin(); itr != surf_ids.end(); ++itr)
                            {
                            	*s_itr = mol_types[ *itr];
                            }

                            for( std::vector< std::string>::const_iterator mol_itr = mol_types.begin(); mol_itr != mol_types.end(); ++mol_itr)
                            {
                            	if( *mol_itr == "all")
                            	{

                            	}
                            	else if( std::find( sforce_types.begin(), sforce_types.end(), *mol_itr) != sforce_types.end())
                            	{
                            		DebugWrite( "sforce is mounted for: " << *mol_itr);
                                    dist_and_id = SURF_GRID->BasicClosestSurfPoint( grid_position, *mol_itr);
                            		DebugWrite( "NSP: " << dist_and_id.second << " dist: " << dist_and_id.first);
                                    closest_connection = SURF_GRID->CalcPositionFromID( dist_and_id.second) - grid_position;
                            		DebugWrite( "closest connection calculated for: " << *mol_itr);
                                    recursive_grid_point.InsertNewKeyAndValue
                                    (
                                    		*mol_itr,
                                    		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
											(
													new store::ConstantGridPoint< mol::Atom, math::Vector3N >
													(
		                                                    closest_connection.SetToLength( SFORCE_LENGTH),
															dist_and_id.second,
															dist_and_id.first * SFORCE_LENGTH
													)
											)
                                    );
                            		DebugWrite( "sforce mounted for: " << *mol_itr);
                            	}
                            	else
                            	{
                            		DebugWrite( "interaction force is mounted to recursive GP for: " << *mol_itr << " mount: " << util::EnumHandler< util::FunctorEnum>().String( recursive_grid_point.GetGridPoint( "all")->GetClassID()));
                            		// copy interaction grid point but recalculate NSP ???

//                            		if( recursive_grid_point.GetGridPoint( "all")->GetClassID() != store::InteractionGridPoint().GetClassID())
//                            		{
//                            			std::cout << "====> incorrect force mounted! should be store::InteractionGridPoint but is: " << recursive_grid_point.GetGridPoint( "all")->GetClassID();
//                            		}

                                    recursive_grid_point.InsertNewKeyAndValue
                                    (
                                    		*mol_itr,
                                    		recursive_grid_point.GetGridPoint( "all")
                                    );
                            	}
                            } // iterate mol types

                            GRID->SetGridPoint
                            (
                                    grid_position,
                            		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
									(
											new store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >( recursive_grid_point)
									)
                            );
                    		DebugWrite( "recursiveTMGP is mounted");

                        } // else value < 0
                    } // iterate grid
            } // iterate grid
            std::cout << std::endl;
        }  // else (multiple mol types)
        DebugWrite( __FUNCTION__ << " min: " << GRID->GetMin() << " max: " << GRID->GetMax());
        return GRID;
    }



    void
    AtomForceGrid::FillSurfaceLayer
    (
            const boost::shared_ptr< geom::SurfGrid> &SURF_GRID,
            const boost::shared_ptr< phys::ForceContainer> &FORCES,
            const float &MIN_FORCE_MAGNITUDE
    )
    {
    	StandardWrite( __FUNCTION__);

        math::Vector3N
            grid_position;

        store::Vector3N< int>
			indices,relative_indices,
			offset;

		std::vector< int>::iterator
			offset_itr = offset.begin(),
			indices_itr;

		std::vector< float>::const_iterator
			grid_itr = m_Grid->GetMin().begin(),
			surf_itr = SURF_GRID->GetMinimum().begin(),
			delta_itr = m_Grid->GetDelta().begin();

		for( ; offset_itr != offset.end(); ++offset_itr, ++surf_itr, ++grid_itr, ++delta_itr)
		{
			*offset_itr = math::Round< int>( ( *grid_itr - *surf_itr) / *delta_itr);
		}

		int
			id;

        std::cout << "offset between force grid and surf grid: " << offset( 0) << " " << offset( 1) << " " << offset( 2) << " (should be zero when not running in parallel)" << std::endl;

        if( SURF_GRID->GetSurfaceIndices().size() == 1)
        {
        	StandardWrite( "no molecule type dependence, i.e. single surf layer");
        	for( std::multimap< short, std::pair< short, short> >::const_iterator itr = SURF_GRID->GetSurfaceIndices()[0].begin(); itr != SURF_GRID->GetSurfaceIndices()[0].end(); ++itr)
        	{
        		indices(0) = itr->first;
        		indices(1) = itr->second.first;
        		indices(2) = itr->second.second;

        		id = SURF_GRID->CalcID( indices);

        		grid_position = SURF_GRID->CalcPositionFromIndices( indices);

				if( grid_position >= m_Grid->GetMin() && grid_position < m_Grid->GetMax())
				{
					m_Grid->SetGridPoint
					(
							grid_position,
							boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
							(
									new store::InteractionGridPoint
									(
											phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, 0.0),
											id
									)
							)
					);
				}
        	}
        }
        else
        {
        	StandardWrite( "molecule type dependent");
        	boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				grid_point;
        	int ii = 0;

        	for( std::vector< std::multimap< short, std::pair< short, short> > >::const_iterator type_itr = SURF_GRID->GetSurfaceIndices().begin(); type_itr != SURF_GRID->GetSurfaceIndices().end(); ++type_itr, ++ii)
        	{
				std::string
					type = SURF_GRID->GetMoleculeTypes()[ ii];

				StandardWrite( "fill for: " << type);
				for( std::multimap< short, std::pair< short, short> >::const_iterator itr = type_itr->begin(); itr != type_itr->end(); ++itr)
				{
					indices = store::Vector3N< int>( itr->first, itr->second.first, itr->second.second);

					id = SURF_GRID->CalcID( indices);

					grid_position = SURF_GRID->CalcPositionFromIndices( indices);

					if( grid_position >= m_Grid->GetMin() && grid_position < m_Grid->GetMax())
					{
					    // indices -= offset
//					    std::transform(indices.begin(), indices.end(), offset.begin(), relative_indices.begin(), std::minus<int>());
//					    std::cout << "surf point is within grid limits (non-trivial for parallel runs)" << std::endl;
					    grid_point = GetGridPoint( grid_position, type);

						if
						(
								grid_point->GetClassID() == util::e_ConstantGridPoint
								|| grid_point->GetClassID() == util::e_TypeMappedGridPoint
						)
						{
							m_Grid->SetGridPoint
							(
									grid_position,
									boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
									(
											new store::InteractionGridPoint
											(
												phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, 0.0),
												id
											)
									)
							);
						}
						else if( grid_point->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
						{
							store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *
								recursive_ptr = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *) grid_point.get();


							boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
								type_mapped_sforce = recursive_ptr->TypeMap()( type);

							if( type_mapped_sforce->GetClassID() == util::e_ConstantGridPoint)
							{
								recursive_ptr->TypeMap()( type) = boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
								(
										new store::InteractionGridPoint
										(
												phys::factory::BuildPotentialForces( FORCES, grid_position, MIN_FORCE_MAGNITUDE, 0.0),
												id
										)
								);

								DebugWrite( "recursive grid point should contain interaction grid point for " << type << ": " << grid_point);

								if( recursive_ptr->TypeMap()( type)->GetClassID() != util::e_InteractionGridPoint)
								{
									std::cout << "grid point was not modified correctly: " << util::EnumHandler< util::FunctorEnum>().String( recursive_ptr->TypeMap()( type)->GetClassID()) << std::endl;
									exit( 1);
								}

							}
						} // if( grid_point->GetClassID() == ...
					}  // if grid position is within grid limits (crucial for subgrids)
				} // TODO: check this // for surf points
        	} // for mol types
        } // if nr layers
    } // end FillSurfLayer



    void
    AtomForceGrid::EmptySurfaceLayer
    (
            const boost::shared_ptr< geom::SurfGrid> &SURF_GRID
    )
    {
    	StandardWrite( __FUNCTION__);

        math::Vector3N
            grid_position;

        store::Vector3N< int>
			indices,
	    relative_indices,
			offset;

        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
        	void_grid_point( new store::SurfGridPoint< mol::Atom, math::Vector3N >());

		std::vector< int>::iterator
			offset_itr = offset.begin(),
			indices_itr;

		std::vector< float>::const_iterator
			grid_itr = m_Grid->GetMin().begin(),
			surf_itr = SURF_GRID->GetMinimum().begin(),
			delta_itr = m_Grid->GetDelta().begin();

		for( ; offset_itr != offset.end(); ++offset_itr, ++surf_itr, ++grid_itr, ++delta_itr)
		{
			*offset_itr = math::Round< int>( ( *grid_itr - *surf_itr) / *delta_itr);
		}

		int
			id;

        std::cout << "offset between force grid and surf grid: " << offset( 0) << " " << offset( 1) << " " << offset( 2) << " (should be zero when not running in parallel)" << std::endl;

        if( SURF_GRID->GetSurfaceIndices().size() == 1)
        {
        	StandardWrite( "no molecule type dependence, i.e. single surf layer");
        	for( std::multimap< short, std::pair< short, short> >::const_iterator itr = SURF_GRID->GetSurfaceIndices()[0].begin(); itr != SURF_GRID->GetSurfaceIndices()[0].end(); ++itr)
        	{
        		indices(0) = itr->first;
        		indices(1) = itr->second.first;
        		indices(2) = itr->second.second;

        		id = SURF_GRID->CalcID( indices);

        		grid_position = SURF_GRID->CalcPositionFromIndices( indices);

				if( grid_position >= m_Grid->GetMin() && grid_position < m_Grid->GetMax())
				{
					m_Grid->SetGridPoint
					(
							grid_position,
							void_grid_point
					);
				}
        	}
        }
        else
        {
        	StandardWrite( "molecule type dependent");
        	boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				grid_point;
        	int ii = 0;

        	for( std::vector< std::multimap< short, std::pair< short, short> > >::const_iterator type_itr = SURF_GRID->GetSurfaceIndices().begin(); type_itr != SURF_GRID->GetSurfaceIndices().end(); ++type_itr, ++ii)
        	{
				std::string
					type = SURF_GRID->GetMoleculeTypes()[ ii];

				StandardWrite( "emptying surf layer for: " << type);
				for( std::multimap< short, std::pair< short, short> >::const_iterator itr = type_itr->begin(); itr != type_itr->end(); ++itr)
				{
					indices = store::Vector3N< int>( itr->first, itr->second.first, itr->second.second);

					id = SURF_GRID->CalcID( indices);

					grid_position = SURF_GRID->CalcPositionFromIndices( indices);
					if( grid_position >= m_Grid->GetMin() && grid_position < m_Grid->GetMax())
					{
					        // indices -= offset
					    std::cout << "set to void: surf point at pos: "; 
					    std::copy( grid_position.begin(), grid_position.end(), std::ostream_iterator<float>( std::cout, " "));
					    std::cout << std::endl;
					    //	        std::transform( indices.begin(), indices.end(), offset.begin(), relative_indices.begin(), std::minus<int>());
//					        std::cout << "surf point is within grid limits (non-trivial for parallel runs)" << std::endl;
						grid_point = GetGridPoint( grid_position, type);

						if
						(
								grid_point->GetClassID() == util::e_ConstantGridPoint
								|| grid_point->GetClassID() == util::e_TypeMappedGridPoint
						)
						{
							m_Grid->SetGridPoint( grid_position, void_grid_point);
						}
						else if( grid_point->GetClassID() == util::e_RecursiveTypeMappedGridPoint)
						{
						    std::cout << "recursive it is" << std::endl;
							store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *
								recursive_ptr = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >*) grid_point.get();

							boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
								type_mapped_sforce = recursive_ptr->TypeMap()( type);

							if( type_mapped_sforce->GetClassID() == util::e_ConstantGridPoint)
							{
								recursive_ptr->TypeMap()( type) = void_grid_point;

								DebugWrite( "recursive grid point should contain interaction grid point for " << type << ": " << grid_point);

								if( recursive_ptr->TypeMap()( type)->GetClassID() != util::e_VoidGridPoint)
								{
									std::cout << "grid point was not modified correctly: " << util::EnumHandler< util::FunctorEnum>().String( recursive_ptr->TypeMap()( type)->GetClassID()) << std::endl;
									exit( 1);
								}
							}
						} // if( grid_point->GetClassID() == ...
						std::cout << "done" << std::endl;
					}  // if grid position is within grid limits (crucial for subgrids)
				} // TODO: check this // for surf points
        	} // for mol types
        } // if nr layers
    } // end EmptySurfLayer




    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*

    const boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > >
    &AtomForceGrid::CalculateGrid
    (
            boost::shared_ptr< store::PositionGrid< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > > &GRID,
            const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &SURF_MAP,
            const boost::shared_ptr< phys::ForceContainer> &FORCES,
            const float &CUTOFF,
            const float &SFORCE_LENGTH
    )
    {
        math::Vector3N
            min( GRID->GetMin()),
            max( GRID->GetMax()),
            delta( GRID->GetDelta());

        // iterate through grid
        float
            x_start( min( 0) + 0.5 * delta( 0)), x( x_start),
            y_start( min( 1) + 0.5 * delta( 1)), y( y_start),
            z_start( min( 2) + 0.5 * delta( 2)), z( z_start),
            probe_radius( 1.4);

        std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator
            surf_itr( SURF_MAP->begin());

        std::map< std::string, std::vector< std::vector< size_t> > >
            neighbor_map; // lists for each object first the object then its neighbors

        for( ; surf_itr != SURF_MAP->end(); ++surf_itr)
        {
            neighbor_map.insert( std::make_pair( surf_itr->first, geom::factory::ListOverlappingObjects( surf_itr->second)));
        }

        surf_itr = SURF_MAP->begin();

        std::map< std::string, std::vector< std::vector< size_t> > >::const_iterator
            neigh_itr( neighbor_map.begin());

        if( SURF_MAP->size() == 1)
        {
            StandardWrite( "  ... uniform grid for all molecule types ...");
            for( ; x < max( 0); x += delta( 0))
                for( y = y_start; y < max( 1); y += delta( 1))
                    for( z = z_start; z < max( 2); z += delta( 2))
                    {
                        math::Vector3N
                            grid_position( new math::Vector3N( x, y, z));

                        GRID->SetGridPoint
                        (
                                grid_position,
                                CalculateGridPoint( grid_position, surf_itr->second, neigh_itr->second, FORCES, CUTOFF, SFORCE_LENGTH, probe_radius)
                        );
                    }
        }
        else
        {
            StandardWrite( "grid depends on molecule type");
            for( ; x < max( 0); x += delta( 0))
                for( y = y_start; y < max( 1); y += delta( 1))
                    for( z = z_start; z < max( 2); z += delta( 2))
                    {
                        math::Vector3N
                            grid_position( new math::Vector3N( x, y, z));

                        std::map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >
                            grid_point_map;

                        bool
                            same_type( true);

                        for( surf_itr = SURF_MAP->begin(), neigh_itr = neighbor_map.begin(); surf_itr != SURF_MAP->end(); ++surf_itr, ++neigh_itr)
                        {
                            grid_point_map.insert( std::make_pair( surf_itr->first, CalculateGridPoint( grid_position, surf_itr->second, neigh_itr->second, FORCES, CUTOFF, SFORCE_LENGTH, probe_radius)));
                            DebugWrite( "inserted: " << surf_itr->first << " " << grid_point_map[ surf_itr->first]);
                        }

                        // are all grid points of the same kind?
                        std::map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >::const_iterator
                            itr( grid_point_map.begin());

                        std::string
                            previous( itr->second->GetClassName());

                        while( itr != grid_point_map.end() && same_type)
                        {
                            if( itr->second->GetClassName() != previous)
                            {
                                DebugWrite( "different grid point types for different molecules");
                                same_type = false;
                            }
                            ++itr;
                        }

                        // if same check content
                        if( same_type)
                        {
                            DebugWrite( "same grid point type for all molecules");
                            bool
                                same_content( true);

                            // is content the same?
                            itr = grid_point_map.begin();

                            math::Vector3N
                                previous_vector( *itr->second->operator()( mol::Atom( new mol::Atom())));

                            if( itr->second->GetClassName() == store::ConstantGridPoint< mol::Atom, math::Vector3N >().GetClassName())
                            {
                                while( same_content && itr != grid_point_map.end())
                                {
                                    if( previous_vector != *itr->second->operator()( mol::Atom( new mol::Atom())))
                                    {
                                        DebugWrite( "different value in grid point");
                                        same_content = false;
                                    }
                                }
                            }

                            // if same content use simple grid point
                            if( same_content)
                            {
                                DebugWrite( "simple grid point");
                                GRID->SetGridPoint
                                (
                                        grid_position,
                                        grid_point_map.begin()->second
                                );
                            }
                            // if content differs use type mapped grid point
                            else
                            {
                                DebugWrite( "type mapped grid point");
                                store::TypeMappedGridPoint< mol::Atom, math::Vector3N >
                                    type_mapped_grid_point;

                                for( itr = grid_point_map.begin(); itr != grid_point_map.end(); ++itr)
                                {
                                    type_mapped_grid_point.InsertNewKeyAndValue( itr->first, itr->second->operator()( mol::Atom( new mol::Atom())));
                                }

                                GRID->SetGridPoint
                                (
                                        grid_position,
                                        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::TypeMappedGridPoint< mol::Atom, math::Vector3N >( type_mapped_grid_point))
                                );
                            }
                        }
                        // if not the same: recursive type mapped grid point
                        else
                        {
                            DebugWrite( "recursive type mapped grid point");
                            store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > recursive_grid_point;

                            for( itr = grid_point_map.begin(); itr != grid_point_map.end(); ++itr)
                            {
                                recursive_grid_point.InsertNewKeyAndValue( itr->first, itr->second);
                            }

                            GRID->SetGridPoint
                            (
                                    grid_position,
                                    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N >( recursive_grid_point))
                            );
                        }
                    }
        }
        return GRID;
    }




    boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
    AtomForceGrid::CalculateGridPoint
    (
            const math::Vector3N &POSITION,
            const store::ShPtrVec< geom::DirectedPointSurfaceObject> &SURF_VEC,
            const std::vector< std::vector< size_t> > &NEIGHBOR_LIST,
            const boost::shared_ptr< phys::ForceContainer> &FORCES,
            const float &CUTOFF,
            const float &SFORCE_LENGTH,
            const float PROBE_RADIUS
    )
    {

        float
            distance,
            smallest_distance( std::numeric_limits< float>::max());

//        bool
//            is_inside_object( false);

        std::vector< std::vector< size_t> >::const_iterator
            neighbor_list_itr( NEIGHBOR_LIST.begin());

        std::vector< size_t>::const_iterator
            neighbor_itr, second_itr;

        math::Vector3N
            surf( new math::Vector3N());

        boost::shared_ptr< store::VoidGridPoint< mol::Atom, math::Vector3N > >
            void_grid_point( new store::VoidGridPoint< mol::Atom, math::Vector3N >());


        for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator object_itr( SURF_VEC.begin()); object_itr != SURF_VEC.end(); ++object_itr, ++neighbor_list_itr)
        {
            distance = ( *object_itr)->Distance( *POSITION);
            if( distance <= 0.0) // point within object => sforce
            {
                DebugWrite ("sforce");
                *surf = ( *object_itr)->ProjectionOnSurface( *POSITION);
                bool
                    is_overlapped = false;

                // check whether surface projection is overlapped by neighboring objects or not
                for( neighbor_itr = neighbor_list_itr->begin(); neighbor_itr != neighbor_list_itr->end() && !is_overlapped; ++neighbor_itr)
                {
                    if( SURF_VEC( *neighbor_itr)->Distance( *surf) <= 0)
                    {
                        is_overlapped = true;
                    }
                }
                if( is_overlapped)
                {
                    DebugWrite( "find closest surf point");
                    // find closest surf point ins surf object
                    boost::shared_ptr< std::pair< boost::shared_ptr< geom::DirectedPointSurfaceObject>, boost::shared_ptr< geom::DirectedSurfacePoint> > >
                        closest( geom::factory::ClosestObjectAndPoint( SURF_VEC, surf));

                    return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                    (
                            new store::ConstantGridPoint< mol::Atom, math::Vector3N >
                            (
                                    // this version creates vectors of constant size (linear repulsive potential), parallel to normal vector
                                    math::Vector3N( new math::Vector3N( ( closest->second->GetPosition() - *POSITION).SetToLength( SFORCE_LENGTH)))
                            )
                    );
                }
                else
                {
                    DebugWrite( "point to surf");
                    return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                    (
                            new store::ConstantGridPoint< mol::Atom, math::Vector3N >
                            (
                                    // this version creates vectors of constant size (linear repulsive potential), parallel to normal vector
                                    math::Vector3N( new math::Vector3N( ( *surf - *POSITION).SetToLength( SFORCE_LENGTH)))
                            )
                    );
                }
//                is_inside_object = true;
            }
            else if( distance <= PROBE_RADIUS) // point potentially in smoothed surface area (=> potentially sforce)
            {
                DebugWrite( "interaction or smooth range");
                // iterate through neighbors to see what's going on
                math::Vector3N all;
                for( neighbor_itr = neighbor_list_itr->begin(); neighbor_itr != neighbor_list_itr->end(); ++neighbor_itr)
                {
                    all += SmoothedSurfaceRepulsion( *object_itr, SURF_VEC( *neighbor_itr), *POSITION, PROBE_RADIUS, SFORCE_LENGTH);
                }

                for( neighbor_itr = neighbor_list_itr->begin(); neighbor_itr != neighbor_list_itr->end(); ++neighbor_itr)
                    for( second_itr = neighbor_itr + 1; second_itr != neighbor_list_itr->end(); ++second_itr)
                    {
                        all += SmoothedSurfaceRepulsion( SURF_VEC( *second_itr), SURF_VEC( *neighbor_itr), *POSITION, PROBE_RADIUS, SFORCE_LENGTH);
                    }

                if( all.IsLengthLargerZero())
                {
                    DebugWrite( "smooth range");
                    return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
                    (
                            new store::ConstantGridPoint< mol::Atom, math::Vector3N >
                            (
                                    math::Vector3N( new math::Vector3N( all.SetToLength( SFORCE_LENGTH)))
                            )
                    );
                }
                smallest_distance = distance;
            }
            else if( distance < smallest_distance)
            {
                smallest_distance = distance;
            }
        }

//        if( !is_inside_object)
//        {
            if( smallest_distance < CUTOFF) // interaction range
            {
                DebugWrite( "interaction range");
                return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( new store::InteractionGridPoint( phys::factory::BuildPotentialForces( FORCES, POSITION, MIN_FORCE_MAGNITUDE, dist_and_id.first)));
            }
//            else // void range
//            {
                DebugWrite( "void range");
                return boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >( void_grid_point);
//            }
//        }
    }

*/



    math::Vector3N
    AtomForceGrid::SmoothedSurfaceRepulsion
    (
            const boost::shared_ptr< geom::DirectedPointSurfaceObject> &FIRST,
            const boost::shared_ptr< geom::DirectedPointSurfaceObject> &SECOND,
            const math::Vector3N &POS,
            const float &PROBE_RADIUS,
            const float &LENGTH
    )
    {
        math::Vector3N result;

        float
            r1p( FIRST->GetRadius() + PROBE_RADIUS),
            r2p( SECOND->GetRadius() + PROBE_RADIUS),
            d1( math::Distance( POS, FIRST->GetPosition())),
            d2( math::Distance( POS, SECOND->GetPosition()));

        if
        (
                d1 < FIRST->GetRadius()
                || d2 < SECOND->GetRadius()
                || d1 > r1p
                || d2 > r2p
        )
        {
            return result;
        }

        float
            d12( math::Distance( FIRST->GetPosition(), SECOND->GetPosition())),
            gamma( math::AngleFromLawOfCosinus( d2, d12, d1)),
            delta( math::AngleFromLawOfCosinus( r2p, d12, r1p)),
            alpha( delta - gamma),
            psi( math::AngleFromLawOfCosinus( r1p, d12, r2p)),
            phi( math::AngleFromLawOfCosinus( d1, d12, d2)),
            beta( psi - phi),
            dp( math::LawOfCosinus( d1, r1p, alpha));

#ifdef DEBUG
        if( !math::IsEqualWithinThreshold( dp, math::LawOfCosinus( d2, r2p, beta), 0.01))
        {
            std::cout << "===> " << dp << " vs " << math::LawOfCosinus( d2, r2p, beta) << " diff: " << fabs( dp - math::LawOfCosinus( d2, r2p, beta)) << std::endl;
        }
#endif

        if( dp > PROBE_RADIUS && alpha > 0 && beta > 0)
        {
            math::Vector3N
                x( ( SECOND->GetPosition() - FIRST->GetPosition()).SetToLength( 1.0)),
                y( ( math::CrossProduct( math::CrossProduct( x, ( POS - FIRST->GetPosition())), x)).SetToLength( 1.0));

            result += r1p * ( cos( delta) * x + sin( delta) * y);
        }


        return result.SetToLength( LENGTH);
    }





    /////////////////////////
    //      Read/Write     //
    /////////////////////////

    std::istream &AtomForceGrid::Read( std::istream &STREAM)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
    	std::string str;
    	STREAM >> str;
    	if( str != mystr::GetClassName( __PRETTY_FUNCTION__))
    	{
    		std::cout << "===> not the correct header <" << str << "> for: " << __FUNCTION__ << std::endl;
    	}

	typedef  boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >  ptr_func;

	STREAM >> str; 

	if( str == "Grid1D")
	{
	    StandardWrite( "mount grid1D");
	    m_Grid = boost::shared_ptr< store::PositionGrid< ptr_func> >( new store::Grid1D< ptr_func>());
	}
	else if( str == "GridVector3D" || str == "store::GridVector3D<t_RETURN>")
	{
	    StandardWrite( "mount grid3D");
	    m_Grid = boost::shared_ptr< store::PositionGrid< ptr_func> >( new store::GridVector3D< ptr_func>());
	}
#ifdef SQLITE
	else if( str == "GridSqlite")
	{
	    StandardWrite( "mount gridSqlite");
	    m_Grid = boost::shared_ptr< store::PositionGrid< ptr_func> >( new store::GridSqlite< ptr_func>());
	}
#endif
	else
	{
	    std::cerr << "undefined grid type: " << str << std::endl;
	    exit( -1);
	}

 	// Memory();

    	return m_Grid->Read( STREAM);
    }

    std::string AtomForceGrid::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }




    std::ostream &AtomForceGrid::Write( std::ostream &STREAM) const
    {
        STREAM << mystr::GetClassName( __PRETTY_FUNCTION__) << std::endl;
        m_Grid->Write( STREAM);
        return STREAM;
    }

    std::ostream &AtomForceGrid::WriteAsPdb( std::ostream &STREAM) const
    {

    	// flush after each char:
    	STREAM.unsetf( std::ios::unitbuf);

        math::Vector3N
            min( m_Grid->GetMin()),
            max( m_Grid->GetMax()),
            delta( m_Grid->GetDelta());


store::VoidGridPoint< mol::Atom, math::Vector3N >();

        // iterate through grid
        float
            x_start = min( 0) + 0.5 * delta( 0),
            x = x_start,
            y_start = min( 1) + 0.5 * delta( 1),
            y = y_start,
            z_start = min( 2) + 0.5 * delta( 2),
            z = z_start;


        int atom_serial = 0;

        float magnitude;

        char atom_name;

        util::FunctorEnum
	    point_type;

        boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
			grid_point;

        mol::Atom atom = mol::CreateNeutralAtom();

        for( ; x < max( 0); x += delta( 0))
        	for( y = y_start; y < max( 1); y += delta( 1))
        		for( z = z_start; z < max( 2); z += delta( 2))
        		{
        			grid_point = m_Grid->GetGridPoint( math::Vector3N( x, y, z));
        			point_type = grid_point->GetClassID();
        			if( point_type == util::e_ConstantGridPoint)
        			{
        				atom_name = 'C';
        			}
        			else if( point_type == util::e_InteractionGridPoint)
        			{
        				atom_name = 'O';
        			}
        			else if( point_type == util::e_VoidGridPoint)
        			{
        				atom_name = 'H';
        			}
        			else
        			{
        				std::cout << "===> " << __FUNCTION__ << " not written yet for the case of added surf objects (grid point type: " << point_type << ")" << std::endl;
        				continue;
        			}
        			magnitude = grid_point->operator()( atom).Length();

        			STREAM << "ATOM  ";
        			STREAM.setf( std::ios::right, std::ios::adjustfield);
        			STREAM.width( 5);
        			STREAM << (++atom_serial)%100000 << " ";
        			STREAM.width( 4);
        			STREAM << atom_name << " ";
        			STREAM << "ALA A   1    ";
        			STREAM.width( 8);
        			STREAM.precision( 3);
        			STREAM << x;
        			STREAM.width( 8);
        			STREAM.precision( 3);
//        			STREAM.fill( '0');
        			STREAM << y;
        			STREAM.width( 8);
        			STREAM.precision( 3);
        			STREAM << z << "      ";
        			STREAM.width( 6);
        			STREAM.precision( 2);
        			STREAM << magnitude;
        			STREAM << std::endl;
        		}
    	return STREAM;
    }

} // end namespace mol



