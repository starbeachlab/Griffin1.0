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


#include "../../include/geometry/volume_functions.h"

#include <time.h>



namespace geom
{

    float VolumeOfCapsOfOverlappingSpheres( const Sphere &A, const Sphere &B)
    {
        float distance( math::Distance( A.GetPosition(), B.GetPosition()));
        return
            math::Pi * pow( A.GetRadius() + B.GetRadius() - distance, 2)
            * ( pow( distance, 2) + 2.0 * distance * B.GetRadius() - 3.0 * pow( B.GetRadius(), 2) + 2.0 * distance * A.GetRadius() + 6.0 * A.GetRadius() * B.GetRadius() - 3.0 * pow( A.GetRadius(), 2))
            / ( 12.0 * distance);
    }


    bool IsSphereInsideCube( const Cube &CUBE, const Sphere &SPHERE)
    {
//        float diagonal( sqrt( 0.75) * CUBE.GetSizes()[ 0]);

        for( size_t i = 0; i < 3; ++i)
        {
            if( SPHERE.GetPosition()( i) < CUBE.GetPosition()( i) - 0.5 * CUBE.GetSizes()[ i] + SPHERE.GetRadius()
                    || SPHERE.GetPosition()( i) > CUBE.GetPosition()( i) + 0.5 * CUBE.GetSizes()[ i] - SPHERE.GetRadius())
            {
                return false;
            }
        }
        return true;
    }

    GeometricRelation IdentifyGeometricRelation( const Cube &CUBE, const Sphere &SPHERE)
    {
        math::Vector3N cube_pos( CUBE.GetPosition());
        store::Vector3N< float> cube_sizes( CUBE.GetSizes());

        DebugWrite( __FUNCTION__);
        DebugWrite( "cube pos: " <<cube_pos);
        DebugWrite( "cube sizes: " << cube_sizes( 0));
        DebugWrite( "sphere pos: " << SPHERE.GetPosition());
        DebugWrite( "sphere radius: " << SPHERE.GetRadius());

//        size_t checksum( 0);

        // check for no overlap
        for( size_t i = 0; i < 3; ++i)
        {
            if( SPHERE.GetPosition()( i) > CUBE.GetPosition()( i) + 0.5 * CUBE.GetSizes()[ i] + SPHERE.GetRadius()
                    ||     SPHERE.GetPosition()( i) < CUBE.GetPosition()( i) - 0.5 * CUBE.GetSizes()[ i] - SPHERE.GetRadius())
            {
                DebugWrite( "no overlap");
                return NO_OVERLAP;
            }
        }

        float diagonal( sqrt( 0.75) * cube_sizes( 0));

        DebugWrite( "diagonal: " << diagonal);

        // check for complete burial of cube in sphere
        if( SPHERE.GetRadius() >= diagonal && math::SquaredDistance( SPHERE.GetPosition(), CUBE.GetPosition()) < pow( SPHERE.GetRadius() - diagonal, 2))
        { // the second condition will be not exact at some areas close to sphere surface
            DebugWrite( "cube buried in sphere");
            return A_BURIED_IN_B;
        }
        // check for complete burial of sphere in cube
        else if( cube_sizes( 0) > 2.0 * SPHERE.GetRadius() && IsSphereInsideCube( CUBE, SPHERE))
        {
            DebugWrite( "sphere inside cube");
            return B_BURIED_IN_A;
        }
        DebugWrite( "partly overlapped");
        return PARTLY_OVERLAPPED;
    }


     store::Limits3D CalcLimits( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE)
    {
        store::Limits3D limits;
        math::Vector3N position;

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator itr( MOLECULE->GetAtoms().begin()); itr != MOLECULE->GetAtoms().end(); ++itr)
        {
            position = ( *itr)->GetPosition();
            if( position( 0) - ( *itr)->GetVanDerWaalsRadius() < limits.GetXMin())
            { limits.SetXMin( position( 0) - ( *itr)->GetVanDerWaalsRadius());}
            if( position( 0) + ( *itr)->GetVanDerWaalsRadius() > limits.GetXMax())
            { limits.SetXMax( position( 0) + ( *itr)->GetVanDerWaalsRadius());}

            if( position( 1)  - ( *itr)->GetVanDerWaalsRadius()< limits.GetYMin())
            { limits.SetYMin( position( 1) - ( *itr)->GetVanDerWaalsRadius());}
            if( position( 1) + ( *itr)->GetVanDerWaalsRadius() > limits.GetYMax())
            { limits.SetYMax( position( 1) + ( *itr)->GetVanDerWaalsRadius());}

            if( position( 2)  - ( *itr)->GetVanDerWaalsRadius()< limits.GetZMin())
            { limits.SetZMin( position( 2) - ( *itr)->GetVanDerWaalsRadius());}
            if( position( 2) + ( *itr)->GetVanDerWaalsRadius() > limits.GetZMax())
            { limits.SetZMax( position( 2) + ( *itr)->GetVanDerWaalsRadius());}
        }
        return limits;
    }

     store::Limits3D CalcLimits( const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLECULES)
    {
        store::Limits3D limits;

        for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator itr( MOLECULES->begin()); itr != MOLECULES->end(); ++itr)
        {
            limits.Merge( CalcLimits( *itr));
        }
        return limits;
    }



//
//    float
//    Volume
//    (
//            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
//            const store::Limits3D &LIMITS,
//            const float &VOXELSIZE,
//            const size_t &NR_LEVELS,
//            float &VOLUME
//    )
//    {
//        bool shall_we_continue;//, sub_grid;
//
//        assert( LIMITS.AreLimitsInAgreementWithVoxelSize( 1e-8));
//
//        float rmax = sqrt( 0.75) * VOXELSIZE;
//        float distance;
//
//        for( float x = LIMITS.GetXMin(); x + VOXELSIZE <= LIMITS.GetXMax(); x += VOXELSIZE)
//            for( float y = LIMITS.GetYMin();y + VOXELSIZE <= LIMITS.GetYMax(); y += VOXELSIZE)
//                for( float z = LIMITS.GetZMin(); z + VOXELSIZE <= LIMITS.GetZMax(); z += VOXELSIZE)
//                {
//                    shall_we_continue = true;
////                    sub_grid = false;
//                    boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > atoms( new mol::SimpleMolecule< mol::Atom>());
//
//                    for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr( MOLECULE->GetAtoms().begin()); atom_itr != MOLECULE->GetAtoms().end() && shall_we_continue; ++atom_itr)
//                    {
//                        distance = math::Distance( ( *itr)->GetPosition(), math::Vector3N( x + 0.5 * VOXELSIZE, y + 0.5 * VOXELSIZE, z + 0.5 * VOXELSIZE));
//                        if( distance + rmax < ( *itr)->GetVanDerWaalsRadius())
//                        {
//                            VOLUME += volume;
//                            shall_we_continue = false;
////                            sub_grid = false;
//                            atoms = new mol::SimpleMolecule< mol::Atom>();
//                        }
//                        else if( distance - rmax > ( *itr)->GetVanDerWaalsRadius())
//                        {
////                            sub_grid = true;
//                            atoms->AddAtom( *itr);
//                        }
//                    }
//
//                    if( atoms->GetAtoms().size() > 0 && NR_LEVELS > 0)
//                    {
//
//                        store::Limits3D new_limits( x, x + VOXELSIZE, y, y + VOXELSIZE, z, z + VOXELSIZE);
//                        Volume( atoms, new_limits, 0.5 * VOXELSIZE, NR_LEVELS - 1, VOLUME);
//                    }
//                }
//
//    }
//
//
//


    float
    Volume
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const store::Limits3D &LIMITS,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS,
            float &VOLUME
    )
    {
        DebugWrite( __FUNCTION__ << " voxelsize: " << VOXELSIZE << " nr-levels: " << NR_LEVELS);

        bool shall_we_continue;//, sub_grid;

        assert( LIMITS.AreLimitsInAgreementWithVoxelSize( VOXELSIZE, 1e-8));

        DebugWrite( "Voxelsize: " << VOXELSIZE << " #levels " << NR_LEVELS << " VOL: " << VOLUME);

        float volume( pow( VOXELSIZE, 3));

        for( float x = LIMITS.GetXMin(); x < LIMITS.GetXMax(); x += VOXELSIZE)
            for( float y = LIMITS.GetYMin();y < LIMITS.GetYMax(); y += VOXELSIZE)
                for( float z = LIMITS.GetZMin(); z < LIMITS.GetZMax(); z += VOXELSIZE)
                {
                    shall_we_continue = true;
                    boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > atoms( new mol::SimpleMolecule< mol::Atom>());

                    DebugWrite( "voxel: " << x << " " << y << " " << z);

                    for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr( MOLECULE->GetAtoms().begin()); atom_itr != MOLECULE->GetAtoms().end() && shall_we_continue; ++atom_itr)
                    {
                        GeometricRelation relation( IdentifyGeometricRelation( Cube( math::Vector3N( x + 0.5 * VOXELSIZE, y + 0.5 * VOXELSIZE, z + 0.5 * VOXELSIZE), VOXELSIZE), Sphere( ( *atom_itr)->GetPosition(), ( *atom_itr)->GetVanDerWaalsRadius())));
                        switch( relation)
                        {
                        case A_BURIED_IN_B:
                            VOLUME += volume;
                            atoms = boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >( new mol::SimpleMolecule< mol::Atom>());
                            shall_we_continue = false;
                            DebugWrite( "cube buried in sphere, add: " << volume);
                            break;
                        case B_BURIED_IN_A:
                          DebugWrite( "sphere buried in cube!?!");
                        case PARTLY_OVERLAPPED:
                          DebugWrite( "add atom");
                            atoms->AddAtom( *atom_itr);
                            break;
                        default: break;
                        }
                    }

                    if( atoms->GetAtoms().size() > 0 && NR_LEVELS > 0)
                    {
                      DebugWrite( "iterate");
                        store::Limits3D new_limits( x, x + VOXELSIZE, y, y + VOXELSIZE, z, z + VOXELSIZE);
                        Volume( atoms, new_limits, 0.5 * VOXELSIZE, NR_LEVELS - 1, VOLUME);
                    }
                    //TODO !!!: performance
//                    else if( atoms->GetAtoms().size() == 1)
//                    {
//                        VOLUME += VolumeOfSphereInCube( );
//                    }
                    else if( atoms->GetAtoms().size() > 0)
                    {
                      DebugWrite( "#atoms: " << atoms->GetAtoms().size());
                      VOLUME +=  volume/ 14.0;
                      DebugWrite( "add " << volume << " leads to: " << VOLUME);
                    }
                }
        return VOLUME;
    }


    float
    Volume
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS
    )
    {
        DebugWrite( __FUNCTION__ << " voxelsize: " << VOXELSIZE << " nr-levels: " << NR_LEVELS);
        store::Limits3D limits( CalcLimits( MOLECULE));
        DebugWrite( "limits: " << limits);
        limits.AdjustToVoxelSize( VOXELSIZE);
        DebugWrite( "adjusted limits: " << limits);
        float volume( 0.0);
        return Volume( MOLECULE, limits, VOXELSIZE, NR_LEVELS, volume);
    }





    std::vector< store::Map< std::string, float> >
    Volume
    (
            const boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURFACES,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS
    )
    {
        DebugWrite( __FUNCTION__ << " voxelsize: " << VOXELSIZE << " nr-levels: " << NR_LEVELS);

        std::vector< store::Map< std::string, float> > sliced_volume( 4);
        for( size_t i = 0; i < 4; ++i)
        {
            DebugWrite( "layer-for-volume: " << i);
            for( std::map< std::string, geom::DirectedPointSurface>::const_iterator itr = SURFACES->begin(); itr != SURFACES->end(); ++itr)
            {
                DebugWrite( "type: " << itr->first);
                sliced_volume[ i].InsertNewKeyAndValue( itr->first, Volume( boost::shared_ptr< geom::DirectedPointSurface>( new DirectedPointSurface( itr->second)), MEMBRANE->FourLayerBoxLimits( i), VOXELSIZE, NR_LEVELS));
            }
        }
        return sliced_volume;
    }



    float
    Volume
    (
            const boost::shared_ptr< geom::DirectedPointSurface> &SURFACE,
            const store::Limits3D &LIMITS,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS
    )
    {
        DebugWrite( __FUNCTION__ << " voxelsize: " << VOXELSIZE << " nr-levels: " << NR_LEVELS);

        float
            volume( 0.0),
            voxel_area( pow( VOXELSIZE, 2)),
            voxel_volume( pow( VOXELSIZE, 3));

        for( float x = LIMITS.GetXMin(); x < LIMITS.GetXMax(); x += VOXELSIZE)
            for( float y = LIMITS.GetYMin();y < LIMITS.GetYMax(); y += VOXELSIZE)
                for( float z = LIMITS.GetZMin(); z < LIMITS.GetZMax(); z += VOXELSIZE)
                {
//                    DebugWrite( x << " " << y << " " << z);

                    Cube
                        voxel( math::Vector3N( x + 0.5 * VOXELSIZE, y + 0.5 * VOXELSIZE, z + 0.5 * VOXELSIZE), VOXELSIZE);
//                    boost::shared_ptr< DirectedSurfacePoint>
//                        closest_surface_point( SURFACE.ClosestSurfacePoint( voxel.GetPosition()));

                    boost::shared_ptr< std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> > >
                        surf_point_map( SURFACE->DistanceSortedSurfacePoints( voxel.GetPosition()));
                    boost::shared_ptr< DirectedSurfacePoint>
                        closest_surface_point( surf_point_map->begin()->second);
                    bool
                        is_surf_point_inside_voxel( voxel.IsPointWithin( closest_surface_point->GetPosition()));
                    math::Vector3N
                        connection( voxel.GetPosition() - closest_surface_point->GetPosition());
                    float
                        angle( math::Angle( connection, closest_surface_point->GetNormalVector()));

                    if( !is_surf_point_inside_voxel)
                    {
                        if( angle > math::HalfPi)
                        {
                            DebugWrite( "voxel is inside volume");
                            volume += voxel_volume;
                        }
//                        else
//                        {
//                            DebugWrite( "voxel is outside volume");
//                        }
                    }
                    // if surface is cutting voxel, either go to next level or estimate volume
                    else
                    {
                        DebugWrite( "voxel is partly inside");
                        // devide voxel into 8
                        if( NR_LEVELS > 0)
                        {
//                            geom::DirectedPointSurface surf;
//                            surf.PushBack( closest_surface_point);

                            boost::shared_ptr< geom::DirectedPointSurface>
                                surf( SurfacePointsWithinCube( voxel, surf_point_map));

                            store::Limits3D new_limits( x, x + VOXELSIZE, y, y + VOXELSIZE, z, z + VOXELSIZE);
                            volume += Volume( surf, new_limits, 0.5 * VOXELSIZE, NR_LEVELS - 1);     ////////////////// COLLECT SURF POINTS AND PASS TO FCT ////////////////////////
                        }
                        // estimate voxel volume
                        else
                        {
                            float distance_to_plane( connection.Length() * cos( angle));
                            DebugWrite( "add volume: " << ( ( 0.5 * VOXELSIZE + distance_to_plane) * voxel_area));
                            volume += ( 0.5 * VOXELSIZE + distance_to_plane) * voxel_area;

//                            if( angle > math::HALF_PI)
//                            {
//                                angle = math::HALF_PI - angle;
//                                float distance_to_plane( connection.Length() * cos( angle));
//                                volume += ( 0.5 * VOXELSIZE + distance_to_plane) * voxel_area;
//                            }
//                            else
//                            {
//                                float distance_to_plane( connection.Length() * cos( angle));
//                                volume += ( 0.5 * VOXELSIZE - distance_to_plane) * voxel_area;
//                            }
                        }
                    }

                }
        return volume;
    }




    boost::shared_ptr< geom::DirectedPointSurface>
    SurfacePointsWithinCube
    (
            const Cube &CUBE,
            const boost::shared_ptr< std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> > > &DISTANCE_MAP
    )
    {
        boost::shared_ptr< DirectedPointSurface>
            surf( new DirectedPointSurface());
        for( std::multimap< float, boost::shared_ptr< DirectedSurfacePoint> >::const_iterator itr = DISTANCE_MAP->begin(); itr != DISTANCE_MAP->end(); ++itr)
        {
            if( CUBE.IsPointWithin( itr->second->GetPosition()))
            {
                surf->PushBack( itr->second);
            }
        }
        return surf;
    }



    bool
    IsMoleculeWithinLimits
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const store::Limits3D &LIMITS
    )
    {
        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr( MOLECULE->GetAtoms().begin()); atom_itr != MOLECULE->GetAtoms().end(); ++atom_itr)
        {
            if( LIMITS.IsWithin( ( *atom_itr)->GetPosition()))
            {
                return true;
            }
        }
        return false;
    }


    size_t
    NumberMoleculesInLimits
    (
            const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
            const store::Limits3D &LIMITS
    )
    {
        size_t count( 0);
        for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( MOLS->begin()); mol_itr != MOLS->end(); ++mol_itr)
            if( IsMoleculeWithinLimits( *mol_itr, LIMITS))
            {
                ++count;
            }
        return count;
    }


    float
    SimpleDensity
    (
            const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
            const store::Limits3D &LIMITS
    )
    {
        size_t count( NumberMoleculesInLimits( MOLS, LIMITS));
        return float( count) / LIMITS.Volume();
    }



//    float
//    SumOverSquaredDepthOfPenetration
//    (
//            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
//            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > MOLECULE
//    )
//    {
//        DebugWrite( __FUNCTION__);
//
//        float
//            score( 0.0),
//            penetration;
////        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( MOLECULE->GetAtoms().begin()); mol_itr != MOLECULE->GetAtoms().end(); ++mol_itr)
////        {
////            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator inserted_itr( INSERTED_MOLECULE->GetAtoms().begin()); inserted_itr != INSERTED_MOLECULE->GetAtoms().end(); ++inserted_itr)
////            {
////                penetration = math::Distance( ( *mol_itr)->GetPosition(), ( *inserted_itr)->GetPosition()) - ( *mol_itr)->GetVanDerWaalsRadius() - ( *inserted_itr)->GetVanDerWaalsRadius();
////                if( penetration < 0)
////                {
////                    score += pow( penetration, 2);
////                }
////            }
////        }
//        return score;
//    }



    std::pair< size_t, float>
    SumOverSquaredDepthOfPenetration
    (
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            geom::DirectedPointSurface SURFACE
    )
    {
        DebugWrite( __FUNCTION__);

        float
            score( 0.0);
        size_t
            buried( 0);


        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( INSERTED_MOLECULE->GetAtoms().begin()); mol_itr != INSERTED_MOLECULE->GetAtoms().end(); ++mol_itr)
        {
            boost::shared_ptr< DirectedSurfacePoint>
                closest_surface_point( SURFACE.ClosestSurfacePoint( ( *mol_itr)->GetPosition()));
            math::Vector3N
                connection( ( *mol_itr)->GetPosition() - closest_surface_point->GetPosition());

            if( connection * closest_surface_point->GetNormalVector() < 0) //  ~ angle > math::HalfPi
            {
                ++buried;
                score += pow( connection.Length(), 2);  // ~ pow( connection.Length() * cos( angle) , 2);
            }
        }

        if( score > 500)
        {
            std::cout << "high scoring = deeply buried, check: " << INSERTED_MOLECULE << std::endl;
            std::cout << "detailed score list: total:" << score << std::endl;
            size_t cc( 0);
            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( INSERTED_MOLECULE->GetAtoms().begin()); mol_itr != INSERTED_MOLECULE->GetAtoms().end(); ++mol_itr, ++cc)
            {
                boost::shared_ptr< DirectedSurfacePoint>
                    closest_surface_point( SURFACE.ClosestSurfacePoint( ( *mol_itr)->GetPosition()));
                math::Vector3N
                    connection( ( *mol_itr)->GetPosition() - closest_surface_point->GetPosition());

                if( connection * closest_surface_point->GetNormalVector() < 0) //  ~ angle > math::HalfPi
                {
                    std::cout << cc << " bury depth: " << connection.Length() << std::endl;
                }
            }

        }

#ifdef DEBUG
        if( score > 0)
        {
            std::cout << INSERTED_MOLECULE->GetType() << " " << INSERTED_MOLECULE->GetAtoms()( 0)->GetResidueID() << " is buried: " << score << std::endl;
        }
#endif

        return std::make_pair( buried, score);
    }




    std::pair< size_t, float>
    SumOverSquaredDepthOfPenetration
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID
    )
    {
//              StandardWrite( __FUNCTION__);
        float
            score( 0.0);
        size_t
            buried( 0);

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( INSERTED_MOLECULE->GetAtoms().begin()); mol_itr != INSERTED_MOLECULE->GetAtoms().end(); ++mol_itr)
        {
            short value = SURF_GRID->operator()( ( *mol_itr)->GetPosition());
            //            std::cout << value << " ";
            //            std::cout.flush();
            if( value < 1)
            {
              //              std::cout << "smaller ";
              //              std::cout.flush();
                if( SURF_GRID->IsInside( value, ( *mol_itr)->GetResidueType()))
                {
                  //                  std::cout << " inside ";
                  //              std::cout.flush();
                    ++buried;
                }
                score += SURF_GRID->BasicClosestSurfPoint( ( *mol_itr)->GetPosition(), ( *mol_itr)->GetResidueType()).first;
                //                std::cout << " scored ";
                //              std::cout.flush();
            }
            //            std::cout << std::endl;
        }
        //        std::cout << __FUNCTION__ << " done" << std::endl;
        return std::make_pair( buried, score);
    }






    bool
    IsHeavyAtomOfMoleculePenetrating
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID
    )
    {
        DebugWrite( __FUNCTION__);

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( INSERTED_MOLECULE->GetAtoms().begin()); mol_itr != INSERTED_MOLECULE->GetAtoms().end(); ++mol_itr)
        {
            DebugWrite( "check for: " << **mol_itr);
            if( ( *mol_itr)->GetMass() > 1.5)
            {
                short value = SURF_GRID->operator()( ( *mol_itr)->GetPosition());
                if( value < 1)
                {
                    if( SURF_GRID->IsInside( value, ( *mol_itr)->GetResidueType()))
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }





    bool
    IsHeavyAtomOfMoleculePenetrating
    (
            boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > INSERTED_MOLECULE,
            geom::DirectedPointSurface SURFACE
    )
    {
        DebugWrite( __FUNCTION__);

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( INSERTED_MOLECULE->GetAtoms().begin()); mol_itr != INSERTED_MOLECULE->GetAtoms().end(); ++mol_itr)
        {
            DebugWrite( "check for: " << **mol_itr);
            if( ( *mol_itr)->GetMass() > 1.5)
            {
                DebugWrite( "heavy atom");
                boost::shared_ptr< DirectedSurfacePoint>
                    closest_surface_point( SURFACE.ClosestSurfacePoint( ( *mol_itr)->GetPosition()));
                math::Vector3N
                    connection( ( *mol_itr)->GetPosition() - closest_surface_point->GetPosition());

                if( connection * closest_surface_point->GetNormalVector() < 0)
                {
                    DebugWrite( "inside");
                    return true;
                }
            }
        }
        return false;
    }


    boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
    RemoveMoleculesInImplicitMembraneVolume
    (
            std::ostream &STREAM,
            boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > SURF_MAP,
            const float &VOXELSIZE,
            const size_t &NR_LEVELS
    )
    {
        StandardWrite( __FUNCTION__);

        time_t
            start,
            end;
        time( &start);

        StandardWrite( "calculate volume ...");
        std::vector< store::Map< std::string, float> >    // vec< map< 'all|other|other', volume> >
            volume_map( Volume( MEMBRANE, SURF_MAP, VOXELSIZE, NR_LEVELS));

        time( &end);
        std::cout << "volume calulation: " << difftime( end, start) << "s" << std::endl;
        time( &start);

        std::cout << "\nvolume map:" << std::endl;
        size_t
            count( 0);

        for( std::vector< store::Map< std::string, float> >::const_iterator layer_itr = volume_map.begin(); layer_itr != volume_map.end(); ++layer_itr, ++count)
        {
            std::cout << "layer: " << count << std::endl;
            for( std::map< std::string, float>::const_iterator itr = layer_itr->begin(); itr != layer_itr->end(); ++itr)
            {
                std::cout << itr->first << ":   " << itr->second << std::endl;
            }
        }
        std::cout << std::endl;

        boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
            remaining( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());

//        std::ofstream write( "output/deleted.pdb");

        // find the ones that are deepest buried in implicit volume for each layer
        for( size_t i = 0; i < 4; ++i)
        {
            StandardWrite( "=== layer: " << i << " ====");
            store::Map< std::string, float> density_map( MEMBRANE->Densities( i));  // molecules per layer volume // map< POPE|TIP3, density>  no 'all'!!

            std::cout << "densities of molecules in layer:" << std::endl;
            for( std::map< std::string, float>::const_iterator itr = density_map.begin(); itr != density_map.end(); ++itr)
            {
                std::cout << itr->first << ": " << itr->second << std::endl;
            }
            std::cout << std::endl;


#ifdef DEBUG
            std::cout << "keys: ";
            std::vector< std::string> keys( MEMBRANE->GetMoleculeLayer( i).GetKeys());
            std::copy(keys.begin(), keys.end(), std::ostream_iterator< std::string>( std::cout, " "));
            std::cout << std::endl << std::endl;
#endif

            // for each layer check each molecule type individually
            for( std::map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > >::const_iterator type_itr = MEMBRANE->GetMoleculeLayer( i).begin(); type_itr != MEMBRANE->GetMoleculeLayer( i).end(); ++type_itr)
            {

                std::string type;
                // PASS to own class
                if( SURF_MAP->IsValidKey( type_itr->first))
                {
                    type = type_itr->first;
                }
                else
                {
                    type = "all";
                }

                StandardWrite( "molecule-type: <" << type_itr->first << "> belongs to surface type: " << type);

                float
                    density( density_map( type_itr->first)),
                    volume( volume_map[ i]( type));
                size_t
                    remaining_lipids( 0),
                    erased_lipids( 0),
                    remaining_waters( 0),
                    erased_waters( 0);

                if( std::find( MEMBRANE->GetLipidKeys().begin(), MEMBRANE->GetLipidKeys().end(), type_itr->first) != MEMBRANE->GetLipidKeys().end())  // lipids
                {
                    int
                        mols_to_be_removed( math::Round< size_t>( volume * density));
                    std::multimap< std::pair< size_t, float>, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >, LargerThanForPairs< size_t, float> >
                        depth_sorted_molecules;

                    StandardWrite( "lipid density: " << density << " volume: " << volume << " #mols-to-be-removed: " << mols_to_be_removed);

                    for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( type_itr->second->begin()); mol_itr != type_itr->second->end(); ++mol_itr)
                    {
                        depth_sorted_molecules.insert( std::make_pair( SumOverSquaredDepthOfPenetration( *mol_itr, ( *SURF_MAP)( type)), *mol_itr));
                    }

                    StandardWrite( "depth sorted molecules (" << depth_sorted_molecules.size() << ")");
                    int count( 0), lc( 0);

                    StandardWrite( "count  deleted/remains  #atoms_buried:  burial score:  mol-ID: ");

                    for( std::multimap< std::pair< size_t, float>, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr = depth_sorted_molecules.begin(); mol_itr != depth_sorted_molecules.end(); ++mol_itr, ++lc)
                    {
                        StandardWriteNoFlush( lc << "  ");
                        if( ++count > mols_to_be_removed)
                        {
                            DebugWrite( "remaining");
                            StandardWriteNoFlush( "  remains  ");
                            remaining->push_back( mol_itr->second);
                            ++remaining_lipids;
                        }
                        else
                        {
                            DebugWrite( "kicked out");
                            StandardWriteNoFlush( "  deleted  ");
                            mol::file::WriteToPdb< mol::Atom>( STREAM, mol_itr->second);
                            ++erased_lipids;
                        }
                        StandardWrite( "  " << mol_itr->first.first << "  " << mol_itr->first.second << "  " << mol_itr->second->GetAtoms()( 0)->GetResidueID());
                    }
                    time( &end);
                    StandardWrite( "remaining lipids: " << remaining_lipids << " erased: " << erased_lipids << " in "<< difftime( end, start) << "s" << std::endl);
                    time( &start);
                }
                else  // water
                {
                    StandardWrite( "removing all waters, whose oxygen atoms are with the surface");
                    for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( type_itr->second->begin()); mol_itr != type_itr->second->end(); ++mol_itr)
                    {
                        if( IsHeavyAtomOfMoleculePenetrating( *mol_itr, ( *SURF_MAP)( type)))
                        {
                            DebugWrite( "kicked out");
                            mol::file::WriteToPdb< mol::Atom>( STREAM, *mol_itr);
                            ++erased_waters;
                        }
                        else
                        {
                            DebugWrite( "remaining");
                            remaining->push_back( *mol_itr);
                            ++remaining_waters;
                        }
                    }
                    time( &end);
                    StandardWrite( " waters: remaining: " << remaining_waters << " erased: " << erased_waters << " in "<< difftime( end, start) << "s" << std::endl);
                    time( &start);

                }

                DebugWrite( "remaining: " << remaining->size() << " atoms");
            }
        }
//        write.close();
//        write.clear();

        return remaining;
    }


//        for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( MOLS->begin()); mol_itr != MOLS->end(); ++mol_itr)
//        {
//            z_coordinate = 0.0;
//            distance = 0.0;
//            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr( ( *mol_itr)->GetAtoms().begin()); atom_itr != ( *mol_itr)->GetAtoms().end(); ++atom_itr)
//            {
//                distance += SumOverSquaredDepthOfPenetration( SURF_MAP, *mol_itr);
//                z_coordinate += ( *atom_itr)->GetPosition()( 2);
//            }
//            if( z_coordinate > 0)
//            {
//                upper_distance_sorted_atoms[ ( *mol_itr)->GetType()].insert( std::make_pair( distance, *mol_itr));
//            }
//            else
//            {
//                lower_distance_sorted_atoms[ ( *mol_itr)->GetType()].insert( std::make_pair( distance, *mol_itr));
//            }
//        }
//
//#ifdef DEBUG
//        std::cout << "write distance sorted atom list: " << std::endl;
//        std::cout << "lower layer: " << std::endl;
//        for( std::map< std::string, std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > > >::const_iterator type_itr( lower_distance_sorted_atoms.begin()); type_itr != lower_distance_sorted_atoms.end(); ++type_itr)
//            for( std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator atom_itr( type_itr->second.begin()); atom_itr != type_itr->second.end(); ++atom_itr)
//            {
//                std::cout << type_itr->first << "  " << atom_itr->first << "  " << ( *atom_itr->second)->GetType() << std::endl;
//            }
//        std::cout << "upper layer: " << std::endl;
//        for( std::map< std::string, std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > > >::const_iterator type_itr( upper_distance_sorted_atoms.begin()); type_itr != upper_distance_sorted_atoms.end(); ++type_itr)
//            for( std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator atom_itr( type_itr->second.begin()); atom_itr != type_itr->second.end(); ++atom_itr)
//            {
//                std::cout << type_itr->first << "  " << atom_itr->first << "  " << ( *atom_itr->second)->GetType() << std::endl;
//            }
//#endif
//
//
//        ////  calc number of atoms to delete  /////
//        float
//            upper_implicit_volume( 0.0),
//            lower_implicit_volume( 0.0);
//
//        // get limits of molecules
//        store::Limits3D
//            limits,
//            upper_limits,
//            lower_limits;
//
//        if( LIMITS.AreDefined())
//        {
//            limits = LIMITS;
//        }
//        else
//        {
//            limits = CalcLimits( MOLS);
//        }
//
//        DebugWrite( "limits before adjust: " << limits);
//        limits.AdjustToVoxelSize( VOXELSIZE);
//        DebugWrite( "limits after adjust: " << limits);
//
//        // slice volume in the middle of lower layer
//        upper_limits = limits;
//        upper_limits.SetZMin( 0.0);
//        lower_limits = limits;
//        lower_limits.SetZMax( 0.0);
//        DebugWrite( "upper_limits: " << upper_limits);
//        DebugWrite( "lower_limits: " << lower_limits);
//
//        if( USER_DEFINED_SURF_OBJECTS->size() > 0)
//        {
//            Volume( INSERTED_MOLECULE, USER_DEFINED_SURF_OBJECTS, MEMBRANE.UpperLipidLayerLimits(), VOXELSIZE, NR_LEVELS, upper_implicit_volume);
//            Volume( INSERTED_MOLECULE, USER_DEFINED_SURF_OBJECTS, MEMBRANE.LowerLipidLayerLimits(), VOXELSIZE, NR_LEVELS, lower_implicit_volume);
//        }
//        else
//        {
//            Volume( INSERTED_MOLECULE, MEMBRANE.UpperLipidLayerLimits(), VOXELSIZE, NR_LEVELS, upper_implicit_volume);
//            Volume( INSERTED_MOLECULE, MEMBRANE.LowerLipidLayerLimits(), VOXELSIZE, NR_LEVELS, lower_implicit_volume);
//        }
//
//
//        size_t
//            upper_nr_atoms_to_delete( math::Round( float( upper_distance_sorted_atoms.size()) * upper_implicit_volume / upper_limits.Volume())),  // or RoundUp // size_t( value + 0.999999999)
//            lower_nr_atoms_to_delete( math::Round( float( lower_distance_sorted_atoms.size()) * lower_implicit_volume / lower_limits.Volume()));
//
//        MOLS->clear();
//
//        // return remaining
//        for( std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator itr( upper_distance_sorted_atoms.begin() + upper_nr_atoms_to_delete); itr != upper_distance_sorted_atoms.end(); ++itr)
//        {
//            MOLS->push_back( *itr);
//        }
//        for( std::multimap< float, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator itr( lower_distance_sorted_atoms.begin() + lower_nr_atoms_to_delete); itr != lower_distance_sorted_atoms.end(); ++itr)
//        {
//            MOLS->push_back( *itr);
//        }
//    }



    std::vector< store::Map< std::string, float> >
    Volume
    (
            const boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID
    )
    {
        DebugWrite( __FUNCTION__);

        std::vector< store::Map< std::string, float> > sliced_volume( 4);

        for( size_t i = 0; i < 4; ++i)
        {
            DebugWrite( "layer-for-volume: " << i);
            store::Limits3D limits = MEMBRANE->FourLayerBoxLimits( i);
            sliced_volume[ i] = SURF_GRID->Volume( limits.GetMin(), limits.GetMax());
        }
        return sliced_volume;
    }








    boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
    RemoveMoleculesInImplicitMembraneVolume
    (
            std::ostream &STREAM,
            boost::shared_ptr< mol::Membrane> &MEMBRANE,
            const boost::shared_ptr< SurfGrid> &SURF_GRID,
            const std::string &NR_LIPIDS_LOWER_LAYER,
            const std::string &NR_LIPIDS_UPPER_LAYER,
            const std::vector< std::pair< std::string, int> > &KEEP_MOLS,
            const std::vector< std::pair< std::string, int> > &REMOVE_MOLS
    )
    {
        DebugWrite( __FUNCTION__);

        time_t
            start,
            end;
        time( &start);


        int
			total_lipids_remaining = 0,
			total_waters_remaining = 0,
			total_lipids_erased = 0,
			total_waters_erased = 0;

//        std::vector< int> nr_lipids_to_remove;
//
//        if( cmd.IsFlagSet( "number_of_lipids_to_carve"))
//        {
//            nr_lipids_to_remove = cmd.GetArgumentsForFlag< int>( "number_of_lipids_to_carve");
//        }
//        else
//        {
//
//        }

        StandardWrite( "calculate volume of bilayer/water box ...");

        std::vector< store::Map< std::string, float> >    // vec: layer : map< type('all|other|other'), volume>
            volume_map = Volume( MEMBRANE, SURF_GRID);

        time( &end);
        std::cout << "volume calulation: " << difftime( end, start) << "s" << std::endl;
        time( &start);

        std::cout << "\nthe limit definitions of the bilayer/water system:" << std::endl;
        std::cout << "layer 0: lower solution layer, from z_lower_total_box_limit to z_lower_lipid_layer_limit\n";
        std::cout << "layer 1: lower bilayer  layer, from z_lower_lipd_layer_limit to z_bilayer_center\n";
        std::cout << "layer 2: upper bilayer  layer, from z_bilayer_center to z_upper_lipid_layer_limit\n";
        std::cout << "layer 3: upper solution layer, from z_upper_lipid_layer_limit to z_upper_total_box_limit" << std::endl;

        size_t
            count( 0);

        std::cout << "\nvolume map of the implicit molecule (as experienced by bilayer molecule types):" << std::endl;
        std::map< std::string, float> volume_sum;//( volume_map[0].size(), 0.0);
        for( std::vector< store::Map< std::string, float> >::const_iterator layer_itr = volume_map.begin(); layer_itr != volume_map.end(); ++layer_itr, ++count)
        {
            std::cout << "layer: " << count << std::endl;
            for( std::map< std::string, float>::const_iterator itr = layer_itr->begin(); itr != layer_itr->end(); ++itr)
            {
                std::cout << itr->first << ":   " << itr->second << " Angstroem^3" << std::endl;
                volume_sum[ itr->first] += itr->second;
            }
        }
        std::cout << std::endl;

        std::cout << "the sum of the volume for each defined molecule type over all 4 layers: " << std::endl;
        for( std::map< std::string, float>::const_iterator map_itr = volume_sum.begin(); map_itr != volume_sum.end(); ++map_itr)
        {
            std::cout << "total-volume " << map_itr->first << ": " << map_itr->second << " Angstroem^3" << std::endl;
        }

        if( volume_sum.size() > 1)
        {
            std::cout << "\nthe difference of the volume for 'all' molecule types and an added, specific one, is the volume that is defined solely for the specific molecule type:" << std::endl;
            for( std::map< std::string, float>::const_iterator map_itr = volume_sum.begin(); map_itr != volume_sum.end(); ++map_itr)
            {
                if( map_itr->first != "all")
                {
                    std::cout << "volume diff " << map_itr->first << " to all: " << fabs( volume_sum[ "all"] - map_itr->second) << " Angstroem^3" << std::endl;
                }
            }
        }
        std::cout << std::endl;



        boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
            remaining( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());


        std::vector< float> z_limits( MEMBRANE->GetZLimits());


//        std::ofstream write( "output/deleted.pdb");

        // find the ones that are deepest buried in implicit volume for each layer
        for( size_t i = 0; i < 4; ++i)
        {
            StandardWrite( "===== layer: " << i << " =====");

            store::Map< std::string, float>
                density_map( MEMBRANE->Densities( i));  // molecules per layer volume // map< POPE|TIP3, density>  no 'all'!!
            float
                bilayer_thickness( z_limits[ i+1] - z_limits[i]);

            std::cout << "molecule types in layer:" << std::endl;
            for( std::map< std::string, float>::const_iterator itr = density_map.begin(); itr != density_map.end(); ++itr)
            {
                std::cout << itr->first << ": " << std::endl;
                std::cout << "\t" << itr->second << " density (molecules per Angstroem^3)" << std::endl;
                std::cout << "\t" << 1.0 / ( itr->second * bilayer_thickness) << " average surface area of the molecule (Angstroem^2 per molecule) (box_x_size * box_y_size / number_molecules)" << std::endl;
            }
            std::cout << std::endl;


#ifdef DEBUG
            std::cout << "keys: ";
            std::vector< std::string> keys( MEMBRANE->GetMoleculeLayer( i).GetKeys());
            std::copy(keys.begin(), keys.end(), std::ostream_iterator< std::string>( std::cout, " "));
            std::cout << std::endl << std::endl;
#endif

            // for each layer check each molecule type individually
            for( std::map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > >::const_iterator type_itr = MEMBRANE->GetMoleculeLayer( i).begin(); type_itr != MEMBRANE->GetMoleculeLayer( i).end(); ++type_itr)
            {

                std::string type;
                // PASS to own class
                if( volume_map[i].IsValidKey( type_itr->first))
                {
                    type = type_itr->first;
                }
                else
                {
                    type = "all";
                }

                StandardWrite( "molecule-type: <" << type_itr->first << "> translates to surface type: " << type);

                float
                    density( density_map( type_itr->first)),
                    volume( volume_map[ i]( type));
                size_t
                    remaining_lipids( 0),
                    erased_lipids( 0),
                    remaining_waters( 0),
                    erased_waters( 0),
                    cc = 0;

                if( std::find( MEMBRANE->GetLipidKeys().begin(), MEMBRANE->GetLipidKeys().end(), type_itr->first) != MEMBRANE->GetLipidKeys().end())  // lipids
                {
                    int
                        mols_to_be_removed( math::Round< size_t>( volume * density));

                    if( cc == 0 && NR_LIPIDS_LOWER_LAYER != "")
                    {
                    	StandardWrite( "lower layer: adjust number of lipids to remove, now: " << mols_to_be_removed);
                    	if( NR_LIPIDS_LOWER_LAYER[0] == '+' || NR_LIPIDS_LOWER_LAYER[0] == '-')
                    	{
                    		mols_to_be_removed += mystr::ConvertStringToNumericalValue< int>( NR_LIPIDS_LOWER_LAYER);
                    	}
                     	else
                    	{
                    		mols_to_be_removed = mystr::ConvertStringToNumericalValue< int>( NR_LIPIDS_LOWER_LAYER);
                    	}
                    	StandardWrite( "adjusted: " << mols_to_be_removed);
                    }
                    else if( cc == 1 && NR_LIPIDS_UPPER_LAYER != "")
                    {
                    	StandardWrite( "upper layer: adjust number of lipids to remove, now: " << mols_to_be_removed);
                    	if( NR_LIPIDS_UPPER_LAYER[0] == '+' || NR_LIPIDS_UPPER_LAYER[0] == '-')
                    	{
                    		mols_to_be_removed += mystr::ConvertStringToNumericalValue< int>( NR_LIPIDS_UPPER_LAYER);
                    	}
                     	else
                    	{
                    		mols_to_be_removed = mystr::ConvertStringToNumericalValue< int>( NR_LIPIDS_UPPER_LAYER);
                    	}
                    	StandardWrite( "adjusted: " << mols_to_be_removed);
                    }

                    std::multimap< std::pair< size_t, float>, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >, LargerThanForPairs< size_t, float> >
                        depth_sorted_molecules;

                    StandardWrite( "\nlipid carving ...");
                    StandardWrite( "lipid density: " << density << " \nvolume of implicit molecule in this layer: " << volume << "\nnumber of lipids to remove: " << mols_to_be_removed);
                    StandardWrite( "average cross section area of the implicit molecule in this layer: " <<  volume / bilayer_thickness);

                    StandardWrite( "\nlipids sorted by number of atoms buried in the implicit molecule and their depth (" << depth_sorted_molecules.size() << " molecules):");

                    StandardWrite( "count  deleted/remains  #atoms_buried:  burial score:  mol-ID: ");

                   int count( 0), lc( 0);


                    for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( type_itr->second->begin()); mol_itr != type_itr->second->end(); ++mol_itr)
                    {
                    	std::pair< std::string, int> pair = std::make_pair( ( *mol_itr)->GetType(), (*mol_itr)->GetID());
                    	if( std::find( KEEP_MOLS.begin(), KEEP_MOLS.end(), pair) != KEEP_MOLS.end())
                    	{
                            DebugWrite( "remaining(list)");
                            StandardWriteNoFlush( "  remains(list)  ");
                            remaining->push_back( *mol_itr);
                            ++remaining_lipids;
                            StandardWrite( "  "  << ( *mol_itr)->GetAtoms()( 0)->GetResidueID()  << "  " << ( *mol_itr)->GetAtoms()( 0)->GetResidueType() << " " << ( *mol_itr)->GetType() << " " << ( *mol_itr)->GetID());
                    	}
                    	else if( std::find( REMOVE_MOLS.begin(), REMOVE_MOLS.end(), pair) != REMOVE_MOLS.end())
                    	{
							DebugWrite( "kicked out(list)");
							StandardWriteNoFlush( "  deleted(list)  ");
							mol::file::WriteToPdb< mol::Atom>( STREAM, *mol_itr);
							++erased_lipids;
							++count;
							StandardWrite( "  " << ( *mol_itr)->GetAtoms()( 0)->GetResidueID()  << "  " << ( *mol_itr)->GetAtoms()( 0)->GetResidueType()<< " " << ( *mol_itr)->GetType() << " " << ( *mol_itr)->GetID());
                    	}
                    	else
                    	{
                    		depth_sorted_molecules.insert( std::make_pair( SumOverSquaredDepthOfPenetration( *mol_itr, SURF_GRID), *mol_itr));
                    	}
                    }

                    for( std::multimap< std::pair< size_t, float>, boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr = depth_sorted_molecules.begin(); mol_itr != depth_sorted_molecules.end(); ++mol_itr, ++lc)
                    {
                        StandardWriteNoFlush( lc << "  ");
                        if( ++count > mols_to_be_removed)
                        {
                            DebugWrite( "remaining");
                            StandardWriteNoFlush( "  remains  ");
                            remaining->push_back( mol_itr->second);
                            ++remaining_lipids;
                        }
                        else
                        {
                            DebugWrite( "kicked out");
                            StandardWriteNoFlush( "  deleted  ");
                            mol::file::WriteToPdb< mol::Atom>( STREAM, mol_itr->second);
                            ++erased_lipids;
                        }
                        StandardWrite( "  " << mol_itr->first.first << "  " << mol_itr->first.second << "  " << mol_itr->second->GetAtoms()( 0)->GetResidueID() << "  " << mol_itr->second->GetAtoms()( 0)->GetResidueType() << " " << mol_itr->second->GetType() << " " << mol_itr->second->GetID());
                    }
                    StandardWrite ( "... lipid carving done");
                    time( &end);
                    total_lipids_remaining += remaining_lipids;
                    total_lipids_erased += erased_lipids;
                    StandardWrite( "remaining lipids: " << remaining_lipids << " erased: " << erased_lipids << " in "<< difftime( end, start) << "s" << std::endl);
                    time( &start);
                    ++cc;
                }
                else  // water++
                {
                    StandardWrite( "\nwater++ carving ...");
                    StandardWrite( "removing all waters, whose oxygen atoms are with the surface");
                    for( std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::const_iterator mol_itr( type_itr->second->begin()); mol_itr != type_itr->second->end(); ++mol_itr)
                    {
                    	std::pair< std::string, int> pair = std::make_pair( ( *mol_itr)->GetType(), (*mol_itr)->GetID());
                    	if( std::find( KEEP_MOLS.begin(), KEEP_MOLS.end(), pair) != KEEP_MOLS.end())
                    	{
                            DebugWrite( "remaining due to list");
                            remaining->push_back( *mol_itr);
                            ++remaining_waters;
                    	}
                    	else
                    	{
							if( IsHeavyAtomOfMoleculePenetrating( *mol_itr, SURF_GRID) || std::find( REMOVE_MOLS.begin(), REMOVE_MOLS.end(), pair) != REMOVE_MOLS.end())
							{
								DebugWrite( "kicked out");
								mol::file::WriteToPdb< mol::Atom>( STREAM, *mol_itr);
								++erased_waters;
							}
							else
							{
								DebugWrite( "remaining");
								remaining->push_back( *mol_itr);
								++remaining_waters;
							}
                    	}
                    }
                    StandardWrite( "... water carving done");
                    time( &end);
                    total_waters_remaining += remaining_waters;
                    total_waters_erased += erased_waters;
                    StandardWrite( " waters: remaining: " << remaining_waters << " erased: " << erased_waters << " in "<< difftime( end, start) << "s" << std::endl);
                    time( &start);

                }

                DebugWrite( "remaining: " << remaining->size() << " atoms");
            }
        }
//        write.close();
//        write.clear();

        StandardWrite( "Total number of lipids remaining: " << total_lipids_remaining << " erased: " << total_lipids_erased);
        StandardWrite( "Total number of non-lipids remaining: " << total_waters_remaining << " erased: " << total_waters_erased);
        StandardWrite( "Total number of erased atoms: " << ( total_lipids_erased + total_waters_erased));
        StandardWrite( "Total number of remaining atoms: " << ( total_lipids_remaining + total_waters_remaining));

        return remaining;
    }







} // end namespace geom


