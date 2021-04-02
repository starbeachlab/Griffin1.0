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


#include "../../include/geometry/surface_factory.h"

namespace geom
{
    namespace factory
    {
        bool DoObjectsOverlap( const boost::shared_ptr< Object> &FIRST, const boost::shared_ptr< Object> &SECOND, const float &PROBE_RADIUS)
        {
//            std::cout << "do they? " << FIRST  << std::endl;
//            std::cout << SECOND << std::endl;
            if( ( FIRST->GetPosition() - SECOND->GetPosition()).Length() > FIRST->MaxRadius() + SECOND->MaxRadius() + PROBE_RADIUS) // exact for 2 spheres, approximated for others
            { return false;}
            return true;
        }

        std::vector< std::vector< size_t> > ListOverlappingObjects( const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS, const float &PROBE_RADIUS)
        {
//            std::cout  << "  list overlapping objects ..." << std::endl;
            std::vector< std::vector< size_t> > list( SURFS.size());
            for( size_t i = 0; i < SURFS.size() - 1; ++i)
                for( size_t k = i + 1; k < SURFS.size(); ++k)
                {
                    if( DoObjectsOverlap( SURFS[ i]->GetObject(), SURFS[ k]->GetObject(), PROBE_RADIUS))
                    {
                        list[ i].push_back( k);
                        list[ k].push_back( i);
                    }
                }
//            std::cout << "  ... done" << std::endl;
            return list;
        }


        boost::shared_ptr< PointSurface> BuildNonOverlappingPointSurface( const store::ShPtrVec< PointSurfaceObject> &SURFS)
        {
            return boost::shared_ptr< PointSurface>(); // TODO!!
        }

        boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> >
        RemoveOverlappingSurfaces( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &MAP)
        {
            boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > new_map( new store::Map< std::string, geom::DirectedPointSurface>());
            for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator itr = MAP->begin(); itr != MAP->end(); ++itr)
            {
                StandardWrite( "remove overlapping surface for: " << itr->first << " ...");
                store::ShPtrVec< DirectedPointSurfaceObject> surfs( itr->second);
                RemoveOverlappingSurface( surfs);
                new_map->InsertNewKeyAndValue( itr->first, DirectedPointSurfaceFromPointSurfaceObjects( surfs));
            }
            return new_map;
        }

        store::ShPtrVec< DirectedPointSurfaceObject> &RemoveOverlappingSurface( store::ShPtrVec< DirectedPointSurfaceObject> &SURFS)
        {
//            std::cout << "  remove overlapping surface ..." << SURFS.size() << std::endl;
            std::vector< std::vector< size_t> > list( ListOverlappingObjects( SURFS));
            if( list.size() != SURFS.size())
            {
                std::cout << "list size: " << list.size() << std::endl;
                std::cout << "objects size: " << SURFS.size() << std::endl;
                std::cout << "exit now" << std::endl;
                exit( -1);
            }
//            std::cout << "list of overlaps: " << std::endl;
//            size_t k( 0);
//            for( std::vector< std::vector< size_t> >::iterator itr( list.begin()); itr != list.end(); ++itr, ++k)
//            {
//                std::cout << "nr points: " << SURFS[k]->GetSurface().GetNrPoints() << std::endl;
//                for( std::vector< size_t>::iterator jtr( itr->begin()); jtr != itr->end(); ++jtr)
//                {
//                    std::cout << *jtr << "  ";
//                }
//                std::cout << std::endl;
//            }
            // iterate though spheres
            size_t i( 0);
            DebugWrite( "iterate through spheres ... ");
            std::vector< std::vector< size_t> >::const_iterator ltr( list.begin());
            for( std::vector< boost::shared_ptr< DirectedPointSurfaceObject> >::const_iterator itr( SURFS.begin()); itr != SURFS.end(); ++itr, ++ltr)
            {
                size_t
//                    cc( 0),
                    deleted(0),
                    kept( 0);
//                    before(( *itr)->GetSurface().GetNrPoints());

#ifdef DEBUG
                if( itr == SURFS.begin())
                {
                    std::cout << "nr points: " <<  ( *itr)->GetSurface().GetNrPoints() << "  " << ( *itr)->DirectedSurfacePoints().size() << std::endl;
                }
#endif

                for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator sptr = ( *itr)->DirectedSurfacePoints().begin(); sptr != ( *itr)->DirectedSurfacePoints().end();)
                {
                    bool found = false;
                    std::vector< size_t>::const_iterator vtr( ltr->begin());
                    while( ( !found) && vtr != ltr->end())
                    {
                        if( SURFS[ *vtr]->IsPointWithin( ( *sptr)->GetPosition()))
                        {
                            ( *itr)->DirectedSurfacePoints().erase( sptr);
                            found = true;
                            ++deleted;
                        }
                        ++vtr;
                    }
                    if( !found)
                    { ++sptr; ++kept;}
                }

//#ifdef DEBUG
//                std::cout << deleted << " deleted + " << kept << " kept = " << deleted + kept << " now nr points: " << ( *itr)->GetSurface().GetNrPoints() << std::endl;
//#endif

                ++i;
            }
//            std::cout << "  ... done" << std::endl;
            return SURFS;
        }

//        boost::shared_ptr< store::Map< std::string, geom::PointSurface> >
//        SmoothSurfaces( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::PointSurfaceObject> > > &MAP, const float &PROBE_RADIUS)
//        {
//            boost::shared_ptr< store::Map< std::string, geom::PointSurface> > new_map( new store::Map< std::string, geom::PointSurface>());
//            for( std::map< std::string, store::ShPtrVec< geom::PointSurfaceObject> >::const_iterator itr = MAP->begin(); itr != MAP->end(); ++itr)
//            {
//                std::cout << " remove overlapping surface for: " << itr->first << " ..." << std::endl;
//                store::ShPtrVec< PointSurfaceObject> surfs( itr->second);
//                SmoothSurface( surfs);
//                std::cout << " ... create point surface from objects ... " << std::endl;
//                new_map->InsertNewKeyAndValue( itr->first, PointSurfaceFromPointSurfaceObjects( surfs));
//                std::cout << " ... done" << std::endl;
//            }
//            return new_map;
//        }


//        store::ShPtrVec< PointSurfaceObject> &SmoothSurface( store::ShPtrVec< PointSurfaceObject> &SURFS, const float &PROBE_RADIUS)
//        {
//            std::vector< std::vector< size_t> > list( ListOverlappingObjects( SURFS, PROBE_RADIUS));
//            for( size_t i = 0; i < list.size(); ++i)
//                for( size_t j = 0; j < list[i].size(); ++j)
//                {
//                    math::Vector3N connection( SURFS( j)->GetPosition() - SURFS( i)->GetPosition());
//                    float distance( connection.Length());
//                    float alpha_max( math::AngleFromLawOfCosinus( SURFS( j)->GetRadius() + PROBE_RADIUS, SURFS( i)->GetRadius() + PROBE_RADIUS, distance));
//
//                    for( std::vector< boost::shared_ptr< SurfacePoint> >::iterator point_itr( SURFS( i)->SurfacePoints().begin()); point_itr != SURFS( i)->SurfacePoints().end(); ++point_itr)
//                    {
//                        math::Vector3N radial( ( *point_itr)->GetPosition() - SURFS( i)->GetPosition());
//                        float angle( math::Angle( radial, connection));
//                        if( angle < alpha_max)
//                        {
//                            math::Vector3N z_axis( math::CrossProduct( connection, radial).Normalize());
//                            math::Vector3N center_connection( );
//                            math::Vector3N probe_center( ( SURFS( i)->GetRadius() + PROBE_RADIUS) * center_connection.NormalizedCopy().RotateAroundAxis( z_axis, beta) + SURFS( i)->GetPosition());
//
//                            if( math::Distance( probe_center,  part of any other object){ erase}
//
//                            float arc( angle * SURFS( i)->GetRadius());
//                            float beta( arc / PROBE_RADIUS);
//                            if( arc < 0.5 * math::Pi * PROBE_RADIUS)
//                            {
//                                ( *point_itr)->SetPosition( ( -PROBE_RADIUS * center_connection.NormalizedCopy()).RotateAroundAxis( z_axis, beta) + probe_center);
//                            }
//                            else
//                            {
//                                SURFS.erase( point_itr);
//                                --point_itr;
//                            }
//                        }
//                    }
//                }
//            return SURFS;
//        }


        PointSurface PointSurfaceFromPointSurfaceObjects( const store::ShPtrVec< geom::PointSurfaceObject> &OBJECTS)
        {
            store::ShPtrVec< geom::SurfacePoint> points;
            for( std::vector< boost::shared_ptr< geom::PointSurfaceObject> >::const_iterator itr( OBJECTS.begin()); itr != OBJECTS.end(); ++itr)
            {
                points += ( *itr)->GetSurfacePoints();
            }
            return PointSurface( points);
        }

        DirectedPointSurface DirectedPointSurfaceFromPointSurfaceObjects( const store::ShPtrVec< geom::DirectedPointSurfaceObject> &OBJECTS)
        {
            store::ShPtrVec< geom::DirectedSurfacePoint> points;
            for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator itr( OBJECTS.begin()); itr != OBJECTS.end(); ++itr)
            {
                points += ( *itr)->GetDirectedSurfacePoints();
            }
            return DirectedPointSurface( points);
        }

//        PointSurface PointSurfaceFromSimpleMolecule( const mol::SimpleMolecule< mol::SurfAtom> &MOL)
//        {
//            store::ShPtrVec< geom::SurfacePoint> points;
//            for( std::vector< boost::shared_ptr< mol::SurfAtom> >::const_iterator itr = MOL.GetAtoms().begin(); itr != MOL.GetAtoms().end(); ++itr)
//            {
//                points += ( *itr)->GetSurf()->GetSurface().GetData();
//            }
//            return PointSurface( points);
//        }

        DirectedPointSurface PointSurfaceFromSimpleMolecule( const mol::SimpleMolecule< mol::SurfAtom> &MOL)
        {
            store::ShPtrVec< geom::DirectedSurfacePoint> points;
            for( std::vector< boost::shared_ptr< mol::SurfAtom> >::const_iterator itr = MOL.GetAtoms().begin(); itr != MOL.GetAtoms().end(); ++itr)
            {
                points += ( *itr)->GetSurf()->GetSurface().GetData();
            }
            return DirectedPointSurface( points);
        }

        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        SurfaceObjectsFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &RESOLUTION,
                const float &FACTOR,
                const float &OFFSET
        )
        {
            boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
                result( new store::ShPtrVec< geom::DirectedPointSurfaceObject>( MOL->GetAtoms().size()));
            std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::iterator
                result_itr( result->begin());

            DebugWrite( __FUNCTION__ << " resolution: " << RESOLUTION << " factor: " << FACTOR << " offset: " << OFFSET);

            for
            (
                    std::vector< boost::shared_ptr< mol::Atom> >::const_iterator itr = MOL->GetAtoms().begin();
                    itr != MOL->GetAtoms().end();
                    ++itr, ++result_itr
            )
            {
                *result_itr = mol::SurfAtom( **itr, RESOLUTION, FACTOR, OFFSET).GetSurf();
            }
            return result;
        }

        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        BuildSmoothSurfacesFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &DELTA,
                const float &FACTOR,
                const float &OFFSET,
                const float &PROBE_RADIUS
        )
        {
            // distance map that is ordered by distance (slick shit)
            std::map< boost::shared_ptr< mol::Atom>,  std::multimap< float, boost::shared_ptr< mol::Atom> > > overlapping_atoms;

            float shift( 2.0 * ( OFFSET + PROBE_RADIUS));
            float distance;
            for
            (
                    std::vector< boost::shared_ptr< mol::Atom> >::const_iterator first = MOL->GetAtoms().begin();
                    first + 1 != MOL->GetAtoms().end(); ++first
            )
            {
                for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator second = first + 1; second != MOL->GetAtoms().end(); ++second)
                {
                    distance = math::Distance( ( *first)->GetPosition(), ( *second)->GetPosition());
                    if( distance < FACTOR * ( ( *first)->GetVanDerWaalsRadius() + ( *second)->GetVanDerWaalsRadius()) + shift)
                    {
                        std::cout << ( *first)->GetAtomID() << "  " << ( *second)->GetAtomID() << " dist: " << distance << std::endl;
                        overlapping_atoms[ *first].insert( std::make_pair( distance, *second));
                        overlapping_atoms[ *second].insert( std::make_pair( distance, *first));
                    }
                }
//                *result_itr = mol::SurfAtom( **itr, RESOLUTION, FACTOR, OFFSET).GetSurf();

            }
            assert( overlapping_atoms.size() == MOL->GetAtoms().size());

            for( std::map< boost::shared_ptr< mol::Atom>,  std::multimap< float, boost::shared_ptr< mol::Atom> > >::const_iterator itr = overlapping_atoms.begin(); itr != overlapping_atoms.end(); ++itr)
            {
                std::cout << itr->first->GetAtomID() << " :  ****" << std::endl;
                for( std::multimap< float, boost::shared_ptr< mol::Atom> >::const_iterator btr = itr->second.begin(); btr != itr->second.end() ; ++btr)
                {
                    std::cout << btr->first << " <=> " << btr->second->GetAtomID() << std::endl;
                }
            }

            // integrate over surface points
            int
                rings( 180 / int( DELTA)),  // Adjust this in case Delta shall be < 1.
                nr_intervalls;
//                count,
            size_t
                snuggle_count;
            float
                theta_1( 0),
                theta_2( 0),
//                s_0( 0),
                delta_phi( 0),
                phi( 0),
                surf( 0),
                delta_radians( DELTA * math::Pi / float( 180.0)),
                diff,
                epsilon( 0.2),
                a,
                small_delta( 5.0 * math::Pi / float( 180.0)),
                alpha;
            bool
                no_overlaps,
                fill;
            std::multimap< float, boost::shared_ptr< mol::Atom> >::const_iterator
                neighbor_itr;
            std::map< boost::shared_ptr< mol::Atom>, std::multimap< float, boost::shared_ptr< mol::Atom> > >::const_iterator
                atom_itr;
            std::map< boost::shared_ptr< mol::Atom>, boost::shared_ptr< geom::DirectedPointSurfaceObject> >
                surf_map;
            math::Vector3N
                factor,
                connection,
                probe_center,
                x_axis,
                y_axis,
                z_axis,
                surf_point,
                surf_conn;
            std::vector< boost::shared_ptr< mol::Atom> >
                snuggle_neighbor;

            for( atom_itr = overlapping_atoms.begin(); atom_itr != overlapping_atoms.end(); ++atom_itr)
            {
                surf_map[ atom_itr->first] =  boost::shared_ptr< DirectedPointSurfaceObject>( new DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Sphere( atom_itr->first->GetPosition(), atom_itr->first->GetVanDerWaalsRadius())), DirectedPointSurface()));
            }
//            DirectedPointSurfaceObject( boost::shared_ptr< Object>( new Sphere( POS, RADIUS)), DirectedPointSurface(
              // the initial surface area
              float
                  surf_0( delta_radians * sin( delta_radians));// * r^2, but r = 1

              for( int i = 0 ; i < rings ; ++i )
              {
                  // calculating delta_theta and delta_phi such that the surface area is as close to surf_0 as possible.
                  theta_1 =  i * delta_radians;
                  theta_2 = (i+1) * delta_radians;
                  delta_phi = surf_0 / ( cos( theta_1 ) - cos( theta_2 )); // /r^2, but r = 1  // could use inverse and spare one division
                  nr_intervalls = int( 2.0 * math::Pi / delta_phi) + 2;
                  delta_phi = 2.0 * math::Pi / float(nr_intervalls);

                  // the actual surface area, assigned to a regarded integration point
                  surf =  delta_phi * ( cos( theta_1) - cos( theta_2) );  // * r^2, but r = 1

                  // testing if a certain point on the surface of a sphere is overlapped by other spheres.
                  for( int j = 0 ; j < nr_intervalls ; j ++ )
                  {
                      phi = j * delta_phi;
                      factor(0) = sin( ( i+0.5) * delta_radians ) * cos( phi);
                      factor(1) = sin( ( i+0.5) * delta_radians ) * sin( phi);
                      factor(2) = cos( ( i+0.5) * delta_radians );
                      //                      s_DefaultSurf->PushBack( boost::shared_ptr< SurfacePoint>( new DirectedSurfacePoint( position, position.NormalizedCopy(), surf)));

                      // iterate through all atoms
                      for( atom_itr = overlapping_atoms.begin(); atom_itr != overlapping_atoms.end(); ++atom_itr)
                      {
                          std::cout << "========================================\n================================== \nID: " << atom_itr->first->GetAtomID() << std::endl;
                          no_overlaps = true;
                          fill = false;
                          neighbor_itr = atom_itr->second.begin();
                          snuggle_count = 0;
                          snuggle_neighbor.clear();

                          while( no_overlaps && neighbor_itr != atom_itr->second.end())
                          {
                              std::cout << "  " << neighbor_itr->second->GetAtomID();
                              probe_center = atom_itr->first->GetPosition() + ( atom_itr->first->GetVanDerWaalsRadius() + PROBE_RADIUS) * factor;
                              diff = math::Distance( probe_center, neighbor_itr->second->GetPosition()) - ( neighbor_itr->second->GetVanDerWaalsRadius() + PROBE_RADIUS);
                              if( diff < -epsilon)
                              {
                                  std::cout << " overlapped ";
                                  no_overlaps = false;
                              }
                              else if( std::abs( diff) < epsilon)
                              {
                                  fill = true;
                                  snuggle_neighbor.push_back( neighbor_itr->second);
                                  ++snuggle_count;
                                  std::cout << "snug";
                              }
                              ++neighbor_itr;
                          }
                          std::cout << std::endl;
                          if( no_overlaps && !fill)
                          {
                              surf_map[ atom_itr->first]->PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>
                              (
                                      new DirectedSurfacePoint
                                      (
                                              atom_itr->first->GetPosition() + atom_itr->first->GetVanDerWaalsRadius() * factor,
                                              factor,
                                              surf * atom_itr->first->GetVanDerWaalsRadius() * atom_itr->first->GetVanDerWaalsRadius()
                                      )
                              ));
                          }
                          else if( no_overlaps && fill)
                          {
                              if( snuggle_count > 1)
                              {
                                  std::cout << "===> more than one neighbor found to snuggle with = to fill with smoothed surface" << std::endl;
                              }
                              for( size_t i = 0; i < snuggle_count; ++i)
                              {
                                  surf_conn = atom_itr->first->GetVanDerWaalsRadius() * factor;
                                  surf_point = atom_itr->first->GetPosition() + surf_conn;
                                  std::cout  << "==============================================\nsurface point: " << surf_point << std::endl;
                                  std::cout << "probe center: " << probe_center << std::endl;
                                  std::cout << "dist probe center/surfacepoint: " << math::Distance( surf_point, probe_center) << std::endl;
                                  std::cout << "dist probe center / snuggling neighbor: " << math::Distance( probe_center, snuggle_neighbor[i]->GetPosition()) << std::endl;
                                  connection = snuggle_neighbor[i]->GetPosition() - atom_itr->first->GetPosition();
                                  std::cout << "dist atoms: " << connection.Length() << " > ? sum radii: " << snuggle_neighbor[i]->GetVanDerWaalsRadius() + atom_itr->first->GetVanDerWaalsRadius() <<  " = " << snuggle_neighbor[i]->GetVanDerWaalsRadius() << " + " << atom_itr->first->GetVanDerWaalsRadius() << std::endl;
                                  alpha = math::Angle( connection, surf_point);
                                  std::cout << "max angle: " << alpha * 180.0 / math::Pi << std::endl;
                                  if( alpha > 0.5 * math::Pi)
                                  {
                                      surf_map[ atom_itr->first]->PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>
                                      (
                                              new DirectedSurfacePoint
                                              (
                                                      atom_itr->first->GetPosition() + surf_conn,
                                                      factor,
                                                      surf * atom_itr->first->GetVanDerWaalsRadius() * atom_itr->first->GetVanDerWaalsRadius()
                                              )
                                      ));
                                      std::cout << "added: " << atom_itr->first->GetPosition() + surf_conn << std::endl;
                                      std::cout << "############  don't smooth here  ###############" << std::endl;
                                      continue;
                                  }
                                  y_axis =  -1.0 * ( surf_conn).Normalize();
                                  z_axis = math::CrossProduct( surf_conn, connection).Normalize();
                                  x_axis = math::CrossProduct( y_axis, z_axis).Normalize();
                                  std::cout << "lengths: " << x_axis.Length() << " " << y_axis.Length() << " " << z_axis.Length() << std::endl;
                                  std::cout << "angles: " << math::Angle( x_axis, y_axis)  * 180.0 / math::Pi << "  " << math::Angle( x_axis, z_axis)  * 180.0 / math::Pi << "  " << math::Angle( y_axis, z_axis)  * 180.0 / math::Pi << std::endl;
                                  for( a = 0.0; a < 0.5*math::Pi - alpha; a += small_delta)
                                  {
                                      math::Vector3N newPoint( probe_center + PROBE_RADIUS * ( sin( a) * x_axis + cos( a) * y_axis));
                                      std::cout << a* 180.0 / math::Pi  << ": " << newPoint << std::endl;
                                      std::cout << "length: " << ( PROBE_RADIUS * ( sin( a) * x_axis + cos( a) * y_axis)).Length() << " dist to surf point: " << math::Distance( surf_point, newPoint) <<  "dist to atom: " << math::Distance( newPoint, atom_itr->first->GetPosition()) << " dist to snuggling: " << math::Distance( newPoint, snuggle_neighbor[i]->GetPosition()) << std::endl;
                                      surf_map[ atom_itr->first]->PushBackDirectedSurfacePoint( boost::shared_ptr< DirectedSurfacePoint>( new DirectedSurfacePoint
                                              (
                                                      probe_center + PROBE_RADIUS * ( sin( a) * x_axis + cos( a) * y_axis),
                                                      -1.0 * y_axis,
                                                      surf * atom_itr->first->GetVanDerWaalsRadius() * atom_itr->first->GetVanDerWaalsRadius()
                                              )
                                      ));
                                  }
                              }
                          }
                      }
                  }
              }

              assert( surf_map.size() == MOL->GetAtoms().size());

              std::map< boost::shared_ptr< mol::Atom>, boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator surf_itr( surf_map.begin());

              boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
                  result( new store::ShPtrVec< geom::DirectedPointSurfaceObject>( MOL->GetAtoms().size()));
              std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::iterator
                  result_itr( result->begin());

              for( ; surf_itr != surf_map.end(); ++surf_itr, ++result_itr)
              {
                  *result_itr = surf_itr->second;
              }


            return result;
        }  // end BuildSmoothSurfaceForSimpleMolecule

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
        IntegrateSmoothSurfacesFromSimpleMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL,
                const float &DELTA,
                const float &FACTOR,
                const float &OFFSET,
                const float &PROBE_RADIUS
        )
        {
            // distance map that is ordered by distance
            std::map< boost::shared_ptr< mol::Atom>,  std::multimap< float, boost::shared_ptr< mol::Atom> > > overlapping_atoms, unique_overlapping_atoms;
//            float shift( 2.0 * ( OFFSET + PROBE_RADIUS));
            float distance, max_distance, r1, r2;

            for
            (
                    std::vector< boost::shared_ptr< mol::Atom> >::const_iterator first = MOL->GetAtoms().begin();
                    first + 1 != MOL->GetAtoms().end();
                    ++first
            )
            {
                for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator second = first + 1; second != MOL->GetAtoms().end(); ++second)
                {
                    distance = math::Distance( ( *first)->GetPosition(), ( *second)->GetPosition());
                    // use distance where smoothed surface reaches a zero thickness as threshold for overlapping criteria
                    r1 = PROBE_RADIUS + ( *first)->GetVanDerWaalsRadius();
                    r2 = PROBE_RADIUS + ( *second)->GetVanDerWaalsRadius();
                    max_distance = std::sqrt( r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cos( math::Pi - asin( PROBE_RADIUS / r1) - asin( PROBE_RADIUS / r2)));
                    if( distance < max_distance)
//                        if( distance < FACTOR * ( ( *first)->GetVanDerWaalsRadius() + ( *second)->GetVanDerWaalsRadius()) + shift)
                    {
                        overlapping_atoms[ *first].insert( std::make_pair( distance, *second));
                        overlapping_atoms[ *second].insert( std::make_pair( distance, *first));
                        unique_overlapping_atoms[ *first].insert( std::make_pair( distance, *second));
                    }
                }

            }
            assert( overlapping_atoms.size() == MOL->GetAtoms().size());

            // integrate over surface points
            int
                rings( 180 / int( DELTA)),  // Adjust this in case Delta shall be < 1.
                nr_intervalls;
//                count,
//                snuggle_count;
            float
                theta_1( 0),
                theta_2( 0),
//                s_0( 0),
                delta_phi( 0),
                phi( 0),
                surf( 0),
                delta_radians( DELTA * math::Pi / float( 180.0)),
//                diff,
                epsilon( 0.2),
//                a,
//                small_delta( 5.0 * math::Pi / float( 180.0)),
//                alpha,
                surf_0,
                square_radius,
                radius,
                two_pi( 2.0 * math::Pi),
                half_pi( 0.5 * math::Pi),
                local_phi,
                dist,
                angle,
                angle_atom,
                angle_neighbor,
                atom_plus_probe_radius,
                neighbor_plus_probe_radius,
//                atom_phi,
                probe_center_diff,
                swing_radius;
            bool
                no_overlaps,
                fill;
            std::multimap< float, boost::shared_ptr< mol::Atom> >::const_iterator
                neighbor_itr;
            std::map< boost::shared_ptr< mol::Atom>, std::multimap< float, boost::shared_ptr< mol::Atom> > >::const_iterator
                atom_itr;
            std::map< boost::shared_ptr< mol::Atom>, boost::shared_ptr< geom::DirectedPointSurfaceObject> >
                surf_map;
            math::Vector3N
                factor,
                connection,
                probe_center,
                x_axis( 1.0, 0.0, 0.0),
                y_axis( 0.0, 1.0, 0.0),
                z_axis( 0.0, 0.0, 1.0),
                surf_point,
                surf_conn,
//                local_x,
//                local_y,
                origin,
                point,
                swing_point,
                swing_axis,
                atom_axis,
                neighbor_axis;
            std::vector< boost::shared_ptr< mol::Atom> >
                snuggle_neighbor;

            // create sphere objects
            for( atom_itr = overlapping_atoms.begin(); atom_itr != overlapping_atoms.end(); ++atom_itr)
            {
                surf_map[ atom_itr->first] =  boost::shared_ptr< DirectedPointSurfaceObject>
                (
                        new DirectedPointSurfaceObject
                            (
                                boost::shared_ptr< Object>
                                (
                                        new Sphere
                                        (
                                                atom_itr->first->GetPosition(),
                                                atom_itr->first->GetVanDerWaalsRadius()
                                        )
                                ),
                                DirectedPointSurface()
                        )
                );
            }

            // build non overlapping point surface for spheres (atoms)
            // iterate through all atoms
            for( atom_itr = overlapping_atoms.begin(); atom_itr != overlapping_atoms.end(); ++atom_itr)
            {
                no_overlaps = true;
                fill = false;
                radius = atom_itr->first->GetVanDerWaalsRadius();
                square_radius = radius * radius;

                if( atom_itr->second.size() > 0)
                {
                    surf_0 = delta_radians * sin( delta_radians) * square_radius;

                    for( int i = 0 ; i < rings ; ++i )
                    {
                        // calculating delta_theta and delta_phi such that the surface area is as close to surf_0 as possible.
                        theta_1 =  i * delta_radians;
                        theta_2 = (i+1) * delta_radians;
                        delta_phi = surf_0 / ( ( cos( theta_1 ) - cos( theta_2 )) * square_radius);
                        nr_intervalls = int( two_pi / delta_phi) + 2;
                        delta_phi = two_pi / float( nr_intervalls);

                        // the actual surface area, assigned to a regarded integration point
                        surf =  delta_phi * ( cos( theta_1) - cos( theta_2)) * square_radius;

                        // testing if a certain point on the surface of a sphere is overlapped by other spheres.
                        for( int j = 0 ; j < nr_intervalls ; ++ j)
                        {
                            phi = j * delta_phi;
                            factor(0) = sin( ( i+0.5) * delta_radians) * cos( phi);
                            factor(1) = sin( ( i+0.5) * delta_radians) * sin( phi);
                            factor(2) = cos( ( i+0.5) * delta_radians);

                            no_overlaps = true;
                            neighbor_itr = atom_itr->second.begin();
                            while( no_overlaps && neighbor_itr != atom_itr->second.end())
                            {
                                probe_center = atom_itr->first->GetPosition() + ( radius + PROBE_RADIUS) * factor;
                                probe_center_diff = math::Distance( probe_center, neighbor_itr->second->GetPosition()) - ( neighbor_itr->second->GetVanDerWaalsRadius() + PROBE_RADIUS);

                                if( probe_center_diff < -epsilon)
                                {
                                    no_overlaps = false;
                                }
                                ++neighbor_itr;
                            }
                            if( no_overlaps)
                            {
                                math::Vector3N local( atom_itr->first->GetPosition() + radius * factor);
//                                          std::cout << "graphics 0 color blue" << std::endl;
//                                          std::cout << "graphics 0 sphere {" << local( 0) << " " << local(1) << " " << local( 2) << "} radius 0.1 resolution 1" << std::endl;
                                surf_map[ atom_itr->first]->PushBackDirectedSurfacePoint
                                (
                                        boost::shared_ptr< DirectedSurfacePoint>
                                        (
                                                new DirectedSurfacePoint
                                                (
                                                        atom_itr->first->GetPosition() + radius * factor,
                                                        factor,
                                                        surf
                                                )
                                        )
                                );
                            } // if( n    o_overlaps)
                        } // for (intervals)
                    } // for( rings)
                } // if( size > 0)
            } // for( atom_itr)


            // torus slice based smoothing of the atom surfaces
            float surf_element( std::numeric_limits< float>::max()); // TODO: try to get this right

            std::map< size_t, std::multimap< size_t, boost::shared_ptr< TorusSlice> > >
                torus_map;
            std::map< size_t, std::multimap< size_t, boost::shared_ptr< TorusSlice> > >::iterator
                torus_map_itr;
            std::map< size_t, boost::shared_ptr< TorusSlice> >::iterator
                torus_itr;
//            std::multimap< size_t, boost::shared_ptr< TorusSlice> >
//                all_tori;

            // loop over unique atom pairs and create torus objects
            for( atom_itr = unique_overlapping_atoms.begin(); atom_itr != unique_overlapping_atoms.end(); ++atom_itr)
            {
                for( neighbor_itr = atom_itr->second.begin(); neighbor_itr != atom_itr->second.end(); ++neighbor_itr)
                {
                    connection = neighbor_itr->second->GetPosition() - atom_itr->first->GetPosition();
                    dist = connection.Length();
                    connection.Normalize();

                    atom_plus_probe_radius = atom_itr->first->GetVanDerWaalsRadius() + PROBE_RADIUS;
                    neighbor_plus_probe_radius = neighbor_itr->second->GetVanDerWaalsRadius() + PROBE_RADIUS;

                    angle_atom = math::AngleFromLawOfCosinus( neighbor_plus_probe_radius, atom_plus_probe_radius, dist);
                    angle_neighbor = math::AngleFromLawOfCosinus( atom_plus_probe_radius, neighbor_plus_probe_radius, dist);

                    origin = atom_itr->first->GetPosition() + ( atom_plus_probe_radius) * cos( angle_atom) * connection;

                    angle_atom = half_pi - angle_atom;
                    angle_neighbor = half_pi - angle_neighbor;


//                    std::cout << "graphics 0 color white" << std::endl;
//                    std::cout << "graphics 0 sphere {" << origin( 0) << " " << origin(1) << " " << origin( 2) << "} radius 0.15 resolution 3" << std::endl;
//                    std::cout << "graphics 0 cylinder {" << neighbor_itr->second->GetPosition()( 0) << " " << neighbor_itr->second->GetPosition()(1) << " " << neighbor_itr->second->GetPosition()( 2) << "} {" << atom_itr->first->GetPosition()( 0) << " " << atom_itr->first->GetPosition()(1) << " " << atom_itr->first->GetPosition()( 2) << "} radius 0.05" << std::endl;

                    atom_axis = ( atom_itr->first->GetPosition() - origin).NormalizedCopy();
                    neighbor_axis = ( neighbor_itr->second->GetPosition() - origin).NormalizedCopy();

                    // distance from connecting axis
                    swing_radius = atom_plus_probe_radius * cos( angle_atom);

                    boost::shared_ptr<  TorusSlice> torus
                    (
                            new TorusSlice
                            (
                                    origin,
                                    atom_axis,
                                    swing_radius,
                                    PROBE_RADIUS,
                                    2.0 * math::Pi - angle_atom,
                                    angle_neighbor,
                                    delta_radians,
                                    delta_radians
                            )
                    );

                    torus_map[ atom_itr->first->GetAtomID()].insert
                    (
                            std::make_pair
                            (
                                    neighbor_itr->second->GetAtomID(),
                                    torus
                            )
                    );

                    torus_map[ neighbor_itr->second->GetAtomID()].insert
                    (
                            std::make_pair
                            (
                                    atom_itr->first->GetAtomID(),
                                    torus
                            )
                    );

//                    all_tori.insert( std::make_pair( atom_itr->first->GetAtomID(), torus_map[ atom_itr->first->GetAtomID()][ neighbor_itr->second->GetAtomID()]));
                }
            }

//            // complete torus map by checking neighbors of neighbors (multimap)
//            for( torus_map_itr = torus_map.begin(); torus_map_itr + 1 != torus_map.end(); ++torus_map_itr)
////                for( torus_map_itr = torus_map_itr + 1; torus_map_itr != torus_map.end(); ++torus_map_itr)
//                {
//                    // check whether there is any overlap between two tori
//                }
//
//            // or iterate once and go from each to neighbors of neighbors
//            for( torus_map_itr = torus_map.begin(); torus_map_itr + 1 != torus_map.end(); ++torus_map_itr)
//            {
//                // the direct neighbors
//                for( std::multimap< size_t, boost::shared_ptr< TorusSlice> >::const_iterator itr = torus_map_itr->second.begin(); itr != torus_map_itr->second.end(); ++itr)
//                {
//                    // the next to next neighbors
//                    for( std::multimap< size_t, boost::shared_ptr< TorusSlice> >::const_iterator next_next = torus_map[ itr->first].begin(); next_next = torus_map[ itr->first].end(); ++next_next)
//                    {
//                        // ensure that they are not the same
//                    }
//                }
//            }

//            std::cout << "THE MAP" << std::endl << torus_map << std::endl;

            // build non overlapping point surfaces for torus
            for( atom_itr = unique_overlapping_atoms.begin(); atom_itr != unique_overlapping_atoms.end(); ++atom_itr)
            {

                torus_map_itr = torus_map.find( atom_itr->first->GetAtomID());

                if( torus_map_itr == torus_map.end())
                {
                    std::cout << "===> torus for atom: " << atom_itr->first->GetAtomID() << " not found" << std::endl;
                }
                else
                {
                    for( neighbor_itr = atom_itr->second.begin(); neighbor_itr != atom_itr->second.end(); ++neighbor_itr)
                    {
                        torus_itr = torus_map_itr->second.find( neighbor_itr->second->GetAtomID());

                        torus_itr->second->CalculateArbitraryXAndYAxis();
//                        torus_itr->second->WriteVmdSurface( std::cout );


// TODO: INTRODUCE SMARTER QUESTIONS FOR SPEED!

                        for( local_phi = 0; local_phi < two_pi; local_phi += torus_itr->second->GetDeltaOuterCircle())
                        {

                            swing_point = torus_itr->second->GetPosition() + torus_itr->second->GetRadius() * ( cos( local_phi) * torus_itr->second->GetArbitraryXAxis() + sin( local_phi) * torus_itr->second->GetArbitraryYAxis());

//                            std::cout << "graphics 0 color white" << std::endl;
//                            std::cout << "graphics 0 sphere {" << swing_point( 0) << " " << swing_point(1) << " " << swing_point( 2) << "} radius 0.05 resolution 1" << std::endl;

                            swing_axis = ( torus_itr->second->GetPosition() - swing_point).NormalizedCopy();

                            for( angle = 0.5 * delta_radians; angle <= 2.0 * math::Pi - torus_itr->second->GetMinimumAngle(); angle += torus_itr->second->GetDeltaInnerRing())
                            {
                                point = swing_point + PROBE_RADIUS * cos( angle) * swing_axis + PROBE_RADIUS * sin( angle) * torus_itr->second->GetNormalVector();
//                                std::cout << "graphics 0 color green" << std::endl;
//                                std::cout << "graphics 0 sphere {" << point( 0) << " " << point(1) << " " << point( 2) << "} radius 0.05 resolution 1" << std::endl;
                                AddToSurfaceMapIfNoOverlapWithOtherAtoms
                                (
                                        PROBE_RADIUS,
                                        point,
                                        ( swing_point - point).NormalizedCopy(),
                                        surf_element,
                                        surf_map[ atom_itr->first],
//                                        all_tori
                                        torus_map_itr->second
                                );
                            }
                            for( angle = 0.5 * delta_radians; angle <= torus_itr->second->GetMaximumAngle(); angle += torus_itr->second->GetDeltaInnerRing())
                            {
                                point = swing_point + PROBE_RADIUS * cos( angle) * swing_axis - PROBE_RADIUS * sin( angle) * torus_itr->second->GetNormalVector();
//                                std::cout << "graphics 0 color yellow" << std::endl;
//                                std::cout << "graphics 0 sphere {" << point( 0) << " " << point(1) << " " << point( 2) << "} radius 0.05 resolution 1" << std::endl;
                                AddToSurfaceMapIfNoOverlapWithOtherAtoms
                                (
                                        PROBE_RADIUS,
                                        point,
                                        ( swing_point - point).NormalizedCopy(),
                                        surf_element,
                                        surf_map[ neighbor_itr->second],
//                                        all_tori
                                        torus_map_itr->second
                                );
                            } // for( inner_circle)
                        } // for( outer_circle)
                    } // for( neigh)
                } // else
            }
              assert( surf_map.size() == MOL->GetAtoms().size());

              std::map< boost::shared_ptr< mol::Atom>, boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator surf_itr( surf_map.begin());

              boost::shared_ptr< store::ShPtrVec< geom::DirectedPointSurfaceObject> >
                  result( new store::ShPtrVec< geom::DirectedPointSurfaceObject>( MOL->GetAtoms().size()));
              std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::iterator
                  result_itr( result->begin());
              for( ; surf_itr != surf_map.end(); ++surf_itr, ++result_itr)
              {
                  *result_itr = surf_itr->second;
              }

            return result;
        }  // end IntegrateSmoothSurfaceForSimpleMolecule

        void AddToSurfaceMapIfNoOverlapWithOtherAtoms
        (
                const float &PROBE_RADIUS,
                const math::Vector3N &POSITION,
                const math::Vector3N &NORMAL,
                const float &SURF_SIZE,
                boost::shared_ptr< geom::DirectedPointSurfaceObject> &SURF,
                const std::multimap< size_t, boost::shared_ptr< TorusSlice> > &TORUS
        )
        {
            bool is_overlapped( false);

            std::multimap< size_t, boost::shared_ptr< TorusSlice> >::const_iterator torus_itr( TORUS.begin());
            std::cout << "nr_tori: " << TORUS.size() << std::endl;
            size_t count( 0);
            while( !is_overlapped && torus_itr != TORUS.end())
            {
                if( IsPointSurroundedByTorusSlice( POSITION, *torus_itr->second, 0.01))
                {
//                    std::cout << "overlapped" << std::endl;
                    is_overlapped = true;
                }
                ++torus_itr;
                ++count;
            }
            std::cout << count << " checked" << std::endl;
            if( !is_overlapped)
            {
                SURF->PushBackDirectedSurfacePoint
                (
                        boost::shared_ptr< DirectedSurfacePoint>
                        (
                                new DirectedSurfacePoint
                                (
                                        POSITION,
                                        NORMAL,
                                        SURF_SIZE
                                )
                        )
                );
            }
        }

        // reads surface objects defined in file into surface map and returns a surface map with only new objects
        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > >
        ReadIntoSurfObjectMap
        (
                std::istream &STREAM,
                boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &SURF
        )
        {
            DebugWrite( __FUNCTION__);
            size_t
                nr_objects;
            STREAM >> nr_objects;
            std::string
                molecule_group,
                surface_object_type;
            boost::shared_ptr< DirectedPointSurfaceObject>
                point_surface_object;
            boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > >
                new_objects( new store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >());
            new_objects->InsertNewKeyAndValue( "all", store::ShPtrVec< geom::DirectedPointSurfaceObject>());

            for( size_t i( 0); i < nr_objects; ++i)
            {
                STREAM >> molecule_group;
                STREAM >> surface_object_type;
                if( surface_object_type == "PointSurfaceSphere")
                {
                    DebugWrite( "PointSurfaceSphere");
                    float x, y, z, r, res;
                    STREAM >> x >> y >> z >> r >> res;
//                    res *= math::Pi / 180.0;
                    point_surface_object = boost::shared_ptr< DirectedPointSurfaceObject>( new PointSurfaceSphere( math::Vector3N( x, y, z), r, res));
                }
                else if( surface_object_type == "PointSurfaceCube")
                {
                    DebugWrite( "PointSurfaceCube found");
                    float x, y, z, xdim, ydim, zdim, res;
                    STREAM >> x >> y >> z >> xdim >> ydim >> zdim >> res;
                    DebugWrite( "PointSurfaceCube values read");
//                    res *= math::Pi / 180.0;
                    point_surface_object = boost::shared_ptr< DirectedPointSurfaceObject>( new PointSurfaceCube( math::Vector3N( x, y, z), xdim, ydim, zdim, res));
                }
                else if( surface_object_type == "PointSurfaceCylinder")
                {
                    DebugWrite( "PointSurfaceCylinder");
                    float x, y, z, radius, height, resolution;
                    STREAM >> x >> y >> z;
                    math::Vector3N pos( x, y, z);
                    STREAM >> x >> y >> z;
                    math::Vector3N direction( x, y, z);
                    STREAM >> radius >> height >> resolution;
                    resolution *= math::Pi / 180.0;
                    point_surface_object
                        = boost::shared_ptr< DirectedPointSurfaceObject>( new PointSurfaceCylinder( pos, direction, radius, height, resolution));
                }

//                if( mystr::AllToLowerCase( molecule_group) == "all")
                if( molecule_group == "all")
                {
                    DebugWrite( "add to all");
                    // push back to all keys of map!
                    std::vector< std::string> keys( SURF->GetKeys());
                    for( std::vector< std::string>::const_iterator itr( keys.begin()); itr != keys.end(); ++itr)
                    {
                        ( *SURF)( *itr).push_back( point_surface_object);
//                        ( *new_objects)( *itr).push_back( point_surface_object);
                    }
                }
                else
                {
                    DebugWrite( "add to " << molecule_group);
                    // if new key: create new key and copy everything from 'all'
                    if( !SURF->IsValidKey( molecule_group))
                    {
                        DebugWrite( "new molecule type");
                        SURF->InsertNewKeyAndValue( molecule_group, ( *SURF)( "all").GetHardCopy()); // HARD COPY !!!!
//                        new_objects->InsertNewKeyAndValue( molecule_group, ( *SURF)( "all").GetHardCopy());
                    }
                    // push back to existing key
                    ( *SURF)( molecule_group).push_back( point_surface_object);
//                    ( *new_objects)( molecule_group).push_back( point_surface_object);
                }
            }
            return new_objects;
        }



        // reads surface objects defined in file into surface map and returns a surface map with only new objects
        boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
        ReadIntoSurfObjectMap
        (
                std::istream &STREAM
        )
        {
            DebugWrite( __FUNCTION__);
            size_t
                nr_objects;
            STREAM >> nr_objects;
            std::string
                molecule_group,
                surface_object_type;
            boost::shared_ptr< Object>
                geom_object;
            boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > >
                new_objects( new store::Map< std::string, store::ShPtrVec< geom::Object> >());
            new_objects->InsertNewKeyAndValue( "all", store::ShPtrVec< geom::Object>());

            for( size_t i( 0); i < nr_objects; ++i)
            {
                STREAM >> molecule_group;
                STREAM >> surface_object_type;
                if( surface_object_type == "Sphere")
                {
                    DebugWrite( "Sphere");
                    float x, y, z, r;
                    STREAM >> x >> y >> z >> r;
//                    res *= math::Pi / 180.0;
                    geom_object = boost::shared_ptr< Object>( new Sphere( math::Vector3N( x, y, z), r));
                }
                else if( surface_object_type == "Cube")
                {
                    DebugWrite( "Cube found");
                    float x, y, z, xdim, ydim, zdim;
                    STREAM >> x >> y >> z >> xdim >> ydim >> zdim;
                    DebugWrite( "Cube values read");
//                    res *= math::Pi / 180.0;
                    geom_object = boost::shared_ptr< Object>( new Cube( math::Vector3N( x, y, z), xdim, ydim, zdim));
                }
                else if( surface_object_type == "Cylinder")
                {
                    DebugWrite( "Cylinder");
                    float x, y, z, radius, height;
                    STREAM >> x >> y >> z;
                    math::Vector3N pos( x, y, z);
                    STREAM >> x >> y >> z;
                    math::Vector3N direction( x, y, z);
                    STREAM >> radius >> height;
                    geom_object
                        = boost::shared_ptr< Object>( new Cylinder( pos, direction, radius, height));
                }

//                if( mystr::AllToLowerCase( molecule_group) == "all")
                if( molecule_group == "all")
                {
                    DebugWrite( "add to all");
                    // push back to all keys of map!
                    std::vector< std::string> keys( new_objects->GetKeys());
                    for( std::vector< std::string>::const_iterator itr( keys.begin()); itr != keys.end(); ++itr)
                    {
                        ( *new_objects)( *itr).push_back( geom_object);
                    }
                }
                else
                {
                    DebugWrite( "add to " << molecule_group);
                    // if new key: create new key and copy everything from 'all'
                    if( !new_objects->IsValidKey( molecule_group))
                    {
                        DebugWrite( "new molecule type");
                        new_objects->InsertNewKeyAndValue( molecule_group, ( *new_objects)( "all").GetHardCopy());
                    }
                    // push back to existing key
                    ( *new_objects)( molecule_group).push_back( geom_object);
                }
            }
            return new_objects;
        }



      bool IsPointInside( const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT, boost::shared_ptr< const math::Vector3N> &POS)
      {
        if( math::Angle( SURF_POINT->GetNormalVector(), SURF_POINT->GetPosition() - *POS) < 0.5 * math::Pi)
          { return false;}
        return true;
      }

      boost::shared_ptr< std::pair< boost::shared_ptr< DirectedPointSurfaceObject>, boost::shared_ptr< DirectedSurfacePoint> > >
      ClosestObjectAndPoint( const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS, const math::Vector3N &POSITION)
      {
          boost::shared_ptr< std::pair< boost::shared_ptr< DirectedPointSurfaceObject>, boost::shared_ptr< DirectedSurfacePoint> > >
          closest_object_and_point( new std::pair< boost::shared_ptr< DirectedPointSurfaceObject>, boost::shared_ptr< DirectedSurfacePoint> >());
          float closest_distance( std::numeric_limits< float>::max());
          float distance;
          for( std::vector< boost::shared_ptr< DirectedPointSurfaceObject> >::const_iterator surf_itr( SURFS.begin()); surf_itr != SURFS.end(); ++surf_itr)
              for( std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator point_itr( ( *surf_itr)->GetDirectedSurfacePoints().begin()); point_itr != ( *surf_itr)->GetDirectedSurfacePoints().end(); ++point_itr)
              {
                  distance = math::Distance( POSITION, ( *point_itr)->GetPosition());
                  if( distance < closest_distance)
                  {
                      closest_distance = distance;
                      closest_object_and_point->second = *point_itr;
                      closest_object_and_point->first = *surf_itr;
                  }
              }
          return closest_object_and_point;
      }


      boost::shared_ptr< std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> > >
      ObjectDistanceMap( const store::ShPtrVec< DirectedPointSurfaceObject> &SURFS, const math::Vector3N &POSITION)
      {
          boost::shared_ptr< std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> > >
              distance_map( new std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> >());
          for( std::vector< boost::shared_ptr< DirectedPointSurfaceObject> >::const_iterator surf_itr( SURFS.begin()); surf_itr != SURFS.end(); ++surf_itr)
          {
              distance_map->insert( std::make_pair( math::Distance( POSITION, ( *surf_itr)->GetPosition()), *surf_itr));
          }
          return distance_map;
      }

      boost::shared_ptr< geom::DirectedSurfacePoint>
      ClosestSurfacePoint
      (
              const boost::shared_ptr< std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> > > MAP,
              const math::Vector3N &POS
      )
      {
          boost::shared_ptr< geom::DirectedSurfacePoint>
              closest;
          float
              dist,
              closest_distance( std::numeric_limits< float>::max()),
              float_max_radius( 4.0),
              max_dist( std::numeric_limits< float>::max());
          bool
              found_first_non_empty_object( false);
          std::multimap< float, boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator
              map_itr( MAP->begin());
          do
          {
              if( map_itr->second->GetDirectedSurfacePoints().size() > 0)
              {
                  if( !found_first_non_empty_object)
                  {
                      found_first_non_empty_object = true;
                      max_dist = map_itr->first + float_max_radius;
                  }
                  for
                  (
                          std::vector< boost::shared_ptr< geom::DirectedSurfacePoint> >::const_iterator surf_itr( map_itr->second->GetDirectedSurfacePoints().begin());
                          surf_itr != map_itr->second->GetDirectedSurfacePoints().end();
                          ++surf_itr
                  )
                  {
                      dist = math::Distance( POS, ( *surf_itr)->GetPosition());
                      if( dist < closest_distance)
                      {
                          closest_distance = dist;
                          closest = *surf_itr;
                      }
                  }
              }
              ++map_itr;
          } while( map_itr->first <= max_dist && map_itr != MAP->end());
          return closest;
      }


      std::vector< std::vector< size_t> >
      FuseClusters
      (
        std::vector< std::vector< size_t> > &CLUSTERS,
        const std::vector< size_t> &TO_FUSE
      )
      {
          size_t count( 0);
          std::vector< std::vector< size_t> > tmp;
          for( std::vector< std::vector< size_t> >::const_iterator itr = CLUSTERS.begin(); itr != CLUSTERS.end(); ++itr, ++count)
          {
              if( IsElement( TO_FUSE, count))
              {
                  if( count == TO_FUSE[ 0])
                  {
                      tmp.push_back( *itr);
                  }
                  else
                  {
                      tmp[ TO_FUSE[ 0]].insert( tmp[ TO_FUSE[ 0]].end(), CLUSTERS[ count].begin(), CLUSTERS[ count].end());
                  }
              }
              else
              {
                  tmp.push_back( *itr);
              }
          }
          return tmp;
      }


      bool
      IsElementOfCluster
      (
              const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT,
              const std::vector< size_t> &CLUSTER_IDS,
              const geom::DirectedPointSurface &SURF,
              const float &CLUSTER_THRESHOLD
      )
      {
          for( std::vector< size_t>::const_reverse_iterator ritr = CLUSTER_IDS.rbegin(); ritr != CLUSTER_IDS.rend(); ++ritr)
          {
//              DebugWrite( __FUNCTION__ << " cluster: " << *ritr << " threshold: " << CLUSTER_THRESHOLD << " surf-size: " << SURF.GetData().size());
              if( math::Distance( SURF_POINT->GetPosition(), SURF.GetData()( *ritr)->GetPosition()) < CLUSTER_THRESHOLD)
              {
                  return true;
              }
          }
          return false;
      }


      std::vector< size_t>
      ElementOfExistingClusters
      (
              const boost::shared_ptr< DirectedSurfacePoint> &SURF_POINT,
              const std::vector< std::vector< size_t> >    &CLUSTER_IDS,
              const geom::DirectedPointSurface &SURF,
              const float &CLUSTER_THRESHOLD
      )
      {
//          DebugWrite( __FUNCTION__);
          std::vector< size_t> result;
          size_t count( 0);
          for( std::vector< std::vector< size_t> >::const_iterator itr = CLUSTER_IDS.begin(); itr != CLUSTER_IDS.end(); ++itr, ++count)
          {
              if( IsElementOfCluster( SURF_POINT, *itr, SURF, CLUSTER_THRESHOLD))
              {
//                  DebugWrite( __FUNCTION__ << ": element of " << count);
                  result.push_back( count);
              }
          }
          return result;
      }


      void
      FilterSurfArtefacts
      (
              boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURF_MAP,
              const size_t &CLUSTER_MIN_SIZE,
              const float &CLUSTER_THRESHOLD
      )
      {
          DebugWrite( __FUNCTION__);
          for( std::map< std::string, geom::DirectedPointSurface>::iterator type_itr( SURF_MAP->begin()); type_itr != SURF_MAP->end(); ++type_itr)
          {
              std::vector< std::vector< size_t> >    clusters;
//                std::vector< std::vector< std::vector< boost::shared_ptr< DirectedSurfacePoint> >::iterator> >    clusters;

              DebugWrite( "mol_type: " << type_itr->first << " elements: " << type_itr->second.GetData().size());

              size_t count( 0);
              for
              (
                      std::vector< boost::shared_ptr< DirectedSurfacePoint> >::const_iterator surf_itr( type_itr->second.GetData().begin());
                      surf_itr != type_itr->second.GetData().end();
                      ++surf_itr, ++count
              )
              {
                  std::vector< size_t> cluster_ids( ElementOfExistingClusters( *surf_itr, clusters, type_itr->second, CLUSTER_THRESHOLD));

                  if( cluster_ids.size() == 1) // add to existing cluster
                  {
//                      DebugWrite( "add to existing cluster");
                      clusters[ cluster_ids[ 0]].push_back( count);
                  }
                  else if( cluster_ids.size() > 1) // fuse cluster  ///////////////////////////////////// add count !!!!!!!!!!!!!!!!!!!!1
                  {
                      clusters[ cluster_ids[ 0]].push_back( count);
                      clusters = FuseClusters( clusters, cluster_ids);
                      DebugWrite( "fuse clusters, now: " << clusters.size() << " clusters");
                  }
                  else // add new cluster and add element to new cluster
                  {
                      clusters.push_back( std::vector< size_t>());
                      clusters.back().push_back( count);
                      DebugWrite( "add new cluster and add element to it, now: " << clusters.size() << " clusters");
                  }
              }
//              DebugWrite( "the clusters: " << clusters);
#ifdef DEBUG
              count = 0;
              for( size_t i = 0; i < clusters.size(); ++i)
              {
                  std::cout << "cluster: " << i << " has " << clusters[i].size() << " elements" << std::endl;
                  count += clusters[i].size();
              }
              std::cout << count << " elements in total" << std::endl;
              count = 0;
              for( std::vector< std::vector< size_t> >::const_iterator cluster_itr = clusters.begin() + 1; cluster_itr != clusters.end(); ++cluster_itr, ++count)
              {
                  //if( cluster_itr->size() < CLUSTER_MIN_SIZE)
                  std::cout << "graphics mol new" << std::endl;
                  std::cout << "graphics " << count << " color red" <<  std::endl;
                  for( std::vector< size_t>::const_iterator id_itr = cluster_itr->begin(); id_itr != cluster_itr->end(); ++id_itr)
                  {
                      // write vmd command;
                      math::Vector3N pos( type_itr->second.GetData()( *id_itr)->GetPosition());
                      std::cout << "graphics " << count << " sphere {" << pos( 0) << " " << pos( 1) << " " << pos( 2) << "} radius 0.05 resolution 3" << std::endl;
                  }
              }
#endif
              geom::DirectedPointSurface reliable;
              std::cout << "adding ";
              count = 0;
              for( std::vector< std::vector< size_t> >::const_iterator cluster_itr = clusters.begin(); cluster_itr != clusters.end(); ++cluster_itr, ++count)
              {
                  if( cluster_itr->size() > CLUSTER_MIN_SIZE)
                      for( std::vector< size_t>::const_iterator id_itr = cluster_itr->begin(); id_itr != cluster_itr->end(); ++id_itr)
                      {
//                          std::cout << *id_itr <<  "  " ;
                          reliable.PushBack( type_itr->second.GetData()( *id_itr));
                          ++count;
                      }
              }
              type_itr->second = reliable;
              std::cout << " total: " << count << std::endl;
          }
      }

        void
        BuildSurface
        (
                const CommandLineManager &CMD,
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE,
                boost::shared_ptr< store::Map< std::string, geom::DirectedPointSurface> > &SURF_MAP,
                boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > > &SURF_OBJECT_MAP
        )
        {
            std::ifstream
                read;
            std::ofstream
                write;
            boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> > >
                user_defined_surface_objects( new store::Map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >());
            float
                surface_factor = 1.0,
                surface_offset = 1.4,
                surface_resolution = 30.0,
                probe_radius = 1.4;


            if( CMD.IsFlagSet( "magnify_radii"))
            {
                std::vector< float> vec( CMD.GetArgumentsForFlag< float>( "magnify_radii"));
                surface_offset = vec[ 0];
                if( vec.size() > 1)
                {
                    surface_factor = vec[ 1];
                }
            }
            StandardWrite( "vdw radii offset: " << surface_offset);
            StandardWrite( "vdw radii factor: " << surface_factor);

            if( CMD.IsFlagSet( "surface_resolution"))
            {
                surface_resolution = CMD.GetArgumentForFlag< float>( "surface_resolution");
            }
            StandardWrite( "surface resolution: " << surface_resolution);

            // read from file
            if( CMD.IsFlagSet( "read_surface_from_file"))
            {
                StandardWrite(  "from file ...");
                Open( read, CMD.GetArgumentStringForFlag( "read_surface_from_file"));
                read >> SURF_OBJECT_MAP;
                Close( read);

            }
            // or build from pdb
            else if( CMD.IsFlagSet( "smooth_surf"))
            {
                StandardWrite(  "from scratch but smooth ...");
                std::vector< float> values( CMD.GetArgumentsForFlag< float>( "smooth_surf"));
                if( values.size() == 1)
                {
                    probe_radius = values[0];
                }
                SURF_OBJECT_MAP->InsertNewKeyAndValue( "all", *geom::factory::IntegrateSmoothSurfacesFromSimpleMolecule( IMPLICIT_MOLECULE, surface_resolution, surface_factor, surface_offset, probe_radius));
            }
            else
            {
                StandardWrite(  "from scratch ...");
                SURF_OBJECT_MAP->InsertNewKeyAndValue( "all", *geom::factory::SurfaceObjectsFromSimpleMolecule( IMPLICIT_MOLECULE, surface_resolution, surface_factor, surface_offset));
            }

            // add user defined surface objects
            if( CMD.IsFlagSet( "add_surface_objects"))
            {
                StandardWrite(  "add surface objects ...");
    #ifdef DEBUG
                for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator surf_itr( SURF_OBJECT_MAP->begin()); surf_itr != SURF_OBJECT_MAP->end(); ++surf_itr)
                {
                    size_t nr_points( 0);
                    for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator obj_itr = surf_itr->second.begin(); obj_itr != surf_itr->second.end(); ++obj_itr)
                    {
                        nr_points += ( *obj_itr)->NumberSurfacePoints();
                    }
                    std::cout << "mol type: " << surf_itr->first << " nr surf points: " << nr_points << std::endl;
                }
    #endif
                // if string != "all": copy surface object for new type, add new object to new surface object
                Open( read, CMD.GetArgumentStringForFlag( "add_surface_objects"));
                user_defined_surface_objects = geom::factory::ReadIntoSurfObjectMap( read, SURF_OBJECT_MAP);
                Close( read);

            }

    #ifdef DEBUG
            for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator surf_itr( SURF_OBJECT_MAP->begin()); surf_itr != SURF_OBJECT_MAP->end(); ++surf_itr)
            {
                size_t nr_points( 0);
                for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator obj_itr = surf_itr->second.begin(); obj_itr != surf_itr->second.end(); ++obj_itr)
                {
                    nr_points += ( *obj_itr)->NumberSurfacePoints();
                }
                std::cout << "mol type: " << surf_itr->first << " nr surf points: " << nr_points << std::endl;
            }
    #endif

//            StandardWrite(  "remove overlapping surface ...");
            SURF_MAP = geom::factory::RemoveOverlappingSurfaces( SURF_OBJECT_MAP);

    #ifdef DEBUG
            for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator surf_itr( SURF_OBJECT_MAP->begin()); surf_itr != SURF_OBJECT_MAP->end(); ++surf_itr)
            {
                size_t nr_points( 0);
                for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator obj_itr = surf_itr->second.begin(); obj_itr != surf_itr->second.end(); ++obj_itr)
                {
                    nr_points += ( *obj_itr)->NumberSurfacePoints();
                }
                std::cout << "mol type: " << surf_itr->first << " nr surf points: " << nr_points << std::endl;
            }
    #endif

            if( CMD.IsFlagSet( "filter_surf_artefacts"))
            {
//                std::time( &now);
//                std::cout << "surface_building_sofar: " << std::difftime( now, begin) << "s" << std::endl;
//                std::time( &begin);

                if( CMD.IsFlagSet( "write_surface"))
                {
                    StandardWrite(  "write surface to file ...");
                    Open( write, CMD.GetArgumentStringForFlag( "write_surface") + ".before");
                    write << SURF_MAP;
                    Close( write);

                }
                std::vector< float> arguments( CMD.GetArgumentsForFlag< float>( "filter_surf_artefacts"));
                size_t cluster_min_size = size_t( arguments[ 0]);
                float threshold( arguments[ 1]);
                StandardWrite( "filter surf_artefacts: threshold: " << threshold << " min size: " << cluster_min_size);
                geom::factory::FilterSurfArtefacts( SURF_MAP, cluster_min_size, threshold);

//                std::time( &now);
//                std::cout << "surface_filtering: " << std::difftime( now, begin) << "s" << std::endl;
//                std::time( &begin);
            }

            if( CMD.IsFlagSet( "write_surface"))
            {
                StandardWrite(  "write surface to file ...");
                Open( write, CMD.GetArgumentStringForFlag( "write_surface"));
                write << SURF_MAP;
                Close( write);

            }
            if( CMD.IsFlagSet( "write_surf_as_pdb"))
              {
                StandardWrite( "write surface as pdb ...");
                Open( write, CMD.GetArgumentStringForFlag( "write_surf_as_pdb"));
                write << "COMMENT molecule surface represented by points written as ATOMs by membrane_carver.exe" << std::endl;
                size_t mol_count( 0);
                for( std::map< std::string, geom::DirectedPointSurface>::const_iterator itr = SURF_MAP->begin(); itr != SURF_MAP->end(); ++itr, ++mol_count)
                  {
                itr->second.WriteAsPdb( write, mol_count);
                  }
                Close( write);
              }
            if( CMD.IsFlagSet( "write_vmd_commands_for_surface") ||  CMD.IsFlagSet( "write_vmd_commands_for_surface_objects"))
            {
                size_t mol_count( 0);
                StandardWrite(  "write vmd commands for surface ...");

                std::vector< std::string> colors;
                colors.push_back( "white");
                colors.push_back( "green");
                colors.push_back( "red");
                colors.push_back( "yellow");
                colors.push_back( "magenta");
                colors.push_back( "blue");
                colors.push_back( "cyan");

                if( CMD.IsFlagSet( "write_vmd_commands_for_surface"))
                {
                    Open( write, CMD.GetArgumentStringForFlag( "write_vmd_commands_for_surface"));
                    for( std::map< std::string, geom::DirectedPointSurface>::const_iterator itr = SURF_MAP->begin(); itr != SURF_MAP->end(); ++itr, ++mol_count)
                    {
                        write << "mol new" << std::endl;
                        write << "graphics " << mol_count << " color " << colors[ mol_count] << std::endl;
                        StandardWriteNoFlush( "write vmd commands for: ");
                        StandardWrite( itr->first);
                        itr->second.WriteVmdCommands( write, mol_count);
                    }

                    Close( write);
                }
                if( CMD.IsFlagSet( "write_vmd_commands_for_surface_objects"))
                {
                    Open( write, CMD.GetArgumentStringForFlag( "write_vmd_commands_for_surface_objects"));
                    for( std::map< std::string, store::ShPtrVec< geom::DirectedPointSurfaceObject> >::const_iterator itr = SURF_OBJECT_MAP->begin(); itr != SURF_OBJECT_MAP->end(); ++itr, ++mol_count)
                    {
                        write << "mol new" << std::endl;
                        write << "graphics " << mol_count << " color " << colors[ mol_count] << std::endl;
                        StandardWriteNoFlush( "write vmd commands for: ");
                        StandardWrite( itr->first);
                        for( std::vector< boost::shared_ptr< geom::DirectedPointSurfaceObject> >::const_iterator obj_itr = itr->second.begin(); obj_itr != itr->second.end(); ++obj_itr)
                        {
                            ( *obj_itr)->GetObject()->WriteVmdCommands( write, mol_count);
                        }
                    }

                    Close( write);
                }
            }
        }



    } // end namespace factory
} // end namespace geom
