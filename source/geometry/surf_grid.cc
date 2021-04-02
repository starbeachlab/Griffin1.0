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


#include "../../include/geometry/surf_grid.h"

extern float g_PDBFactor;

namespace geom
{

    void SurfGrid::AnalyzeObjectMap(  const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS) const
    {
        DebugWrite( __FUNCTION__);
        for( std::map< std::string, store::ShPtrVec< geom::Object> >::const_iterator type_itr = OBJECTS->begin(); type_itr != OBJECTS->end(); ++type_itr)
        {
            if( type_itr->first != "all")
            {
                m_MolTypes.push_back( type_itr->first);
                //m_SurfIndices.push_back( std::vector< store::Vector3N< short> >());
                m_SurfIndices.push_back( std::multimap< short, std::pair< short, short> >());
            }
        }
    }

    void SurfGrid::AdjustLimitsForSurfaceObjects( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS)
    {
        AnalyzeObjectMap( OBJECTS);

        for( std::map< std::string, store::ShPtrVec< geom::Object> >::const_iterator type_itr = OBJECTS->begin(); type_itr != OBJECTS->end(); ++type_itr)
        {
            for( std::vector< boost::shared_ptr< geom::Object> >::const_iterator obj_itr = type_itr->second.begin(); obj_itr != type_itr->second.end(); ++obj_itr)
            {
                AdjustLimitsForSurfaceObject( *obj_itr);
            }
        }
    }


    void SurfGrid::AddSurfaceObjects( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS)
    {
      //        AnalyzeObjectMap( OBJECTS);

        int mol_type_id;

        for( std::map< std::string, store::ShPtrVec< geom::Object> >::const_iterator type_itr = OBJECTS->begin(); type_itr != OBJECTS->end(); ++type_itr)
        {
            mol_type_id = std::find( m_MolTypes.begin(), m_MolTypes.end(), type_itr->first) - m_MolTypes.begin();
            for( std::vector< boost::shared_ptr< geom::Object> >::const_iterator obj_itr = type_itr->second.begin(); obj_itr != type_itr->second.end(); ++obj_itr)
            {
                AddSurfaceObject( *obj_itr, IndexToKey( mol_type_id));
            }
        }
    }



  void SurfGrid::AdjustLimitsForSurfaceObject( const boost::shared_ptr< geom::Object> &OBJECT)
    {
        StandardWrite( __FUNCTION__);
        store::Limits3D
            limits = OBJECT->BoxLimits();

        std::vector< float> min( limits.GetMin()), max( limits.GetMax());

        int new_index;

        if
          (
           min[0] - 0.5 < m_Min[0]
           )
          {
            m_NrElements[0] += int( ( m_Min[0] - min[0] + 0.5) / m_Delta[0]) + 1;
            StandardWrite( "=> minimum value of grid in x changed from " << m_Min[0] << " to " << min[0] - 0.5 << std::endl);
            m_Min[0] = min[0] - 0.5;
          }
        if
          ( 
           min[1] - 0.5 < m_Min[1]
           )
          {
            m_NrElements[1] += int( ( m_Min[1] - min[1] + 0.5) / m_Delta[1]) + 1;
            StandardWrite( "=> minimum value of grid in y changed from " << m_Min[1] << " to " << min[1] - 0.5 << std::endl);
            m_Min[1] = min[1] - 0.5;
          }
        if
          ( 
           min[2] - 0.5 < m_Min[2]
           )
          {
            m_NrElements[2] += int( ( m_Min[2] - min[2] + 0.5) / m_Delta[2]) + 1;
            StandardWrite( "=> minimum value of grid in z changed from " << m_Min[2] << " to " << min[2] - 0.5 << std::endl);
            m_Min[2] = min[2] - 0.5;
          }
        if
          (
           max[0] + 0.5 > m_Min[0] + m_NrElements[0] * m_Delta[0]
           )
          {
            new_index = int( (max[0] + 0.5 - m_Min[0]) / m_Delta[0]) + 1;
            StandardWrite( "=> maximum value of grid in x changed from " << m_Min[0] + m_NrElements[0] * m_Delta[0] << " to " << m_Min[0] + new_index * m_Delta[0] << " (" << new_index << ")(" << m_NrElements[0] << ") to match: " << max[0] + 0.5 << " value of added object");
            m_NrElements[0] = new_index;
            //            max[0] = m_Min[0] + new_index * m_Delta[0];
          }
        if
          (
           max[1] + 0.5 >  m_Min[1] + m_NrElements[1] * m_Delta[1]
           )
          {
            new_index = int( (max[1] + 0.5 - m_Min[1]) / m_Delta[1]) + 1;
            StandardWrite( "=> maximum value of grid in y changed from " << m_Min[1] + m_NrElements[1] * m_Delta[1] << " to " << m_Min[1] + new_index * m_Delta[1] << " (" << new_index << ")(" << m_NrElements[1] << ") to match: " << max[1] + 0.5 << " value of added object");
            m_NrElements[1] = new_index;
          }
        if
          (
           max[2] + 0.5 > m_Min[2] + m_NrElements[2] * m_Delta[2]
           )
          {
            new_index = int( (max[2] + 0.5 - m_Min[2]) / m_Delta[2]) + 1;
            StandardWrite( "=> maximum value of grid in z changed from " << m_Min[2] + m_NrElements[2] * m_Delta[2] << " to " << m_Min[2] + new_index * m_Delta[2] << " (" << new_index << ")(" << m_NrElements[2] << ") to match: " << max[2] + 0.5 << " value of added object");
            m_NrElements[2] = new_index;
          }
    }


    void SurfGrid::AddSurfaceObject( const boost::shared_ptr< geom::Object> &OBJECT, const short &MOL_TYPE_KEY)
    {
        StandardWrite( __FUNCTION__);
        DebugWrite( "mol_type_key: " << MOL_TYPE_KEY << " object: " << OBJECT);
        store::Limits3D
            limits = OBJECT->BoxLimits();
        store::Vector3N< int>
            min( math::TupleXD< short, 3>::CalcIndices( limits.GetMin())),
            max( math::TupleXD< short, 3>::CalcIndices( limits.GetMax()));
        DebugWrite( "object_limits: " << limits);
        DebugWrite( "min_indices: " << min[0] << " " << min[1] << " " << min[2]);
        DebugWrite( "max_indices: " << max[0] << " " << max[1] << " " << max[2]);

        if
        (
        		min[0] < 0
        		|| min[1] < 0
        		|| min[2] < 0
        		|| max[0] > m_NrElements[0] - 1
        		|| max[1] > m_NrElements[1] - 1
        		|| max[2] > m_NrElements[2] - 1
        )
        {
            std::cout << "===> object out of limits of surf grid - not added!\nIndices:" << std::endl;
            std::cout << min[0] << " should not be smaller than 0" << std::endl;
            std::cout <<  min[1] << " should not be smaller than 0" << std::endl;
            std::cout <<  min[2] << " should not be smaller than 0" << std::endl;
            std::cout <<  max[0] << " should not be larger than  " << m_NrElements[0] - 1 << std::endl;
            std::cout <<  max[1] << " should not be larger than  " << m_NrElements[1] - 1 << std::endl;
            std::cout <<  max[2] << " should not be larger than  " << m_NrElements[2] - 1 << std::endl;
            std::cout << "mol-type-key: " << MOL_TYPE_KEY << std::endl;
            std::cout << "object: " << *OBJECT;
            return;
        }

        math::Vector3N
            pos;
        short
            value;
        std::vector< int> keys;

        for( int i = min[0]; i <= int( max[0]); ++i)
            for( int j = min[1]; j <= int( max[1]); ++j)
                for( int k = min[2]; k <= int( max[2]); ++k)
                {
                    pos = m_Min;
                    pos[0] += (i + 0.5) * m_Delta[ 0];
                    pos[1] += (j + 0.5) * m_Delta[ 1];
                    pos[2] += (k + 0.5) * m_Delta[ 2];
                    DebugWrite( "pos: " << pos);
                    if( OBJECT->IsPointWithin( pos))
                    {
                        value = operator()( store::Vector3N< int>( i, j, k));
                        DebugWrite( "inside, value: " << value);
                        // if point is already molecule dependently inside add this key
                        if( value < 0)
                        {
                          keys = KeyToIndices( value);
                          DebugWrite( "keys: " << keys[0] << " " << keys[1]	<< " " << keys[2]);
                          if( find( keys.begin(), keys.end(), MOL_TYPE_KEY) == keys.end())
                            {
                        	  DebugWrite( "found " << MOL_TYPE_KEY);
                              operator()( store::Vector3N< int>( i, j, k))  +=  MOL_TYPE_KEY;
                            }
                        }
                        // if outside then set to key of molecule
                        else if( value > 0)
                        {
                        	DebugWrite( "fresh");
                          operator()( store::Vector3N< int>( i, j, k))    =  MOL_TYPE_KEY;
                        }
                        // nothing changes for value == 0, means inside for 'all'
                    }
                }
    }


    void SurfGrid::BuildFromMolecule
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            const bool &SMOOTH,
            const float &PROBE_RADIUS
    )
    {
///        BuildSurfGrid
//        (
//                MOLECULE,
//                m_Data,
//                m_Min,
//                m_NrElements,
//                m_Delta[0],
//                SMOOTH,
//                PROBE_RADIUS
//        );
    }
  
    void
    SurfGrid::Prepare
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
            float *X_COO,
            float *Y_COO,
            float *Z_COO,
            float *VDW_RADII,
            std::vector< float> &CENTER,
            const float &FACTOR,
            const float &OFFSET,
	    const float &MARGIN
    )
    {
        StandardWrite( __FUNCTION__);
        std::vector< float> 
	    max( 3, -std::numeric_limits< float>::max());
        //    std::cout << "max: ";
        //    std::copy( max.begin(), max.end(), std::ostream_iterator< float>( std::cout, "  "));
        //    std::cout << std::endl;
        //    std::cout << std::endl;

        PrepareSurfGrid
        (
                MOLECULE,
                m_Min,
                max,
                m_NrElements,
                m_Delta[0],
                X_COO,
                Y_COO,
                Z_COO,
                VDW_RADII,
                CENTER,
                FACTOR,
                OFFSET, 
		MARGIN
        );
        for( int i = 0; i < 3; ++i)
        {
            m_NrElements[i] = int( (max[i] - m_Min[i]) / m_Delta[i]) + 1;
        }
    }

        void
        SurfGrid::Build
        (
                const bool &SMOOTH,
                const float &PROBE_RADIUS,
                const size_t &SIZE,
                const float *X_COO,
               const float *Y_COO,
               const float *Z_COO,
               const float *VDW_RADII,
               const std::vector< float> &CENTER
        )
        {

            BuildSurfGrid
            (
                    m_Data,
                    m_Min,
                    m_NrElements,
                    m_Delta[0],
                    SMOOTH,
                    PROBE_RADIUS,
                    SIZE,
                    X_COO,
                    Y_COO,
                    Z_COO,
                    VDW_RADII,
                    CENTER
            );
        }



    int SurfGrid::IndexToKey( const int &ID) const
    {
        switch( ID)
        {
        case 0:  return  0;  // all:
        case 1:  return -1;  // moltype 1
        case 2:  return -2;  // moltype 2
        case 3:  return -4;  // moltype 2
        default:
            std::cerr << __PRETTY_FUNCTION__ << " only three molecules types possible at this point (beside 'all')" << std::endl;
        }
        return std::numeric_limits< int>::max();
    }

    //! returns the vector indices of a grid value below zero
    std::vector< int> SurfGrid::KeyToIndices( const int &KEY) const
    {
        // chmod like key combinatorics, disected into vector
        std::vector< int> ids;
        switch( KEY)
        {
        case 1:
        	std::cout << "====> CHECK case 1: " << __FILE__ << " " << __LINE__ << std::endl;
//            ids.push_back( 1);
            break;
        case 0:
        	std::cout << "====> CHECK case 0: " << __FILE__ << " " << __LINE__ << std::endl;
//            ids.push_back( 0);
            break;
        case -1:
            ids.push_back( 1);
            break;
        case -2:
            ids.push_back( 2);
            break;
        case -3:
            ids.push_back( 1);
            ids.push_back( 2);
            break;
        case -4:
            ids.push_back( 3);
            break;
        case -5:
            ids.push_back( 1);
            ids.push_back( 3);
            break;
        case -6:
            ids.push_back( 2);
            ids.push_back( 3);
            break;
        case -7:
            ids.push_back( 1);
            ids.push_back( 2);
            ids.push_back( 3);
            break;
        case -8:
            ids.push_back( 4);
            break;
        case -9:
            ids.push_back( 1);
            ids.push_back( 4);
            break;
        case -10:
            ids.push_back( 2);
            ids.push_back( 4);
            break;
        case -11:
            ids.push_back( 1);
            ids.push_back( 2);
            ids.push_back( 4);
            break;
        case -12:
            ids.push_back( 3);
            ids.push_back( 4);
            break;
        case -13:
            ids.push_back( 1);
            ids.push_back( 3);
            ids.push_back( 4);
            break;
        case -14:
            ids.push_back( 2);
            ids.push_back( 3);
            ids.push_back( 4);
            break;
        case -15:
            ids.push_back( 1);
            ids.push_back( 2);
            ids.push_back( 3);
            ids.push_back( 4);
            break;
        default:
            std::cerr << "====> " << __PRETTY_FUNCTION__ << " only three molecules types possible at this point (beside 'all')" << std::endl;
        }
        return ids;
    }


    //std::vector< store::Vector3N< short> >
    void
    SurfGrid::IndicesOfSurfaceGridPoints() const
    {
        StandardWrite( __FUNCTION__);

        store::Vector3N< int> sizes( GetNrElements());

        int
            shell = 1,
            yz_size = sizes[ 1] * sizes[ 2],
            imin, imax,
            jmin, jmax,
            kmin, kmax;
        short
            neighbor_value,
            value;
        bool
            surf;

        // iterate through grid
        for( int i = 0; i < int( sizes[0]); ++i)
            for( int j = 0; j < int( sizes[1]); ++j)
                for( int k = 0; k < int( sizes[2]); ++k)
                {
//                      std::cout << i << "  " << j << "  " << k << "  \n";
//                      std::cout . flush();

                      // if point is inside (value < 1) investigate one shell of points
                    value = math::TupleXD< short, 3>::m_Data[ i * yz_size + j * sizes[ 2] + k];

                    // if point is outside
                    if( value > short( 0))
                    {
                        continue;
                    }


                  // adjust limits for points that are at the surf of the grid
                  if( i == 0)
                  {
                      imin = 1;
                  }
                  else
                  {
                      imin = i - shell;
                  }
                  if( j == 0)
                  {
                      jmin = 1;
                  }
                  else
                  {
                      jmin = j - shell;
                  }
                  if( k == 0)
                  {
                      kmin = 1;
                  }
                  else
                  {
                      kmin = k - shell;
                  }
                  if( i == int( sizes[0]) - 1)
                  {
                      imax = int( sizes[0]) - 2;
                  }
                  else
                  {
                      imax =  i + shell;
                  }
                  if( j == int( sizes[1]) - 1)
                  {
                      jmax = int( sizes[1]) - 2;
                  }
                  else
                  {
                      jmax =  j + shell;
                  }
                  if( k == int( sizes[2]) - 1)
                  {
                      kmax = int( sizes[2]) - 2;
                  }
                  else
                  {
                      kmax =  k + shell;
                  }

                  // iterate 'cross-wise' along the xyz axes
                  surf = false;
//                      for( int ii = imin; ii <= imax; ii += 2)
//                      {
//                          neighbor_value = math::TupleXD< short, 3>::m_Data[ ii * yz_size + j * sizes[ 2] + k];
//                          if
//                          (
//                                  (value == 0 && neighbor_value != 0)
//                                  ||
//                                  (value < 0 && neighbor_value == 1)
//                          )
//                          {
//                              surf = true;
//                              break;
//                          }
//                      }

                  // collect all type-keys found in the first shell (unique, sorted)
                  std::set< short> neighbor_grid_values;

                  // MAYBE: introduce flag for using either cross-wise iteration or complete box
//                  bool use_complete_first_shell = false;
//                  if( use_complete_first_shell)
//                  {
//                      for( int ii = imin; ii <= imax; ++ii)
//                          for( int jj = jmin; jj <= jmax; ++jj)
//                              for( int kk = kmin; kk <= kmax; ++kk)
//                              {
//                                  neighbor_value = math::TupleXD< short, 3>::m_Data[ ii * yz_size + jj * sizes[ 2] + kk];
//                                  if( !( ii == i && jj == j && kk == k) && value != neighbor_value)
//                                  {
//                                      neighbor_grid_values.insert( neighbor_value);
//                                  }
//                              }
//                  }
//                  else
//                  {
                      for( int ii = imin; ii <= imax; ++ii)
                      {
                          neighbor_value = math::TupleXD< short, 3>::m_Data[ ii * yz_size + j * sizes[ 2] + k];
                          if( ii != i && value != neighbor_value)
                          {
                              neighbor_grid_values.insert( neighbor_value);
                          }
                      }
                      for( int jj = jmin; jj <= jmax; ++jj)
                      {
                          neighbor_value = math::TupleXD< short, 3>::m_Data[ i * yz_size + jj * sizes[ 2] + k];
                          if( jj != j && value != neighbor_value)
                          {
                              neighbor_grid_values.insert( neighbor_value);
                          }
                      }
                      for( int kk = kmin; kk <= kmax; ++kk)
                      {
                          neighbor_value = math::TupleXD< short, 3>::m_Data[ i * yz_size + j * sizes[ 2] + kk];
                          if( kk != k && value != neighbor_value)
                          {
                              neighbor_grid_values.insert( neighbor_value);
                          }
                      }
//                  }





//                  std::copy( neighbor_grid_values.begin(), neighbor_grid_values.end(), std::ostream_iterator< short>( std::cout, " "));
//                  std::cout << std::endl;

                  // all neighbors are of the same type than 'this' grid point
                  if( neighbor_grid_values.size() == 0)
                  {
                      continue;
                  }

//                  store::Vector3N< short> grid_point_indices( i, j, k);

                  std::pair< short, std::pair< short, short> >
                    grid_point_indices( std::make_pair( i, std::make_pair( j, k)));

                  std::vector< int> mol_type_ids_neighbor;

                  std::vector< std::vector< int> >
                      neighbor_collector;


                  if( value == 0)
                  {
//                      std::cout << "inside vol" << std::endl;
                      std::set< short> next_neighbors;

                      //                      m_SurfIndices[0].push_back( grid_point_indices);
                      
                      m_SurfIndices[0].insert( grid_point_indices);

                      for( std::set< short>::iterator itr = neighbor_grid_values.begin(); itr != neighbor_grid_values.end(); ++itr)
                          if( *itr < 0)
                          {
                              mol_type_ids_neighbor = KeyToIndices( *itr);
                              neighbor_collector.push_back( mol_type_ids_neighbor);
                              for( std::vector< int>::const_iterator tid_itr = mol_type_ids_neighbor.begin(); tid_itr != mol_type_ids_neighbor.end(); ++tid_itr)
                              {
                                  next_neighbors.insert( *tid_itr);
                              }
                          }

                      // add this point to the surface of all molecule types that are not in the next neighbor list
                      for( short ii = 1; ii < short( m_SurfIndices.size()); ++ii)
                      {
                    	  if( next_neighbors.find( ii) == next_neighbors.end())
                    	  {
                    		  //                              m_SurfIndices[ ii].push_back( grid_point_indices);
                    		  m_SurfIndices[ii].insert( grid_point_indices);
                    	  }
                      }

                      // add molecule types that have interfacing surfaces at this point
                      // add molecule types that are contained in one neigbor but not the others
                      bool found;
                      if( neighbor_collector.size() > 1)
                      {
                          for( int ii = 0; ii < int( neighbor_collector.size()); ++ii)
                              for( std::vector< int>::const_iterator itr = neighbor_collector[ii].begin(); itr != neighbor_collector[ii].end(); ++itr)
                              {
                                  found = false;
                                  for( int jj = 0; !found && jj < int( neighbor_collector.size()); ++jj)
                                  {
                                      if( ii == jj)
                                      {
                                    	  continue;
                                      }
                                      for( std::vector< int>::const_iterator sec = neighbor_collector[jj].begin(); !found && sec != neighbor_collector[jj].end(); ++sec)
                                      {
                                    	  if( *itr == *sec)
                                    	  {
                                    		  found = true;
                                    	  }
                                      }
                                  }
                                  if( !found)
                                  {
                                	  m_SurfIndices[ii].insert( grid_point_indices);
                                  }
                              }
                      }
                  }
                  else
                  {

                      // now: value <= 0 (should it be only < ? )
                      std::vector< int> molecule_type_ids_of_this_grid_point = KeyToIndices( value);
//                      std::cout << "value: " << value << " <= 0, mol types: ";
//                      std::copy( molecule_type_ids_of_this_grid_point.begin(), molecule_type_ids_of_this_grid_point.end(), std::ostream_iterator< int>( std::cout, " "));
//                      std::cout << std::endl;

                      int last = *neighbor_grid_values.rbegin(); // try largest value (last, sorted)
                      if( last > 0) // exposed to 'outer' world, outer for all molecule types
                      {
                    	  //                      std::cout << "exposed to outer world" << std::endl;
                    	  // the grid point is surface point for all molecule types defined here
                          for( std::vector< int>::const_iterator itr = molecule_type_ids_of_this_grid_point.begin(); itr != molecule_type_ids_of_this_grid_point.end(); ++itr)
                          {
//                              std::cout << *itr << " inserts: " << grid_point_indices << std::endl;
//                              std::copy( grid_point_indices.begin(), grid_point_indices.end(), std::ostream_iterator< int>( std::cout , " "));
//                              std::cout << std::endl;
                            //                          m_SurfIndices[*itr].push_back( grid_point_indices);
                        	  m_SurfIndices[*itr].insert( grid_point_indices);
                          }
                          continue;
                      }

    //                  std::cout << "multiple mol type dep surfs overlap" << std::endl;
                      // now we are in the region where multiple molecule type dependent surfaces overlap
                      for( std::set< short>::iterator itr = neighbor_grid_values.begin(); itr != neighbor_grid_values.end(); ++itr)
                          if( *itr < 0)
                          {
                              int diff = abs( value - *itr);
                              // m_SurfIndices[ diff].push_back( grid_point_indices);
                              m_SurfIndices[ diff].insert( grid_point_indices);
                          }

                  }


//                  return;


//                  if( neighbor_grid_values.size() == 1)
//                  {
//                      if( value == 0)
//                      {
//                          if( *neighbor_grid_values.begin() > 0)
//                          {
//                              // insert into surf of all molecules
//                              for( int ii = 0; ii < int( m_SurfIndices.size()); ++ii)
//                              {
//                                  m_SurfIndices[ii].push_back( grid_point_indices);
//                              }
//                          }
//                          else
//                          {
//                              std::vector< int> ids = KeyToIndices( *neighbor_grid_values.begin());
//                              // insert into surf of all molecules that are not in 'ids'
//                              for( int ii = 0; ii < int( m_SurfIndices.size()); ++ii)
//                              {
//                                  if( find( ids.begin(), ids.end(), ii) == ids.end())
//                                  {
//                                      continue;
//                                  }
//                                  m_SurfIndices[ii].push_back( grid_point_indices);
//                              }
//                          }
//                      }
//                      // being in molecule dependent volume
//                      else
//                      {
//                          if( *neighbor_grid_values.begin() > 0) // exposed to outside world
//                          {
//                              std::vector< int> ids = KeyToIndices( *neighbor_grid_values.begin());
//                              for( std::vector< int>::const_iterator itr = ids.begin(); itr != ids.end(); ++itr)
//                              {
//                                  m_SurfIndices[*itr].push_back( grid_point_indices);
//                              }
//                          }
//                          // noting to do if neighbor key is 0
//                      }
//                  }






//                      // if surf point
//                      if( surf)
//                      {
//                          store::Vector3N< short>
//                              indices( i, j, k);
//                          if( value == 0)
//                          {
//                              // add to all surfaces
//                              for( int ii = 0; ii < m_SurfIndices.size(); ++ii)
//                              {
//                                  m_SurfIndices[ii].push_back( indices);
//                              }
//                          }
//                          else
//                          {
//                              // add to specific surface
//                              m_SurfIndices[ abs( value)].push_back( indices);
//                          }
//                      }

              }
    }












//    std::vector< store::Vector3N< short> >
//    SurfGrid::IndicesOfSurfaceGridPoints() const
//    {
//        store::Vector3N< size_t> sizes( GetNrElements());
//
//        int
//            shell = 1,
//            yz_size = sizes[ 1] * sizes[ 2],
//            imin, imax,
//            jmin, jmax,
//            kmin, kmax,
//            id;
//        short
//            value;
//        bool
//            run;
//
//        // iterate through grid
//          for( int i = 0; i < sizes[0]; ++i)
//            for( int j = 0; j < sizes[1]; ++j)
//              for( int k = 0; k < sizes[2]; ++k)
//              {
//                  value = math::TupleXD< short, 3>::m_Data[ i * yz_size + j * sizes[ 2] + k];
//                  // if point is inside (value < 1) investigate one shell of points
//                  if( value < short( 1))
//                  {
//                      imin = std::max( 0, i - shell);
//                      jmin = std::max( 0, j - shell);
//                      kmin = std::max( 0, k - shell);
//                      imax = std::min( int( sizes[0]) - 1, i + shell);
//                      jmax = std::min( int( sizes[1]) - 1, j + shell);
//                      kmax = std::min( int( sizes[2]) - 1, k + shell);
//
//                      run = true;
//                      for( int ii = imin; ii <= imax && run; ++ii)
//                        for( int jj = jmin; jj <= jmax && run; ++jj)
//                          for( int kk = kmin; kk <= kmax; ++kk)
//                          {
//                              // if one of the neighboring grid points is outside (value > 0) then add indices of center point to surface list
//                              if( math::TupleXD< short, 3>::m_Data[ ii * yz_size + jj * sizes[ 2] + kk] > 0)
//                              {
//                                  m_SurfIndices[ abs( value)].push_back( store::Vector3N< short>( i, j, k));
//                                  run = false;
//                                  break;
//                              }
//                          }
//                  }
//              }
//    }

    size_t
    SurfGrid::CountAtomsInside( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE) const
    {
        size_t counter( 0);
        int value;
        std::vector< int> indices;

        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator itr = MOLECULE->GetAtoms().begin(); itr != MOLECULE->GetAtoms().end(); ++itr)
        {
            value = math::TupleXD< short, 3>::operator()( ( *itr)->GetPosition());
            if( value == 0)  // inside for all molecule types
            {
                ++counter;
            }
            else if( value < 0)  // inside for specific molecule types
            {
                // check whether the residue type of this molecule is given in the member list of molecule types for which specific volumes are defined
                std::vector< std::string>::iterator str_itr = find( m_MolTypes.begin(), m_MolTypes.end(), ( *itr)->GetResidueType());
                // if it is defined
                if( str_itr != m_MolTypes.end())
                {
                    // get id of molecule type
                    int id = str_itr - m_MolTypes.begin();
                    // dissect the key value for this grid point (chmod-like combinatoric)
                    indices = KeyToIndices( value);
                    // if molecule type is in agreement with
                    if( find( indices.begin(), indices.end(), id) != indices.end())
                    {
                        ++counter;
                    }
                }
            }
        }
        return counter;
    }



    float SurfGrid::CalcVoxelVolume() const
    {
        return  m_Delta[0] * m_Delta[1] * m_Delta[2];
    }



    store::Vector3N< int>
    SurfGrid::AdjustMaxIndices( const store::Vector3N< int> &MAX) const
    {
        store::Vector3N< int> max( MAX);
        for( int i = 0; i < 3; ++i)
        {
            if( max[i] >= math::TupleXD< short, 3>::m_NrElements[i])
            {
                std::cout << "=> maximum value of surf grid is adjusted to match definitions (" << i << ")" << std::endl;
                max[i] = math::TupleXD< short, 3>::m_NrElements[i] - 1;
            }
        }
        return max;
    }

    store::Vector3N< int>
    SurfGrid::AdjustMinIndices( const store::Vector3N< int> &MIN) const
    {
        store::Vector3N< int> min( MIN);
        for( int i = 0; i < 3; ++i)
        {
            if( min[i] < 0)
            {
                std::cout << "=> minimum value of surf grid is adjusted to match definitions (" << i << ")" << std::endl;
                min[i] = 0;
            }
        }
        return min;
    }


    //! returns a volume value for each molecule type defined by surface objects
    std::vector< float>
    SurfGrid::Volume( const store::Vector3N< int> &MIN, const store::Vector3N< int> &MAX) const
    {
    	StandardWrite( __PRETTY_FUNCTION__ << " from: " << MIN[0] << " " << MIN[1] << " " << MIN[2] << " to: " << MAX[0] << " " << MAX[1] << " " << MAX[2]);

    	float voxel_volume( CalcVoxelVolume());

	StandardWrite( "voxel-volume: " << voxel_volume);

        std::vector< float> volume( m_MolTypes.size());  // volume[0]: all, volume[1] e.g.: POPE
        short value;
        std::vector< int> indices;

        store::Vector3N< int> max = AdjustMaxIndices( MAX);

        store::Vector3N< int> min = AdjustMinIndices( MIN);

        StandardWrite( "adjusted: " << " from: " << min[0] << " " << min[1] << " " << min[2] << " to: " << max[0] << " " << max[1] << " " << max[2]);

        for( int i = min[0]; i <= max[0]; ++i)
            for( int j = min[1]; j <= max[1]; ++j)
                for( int k = min[2]; k <= max[2]; ++k)
                {
                    value = m_Data[ CalcID( store::Vector3N< int> ( i, j, k))];
                    if( value == 0)
                    {
                        for( std::vector< float>::iterator itr = volume.begin(); itr != volume.end(); ++itr)
                        {
                            *itr += voxel_volume;
                        }
                    }
                    else if( value < 0)
                    {
                        indices = KeyToIndices( value);
//                        std::cout << __FUNCTION__ << " value: " << value << " indices: ";
//                        std::copy( indices.begin(), indices.end(), std::ostream_iterator< int>( std::cout , " "));
//                        std::cout << " size: " << volume.size() << std::endl;
                        for( std::vector< int>::const_iterator itr = indices.begin(); itr != indices.end(); ++itr)
                        {
                            assert( *itr < int( volume.size()));
                            volume[ *itr] += voxel_volume;
                        }
                    }
                }
        return volume;
    }


    store::Map< std::string, float>
    SurfGrid::Volume
    (
            const math::Vector3N &MIN,
            const math::Vector3N &MAX
    ) const
    {
    	StandardWrite( __PRETTY_FUNCTION__ << " from: " << MIN[0] << " " << MIN[1] << " " << MIN[2] << " to: " << MAX[0] << " " << MAX[1] << " " << MAX[2]);

    	store::Map< std::string, float>
            map;
        store::Vector3N< int>
            min( math::TupleXD< short, 3>::CalcIndices( MIN)),
            max( math::TupleXD< short, 3>::CalcIndices( MAX));

        std::vector< float> volume( Volume( min, max));
        std::vector< float>::const_iterator vol_itr = volume.begin();
        std::vector< std::string>::const_iterator mol_itr = m_MolTypes.begin();
        for( ; mol_itr != m_MolTypes.end() && vol_itr != volume.end(); ++mol_itr, ++vol_itr)
        {
            map.InsertNewKeyAndValue( *mol_itr, *vol_itr);
        }
        return map;
    }




    std::pair< float, int>
    SurfGrid::ClosestSurfaceGridPoint( const math::Vector3N &POSITION, const std::string &TYPE) const // TODO: write a general algorithm to be called here, fast search of closest points
    {
        std::multimap< short, std::pair< short, short> >::iterator
          itr;

        std::pair< std::multimap< short, std::pair< short, short> >::iterator, std::multimap< short, std::pair< short, short> >::iterator>
          itr_pair;

        std::vector< int>
          indices = ProjectOnGrid( CalcIndices( POSITION));

        int
            best_dist( std::numeric_limits<int>::max()),
            d,
            ipos = indices[0],
            iup = ipos,
                idown = ipos - 1,
            cc = 0,
            type;

        store::Vector3N< int>
          surf_indices,
          best;

        std::vector< std::string>::const_iterator
          type_itr = find( m_MolTypes.begin(), m_MolTypes.end(), TYPE);
        
        // specified mol type
        if( type_itr != m_MolTypes.end())
        {
            // translate iterator to id to access type dependent vector
            type = type_itr - m_MolTypes.begin();
        }
        else
        {
            // 'all'
            type = 0;
        }

        while( ( std::pow( float( ipos - iup), 2) < best_dist || std::pow( float( ipos - idown),2) < best_dist) && ( iup < m_NrElements[0] || idown >= 0))
        {
          //          std::cout << ipos << " " << iup << " " << idown << " " << sqrt( best_dist) << " " << m_NrElements[0] << std::endl;

          if( idown >= 0)
            {
                itr_pair = m_SurfIndices[ type].equal_range( idown);
                //                std::cout << "down: " << idown << " " << idown + 1 << " " << *itr_pair.first << " " << *itr_pair.second << std::endl;
                for( itr = itr_pair.first; itr != itr_pair.second; ++itr)
                {
                  //                  std::cout << "down: " << *itr << std::endl;
                    d = int( std::pow( float( ipos - (*itr).first),2) + std::pow( float( indices[1] - (*itr).second.first),2) + std::pow( float( indices[2] - (*itr).second.second),2));
                    if( d < best_dist)
                    {
                      //                      std::cout << "improved down" << std::endl;
                        best_dist = d;
                        best = GetIndicesOfSurfacePoint( itr);
                    }
                }
                --idown;
            }
          if( iup < m_NrElements[0])
            {
                itr_pair = m_SurfIndices[ type].equal_range( iup);
                //                std::cout << "up: " << iup << " " << iup + 1 << " " << *itr_pair.first << " " << *itr_pair.second << std::endl;
                for( itr = itr_pair.first; itr != itr_pair.second; ++itr)
                {
                  //                  std::cout << "up: " << *itr << std::endl;
                    d = int( std::pow( float( ipos - (*itr).first),2) + std::pow( float( indices[1] - (*itr).second.first),2) + std::pow( float( indices[2] - (*itr).second.second),2));
                    if( d < best_dist)
                    {
                      //                      std::cout << "improved up" << std::endl;
                        best_dist = d;
                        best = GetIndicesOfSurfacePoint( itr);
                    }
                }
                ++iup;
            }
            ++cc;
        }

        std::pair< float, int> result = std::make_pair( math::Distance( POSITION, CalcPositionFromIndices( best)), CalcID( best));

        // check for consistency:
        std::pair< float, int> check_pair = BasicClosestSurfPoint( POSITION, TYPE);
        if( result != check_pair)
        {
        	std::cout << "====> incorrect closest point found in " << __FUNCTION__
        			<< ", should be: (dist: " << check_pair.first << ", id: " << check_pair.second
					<< ") found: (dist: " << result.first << ", id: " << result.second << ")" << std::endl;
        }
        else
        {
			std::cout << "=> correct" << std::endl;
        }

        return result;
    }


    std::pair< float, int>
    SurfGrid::BasicClosestSurfPoint( const math::Vector3N &POSITION, const std::string &TYPE) const
    {
    	float
			best_dist( std::numeric_limits< float>::max()),
			dist;
		int
			type;

        std::vector< std::string>::const_iterator
			type_itr = find( m_MolTypes.begin(), m_MolTypes.end(), TYPE);

        store::Vector3N< int>
			best;
//        std::vector< int>
//           indices = ProjectOnGrid( CalcIndices( POSITION));

        // specified mol type
        if( type_itr != m_MolTypes.end())
        {
            // translate iterator to id to access type dependent vector
            type = type_itr - m_MolTypes.begin();
        }
        else
        {
            // 'all'
            type = 0;
        }


        for( std::multimap< short, std::pair< short, short> >::const_iterator itr = m_SurfIndices[ type].begin(); itr != m_SurfIndices[ type].end(); ++itr)
        {
        	dist = math::SquaredDistance( POSITION, CalcPositionFromIndices( store::Vector3N< int>( itr->first, itr->second.first, itr->second.second)));
        	if( dist < best_dist)
        	{
        		best_dist = dist;
                best = store::Vector3N< int>( itr->first, itr->second.first, itr->second.second);
        	}
        }

        return std::make_pair( sqrt( best_dist), CalcID( best));
    }



    bool
    SurfGrid::IsInside( const short &VALUE, const std::string &NAME) const
    {
      //      StandardWrite( __FUNCTION__ << " value: " << VALUE << " name: " << NAME);
        if( VALUE == 0)
        {
            return true;
        }
        std::vector< int> mol_ids =  KeyToIndices( int( VALUE));
        for( std::vector< int>::const_iterator itr = mol_ids.begin(); itr != mol_ids.end(); ++itr)
        {
          //          std::cout << *itr << " <=> " << m_MolTypes[ *itr] << std::endl;
            if( m_MolTypes[ *itr] == NAME)
            {
              //              std::cout << "found" << std::endl;
                return true;
            }
        }
        //        std::cout << "outside" << std::endl;
        return false;
    }

    std::ofstream &
    SurfGrid::WriteSurfAsPdb( std::ofstream &STREAM) const
    {
	const char *array[] = {"H", "O", "C", "N","S","Cl","Na"};
        std::vector< std::string> atoms(7);
        std::copy( array, array + 7, atoms.begin());
        STREAM << "REMARK   ";
        std::copy( atoms.begin(), atoms.end(), std::ostream_iterator< std::string>( STREAM, "  "));
        STREAM << std::endl;
        std::vector< int> indices( 3);
        size_t
            atom_id = 0,
            mol_id = 0;
        math::Vector3N pos;

        STREAM << "REMARK  visualization for surface as given by GRIFFINs surf grid, " << m_SurfIndices.size() << " molecule type defined" << std::endl;
        int count = 0;
        for( size_t i = 0; i < m_SurfIndices.size(); ++i)
        {
            STREAM << "REMARK  " << m_SurfIndices[i].size() << " surf points for: " << atoms[i] <<  std::endl;
            count += m_SurfIndices[i].size();
        }
        STREAM << "REMARK  " << count << " surf points in total" << std::endl;

        for( size_t i = 0; i < m_SurfIndices.size(); ++i)
        {
            for( std::multimap< short, std::pair< short, short> >::const_iterator itr = m_SurfIndices[i].begin(); itr != m_SurfIndices[i].end(); ++itr, ++atom_id)
            {
                indices = GetIndicesOfSurfacePoint( itr);
                pos = CalcPositionFromIndices( indices);

                atom_id = atom_id % 100000;

                if( atom_id == 0)
                {
                    ++mol_id;
                }
                STREAM << "ATOM  ";
                STREAM.width( 5);
                STREAM <<  atom_id << "  ";
                STREAM.width( 3);
                STREAM << atoms[i] << " SUR ";
                STREAM.width(5);
                STREAM << mol_id % 100000 << "    ";
                STREAM.width( 8);
                STREAM.precision( 3);
                STREAM << g_PDBFactor * pos( 0);
                STREAM.width( 8);
                STREAM.precision( 3);
                STREAM << g_PDBFactor * pos( 1);
                STREAM.width( 8);
                STREAM.precision( 3);
                STREAM << g_PDBFactor * pos( 2);
                STREAM.width( 15);
                STREAM << " ";
                STREAM << std::endl;
            }
          }
        return STREAM;
    }


    std::ostream &
    SurfGrid::Write( std::ostream &STREAM) const
    {
        STREAM << "SurfGrid" << std::endl;
        STREAM << "moltypes:    " << m_MolTypes.size() << std::endl;
        for( std::vector< std::string>::const_iterator itr = m_MolTypes.begin(); itr != m_MolTypes.end(); ++itr)
        {
            STREAM << *itr << std::endl;
        }
        STREAM << "surfpoints:  " << m_SurfIndices.size() << std::endl;
        for( size_t i = 0; i < m_SurfIndices.size(); ++i)
        {
            STREAM << m_SurfIndices[i].size() << std::endl;
            for( std::multimap< short, std::pair< short, short> >::const_iterator itr = m_SurfIndices[i].begin(); itr != m_SurfIndices[i].end(); ++itr)
            {
                STREAM << itr->first << " " << itr->second.first << " " << itr->second.second << std::endl;
            }
        }

        math::TupleXD< short, 3>::Write( STREAM);

        return STREAM;
    }

    std::istream &
    SurfGrid::Read( std::istream& STREAM)
    {
    	DebugWrite( __PRETTY_FUNCTION__);
    	if( !STREAM)
    	{
    		std::cout << "===> stream closed in " << __PRETTY_FUNCTION__ << std::endl;
    		exit( -1);
    	}
        std::string str;
        int size;
        short x, y, z;

        STREAM >> str;
        if( str != "SurfGrid")
	{
	    std::cout << __PRETTY_FUNCTION__ << " in " << __FILE__ << " line " << __LINE__ << std::endl;
	    std::cout << "wrong class id string: <" << str << ">, expected \'SurfGrid\'" << std::endl;
	    exit( -1);
	}

        STREAM >> str;
        if( str != "moltypes:")
	{
	    std::cout << "string should be \'moltypes\' but is <" << str << ">" << std::endl;
	    exit( -1);
	}
        STREAM >> size;
        m_MolTypes = std::vector< std::string>( size);
        for( std::vector< std::string>::iterator itr = m_MolTypes.begin(); itr != m_MolTypes.end(); ++itr)
        {
            STREAM >> *itr;
        }

        STREAM >> str;
        assert( str == "surfpoints:");
        STREAM >> size;
        m_SurfIndices = std::vector< std::multimap< short, std::pair< short, short> > >( size);
        for( std::vector< std::multimap< short, std::pair< short, short> > >::iterator type_itr = m_SurfIndices.begin(); type_itr != m_SurfIndices.end(); ++type_itr)
        {
            STREAM >> size;
            for( int i = 0; i < size; ++i)
            {
                STREAM >> x >> y >> z;
                type_itr->insert( std::make_pair( x, std::make_pair( y, z)));
            }
        }

        math::TupleXD< short, 3>::Read( STREAM);

        return STREAM;
    }

    std::istream & operator >> ( std::istream &STREAM, SurfGrid &GRID)
    {
        return GRID.Read( STREAM);
    }

    std::ostream & operator << ( std::ostream &STREAM, const SurfGrid &GRID)
    {
        return GRID.Write( STREAM);
    }



    //! Generate random point on the sphere surface (used in MAYER)
    void GenerateRandomPointOnUnitSphereSurface
    (
            float &XP,
            float &YP,
            float &ZP
    )
    {
//        StandardWrite( __FUNCTION__);
        float rp;
        do
        {
            XP = math::BasicRandom( -1.0, 1.0); // equal behaviour of random generators???
            YP = math::BasicRandom( -1.0, 1.0);
            ZP = math::BasicRandom( -1.0, 1.0);
            rp = XP*XP + YP*YP + ZP*ZP;
//            StandardWrite( "rp: " << rp);
        }
        while( rp == 0 || rp > 1.0);
        rp = sqrt( rp);
        XP=XP/rp;
        YP=YP/rp;
        ZP=ZP/rp;
    }





    //  Segment part formula (used in MAYER)
    float
    SegmentPartFormula
    (
            const float &FIRST_VDW_PLUS_PROBE_RADIUS,
            const float &SECOND_VDW_PLUS_PROBE_RADIUS,
            const float &ATOM_DISTANCE
    )
    {
        float v1, v2, segpar;

        v1 = FIRST_VDW_PLUS_PROBE_RADIUS + SECOND_VDW_PLUS_PROBE_RADIUS;
        v2 = fabs( FIRST_VDW_PLUS_PROBE_RADIUS - SECOND_VDW_PLUS_PROBE_RADIUS);
        if( ATOM_DISTANCE > v1)
        {
            segpar=0.;
        }
        else if( ATOM_DISTANCE < v2)
        {
            if( FIRST_VDW_PLUS_PROBE_RADIUS > SECOND_VDW_PLUS_PROBE_RADIUS)
            {
                segpar=0.;
            }
            else
            {
                segpar=1.;
            }
        }
        else
        {
            segpar = ( SECOND_VDW_PLUS_PROBE_RADIUS * SECOND_VDW_PLUS_PROBE_RADIUS - pow( ( FIRST_VDW_PLUS_PROBE_RADIUS - ATOM_DISTANCE), 2)) / (4 * ATOM_DISTANCE * FIRST_VDW_PLUS_PROBE_RADIUS);
        }
        return segpar;
    }

    //TRANX = HALF*(NCLX-1)*DCEL;


    void
    MayerContactAndReentrance
    (
     const  int     &NTPRP,  /// nr atoms
     const  float  *X, // atom pos of protein
     const  float  *Y,
     const  float  *Z,
     const  float  *PBRAD, // vdw radii
     const  float  &RADW, // probe radius
     const  int     &NCLX,  // number of cells of grid
     const  int     &NCLY,
     const  int     &NCLZ,
     const  float  &DCEL,   // cell size
     const  float  &TRANX,  //
     const  float  &TRANY,
     const  float  &TRANZ,
      const  float  &XBCEN,  // center of the grid
     const  float  &YBCEN,
     const  float  &ZBCEN,
     std::vector< short>  &RMAY,
            float  *POX,     // sphere surface positions
            float  *POY,
            float  *POZ,
     const  int     &MAPT   // number of spheres points
    )
    {

        StandardWrite( __FUNCTION__);

        int  iseed, jl, i, j, l, m, nl, ml, kl, ll, nfil, ix, iy, iz;
        int  jx1, jx2, jy1, jy2, jz1, jz2, k, ipx, ipy, ipz, ncyz, cc;
        float  wrsf, raneti, ranetj, xi, yi, zi, xj, yj, zj, aij, fistr0=0.0;
        float  xp, yp, zp, bisq, xc, yc, zc, xsq, ysq, zsq, dsq;
//        float alpha, beta, oxp, oyp, ozp;
        //  float  segpar;

        float rsmall = 1e-6; // as of squantm/sqnt_qm2_mopac.src90

        iseed = 314159265;
        ncyz = NCLY * NCLZ;

//#ifdef NEWRNG
        srand( iseed);
        //iseed = 1;
        //Cmh050712      ISEED=MYNODP              !##PARALLEL
//#endif

        // generate points on the sphere surface and adjust parameter
        wrsf = 0.0;
        for( i = 0; i < NTPRP; ++i)
        {
            raneti = PBRAD[i] + RADW;
            if( wrsf < raneti) // get max atom radius plus probe radius
            {
                wrsf = raneti;
            }
        }
        wrsf = MAPT / ( wrsf * wrsf) + 0.000001; // prestep for subsequent radius adjustment of MAPT

        StandardWrite( "sphere surface point creation ...");
        for( i = 0; i < MAPT; ++i)
        {
            GenerateRandomPointOnUnitSphereSurface( POX[ i], POY[ i], POZ[ i]);
        }

        StandardWrite( "main loop ...");

        // main loop by atoms
        for( i = 0; i < NTPRP; ++i)
        {
            nl = 0;
            std::multimap< float, int> neighbor_map;
//            std::map< size_t, Triplet< float, float, float> > angle_dist_map;

            xi = X[ i];
            yi = Y[ i];
            zi = Z[ i];

            raneti = PBRAD[ i] + RADW;

//            StandardWrite( "atom " << i << " sort closest neighbors");

            // create the list of closest neighbors for the atom from the main loop
            // according to their screening ability of the neighbors
            // 0<Segpar<1
            // Segpar=0  if the neighbor does not reduce the accessible surface of the atom
            //           and it means that it is not a neighbor at all.
            // Segpar=1  if the neighbor buries the atom totally.
            for( j = 0; j < NTPRP; ++j)
            {
                if( i == j)
                {
                    continue;
                }

                xj = X[ j];
                yj = Y[ j];
                zj = Z[ j];

                ranetj = PBRAD[ j] + RADW;

                // atom distance
                aij = sqrt( ( xi - xj) * ( xi - xj) + ( yi - yj) * ( yi - yj) + ( zi - zj) * ( zi - zj));

                // surf influence factor
                fistr0 = SegmentPartFormula( raneti, ranetj, aij);

                //                std::cout << i << ": " << raneti << ", " << j << ": " << ranetj << " dist: " << aij << " segpar: " << fistr0 << std::endl;
                if( fistr0 == 0.)
                {
//                    std::cout << "this atom cannot have any effect on surface, check next atom in second loop" << std::endl;
                    continue;
                }
                if( fistr0 == 1.)
                {
                    DebugWrite( "atom " << j << " buries atom " << i << ", go to next atom in main loop ...");
                    break;
                }

                ++nl;

                neighbor_map.insert( std::make_pair( fistr0, j));

//                alpha = math::AngleFromLawOfCosinus( ranetj, aij, raneti);
//                beta = math::AngleFromLawOfCosinus( raneti, aij, ranetj);
//                angle_dist_map.insert( std::make_pair( j, Triplet( alpha, beta, aij)));
            }

            if( fistr0 == 1.)
            {
//                std::cout << "... go to next atom in main loop ..." << std::endl;
                continue;
            }

//            std::cout << "neigbours of atom " << i << ": " << std::endl;

//            for( std::multimap< float, int>::const_reverse_iterator itr = neighbor_map.rbegin(); itr != neighbor_map.rend(); ++itr)
//            {
//                std::cout <<  itr->first << " : " << itr->second << std::endl;
//            }

            ml = int( wrsf * raneti * raneti);  // adjust nr surf points to radius: ml = min( MAPT, MAPT * (r_i+r_p)^2 / (r_max+r_p)^2)
            if( ml > MAPT)
            {
                ml = MAPT;
            }

//            StandardWrite( "loop over sphere surface");

            // loop over points on the sphere surface
            for( jl = 0; jl < ml; ++jl)
            {
                // surf point of atom i
                xp = xi + raneti * POX[ jl];
                yp = yi + raneti * POY[ jl];
                zp = zi + raneti * POZ[ jl];

//                std::cout << "check for unity: " << POX[ jl] * POX[ jl] + POY[ jl] * POY[ jl] + POZ[ jl] * POZ[ jl] << std::endl;

                // if atom has no neighbor surf point cannot be within smooth volume
                if( nl > 0)
                {
//                    std::cout << "enter (" << nl << " neighbors)" << std::endl;
                    // check sorted neighbors
                    std::multimap< float, int>::const_reverse_iterator m_itr = neighbor_map.rbegin();
                    bool leave_iteration( false);
                    for( kl = 0; kl < nl; ++kl, ++m_itr)
                    {
                        ll = m_itr->second;
                        aij = pow( ( xp - X[ ll]), 2) + pow( ( yp - Y[ ll]), 2) + pow( ( zp - Z[ ll]), 2);  // distance of surface point to neighbor atom
                        // if surf point within neighbor volume everything should be correct = nothing further to do
                        if( aij < pow( PBRAD[ll] + RADW, 2))
                        {
//                            std::cout << "surf point within neighbor volume: " << jl << " neigb: " << kl << " id: " << ll << " dist: " << ( sqrt( aij)) << " < " << PBRAD[ll] << " + " <<  RADW  << std::endl;
                            leave_iteration = true;
                            break;
                        }
                    }
                    if( leave_iteration)
                    {
//                        std::cout << "... sneak out of this iteration" << std::endl;
                        continue;
                    }
                }
                // else reset all guys to 1


//                std::cout << "continue" << std::endl;

                // reset grid points inside a probe sphere around any
                // of random accessible points as a dielectric media.
                bisq = RADW;
                // nr of surrounding gridpoints that are within probe sized test shell
                nfil = int( bisq / DCEL) + 2;
                // relative sphere point position in grid
                xp = xp + TRANX - XBCEN;
                yp = yp + TRANY - YBCEN;
                zp = zp + TRANZ - ZBCEN;
                bisq = bisq * bisq + rsmall;

                // grid indices of surf point
                ix = int( xp / DCEL);
                iy = int( yp / DCEL);
                iz = int( zp / DCEL);

                // subgrid around surf point, adjust if it hits limits of main grid
                jx1 = ix - nfil;
                if( jx1 < 0){ jx1 = 0;}
                jx2 = ix + nfil;
                if( jx2 > NCLX){ jx2 = NCLX;}
                jy1 = iy - nfil;
                if( jy1 < 0){ jy1 = 0;}
                jy2 = iy + nfil;
                if( jy2 > NCLY){ jy2 = NCLY;}
                jz1 = iz - nfil;
                if( jz1 < 0){ jz1 = 0;}
                jz2 = iz + nfil;
                if( jz2 > NCLZ){ jz2 = NCLZ;}

                for( k = jx1; k < jx2; ++k)
                {
                    ipx = k * ncyz;
                    xc  = k * DCEL;
                    xsq = (xc - xp) * (xc - xp);
                    for( l = jy1; l < jy2; ++l)
                    {
                        ipy = l * NCLZ;
                        yc  = l * DCEL;
                        ysq = (yc - yp) * (yc - yp);
                        for( m = jz1; m < jz2; ++m)
                        {
                            ipz = m + ipy + ipx;
                            zc = m * DCEL;
                            zsq = (zc - zp) * (zc - zp);
                            // distance of subgrid point with original surfpoint
                            dsq = xsq + ysq + zsq;

                            if( RMAY[ ipz] < 0.0)
                            {
                                if( dsq < bisq)
                                {
//                                    std::cout << "grid point is outside " << ( sqrt( dsq)) << " < " << ( sqrt( bisq)) << " | " << k << " " << l << " " << m << " " << ipz << std::endl;
                                    assert( ipz < int( RMAY.size()));
                                    RMAY[ ipz] = -RMAY[ ipz]; // this defines points in the probe shell to be outside of smoothed surface // ! bulk kappa restored
                                }
                            }
                        }
                    }
                }
            }
        }
        StandardWrite( "set remaining negative values in rmay to zero = inside smoothed surf");
        cc =0;
        for( k = 0; k < NCLX; k++)
        {
            ipx = k * ncyz;
            for( l = 0; l < NCLY; l++)
            {
                ipy = l * NCLZ;
                for( m = 0; m < NCLZ; ++m, ++cc)
                {
                    ipz = m + ipy + ipx;
//                    std::cout << k << " " << l << " " << m << " " << ipz << " " << cc << std::endl;
                    if(  RMAY[ipz] < 0.)
                    {
//                        std::cout << "generously set to zero" << std::endl;
                        RMAY[ipz] = 0;
                    }
                }
            }
        }
        std::cout << RMAY.size() << std::endl;
    } // end MayerContactAndReentrance



    //  when the dielectric boundary is defined by the van der Waals surface
    //  M = 0 : inside solute
    //      1 : outside solute
    void
    MayerMStep
    (
            const  int    &NTPRP,  // number of atoms
            const  float *X,
            const  float *Y,
            const  float *Z,
            const  float *PBRAD,  // radii
            const  float &RADW,  // water radius
            const  int    &NCLX,  // number of cells in x
            const  int    &NCLY,
            const  int    &NCLZ,
            const  float &DCEL,  // cell size
            const  float &TRANX,
            const  float &TRANY,
            const  float &TRANZ,
            const  float &XBCEN, // box center x coordin
            const  float &YBCEN,
            const  float &ZBCEN,
            std::vector< short> &RMAY
    )
    {
        StandardWrite( __FUNCTION__);

        int      i, k, l, m, ix, iy, iz, jx1, jx2, jy1, jy2, jz1, jz2;
        int      ipx, ipy, ipz;
        int      nfil, ncyz;
        float   xi, yi, zi, xc, yc, zc, dsq, xsq, ysq, zsq;
        float   sqr, sqrw;


        //cwi      ncyz = NCLX * NCLY
        ncyz = NCLY * NCLZ;
        for( i = 0; i < NTPRP; ++i)
        {
            sqr  = PBRAD[i];
            if(sqr <= 0.0)
            {
                break;
            }
            xi = X[i] + TRANX - XBCEN;
            yi = Y[i] + TRANY - YBCEN;
            zi = Z[i] + TRANZ - ZBCEN;
            sqrw = sqr + RADW;
            nfil = int( sqrw / DCEL) + 2;
            sqr  = sqr  * sqr;
            sqrw = sqrw * sqrw;
            ix = int(xi / DCEL);
            iy = int(yi / DCEL);
            iz = int(zi / DCEL);

//          jx1 = ix - nfil + 1;
//          if(jx1 < 1){ jx1 = 1;}
//          jx2 = ix + nfil - 1;
//          if(jx2 > NCLX){ jx2 = NCLX;}
//          jy1 = iy - nfil + 1;
//          if(jy1 < 1){ jy1 = 1;}
//          jy2 = iy + nfil - 1;
//          if(jy2 > NCLY){ jy2 = NCLY;}
//          jz1 = iz - nfil + 1;
//          if(jz1 < 1){ jz1 = 1;}
//          jz2 = iz + nfil - 1;
//          if(jz2 > NCLZ){ jz2 = NCLZ;}

            // subgrid around surf point, adjust if it hits limits of main grid
            jx1 = ix - nfil;
            if( jx1 < 0){ jx1 = 0;}
            jx2 = ix + nfil;
            if( jx2 > NCLX){ jx2 = NCLX;}
            jy1 = iy - nfil;
            if( jy1 < 0){ jy1 = 0;}
            jy2 = iy + nfil;
            if( jy2 > NCLY){ jy2 = NCLY;}
            jz1 = iz - nfil;
            if( jz1 < 0){ jz1 = 0;}
            jz2 = iz + nfil;
            if( jz2 > NCLZ){ jz2 = NCLZ;}

            for( k = jx1; k < jx2; ++k)
            {
                ipx = k * ncyz;
                xc = k * DCEL;
                for( l = jy1; l < jy2; ++l)
                {
                    ipy = l * NCLZ;
                    yc = l * DCEL;
                    for( m = jz1; m < jz2; ++m)
                    {
                        ipz = m + ipy + ipx;
                        zc = m * DCEL;
                        xsq = (xc - xi) * (xc - xi);
                        ysq = (yc - yi) * (yc - yi);
                        zsq = (zc - zi) * (zc - zi);

                        if(RMAY[ipz] != 0)
                        {
                            dsq = xsq + ysq + zsq;  // distance squared to neares atom center
                            if(dsq <= sqr)   // if square of distance is smaller than square of the radius => zero = inside of the atom
                            {
                                // Zero the Debye-Huckel factor inside the solute
//                                          std::cout << "// Zero the Debye-Huckel factor inside the solute" << std::endl;
                                RMAY[ipz] = 0;
                            }
                            else if(dsq > sqr && dsq <= sqrw)
                            {
                                if(RMAY[ipz] > 0)
                                {
//                                  std::cout << "invert rmay" << std::endl;
                                    RMAY[ipz] = -RMAY[ipz];
                                }
                            }
                        }
                    }
                }
            }
        }
    }// end MayerMStep


    void
    PrepareSurfGrid
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            //            std::vector< short> &RMAY,
            std::vector< float> &MIN,
            std::vector< float> &MAX,
            std::vector< int> &NR_BINS,
            const float &DELTA,
            float *X_COO,
            float *Y_COO,
            float *Z_COO,
            float *VDW_RADII,
            std::vector< float> &CENTER,
            const float &FACTOR,
            const float &OFFSET, 
	    const float &MARGIN
    )
    {
        StandardWrite( __FUNCTION__);

        int size = MOLECULE->GetAtoms().size();

        math::Vector3N 
          pos, 
          total_delta;

        for( int i = 0; i < size; ++i)
        {
            pos = MOLECULE->GetAtoms()( i)->GetPosition();
            X_COO[ i] = pos( 0);
            Y_COO[ i] = pos( 1);
            Z_COO[ i] = pos( 2);

            for( int j = 0; j < 3; ++j)
            {
                if( pos( j) - MARGIN < MIN[j])
                {
                  //                  std::cout << "=> min " << j << " changed from: " << MIN[j] << " to : " << pos(j) - 3.0 << std::endl;
                    MIN[j] = pos( j) - MARGIN;
                }
                if( pos( j) + MARGIN > MAX[j])
                {
                  //                  std::cout << "=> max " << j << " changed from: " << MAX[j] << " to : " << pos(j) + 3.0 << std::endl;
                  MAX[j] = pos( j) + MARGIN;
                }
            }

            VDW_RADII[ i] = FACTOR * MOLECULE->GetAtoms()( i)->GetVanDerWaalsRadius() + OFFSET;
//            std::cout <<  X_COO[i] << "  " << Y_COO[ i] << "  " << Z_COO[ i] << "  " << VDW_RADII[i] << std::endl;
        }

        for( int i = 0; i < 3; ++i)
        {
            NR_BINS[ i] = int( (MAX[i] - MIN[i]) / DELTA) + 1;
            total_delta(i) = NR_BINS[i] * DELTA;
            CENTER[i] = 0.5 * ( MAX[i] + MIN[i]);
        }
    }



    void
    BuildSurfGrid
    (
     //            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            std::vector< short> &RMAY,
            std::vector< float> &MIN,
            std::vector< int> &NR_BINS,
            const float &DELTA,
            const float &SMOOTH, 
            const float &PROBE_RADIUS,
            const size_t &SIZE,
            const float *X_COO,
            const float *Y_COO,
            const float *Z_COO,
            const float *VDW_RADII,
            const std::vector< float> &CENTER
    )
    {
        StandardWrite( __FUNCTION__);

        int total_array_length = NR_BINS[0] * NR_BINS[1] * NR_BINS[2];
        StandardWrite( "total grid-array length: " << total_array_length);
        std::cout << "nr grid-points in x y z:  ";
        std::copy( NR_BINS.begin(), NR_BINS.end(), std::ostream_iterator< size_t>( std::cout , " "));
        std::cout << "delta: " << DELTA << std::endl;

        StandardWrite( "initially set all grid values to 1");
        RMAY = std::vector< short>( total_array_length, short( 1));

        StandardWrite( "calc vdW surf");
        MayerMStep
        (
                SIZE,
                X_COO,
                Y_COO,
                Z_COO,
                VDW_RADII,
                PROBE_RADIUS,
                NR_BINS[0],
                NR_BINS[1],
                NR_BINS[2],
                DELTA,
                0.5 * ( NR_BINS[ 0] - 1) * DELTA,
                0.5 * ( NR_BINS[ 1] - 1) * DELTA,
                0.5 * ( NR_BINS[ 2] - 1) * DELTA,
                CENTER[0],
                CENTER[1],
                CENTER[2],
                RMAY
        );

//        std::ofstream write;
//        Open( write, "rmay_grid_mstep.txt");
//        std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( write , " "));
//        write << std::endl;
//        Close( write);

        if( SMOOTH)
        {
            StandardWrite( "surf smoothing");
            int number_surface_points = 10000;
            float surf_point_x[ number_surface_points];
            float surf_point_y[ number_surface_points];
            float surf_point_z[ number_surface_points];


            MayerContactAndReentrance
            (
                    SIZE,
                    X_COO,
                    Y_COO,
                    Z_COO,
                    VDW_RADII,
                    PROBE_RADIUS,
                    NR_BINS[0],
                    NR_BINS[1],
                    NR_BINS[2],
                    DELTA,
                    0.5 * ( NR_BINS[ 0] - 1) * DELTA,
                    0.5 * ( NR_BINS[ 1] - 1) * DELTA,
                    0.5 * ( NR_BINS[ 2] - 1) * DELTA,
                    CENTER[0],
                    CENTER[1],
                    CENTER[2],
                    RMAY,
                    surf_point_x,     // sphere surface positions
                    surf_point_y,
                    surf_point_z,
                    number_surface_points
            );
        }
        else // cleanup grid if no smoothing was performed
        {
            std::replace( RMAY.begin(), RMAY.end(), -1, 1);
        }


//        Open( write, "rmay_grid_nreen.txt");
//        std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( write , " "));
//        write << std::endl;
//        Close( write);

//        int changed = 0;
//        for( std::vector< float>::const_iterator r_itr = RMAY.begin(), c_itr = copy.begin(); r_itr != RMAY.end() && c_itr != copy.end(); ++r_itr, ++c_itr)
//        {
//            if( *r_itr != *c_itr)
//            {
//                ++changed;
//            }
//        }
//        std::cout << changed << " positions in array have changed due to smoothing" << std::endl;


//        std::cout << "rmay: " << std::endl;
//        std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( std::cout , " "));
//        std::cout << std::endl;

    }


    float
    GetGridPoint
    (
            const math::Vector3N &POSITION,
            const std::vector< float> &RMAY,
            const math::Vector3N &MIN,
            const store::Vector3N< int> &NR_BINS,
            const float &DELTA
    )
    {
        std::vector< int> ids( 3);
        math::Vector3N tmp = (POSITION - MIN) / DELTA;
        for( int i = 0; i < 3; ++i)
        {
            ids[i] = size_t( tmp( i));
        }
        int id = math::TupleXD< float, 3>( NR_BINS).CalcID( ids);

        return RMAY[ id];
    }


    void
    CalcSurfGrid
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            const std::string &NAME,
            const size_t &NR_POINTS
    )
    {
        float value;
        std::vector< short>  rmay;
        math::Vector3N pos, min, max;
        store::Vector3N< int> nr_cells;
        float delta = 0.2;
        std::cout << __FUNCTION__ << ": build surf grid" << std::endl;

        std::string name = NAME.substr( 0, NAME.size() - 4);

        std::ofstream
            write;

        Open( write, name + "_grid.pdb");

        //        BuildSurfGrid( MOLECULE, rmay, min, nr_cells, delta);

        max[0] = min[0] + delta * float( nr_cells[0]);
        max[1] = min[1] + delta * float( nr_cells[1]);
        max[2] = min[2] + delta * float( nr_cells[2]);
        int cc = 0;
        for( int i = 0; i < int( nr_cells[0]); ++i)
            for( int j = 0; j < int( nr_cells[1]); ++j)
                for( int k = 0; k < int( nr_cells[2]); ++k, ++cc)
                {
                    int l = k + j * nr_cells[2] + i * nr_cells[2] * nr_cells[1];
//                    std::cout.width(4);
//                    std::cout << i << "  ";
//                    std::cout.width(4);
//                    std::cout << j << "  " ;
//                    std::cout.width(4);
//                    std::cout << k << "  ";
//                    std::cout.width(8);
//                    std::cout << l << "   ";
//                    std::cout.width(8);
//                    std::cout << cc << std::endl;
                    value = rmay[ l];
                    if( value == 0)
                    {
                        pos = min;
                        pos( 0) += (i + 0.5) * delta;
                        pos( 1) += (j + 0.5) * delta;
                        pos( 2) += (k + 0.5) * delta;
//                        geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(cc) / 1e5), cc);
                    }
                    if( value == -1)
                    {
                        pos = min;
                        pos( 0) += (i + 0.5) * delta;
                        pos( 1) += (j + 0.5) * delta;
                        pos( 2) += (k + 0.5) * delta;
//                        geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(cc) / 1e5), cc, "O");
                    }
                }


//        for( int i = 0; i < NR_POINTS; ++i)
//        {
//            pos.Randomize( min, max);
//            value = GetGridPoint( pos, rmay, min, nr_cells, delta);
//            if( value == 0)
//            {
//                geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(i) / 1e5), i);
//            }
//        }

        Close( write);
    }



} // end namespace geom




