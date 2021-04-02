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


#include "../../include/utilities/parallelization.h"

namespace util
{

//    bool ParallelizeProcess
//    (
//              math::Vector &MIN,
//              math::Vector &MAX,
//              const math::Vector &DELTA,
//              const size_t &NR_PROCESSES,
//              const size_t &PROCESS_ID
//    )
//    {
//        std::vector< std::vector< size_t> > possible_disections(  PossibleDevisors( MIN, MAX, DELTA));
//        std::vector< std::vector< size_t> > matching_disections(  MatchingDevisors( possible_disections, NR_PROCESSES));
//        if( matching_disections.size() == 0)
//        {
//            std::cout << "given limits are not disectible into " << NR_PROCESSES << " subprocesses" << std::endl;
//            std::cout << "possible disections are: " << std::endl;
//            size_t i( 0);
//            for( std::vector< std::vector< size_t> >::const_iterator itr( possible_disections.begin()); itr != possible_disections.end(); ++itr, ++i)
//            {
//                std::cout << "dimension: " << i << std::endl;
//                for( std::vector< size_t>::const_iterator m_itr( itr->begin()); m_itr != itr->end(); ++m_itr)
//                {
//                    std::cout << *m_itr << "  ";  // math::Vector( *itr).ProductOfElements() << std::endl;
//                }
//                std::cout << std::endl;
//            }
//            return false;
//        }
//        else
//        {
//            std::vector< size_t> balanced( SelectMostBalancedDevisors( matching_disections));
//            AdjustLimits( MIN, MAX, balanced, PROCESS_ID);
//            return true;
//        }
//    }

    bool ParallelizeProcess
    (
            math::Vector3N &MIN,
            math::Vector3N &MAX,
            const math::Vector3N &DELTA,
            const size_t &NR_PROCESSES,
            const size_t &PROCESS_ID
    )
    {
    	math::Vector3N
//			min = MIN,
			max = MAX;
        std::vector< std::vector< size_t> >    big_bin( 3);
        size_t count( 0);
        for( size_t i = 0; i < 3; ++i)
        {
            size_t nr_bins =  math::Round< size_t>( ( MAX[i] - MIN[i]) / DELTA[i]);
            if( fabs( float( nr_bins) * DELTA[i] - ( MAX[i] - MIN[i])) > 1e-3)
            {
            	std::cout << "===> limits not well adjusted: " << float( nr_bins) * DELTA[i] << " vs " << ( MAX[i] - MIN[i]) << std::endl;
            }
            float bins_per_process = float( nr_bins) / float( NR_PROCESSES);
            int n_extra = math::Round< int>( bins_per_process - int( bins_per_process)) * NR_PROCESSES;

            DebugWrite( i << " nr_bins: " << nr_bins << " disect: " << bins_per_process << " n_extra: " << n_extra);

            std::vector< size_t> bins( NR_PROCESSES);

            for( size_t j = 0; j < NR_PROCESSES; ++j)
            {
                if( (int) j < n_extra)
                {
                    bins[ j] = int( bins_per_process) + 1;
                }
                else
                {
                    bins[ j] = int( bins_per_process);
                }
                DebugWriteNoFlush( bins[j]);
            }
            DebugWrite( "");

#ifdef DEBUG
            float mini = MIN[i];
            for( size_t j = 0; j < NR_PROCESSES; ++j)
            {
                std::cout << count++ << ":  " << mini << " ";
                mini += bins[j] * DELTA[ i];
                std::cout << mini << std::endl;
            }
#endif

            big_bin[i] = bins;
        }

        count = 0;
        math::Vector3N vals;
        vals( 2) = MIN[2];
        for( size_t i = 0; i < NR_PROCESSES; ++i)
        {
            vals( 1) = MIN[1];
            for( size_t j = 0; j < NR_PROCESSES; ++j)
            {
                vals( 0) = MIN[0];
                for( size_t k = 0; k < NR_PROCESSES; ++k)
                {
                    if( count == PROCESS_ID)
                    {
                        for( size_t l = 0; l < 3; ++l)
                        {
                            MIN( l) = vals( l);
                            MAX( l) = std::min( max( l), vals( l) + big_bin[ l][k] * DELTA[l]);
                        }
                        return true;
                    }
#ifdef DEBUG
                    std::cout << "==> " << count << std::endl;
                    std::cout << vals( 0) << "  " << vals( 0) + big_bin[ 0][k] * DELTA[ 0] << std::endl;
                    std::cout << vals( 1) << "  " << vals( 1) + big_bin[ 1][j] * DELTA[ 1] << std::endl;
                    std::cout << vals( 2) << "  " << vals( 2) + big_bin[ 2][i] * DELTA[ 2] << std::endl;
#endif
                    ++count;
                    vals( 0) += big_bin[ 0][k] * DELTA[0];
                }
                vals( 1) += big_bin[ 1][j] * DELTA[1];
            }
            vals( 2) += big_bin[ 2][ i] * DELTA[ 2];
        }
        return false;
    }


  std::vector< std::vector< size_t> > PossibleDevisors( math::Vector &MIN, math::Vector &MAX, const math::Vector &DELTA)
  {
    std::vector< std::vector< size_t> > result( MIN.size());
    size_t nr;
    for( size_t i( 0); i < MAX.size(); ++i)
      {
    nr = math::Round< size_t>( ( MAX(i) - MIN(i)) / DELTA( i)); // check value??
    //    std::cout << "total number of bins: " << nr << std::endl;
    std::vector< size_t> devisors( math::Devisors<size_t>( nr));
    //    std::cout << devisors << std::endl;
    result[ i] = devisors;
      }
    return result;
  }


  std::vector< std::vector< size_t> > MatchingDevisors( const std::vector< std::vector< size_t> > &SUB_BLOCKS, const size_t &NR_PROCESSES)
  {
    std::vector< std::vector< size_t> > combinations;
    for( std::vector< size_t>::const_iterator a_itr( SUB_BLOCKS[0].begin()); a_itr != SUB_BLOCKS[0].end(); ++a_itr)
      for( std::vector< size_t>::const_iterator b_itr( SUB_BLOCKS[1].begin()); b_itr != SUB_BLOCKS[1].end(); ++b_itr)
    for( std::vector< size_t>::const_iterator c_itr( SUB_BLOCKS[2].begin()); c_itr != SUB_BLOCKS[2].end(); ++c_itr)
      if( *a_itr * *b_itr * *c_itr == NR_PROCESSES)
        {
          std::vector< size_t> vec( 3);
          vec[0] = *a_itr;
          vec[1] = *b_itr;
          vec[2] = *c_itr;
          combinations.push_back( vec);
          //          std::cout << "match: " << vec;
        }
    return combinations;
  }


  std::vector< size_t> SelectMostBalancedDevisors( const std::vector< std::vector< size_t> > &SUB_BLOCKS)
  {
    float best( std::numeric_limits< float>::max());
    std::vector< size_t> favorit;
    float tmp;
    for( std::vector< std::vector< size_t> >::const_iterator itr = SUB_BLOCKS.begin(); itr != SUB_BLOCKS.end(); ++itr)
      {
    tmp = math::Vector( *itr).RmsdOfElements().second;
    if( tmp < best)
      {
        best = tmp;
        favorit = *itr;
      }
      }
    //    std::cout << "mr balance: " << favorit;
    return favorit;
  }


  void AdjustLimits( math::Vector &MIN, math::Vector &MAX, const std::vector< size_t> &SUB_BLOCK, const size_t &PROCESS_ID)
  {
      float delta;
      std::vector< size_t> ids( NrToIDs( SUB_BLOCK, PROCESS_ID));
      for( size_t i( 0); i < MIN.size(); ++i)
      {
          delta = ( MAX(i) - MIN(i)) / float( SUB_BLOCK[ i]);
          MIN(i) += ids[ i] * delta;
        MAX(i) = MIN(i) + delta;
      }
  }

  std::vector< size_t> NrToIDs( const std::vector< size_t> &LIMITS, const size_t &NR)
  {
      size_t value( NR);
      std::vector< size_t> ids( LIMITS.size());
      for( size_t i = 0; i < LIMITS.size(); ++i)
      {
          size_t size( (size_t) math::Vector( LIMITS).SubProductOfElements( i + 1, LIMITS.size() - i - 1));
          size_t id( value / size);
          value -= id * size;
          //    std::cout << "size: " << size << " id: " << id << " value: " << value << std::endl;
          ids[i] = id;
      }
      return ids;
  }
    // TEST
//    size_t count( 0);
//    for( size_t i = 0; i < LIMITS[0]; ++i)
//      for( size_t j = 0; j < LIMITS[1]; ++j)
//    for( size_t k = 0; k < LIMITS[2]; ++k)
//      {
//        if( count == NR)
//          {
//        std::vector< size_t> vec( 3);
//        vec[0] = i;
//        vec[1] = j;
//        vec[2] = k;
//        std::cout << "score! " << count << "  " << vec;
//        return vec;
//
//          }
//        ++count;
//      }
//    return std::vector< size_t>();
//  }


//    size_t remain( PROCESS_ID);
//    size_t value;
//    float delta;
//    value = size_t( math::Round( float( remain) /  math::Vector( SUB_BLOCK).SubProductOfElements( 0, MIN.size() - i - 1)));
//    remain -= value * math::Vector( SUB_BLOCK).SubProductOfElements( 0, MIN.size() - i - 1);
//    std::cout << "delta: " << delta << " value: " << value <<  " remain: " << remain << " subprod: " << math::Vector( SUB_BLOCK).SubProductOfElements( 0, MIN.size() - i - 1) <<  std::endl;

} // namespace util
