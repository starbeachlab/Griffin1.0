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


#ifdef MPI_PARALLEL
#ifndef MPI_MOLECULE_ITERATOR_H_
#define MPI_MOLECULE_ITERATOR_H_

#include <mpi.h>
#include <boost/format.hpp>
#include "periodic_molecule_iterator.h"
//#include "../math/matrix3x3N.h"

namespace mol
{


	class MPIMoleculeIterator
	: public MoleculeIterator
	{

	private:
		MoleculeIterator		       m_MolIter;
		int 					       m_NrProcesses;
		int 					       m_ThisProcess;
		int 						   m_TotalSize;
		std::vector< int>			   m_ChunkSizes;
		std::vector< int>			   m_ChunkStart;
		std::vector< float>		       m_Forces;                     //!< linearized forces of all atoms
		std::vector< float>::iterator  m_ForceItr;
#ifdef TIMERS
		ContinuousTimer                                m_Timer;
#endif

	public:
		MPIMoleculeIterator( const MoleculeIterator &ITER) // ,  int ARGC, char *ARGV[])
		: MoleculeIterator(),
		m_MolIter( ITER),
		m_NrProcesses(),
		m_ThisProcess(),
		m_TotalSize( 0),
		m_ChunkSizes(),
		m_ChunkStart(),
		m_Forces(),
		m_ForceItr()
#ifdef TIMERS
		    , m_Timer( "MPIMoleculeIterator:inside-master")
#endif
		{
			StandardWrite( "\nopen MPI connection");

//			MPI_Init( &ARGC, &ARGV);   // mv to main
			MPI_Comm_size( MPI_COMM_WORLD, &m_NrProcesses);
			MPI_Comm_rank( MPI_COMM_WORLD, &m_ThisProcess);

			m_ChunkSizes.resize( m_NrProcesses, 0);
			m_ChunkStart.resize( m_NrProcesses - 1, 0);

			StandardWrite( "total nr processes: " << m_NrProcesses << " ID this process: " << m_ThisProcess);
		}


		virtual ~MPIMoleculeIterator()
		{
//			MPI_Finalize();   // mv to main
//			StandardWrite( "\nMPI connection closed");
		}

		virtual std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator
		    Begin()
		    {
			return m_MolIter.begin();
		    }

		
		virtual std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator
		    End()
		    {
			return m_MolIter.end();
		    }


		virtual void
		Reset
		(
#ifdef CPPIO
    			std::ostream &WRITE,
#else
    			FILE *WRITE,
#endif
			std::ostream &LOG = std::cout
		)
		{
			std::fill( m_Forces.begin(), m_Forces.end(), 0.0);
			m_ForceItr = m_Forces.begin();
			m_MolIter.Reset( WRITE, LOG);
#ifndef FULL_PARALLEL
			m_MolIter.SetAtomCount( Offset() + 1);
#endif
		        DebugPrint( LOG, __FUNCTION__ << " atomcount: " << m_MolIter.GetAtomCount());
		}


		virtual bool
		IterateMol( std::ostream &STREAM = std::cout)
		{
			return m_MolIter.IterateMol( STREAM);
		}


		virtual bool
		IterateAtom( std::ostream &STREAM = std::cout)
		{
			return m_MolIter.IterateAtom( STREAM);
		}


		virtual
		boost::shared_ptr< SimpleMolecule< Atom> > &
		GetThisMolecule()
		{
			return m_MolIter.GetThisMolecule();
		}


		virtual
		boost::shared_ptr< Atom> &
		GetThisAtom()
		{
			return m_MolIter.GetThisAtom();
		}


		virtual void
		SetMols( const store::ShPtrVec< SimpleMolecule< Atom> > &MOLS)
		{
			// pass all molecules to base of this class
			store::ShPtrVec< SimpleMolecule< Atom> >::Data() = MOLS;
			// pass parallel subset of molecules to member molitr, calc all needed values
			ParallelChunk( MOLS);
		}


		void ParallelChunk( const  store::ShPtrVec< SimpleMolecule< Atom> >&MOLS)
		{
			StandardWrite( __FUNCTION__ << " MPI setup process: " << m_ThisProcess);
			store::ShPtrVec< SimpleMolecule< Atom> >
				mols;

			std::vector< int>
				limits( m_NrProcesses),
				mol_sizes( MOLS.size(), 0);

			 std::vector< int>::iterator
				 limit_itr = m_ChunkStart.begin(),
				 size_itr = mol_sizes.begin();

			 int
				 mol_size = 0;

			 m_TotalSize = 0; // ??

			 for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = MOLS.begin(); mol_itr != MOLS.end(); ++mol_itr,  ++size_itr)
			 {
				 mol_size = ( *mol_itr)->GetAtoms().size();
				 m_TotalSize += mol_size;
				 *size_itr = mol_size;
			 }

			 // initial distribution: half the processes get a lower (rounded off), half the processes a rounded up number of mols

			 float
				 average =  float( m_TotalSize) / float( m_NrProcesses),
				 diff = 0.0,
			         prev_diff = 0.0;
			 int
				 process = 0,
				 total_size = 0,
				 prev_process_size = 0,
			         process_size = 0,
				 nr_mols,
			         count = 0,
			         sum = 0;

			 std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator
			     first_itr,
			     second_itr;


			 size_itr = mol_sizes.begin();

			 while( size_itr != mol_sizes.end() && process < m_NrProcesses)
			 {
				 total_size += *size_itr;
				 process_size += *size_itr;

				 diff = float( process_size) - average;

				 if( diff >= 0.0)
				 {
				     if( diff > -prev_diff)  // instabel, check molsize only little above average
				     {
						 m_ChunkSizes[ process] = prev_process_size;  // nr atoms of this chunk
						 process_size = *size_itr;                    // nr atoms for next chunk
						 *limit_itr++ = count;                        // mol id of chunk
				     }
				     else
				     {
						 m_ChunkSizes[ process] = process_size;
						 process_size = 0;
						 *limit_itr++ = count + 1;
				     }
				     ++process;
				 }
				 prev_diff = diff;
				 prev_process_size = process_size;
				 ++size_itr;
				 ++count;
			 }

			 // collect missing atoms (if m_NrProcesses condition was exceeded while there were still atoms)
			 while( size_itr != mol_sizes.end())
			 {
				 total_size += *size_itr;
				 process_size += *size_itr;
				 ++size_itr;
			 }
			 
			 // add missing atoms to last chunk
			 if( process_size > 0)
			 {
			     std::cout << "adjust last chunk from " << m_ChunkSizes.back() << " by " << process_size;
			     m_ChunkSizes.back() += process_size;
			     process_size = 0;
			     std::cout << " to " << m_ChunkSizes.back() << std::endl;
			 }

			 nr_mols = m_ChunkSizes[ m_ThisProcess];

			 StandardWrite( "total number of atoms: " << m_TotalSize << " (" << total_size << ")");
			 StandardWrite( "average number of atoms per chunk: " << average);
			 StandardWrite( "number of molecules: " << MOLS.size());
			 StandardWrite( "this chunk size: " << nr_mols);

			 std::cout << "distribution of atoms (chunk sizes): ";
			 std::copy( m_ChunkSizes.begin(), m_ChunkSizes.end(), std::ostream_iterator< int>( std::cout, "  "));
			 std::cout << std::endl;

			 std::cout << "chunk limits mol-ids: ";
			 std::copy( m_ChunkStart.begin(), m_ChunkStart.end(), std::ostream_iterator< int>( std::cout, "  "));
			 std::cout << std::endl;

			 // set m_Force to correct size;
			 m_Forces.resize( 3 * nr_mols, 0.0);

//			 mols.resize( nr_mols);

			 if( m_ThisProcess == 0)
			 {
			     StandardWrite( "first mol-id of this chunk: 0");
			     first_itr  = MOLS.begin();
			     second_itr = MOLS.begin() + m_ChunkStart[ 0];
			     mols.resize( m_ChunkStart[0]);
			 }
			 else if( m_ThisProcess == m_NrProcesses - 1)
			 {
			     StandardWrite( "first mol-id of this chunk: " << m_ChunkStart[ m_ThisProcess - 1]);
			     first_itr  = MOLS.begin() + m_ChunkStart[ m_ThisProcess - 1];
			     second_itr = MOLS.end();
			     mols.resize( MOLS.size() - m_ChunkStart[ m_ThisProcess - 1]);
			 }
			 else
			 {
			     StandardWrite( "first mol-id of this chunk: " << m_ChunkStart[ m_ThisProcess - 1]);
			     first_itr  = MOLS.begin() + m_ChunkStart[ m_ThisProcess - 1];
			     second_itr = MOLS.begin() + m_ChunkStart[ m_ThisProcess];
			     mols.resize( m_ChunkStart[ m_ThisProcess] - m_ChunkStart[ m_ThisProcess - 1]);
			 }

			 std::cout << "copy mols from " << first_itr - MOLS.begin() << " to " << second_itr - MOLS.begin()  << std::endl;

			 std::copy( first_itr, second_itr, mols.begin());

			 m_MolIter.SetMols( mols);

			 std::cout << mols.size() << " mols read" << std::endl;
			 for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = mols.begin(); mol_itr != mols.end(); ++mol_itr,  ++size_itr)
			 {
			     sum += ( *mol_itr)->GetAtoms().size();
			 }
			 std::cout << sum << " atoms read" << std::endl;
			 StandardWrite( "MPI setup complete\n");
		}



		virtual void CheckAndWriteForce
		(
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
				const math::Vector3N &FORCE,
				const util::FunctorEnum &TYPE = util::e_UNDEFINED_Functor
		)
		{

#ifdef FULL_PARALLEL		  
#ifdef TIMERS
			m_MolIter.m_WriteTimer.Start();
#endif
		    WriteForce( WRITE, FORCE);
#ifdef TIMERS
		    m_MolIter.m_WriteTimer.Stop();
#endif

#else
		    return m_MolIter.CheckAndWriteForce( STREAM, FORCE, TYPE);
#endif
		}

		void WriteForce
		(
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
				const math::Vector3N &FORCE
		)
		{
#ifdef FULL_PARALLEL		  
//			StandardWrite( __PRETTY_FUNCTION__);
//		    fprintf( WRITE, "force %f %f %f\n",FORCE[0],FORCE[1],FORCE[1]); 
		    *m_ForceItr++ = FORCE[0];
		    *m_ForceItr++ = FORCE[1];
		    *m_ForceItr++ = FORCE[2];
#else
			m_MolIter.WriteForce( WRITE, FORCE);
#endif
		}

		virtual void WriteZeroForces
		(
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
				const int &N
		)
		{
#ifdef FULL_PARALLEL		  
			StandardWrite( __PRETTY_FUNCTION__);
			// move iterator by 3*N positions
			m_ForceItr += 3 * N;
#endif
		}

		virtual store::ShPtrVec< SimpleMolecule< Atom> > &Mols()
		{
			StandardWrite( "===> " << __FUNCTION__ << " should not be called from " << __FILE__ << std::endl);
			return m_MolIter.Mols();
		}


//		std::istream &Read( std::istream &STREAM)
//		{
//			// move stream to correct position
//			// read subset
//
//			return STREAM;
//		}

		void WriteForceArray
		(
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
				const std::vector< float> &FORCES,
				int &COUNT
		)
		{
			std::vector< float>::const_iterator itr = FORCES.begin();
			while( itr != FORCES.end())
			{
#ifdef CPPIO
				WRITE << COUNT++ << "  0  " << *itr++ << "  ";
				WRITE << *itr++ << "  ";
				WRITE << *itr++ << std::endl;
#else
//				fprintf( WRITE, "%d 0 %f %f %f\n", COUNT++, *itr++, *itr++, *itr++);  ///  <=== THIS reverses order of output!!!!
				fprintf( WRITE, "%d 0 %f %f %f\n", COUNT++, *itr, *(itr+1), *(itr+2));
				itr += 3;
#endif
			}
		}

		void MergeDataAndWrite
		(
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
				std::ostream &LOG,
				const math::Matrix3x3N &VIRIAL,
				const float &ENERGY,
				const int &SFORCE_COUNT,
				const float &MAX_DIST,
				const int &REDIRECTED_COUNT,
				const float &MAX_REDIRECTED_DIST
		)
		{
#ifdef TIMERS
			m_Timer.Reset();
			m_Timer.Start();
#endif
		    MPI_Status status;
#ifdef FULL_PARALLEL		  
//		    LOG << __FUNCTION__ << " full parallel" << std::endl;
		    // daddy collects data and prints all 
		    if( m_ThisProcess == 0)
		    {

//			LOG << "head writes its stuff" << std::endl;

				math::Matrix3x3N
					virial,
					total_virial = VIRIAL;
				float
					energy,
					max_dist,
					max_redirected_dist,
					total_energy = ENERGY,
					total_max_dist = MAX_DIST,
					total_max_redirected_dist = MAX_REDIRECTED_DIST;
				int
					atom_count = 1,
					sforce_count,
					redirected_count,
				        size,
					total_sforce_count = SFORCE_COUNT,
					total_redirected_count = REDIRECTED_COUNT;

				std::vector< float>
//					kid_array,
					kid_force( 14);


				WriteForceArray( WRITE, m_Forces, atom_count);

//				LOG << "daddy wrote its stuff, count: " << atom_count << std::endl;

#ifdef TIMERS
				m_Timer.Intermediate();
				m_Timer.WriteStatus( LOG);
#endif

				for( int i = 1; i < m_NrProcesses; ++i)
				{
				    size = 3 * m_ChunkSizes[i] + 14;

				    float kid_array [ size];
//				    LOG << "head expects stuff from kid " << i << " datasize: " << size << std::endl;
					// receive array from kids and glue to total array
				    //	kid_array.resize( size);
				    MPI_Recv( &kid_array[0], size, MPI_FLOAT, i, 1000, MPI_COMM_WORLD, &status);

#ifdef TIMERS
				    m_Timer.Intermediate();
				    m_Timer.WriteStatus( LOG);
#endif

//				    LOG << "daddy got it" << std::endl;
//
//					int zeros = 0;
//					for( int i = 0; i < size; ++i)
//					{
//					    if( kid_array[i] == 0)
//					    {
//						++zeros;
//					    }
//					}
//					LOG << "zeros: " << zeros << " out of " << size << std::endl;


					kid_force.resize( size - 14);
					ArrayToValues( kid_array, kid_force, virial, energy, sforce_count, max_dist, redirected_count, max_redirected_dist);

					total_virial += virial;
//					LOG << "total energy before: " << total_energy << " plus " << energy;
					total_energy += energy;
//					LOG << " is " << total_energy << std::endl;
					total_sforce_count += sforce_count;
					if( max_dist > total_max_dist)
					{
						total_max_dist = max_dist;
					}
					total_redirected_count += redirected_count;
					if( max_redirected_dist > total_max_redirected_dist)
					{
						total_max_redirected_dist = max_redirected_dist;
					}
//					LOG << "daddy writes kids stuff, count: " << atom_count << std::endl;
					WriteForceArray( WRITE, kid_force, atom_count);
//					LOG << "daddy wrote kids stuff, count: " << atom_count << std::endl;
#ifdef TIMERS
					LOG << "chunk " << i << " written: ";
					m_Timer.Intermediate();
					m_Timer.WriteStatus( LOG);
#endif
				}

#ifdef CPPIO
				WRITE << boost::format( "%10.5f\n") %  total_energy;
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 0, 0) % total_virial( 0, 1) % total_virial( 0, 2);
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 1, 0) % total_virial( 1, 1) % total_virial( 1, 2);
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 2, 0) % total_virial( 2, 1) % total_virial( 2, 2);
#else
				fprintf( WRITE, "%f\n%f %f %f\n%f %f %f\n%f %f %f\n", total_energy, total_virial( 0, 0), total_virial( 0, 1), total_virial( 0, 2), total_virial( 1 , 0), total_virial( 1, 1), total_virial( 1, 2), total_virial( 2, 0), total_virial( 2, 1), total_virial( 2, 2));
#endif

		    	LOG << "#### buried-non-redirected: " << total_sforce_count
					<< " max-dist(surf): " << 	 total_max_dist
					<< " buried-redirected: " << total_redirected_count
					<< " max-dist(mol-geom-center,redirected): " << total_max_redirected_dist << std::endl;


			}
			// kids send their stuff in an array to daddy
			else
			{
//			    LOG << "kid " << m_ThisProcess << " sends its stuff" << std::endl;

			    int size = 3 * m_ChunkSizes[ m_ThisProcess] + 14;

			    float array[ size];

			    ValuesToArray( &array[0], VIRIAL, ENERGY, SFORCE_COUNT, MAX_DIST, REDIRECTED_COUNT, MAX_REDIRECTED_DIST);

//				int zeros = 0;
//				for( int i = 0; i < size; ++i)
//				{
//				    if( array[i] == 0)
//				    {
//					++zeros;
//				    }
//				}
//				LOG << "zeros: " << zeros << " out of " << size << std::endl;
				

				MPI_Send( &array[0], size, MPI_FLOAT, 0, 1000, MPI_COMM_WORLD);
//				LOG << "kid sent" << std::endl;				    
			}
#ifdef TIMERS
		    LOG << "done: ";
		    m_Timer.Intermediate();
		    m_Timer.WriteStatus( LOG);
#endif

#else    //  ifdef FULL_PARALLEL		  

		     LOG << __FUNCTION__ << " flat parallel" << std::endl;
		     // last process collects data and prints virial etc 
		     if( m_ThisProcess == m_NrProcesses - 1)
		     {
			 LOG << "last process writes virial, energy and stuff" << std::endl;

				math::Matrix3x3N
					virial,
					total_virial = VIRIAL;
				float
					energy,
					max_dist,
					max_redirected_dist,
					total_energy = ENERGY,
					total_max_dist = MAX_DIST,
					total_max_redirected_dist = MAX_REDIRECTED_DIST;
				int
					atom_count = 1,
					sforce_count,
					redirected_count,
				        size,
					total_sforce_count = SFORCE_COUNT,
					total_redirected_count = REDIRECTED_COUNT;

				std::vector< float>
					kid_force( 14);


				for( int i = 0; i < m_NrProcesses - 1; ++i)
				{
				    size = 14;

				    float kid_array [ size];
				    LOG << "last process expects stuff from process " << i << " datasize: " << size << std::endl;
					// receive array from kids and glue to total array
				    //	kid_array.resize( size);
				    MPI_Recv( &kid_array[0], size, MPI_FLOAT, i, 1000, MPI_COMM_WORLD, &status);

				    ArrayToValues( kid_array, virial, energy, sforce_count, max_dist, redirected_count, max_redirected_dist);

					total_virial += virial;
					LOG << "total energy before: " << total_energy << " plus " << energy;
					total_energy += energy;
					LOG << " is " << total_energy << std::endl;
					total_sforce_count += sforce_count;
					if( max_dist > total_max_dist)
					{
						total_max_dist = max_dist;
					}
					total_redirected_count += redirected_count;
					if( max_redirected_dist > total_max_redirected_dist)
					{
						total_max_redirected_dist = max_redirected_dist;
					}
				}


				WRITE << boost::format( "%10.5f\n") %  total_energy;
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 0, 0) % total_virial( 0, 1) % total_virial( 0, 2);
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 1, 0) % total_virial( 1, 1) % total_virial( 1, 2);
				WRITE << boost::format(  "%10.5f %10.5f %10.5f\n") % total_virial( 2, 0) % total_virial( 2, 1) % total_virial( 2, 2);


				LOG << "#### buried-non-redirected: " << total_sforce_count
					<< " max-dist(surf): " << 	 total_max_dist
					<< " buried-redirected: " << total_redirected_count
					<< " max-dist(mol-geom-center,redirected): " << total_max_redirected_dist << std::endl;

			}
			// kids send their stuff in an array to daddy
			else
			{
			    LOG << "process " << m_ThisProcess << " sends its stuff" << std::endl;

			    int size = 14;

			    float array[ size];

			    ValuesToArray( &array[0], VIRIAL, ENERGY, SFORCE_COUNT, MAX_DIST, REDIRECTED_COUNT, MAX_REDIRECTED_DIST);

			    MPI_Send( &array[0], size, MPI_FLOAT, m_NrProcesses - 1, 1000, MPI_COMM_WORLD);

			}
#endif   //  FULL_PARALLEL
		}

		void
		ArrayToValues
		(
				const float * KID_ARRAY,
#ifdef FULL_PARALLEL
				std::vector< float> &FORCE,
#endif
				math::Matrix3x3N &VIRIAL,
				float &ENERGY,
				int &SFORCE_COUNT,
				float &MAX_DIST,
				int &REDIRECTED_COUNT,
				float &MAX_REDIRECTED_DIST
		)
		{
//		    std::cout << __FUNCTION__ << std::endl;
			float * itr = (float *) &KID_ARRAY[0];
			VIRIAL = math::Matrix3x3N( itr);
			ENERGY = *itr++;
			SFORCE_COUNT = math::Round< int>( *itr++);
			MAX_DIST = *itr++;
			REDIRECTED_COUNT = math::Round< int>( *itr++);
			MAX_REDIRECTED_DIST = *itr++;

#ifdef FULL_PARALLEL
//			FORCE.resize( KID_ARRAY.size() - 14);
//			std::copy( itr, KID_ARRAY.end(), FORCE.begin());			
			for( std::vector< float>::iterator ftr = FORCE.begin(); ftr != FORCE.end(); ++ftr, ++itr)  /// todo: dangerous !
			{
			    *ftr = *itr;
			}
#endif

//		    std::cout  << __FUNCTION__ << " done" << std::endl;
		}

		void
		ValuesToArray
		(
		    float *ARRAY,
				const math::Matrix3x3N &VIRIAL,
				const float &ENERGY,
				const int &SFORCE_COUNT,
				const float &MAX_DIST,
				const int &REDIRECTED_COUNT,
				const float &MAX_REDIRECTED_DIST
		)
		{
//		        std::cout  << __FUNCTION__ << std::endl;
//			int size = 3 * m_ChunkSizes[ m_ThisProcess] + 14;
//			std::vector< float> array( size, 0.0);
			float * itr = &ARRAY[0];
			math::Linearize( VIRIAL, itr);
			*itr++ = ENERGY;
			*itr++ = float( SFORCE_COUNT);
			*itr++ = MAX_DIST;
			*itr++ = float( REDIRECTED_COUNT);
			*itr++ = MAX_REDIRECTED_DIST;

#ifdef FULL_PARALLEL
//				int zeros = 0;
				for( std::vector< float>::const_iterator ctr = m_Forces.begin(); ctr != m_Forces.end(); ++ctr, ++itr)
				{
				    *itr = *ctr;
//				    if( *ctr == 0)
//				    {
//					++zeros;
//				    }
				}
//				std::cout << m_ThisProcess << ": forces: zeros: " << zeros << " out of " << m_Forces.size() << std::endl;
				
//			std::copy( m_Forces.begin(), m_Forces.end(), itr);
//			std::cout << __FUNCTION__ << " done" << std::endl;
#endif
		}


		int Offset() const
		{
		    int size = 0;

		    for( int i = 0; i < m_ThisProcess; ++i)
		    {
			size += m_ChunkSizes[i];
		    }

		    return size;
		}
#ifdef CPPIO
		void PrepareStream( std::istream &STREAM) const
		{
		    std::string tmp;
		    if( m_ThisProcess == 0) // not head but first process !
		    {
			return;
		    }

		    int size = Offset();

//		    std::cout << "ignore first " << size << " lines in input file for kid " << m_ThisProcess << std::endl;

		    for( int i = 0; i < size; ++i)
		    {
			std::getline( STREAM, tmp);
		    }
//		    std::cout << "the last line ignored: " << std::endl << tmp << std::endl;
		}
#else
		void PrepareStream( FILE *STREAM) const
		{
			char buf[ 120];
			if( m_ThisProcess == 0)
			{
				return;
			}
			for( int i = 0; i < Offset(); ++i)
			{
				fgets( buf, 120, STREAM);
			}
		}
#endif


		virtual std::string GetClassName() const
		{
			return "MPIMoleculeIterator";
		}

#ifdef TIMERS

		virtual void WriteStatus( std::ostream &LOG = std::cout)
		{
			m_MolIter.WriteStatus( LOG);
			m_Timer.WriteStatus( LOG);
		}


		virtual void ResetTimers()
		{
			m_MolIter.ResetTimers();
			m_MolIter.ResetTimers();
		}
#endif

	}; // end class MPIMoleculeIterator




} // end namespace mol

#endif /* MPI_MOLECULE_ITERATOR_H_ */
#endif // MPI_PARALLEL
