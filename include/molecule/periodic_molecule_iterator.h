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


#ifndef PERIODIC_MOLECULE_ITERATOR_H_
#define PERIODIC_MOLECULE_ITERATOR_H_

#include "molecule_iterator.h"

namespace mol
{

	class PeriodicMoleculeIterator
	: public MoleculeIterator
	{
	private:
		store::ShPtrVec< SimpleMolecule< Atom> >                 m_NonZeroForceMols;
		std::vector< int>                                        m_ZeroCounter;
		std::vector< int>::iterator                     		 m_CountIter;
		int 										             m_Frequency;
		int                                                      m_Round;
		bool 													 m_Update;
		boost::shared_ptr< SimpleMolecule< Atom> >          	 m_ThisMol;

	public:
		PeriodicMoleculeIterator( store::ShPtrVec< SimpleMolecule< Atom> > &MOL, const int FREQUENCY)
		: MoleculeIterator( MOL),
		  m_NonZeroForceMols(),  // omit zero force mols
		  m_ZeroCounter( 1),
		  m_CountIter(),
		  m_Frequency( FREQUENCY),
		  m_Round( 0),
		  m_Update( true),
		  m_ThisMol()
		  {
			m_CountIter = m_ZeroCounter.begin();
			StandardWrite( __PRETTY_FUNCTION__);
		  }

		PeriodicMoleculeIterator( const int FREQUENCY)
		: MoleculeIterator(),
		  m_NonZeroForceMols(),
		  m_ZeroCounter( 1),
		  m_CountIter(),
		  m_Frequency( FREQUENCY),
		  m_Round( 0),
		  m_Update( true),
		  m_ThisMol()
		  {
			m_CountIter = m_ZeroCounter.begin();
		  }


		virtual ~PeriodicMoleculeIterator(){}

		virtual std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator
		Begin()
		{
			return begin();
		}

		
		virtual std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator
		End()
		{
			return end();
		}



		virtual void Reset
		(
#ifdef CPPIO
    		std::ostream &WRITE,
#else
    		FILE *WRITE,
#endif
				std::ostream &LOG = std::cout
		)
		{
		    //  std::cout << __FUNCTION__ << std::endl;
		    //  std::cout << "all: " << MoleculeIterator::size() << std::endl;
		    //  std::cout << " nonzeromols: " << m_NonZeroForceMols.size() << std::endl;
			if( m_Round % m_Frequency == 0)
			{
				// recalculate all forces
				StandardPrint( LOG, __FUNCTION__ << ": recalculate all forces");
				m_Update = true;
				m_NonZeroForceMols.clear();
				MoleculeIterator::Reset( WRITE, LOG);
				m_ZeroCounter = std::vector< int>( 1, 0);
				m_CountIter = m_ZeroCounter.begin();
				//  std::cout << "start diff: " << ( m_CountIter - m_ZeroCounter.begin()) << std::endl;
				m_ThisMol = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>());
			}
			else
			{
				if( m_Round % m_Frequency == 1)
				{
					int nr_non_zero = 0;
					for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator itr = m_NonZeroForceMols.begin(); itr != m_NonZeroForceMols.end(); ++itr)
					{
						nr_non_zero += ( *itr)->GetAtoms().size();
					}
					int zero = 0;
					for( std::vector< int>::const_iterator itr = m_ZeroCounter.begin(); itr != m_ZeroCounter.end(); ++itr)
					{
						zero += *itr;
					}
					StandardPrint( LOG, ">>>> atoms after update: non zero force: " << nr_non_zero << " zero force: " << zero << " total: " << zero + nr_non_zero);
				}
				StandardPrint( LOG, __FUNCTION__ << ": recalculate only non-zero forces");
				m_Update = false;
				m_MolIter = m_NonZeroForceMols.begin();
				m_AtomIter = ( *m_MolIter)->Atoms().begin();
				m_CountIter = m_ZeroCounter.begin();
				m_AtomCount = 0;
				WriteZeroForces( WRITE, *m_CountIter);
				++m_AtomCount;
				++m_CountIter;
			}
			++m_Round;
		}


		virtual bool IterateMol( std::ostream &STREAM = std::cout)  // todo: REORDER !!!
		{
			//	    static int jo = 0;
			//  	std::cout << __FUNCTION__ << " b: " << Mols().size() << " call: " << jo << std::endl;
#ifdef TIMERS
			m_Timer.Start();
#endif
			++m_MolIter;
			if( !m_Update && m_MolIter == m_NonZeroForceMols.end())
			{
			    //  	std::cout << __FUNCTION__ << " ar: " << Mols().size() << std::endl;
#ifdef TIMERS
				m_Timer.Stop();
#endif
				return false;
			}
			else if( m_Update)
			{
				if( m_ThisMol->GetAtoms().size() > 0)
				{
					m_NonZeroForceMols.push_back( m_ThisMol);
					m_ThisMol = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>());
					//  	std::cout << __FUNCTION__ << " au: " << Mols().size() << std::endl;
				}
				if( m_MolIter == store::ShPtrVec< SimpleMolecule< Atom> >::end())
				{
				    //  	std::cout << __FUNCTION__ << " af: " << Mols().size() << std::endl;
#ifdef TIMERS
					m_Timer.Stop();
#endif
					return false;
				}
			}
			m_AtomIter = ( *m_MolIter)->Atoms().begin();
			// 	 	std::cout << __FUNCTION__ << " a: " << Mols().size() << " call: " << jo << std::endl;
#ifdef TIMERS
			m_Timer.Stop();
#endif
			return true;
		}



		virtual bool IterateAtom()
		{
#ifdef TIMERS
			m_Timer.Start();
#endif
			++m_AtomIter;
			++m_AtomCount;

//			if( !m_Update)
//			{
//				WriteZeroForces( STREAM, *m_CountIter);
//				++m_CountIter;
//			}
			if( m_AtomIter == ( *m_MolIter)->GetAtoms().end())
			{
#ifdef TIMERS
				m_Timer.Stop();
#endif
				return false;
			}
#ifdef TIMERS
			m_Timer.Stop();
#endif
			return true;
		}



//		virtual
//		boost::shared_ptr< SimpleMolecule< Atom> > &
//		GetThisMolecule()
//		{
//			return *m_MolIter;
//		}



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
#ifdef TIMERS
			m_WriteTimer.Start();
#endif
//			std::cout << __FUNCTION__ << ": " << FORCE[0] << " " << FORCE[1] << "  " << FORCE[2] << std::endl;
			WriteForce( WRITE, FORCE);
			//  std::cout << __FUNCTION__ << " b: " << Mols().size() << " " << FORCE[0] << " " << FORCE[1] << " " << FORCE[2] << std::endl;
			if( m_Update)
			{
				if( !FORCE.AreAllElements( 0.0) || TYPE == util::e_ConstantGridPoint) // add non-zero forces and buried atoms
				{
				    //	  std::cout << __FUNCTION__ << " 0: " << Mols().size() << " atoms: " << m_ThisMol->Atoms().size()  << std::endl;
					m_ThisMol->Atoms().push_back( *m_AtomIter);
				    // 	 std::cout << __FUNCTION__ << " 1: " << Mols().size() << " atoms: " << m_ThisMol->Atoms().size()  << " size: " << m_ZeroCounter.size() << std::endl;
					m_ZeroCounter.push_back( 0);
				    // 	 std::cout << __FUNCTION__ << " 2: " << Mols().size() << " size: " << m_ZeroCounter.size()  << " diff: " << ( m_CountIter -  m_ZeroCounter.begin()) << std::endl;
					m_CountIter = m_ZeroCounter.begin() + ( m_ZeroCounter.size() - 1);
					// 	 std::cout << __FUNCTION__ << " 3: " << Mols().size() << " diff: " << ( m_CountIter -  m_ZeroCounter.begin()) << std::endl;
				}
				else
				{
					//				    std::cout << __FUNCTION__ << " 4: " << Mols().size() << " " << *m_CountIter << " diff: " << ( m_CountIter -  m_ZeroCounter.begin()) << std::endl;
					++( *m_CountIter);
					//				    std::cout << __FUNCTION__ << " 5: " << Mols().size() << " " << *m_CountIter << " diff: " << ( m_CountIter -  m_ZeroCounter.begin()) << std::endl;
				}
			}
			else
			{
				WriteZeroForces( WRITE, *m_CountIter);
				++m_CountIter;
			}
			//  std::cout << __FUNCTION__ << " a: " << Mols().size() << std::endl;

#ifdef TIMERS
			m_WriteTimer.Stop();
#endif
		}


#ifdef CPPIO
		virtual void WriteZeroForces
		(
        		std::ostream &WRITE,
				const int &N
		)
		{
			for( int i = 0; i < N; ++i)
			{
				WRITE << ++m_AtomCount << " 0  0.0  0.0  0.0" << std::endl;
			}
		}
#else
		virtual void WriteZeroForces
		(
        		FILE *WRITE,
				const int &N
		)
		{
			for( int i = 0; i < N; ++i)
			{
				fprintf( WRITE, "%d 0  0.0  0.0  0.0\n", ++m_AtomCount);
			}
		}
#endif

		virtual const store::ShPtrVec< SimpleMolecule< Atom> > &
		GetNonZeroForceMols() const
		{
			return m_NonZeroForceMols;
		}

		const std::vector< int> &GetZeroCounter() const
		{
			return m_ZeroCounter;
		}

		virtual std::string GetClassName() const
		{
			return "PeriodicMoleculeIterator";
		}

	}; // end class MoleculeIterator


} // end namespace mol




#endif /* PERIODIC_MOLECULE_ITERATOR_H_ */
