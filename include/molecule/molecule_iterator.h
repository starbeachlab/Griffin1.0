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


#ifndef MOLECULE_ITERATOR_H_
#define MOLECULE_ITERATOR_H_

//#include "../storage/shared_pointer_vector.t.h"
//#include "simple_molecule.t.h"

namespace mol
{

	class MoleculeIterator
	: public store::ShPtrVec< SimpleMolecule< Atom> >
	{
	protected:
		store::ShPtrVec< SimpleMolecule< Atom> >::iterator   m_MolIter;
		store::ShPtrVec< Atom>::iterator                     m_AtomIter;
		int                                                  m_AtomCount;
	public:
#ifdef TIMERS
		SumTimer											 m_Timer;
		SumTimer                                             m_WriteTimer;
#endif

	public:
		MoleculeIterator()
		: store::ShPtrVec< SimpleMolecule< Atom> >(),
		  m_MolIter(),
		  m_AtomIter(),
		  m_AtomCount( 1)
#ifdef TIMERS
		  , m_Timer( "MoleculeIterator")
		  , m_WriteTimer( "MoleculeIterator:write-force-timer")
#endif
	    {
	    	DebugWrite( __PRETTY_FUNCTION__);
			  m_MolIter = this->store::ShPtrVec< SimpleMolecule< Atom> >::begin();
//			  m_AtomIter = ( *store::ShPtrVec< SimpleMolecule< Atom> >::begin())->Atoms().begin();
	    }

		MoleculeIterator( store::ShPtrVec< SimpleMolecule< Atom> > &MOL)
		: store::ShPtrVec< SimpleMolecule< Atom> >( MOL),
		  m_MolIter( MOL.begin()),
		  m_AtomIter( ( *MOL.begin())->Atoms().begin()),
		  m_AtomCount( 1)
#ifdef TIMERS
		  , m_Timer( "MoleculeIterator")
		  , m_WriteTimer( "MoleculeIterator:write-force-timer")
#endif
	    {
			DebugWrite( __PRETTY_FUNCTION__);
	    }


		virtual ~MoleculeIterator()
		{
#ifdef TIMERS
			m_Timer.WriteStatus();
			m_WriteTimer.WriteStatus();
#endif
		}


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



		virtual store::ShPtrVec< SimpleMolecule< Atom> > &Mols()
		{
			return *this;
		}

		void SetAtomCount( const int &COUNT)
		{
		    m_AtomCount = COUNT;
		}


		int GetAtomCount() const
		{
		    return m_AtomCount;
		}


		virtual void SetMols( const store::ShPtrVec< SimpleMolecule< Atom> > &MOLS)
		{
			store::ShPtrVec< SimpleMolecule< Atom> >::Data() = MOLS;
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
			DebugPrint( LOG, "plain calculation of all forces");
			m_MolIter = this->store::ShPtrVec< SimpleMolecule< Atom> >::begin();
			m_AtomIter = ( *m_MolIter)->Atoms().begin();
			m_AtomCount = 1;
		}


		virtual bool IterateMol( std::ostream &STREAM = std::cout)
		{
		    DebugPrint( STREAM, __FUNCTION__);
#ifdef TIMERS
			m_Timer.Start();
#endif
			++m_MolIter;
			if( m_MolIter == store::ShPtrVec< SimpleMolecule< Atom> >::end())
			{
//			    STREAM << __FUNCTION__ << " final mol" << std::endl;
				return false;
			}
			m_AtomIter = ( *m_MolIter)->Atoms().begin();
#ifdef TIMERS
			m_Timer.Stop();
#endif
			return true;
		}


		virtual bool IterateAtom( std::ostream &STREAM = std::cout)
		{
		    DebugPrint( STREAM, __FUNCTION__);
#ifdef TIMERS
			m_Timer.Start();
#endif
			++m_AtomIter;
			++m_AtomCount;
			if( m_AtomIter == ( *m_MolIter)->GetAtoms().end())
			{
#ifdef TIMERS
				m_Timer.Stop();
#endif
//				STREAM << __FUNCTION__  << " isatendofmol" << std::endl;
				return false;
			}
#ifdef TIMERS
			m_Timer.Stop();
#endif
			return true;
		}


#ifdef TIMERS

		virtual void WriteStatus( std::ostream &LOG = std::cout)
		{
			m_Timer.WriteStatus( LOG);
			m_WriteTimer.WriteStatus( LOG);
		}


		virtual void ResetTimers()
		{
			m_Timer.Reset();
			m_WriteTimer.Reset();
		}
#endif


		virtual
		boost::shared_ptr< SimpleMolecule< Atom> > &
		GetThisMolecule()
		{
			return *m_MolIter;
		}

		virtual
		boost::shared_ptr< Atom> &
		GetThisAtom()
		{
			return *m_AtomIter;
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
#ifdef TIMERS
			m_WriteTimer.Start();
#endif
//			std::cout << __FUNCTION__ << ": " << FORCE[0] << " " << FORCE[1] << "  " << FORCE[2] << std::endl;
			WriteForce( WRITE, FORCE);
#ifdef TIMERS
			m_WriteTimer.Stop();
#endif
		}


#ifdef CPPIO
		void WriteForce
		(
        		std::ostream &WRITE,
				const math::Vector3N &FORCE
		)
		{
			WRITE << m_AtomCount << "  0  ";

			for( std::vector< float>::const_iterator force_itr = FORCE.begin(); force_itr != FORCE.end(); ++force_itr)
			{
				WRITE << *force_itr << " ";
			}
			WRITE << std::endl;

//			++m_AtomCount;
		}
#else
		void WriteForce
		(
				FILE *WRITE,
				const math::Vector3N &FORCE
		)
		{
//			if( !WRITE)
//			{
//			    std::cout << "====> FILESTREAM closed !!!" << std::endl;
//			}
//			std::cout << __FUNCTION__ << ": " << FORCE[0] << " " << FORCE[1] << "  " << FORCE[2] << std::endl;
			fprintf( WRITE, "%d 0 %f %f %f\n", m_AtomCount, FORCE[0], FORCE[1], FORCE[2]);
		}
#endif

//		virtual int GetAndIncrementAtomCount()
//		{
//			return m_AtomCount++;
//		}

		virtual std::string GetClassName() const
		    {
			return "MoleculeIterator";
		    }

	};


} // end namespace mol



#endif /* ITERATION_MOLECULE_H_ */
