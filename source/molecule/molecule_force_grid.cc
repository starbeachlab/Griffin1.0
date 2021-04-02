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


#include "../../include/molecule/molecule_force_grid.h"


namespace mol
{
    //! copy constructor
    MoleculeForceGrid *MoleculeForceGrid::Clone() const{ return new MoleculeForceGrid( *this);}

    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    //! returns a vector of forces for molecule
    boost::shared_ptr< std::vector< math::Vector3N> > & MoleculeForceGrid::operator()( const boost::shared_ptr< SimpleMolecule< Atom> > &MOL) const
    {
        boost::shared_ptr< std::vector< math::Vector3N> > result( new std::vector< math::Vector3N>( MOL->GetAtoms().size()));
        std::vector< math::Vector3N >::iterator result_itr( result->begin());
        std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr( MOL->GetAtoms().begin());
        for( ; result_itr != result->end() && atom_itr != MOL->GetAtoms().end(); ++atom_itr, ++result_itr)
        {
            *result_itr = m_Grid->operator()( **atom_itr);
        }
        return math::Function< boost::shared_ptr< SimpleMolecule< Atom> >, boost::shared_ptr< std::vector< math::Vector3N> > >::s_Tmp = result;
    }

    // THIS SHOULD NOT BE HERE, THIS SHOULD NOT BE ABOUT
    //! returns a vector of forces for shared pointer vector of molecules
    boost::shared_ptr< std::vector< math::Vector3N> > MoleculeForceGrid::operator()( const boost::shared_ptr< store::ShPtrVec< SimpleMolecule< Atom> > > &MOL) const
    {
        boost::shared_ptr< std::vector< math::Vector3N> > result;
        for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator itr = MOL->begin(); itr != MOL->end(); ++itr)
        {
            boost::shared_ptr< std::vector< math::Vector3N> > tmp( this->operator()( *itr));
            result->insert( result->end(), tmp->begin(), tmp->end());
        }
        return result;
    }


    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////



    int MoleculeForceGrid::Execute
    (
            const std::string &INPUT_FILE,
            const std::string &OUTPUT_FILE,
            boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
            std::ostream &LOG
    )
    {
        DebugWrite( __PRETTY_FUNCTION__);
    	Timer timer( __FUNCTION__, LOG);

//        LOG << __FUNCTION__ << std::endl;

#ifdef CPPIO
        std::ifstream read;
        std::ofstream write;
        read.open( INPUT_FILE.c_str());
#else
        FILE *read = fopen( INPUT_FILE.c_str(), "r");
        FILE *write = NULL;
#endif
        if( read)
        {

#ifdef MPI_PARALLEL
	    mol::MPIMoleculeIterator * mpi_itr = ( mol::MPIMoleculeIterator *) MOL_ITR.get();
	    mpi_itr->PrepareStream( read);
#endif  // MPI_PARALLEL
	    
            if( m_ForceFieldType ==  "gromacs")
            {
            	Timer read_timer( "read new coor gromacs", LOG);
                LOG << "new coordinates Gromacs" << std::endl;

                for( std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator itr = MOL_ITR->Begin(); itr != MOL_ITR->End(); ++itr)
                {
                    file::ReadCoordinatesFromXyz( read, *itr, LOG);
                }
    //                file::ReadNewCoordinatesFromGmxDump( read, MOL, LOG);
            }
            else
            {
//                LOG << "new coordinates charmm/namd" << std::endl;
                Timer read_timer( "read new coordinates in namd format", LOG);
                for( std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator itr = MOL_ITR->Begin(); itr != MOL_ITR->End(); ++itr)
                {
                    file::ReadNewPositionsFromNamdFile( read, *itr);
                }
                //                    file::ReadNewCoordinatesFromPdb( read, MOL, LOG);
            }


#ifdef MPI_PARALLEL

//	    LOG << "position of first atom: " << ( *MOL_ITR->Begin())->GetAtoms()( 0)->GetPosition() << std::endl;
//	    LOG << "first atom: " << ( *MOL_ITR->Begin())->GetAtoms()( 0) << std::endl;
	    int process;
	    MPI_Comm_rank( MPI_COMM_WORLD, &process);
#endif


#ifdef FULL_PARALLEL
	    if( process == 0)
	    {
#endif


#ifdef CPPIO
	    	Open( write, OUTPUT_FILE, LOG);
#else
	    	write = fopen( OUTPUT_FILE.c_str(), "w");
#endif


#ifdef FULL_PARALLEL
	    }
#endif


//                LOG << "calculate forces, energy, virial ..." << std::endl;
	    	CalculateAndWriteForcesEnergyVirial( write, MOL_ITR, LOG);


#ifdef FULL_PARALLEL
		if( process == 0)
 		{
#endif
#ifdef CPPIO
          	Close( write);
#else
		fclose( write);
#endif
#ifdef FULL_PARALLEL
		}
#endif


//		LOG << "write closed" << std::endl;
        }
        else
        {
            LOG << "===> input <" << INPUT_FILE + "> could not be opened!" << std::endl;
        }
#ifdef CPPIO
        Close( read);
#else
	fclose( read);
#endif
        return 0;
    }




    void MoleculeForceGrid::CalculateAndWriteForcesEnergyVirial
    (
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
          boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
  		  std::ostream &LOG
    )
    {
    	DebugPrint( LOG, __FUNCTION__);
    	Timer timer( (std::string) __FUNCTION__ + " all molecules and merge (if parallel)", LOG);

#ifdef TIMERS
    	ResetTimers();
    	MOL_ITR->ResetTimers();
#endif

    	math::Matrix3x3N virial;
    	float energy = 0.0;
    	float max_squared_dist = 0.0;
    	float max_squared_redirected_dist = 0.0;
    	int sforce_count = 0;
    	int redirected = 0;

#ifdef ANALYZE_BURIAL
    	static int
			s_Count = 0;
    	math::Histogram< float>
			surf_distances( 0.0, 0.5, 60);
    	math::Histogram< float>
			redirected_distances( 0.0, 0.5, 40);
#endif   // ANALYZE_BURIAL

    	MOL_ITR->Reset( WRITE, LOG);  // todo: useless first stream

    	do
    	{
    		CalculateAndWriteForcesEnergyVirial
    		(
    			WRITE,
    			MOL_ITR,
				virial,
				energy,
				max_squared_dist,
				max_squared_redirected_dist,
				sforce_count,
				redirected,
				LOG
#ifdef ANALYZE_BURIAL
				, redirected_distances,
				surf_distances
#endif
		    );

        } while( MOL_ITR->IterateMol( LOG));

//	std::cout << "now: all: " << MOL_ITR->Mols().size() << std::endl;
//	if( MOL_ITR->GetClassName() == mol::PeriodicMoleculeIterator( 0).GetClassName())
//	{
//	    mol::PeriodicMoleculeIterator * per = ( mol::PeriodicMoleculeIterator *) MOL_ITR.get();
//	    std::cout << "nonzeros: " << per->GetNonZeroForceMols().size() << std::endl;
//	    int sum = 0;
//	    for( std::vector< int>::iterator itr = per->GetZeroCounter().begin(); itr !=per->GetZeroCounter().end(); ++itr)
//		{
//		    sum += *itr;
//		}
//	    std::cout << "zeros: " << sum  << "/" << per->GetZeroCounter().size() << std::endl;
//	}


#ifdef MPI_PARALLEL  /// todo: this section should be merged with the mpi-mol-iter->CalcAndWriteForces() 

#ifdef TIMERS
    	{
        	Timer mpi_timer( (std::string) __FUNCTION__ + ": mpi merge: ", LOG);
#endif

		mol::MPIMoleculeIterator * mpi_itr = ( mol::MPIMoleculeIterator *) MOL_ITR.get();
		mpi_itr->MergeDataAndWrite( WRITE, LOG, virial, energy, sforce_count, sqrt( max_squared_dist), redirected, sqrt( max_squared_redirected_dist));

#ifdef TIMERS
    	}
#endif

#else  // MPI_PARALLEL

#ifdef TIMERS
    	{
        	Timer mpi_timer( (std::string) __FUNCTION__ + ": write energy and virial: ", LOG);
#endif

#ifdef CPPIO
    	char buffer[50];
    	sprintf( buffer, "%10.5f", energy);
    	WRITE << buffer << std::endl;
    	sprintf( buffer, "%10.5f %10.5f %10.5f", virial( 0, 0), virial( 0, 1), virial( 0, 2));
    	WRITE << buffer << std::endl;
    	sprintf( buffer, "%10.5f %10.5f %10.5f", virial( 1, 0), virial( 1, 1), virial( 1, 2));
    	WRITE << buffer << std::endl;
    	sprintf( buffer, "%10.5f %10.5f %10.5f", virial( 2, 0), virial( 2, 1), virial( 2, 2));
    	WRITE << buffer << std::endl;
#else   // CPPIO
		fprintf( WRITE, "%f\n%f %f %f\n%f %f %f\n%f %f %f\n", energy, virial( 0, 0), virial( 0, 1), virial( 0, 2), virial( 1 , 0), virial( 1, 1), virial( 1, 2), virial( 2, 0), virial( 2, 1), virial( 2, 2));
#endif  // CPPIO


    	LOG << "#### buried-non-redirected: " << sforce_count
			<< " max-dist(surf): " <<  sqrt( max_squared_dist)
			<< " buried-redirected: " << redirected
			<< " max-dist(mol-geom-center,redirected): " << sqrt( max_squared_redirected_dist) << std::endl;
#ifdef TIMERS
    	}
#endif

#endif  // MPI_PARALLEL


#ifdef TIMERS
    	WriteTimerStatus( LOG);
    	MOL_ITR->WriteStatus( LOG);
#endif


//    	++store::TransientInteractionGridPoint::s_pub_Step;

#ifdef ANALYZE_BURIAL
    	m_AngleChecker.Reset();
#endif

#ifdef ANALYZE_BURIAL
    	std::stringstream strstr;
    	strstr << s_Count++;
    	std::string str;
		strstr >> str;
    	std::ofstream out;
    	out.open( ("hist/surf_histogram" + str + ".txt").c_str());
    	out << surf_distances;
    	out.close();
    	out.clear();
    	out.open( ("hist/redirected_histogram" + str + ".txt").c_str());
    	if( !out)
    	{
    		std::cout << "=====> histogram could not be opened! make sure that hist/ dir exists!" << std::endl;
    	}
    	out << redirected_distances;
    	out.close();
    	out.clear();
    	surf_distances.Clear();
    	redirected_distances.Clear();
#endif

    } // end CalculateAndWriteForcesEnergyVirial( MOLS)



    void MoleculeForceGrid::CalculateAndWriteForcesEnergyVirial
    (
#ifdef CPPIO
    		std::ostream &WRITE,
#else
    		FILE *WRITE,
#endif
          boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
  		  math::Matrix3x3N &VIRIAL,
  		  float &ENERGY,
  		  float &MAX_SQUARED_DIST,
  		  float &MAX_SQUARED_REDIRECTED_DIST,
  		  int &SFORCE_COUNT,
  		  int &REDIRECTED_COUNT,
  		  std::ostream &LOG
#ifdef ANALYZE_BURIAL 
		  , math::Histogram< float> & REDIRECTED_DISTANCES,
		  math::Histogram< float> & SURF_DISTANCES
#endif
    )
    {
#ifdef TIMERS
    	m_CalcForceVirialEnergyTimer.Start();
#endif
    	DebugPrint( LOG, __FUNCTION__);

    	std::string
			mol_type = MOL_ITR->GetThisMolecule()->GetAtoms()( 0)->GetResidueType();

//	LOG << "moltype: " << mol_type << std::endl;

    	if( std::find( m_LipidNames.begin(), m_LipidNames.end(), mol_type) != m_LipidNames.end())
    	{

    		float
				squared_dist;
        	math::Vector3N
				zero,
    			force,
    			connect;

#ifdef TIMERS
        	m_GetGPTimer.Start();
#endif
        	math::Vector3N
				geometric_center = MOL_ITR->GetThisMolecule()->GeometricCenter();

        	util::FunctorEnum
				gridpoint_type = m_Grid->GetGridPointType( geometric_center, mol_type);
#ifdef TIMERS
        	m_GetGPTimer.Stop();
#endif



    		boost::shared_ptr< Atom>
				atom_ptr;


#ifdef DEBUG
    		LOG << "MOL-GMC: " << geometric_center( 0) << " " << geometric_center( 1) << " " << geometric_center( 2) << " type: " << util::EnumHandler< util::FunctorEnum>().String( gridpoint_type) << std::endl;
#endif

    		if( gridpoint_type == util::e_VoidGridPoint)
    		{
    			SimpleCalcAndWriteForcesEnergyVirial( WRITE, MOL_ITR, VIRIAL, ENERGY, SFORCE_COUNT, MAX_SQUARED_DIST, LOG);
    		}
    		else
    		{

    			do
    			{

#ifdef TIMERS
    				m_GetGPTimer.Start();
#endif

    				Memory();

    				atom_ptr = MOL_ITR->GetThisAtom();

    				// exclude hydrogens before interpolation takes place
    				if
    				(
						m_IgnoreHydrogens
						&& atom_ptr->GetMass() < 2.0
    				)
    				{
    					util::FunctorEnum
							local_type = m_Grid->GetGridPointType( atom_ptr->GetPosition(), mol_type);

    					if( local_type == util::e_ConstantGridPoint)
    					{
							++SFORCE_COUNT;
							MOL_ITR->CheckAndWriteForce( WRITE, zero, local_type);
							DebugPrint( LOG, "hydrogen excluded: " << atom_ptr->GetMass() << " gridpoint type: " << util::EnumHandler< util::FunctorEnum>().String( local_type));
#ifdef TIMERS
    							m_GetGPTimer.Stop();
#endif
							continue;
    					}
    				}
#ifdef FORCE_INDICES
    				// builds sumgridpoint if interpolation is switched on
    				boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
				    atom_grid_point = m_Grid->GetGridPoint( atom_ptr->GetPosition(), mol_type, WRITE);
#else
    				boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
						atom_grid_point = m_Grid->GetGridPoint( atom_ptr->GetPosition(), mol_type);
#endif 
#ifdef TIMERS
    				m_GetGPTimer.Stop();
#endif
				DebugPrint( LOG, "atom GP: " << atom_ptr->GetPosition()( 0) << " " << atom_ptr->GetPosition()( 1) << " " << atom_ptr->GetPosition()( 2) << "  " << util::EnumHandler< util::FunctorEnum>().String( atom_grid_point->GetClassID()));

#ifdef FORCE_COORDINATES
				for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
    				{
#ifdef CPPIO
    					WRITE << *itr << "  ";
#else
//    					LOG << *itr << "  ";
					fprintf( WRITE, "%f  ", *itr);
#endif
    				}
#endif
#ifdef NO_INTERACTION_FORCES
    				if(  atom_grid_point->GetClassID() == util::e_InteractionGridPoint)
    				{
    					MOL_ITR->CheckAndWriteForce( STREAM, zero);
    					continue;
    				}
#endif
#ifdef SPLIT_FORCES
    				if(   atom_grid_point->GetClassID() == util::e_InteractionGridPoint)
    				{
    					store::InteractionGridPoint * iptr = ( store::InteractionGridPoint *) atom_grid_point.get();
    					boost::shared_ptr< phys::PotentialForceContainer> po_cont = iptr->GetSetPotentialForcesContainer();
    					for( std::map< std::string, boost::shared_ptr< phys::PotentialForce> >::const_iterator po_itr = po_cont->begin(); po_itr != po_cont->end(); ++po_itr)
    					{
    						STREAM << po_itr->second->GetClassName() << ": ";
    						force = po_itr->second->operator()( atom_ptr);
    						for( std::vector< float>::const_iterator force_itr = force->begin(); force_itr != force->end(); ++force_itr)
    						{
    							STREAM << *force_itr << " ";
    						}
    					}
    					force = atom_grid_point->operator ()( atom_ptr);
    					STREAM << " total: ";
    					for( std::vector< float>::const_iterator force_itr = force->begin(); force_itr != force->end(); ++force_itr)
    					{
    						STREAM << *force_itr << " ";
    					}
    					STREAM << std::endl;
    				}
#endif
#ifdef TIMERS
    				m_ForceTimer.Start();
#endif


    				force = atom_grid_point->operator ()( *atom_ptr);


#ifdef TIMERS
    				m_ForceTimer.Stop();
    				m_RedirectTimer.Start();
#endif
					 // if atom is within implicit volume
					 if( atom_grid_point->GetClassID() == util::e_ConstantGridPoint)
					 {
//							 ++SFORCE_COUNT;
					     DebugPrint( LOG, "atom within implicit volume");
					     connect = geometric_center - ( atom_ptr)->GetPosition();
					     // if force points away from NSP of geometric center of molecule
					     if( force.AreAllElements( 0.0))
					     {
						 std::cout << "====> force is zero!!" << std::endl;
						 std::cout << "atom: " << *atom_ptr << std::endl;
						 std::cout << atom_grid_point << std::endl;
					     }

					     if( m_AngleChecker.IsAngleLargerThreshold( connect, force))
					     {
							 if( force.AreAllElements( 0.0))
							 {
								 std::cout << "=====> WHY???" << std::endl;
							 }

							 ++REDIRECTED_COUNT;

							 squared_dist = connect.SquaredLength();

							 if( squared_dist > MAX_SQUARED_REDIRECTED_DIST)
							 {
								 MAX_SQUARED_REDIRECTED_DIST = squared_dist;
							 }


#ifdef ANALYZE_BURIAL
							 REDIRECTED_DISTANCES.InsertValue( squared_dist);
#endif
#ifdef DEBUG
							 LOG << "=> angle checker strikes for ";
							 LOG << " res_id: " << atom_ptr->GetResidueID();
							 LOG << " res_type: " << atom_ptr->GetResidueType();
							 LOG << " atom_id: " << atom_ptr->GetAtomID();
							 LOG << " atom_type: " << atom_ptr->GetType();
							 LOG << " dist: " << squared_dist << std::endl;
#endif
#ifdef VMD_OUTPUT
			if( force.Length() > 0.05)
			{
							 // vmd visualization of force before redirection
//							 LOG << "graphics 0 color red" << std::endl;
							 LOG << "graphics 0 cone {";
							 for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
							 {
								 LOG << *itr << "  ";
							 }
							 LOG << "} {";
							 for( std::vector< float>::const_iterator force_itr = force.begin(), itr = atom_ptr->GetPosition().begin(); force_itr != force.end(); ++force_itr, ++itr)
							 {
								 LOG << ( *force_itr + *itr) << " ";
							 }
							 LOG << "} radius 0.1" << std::endl;
			}
#endif


							 // redirect force, pointing to geometric center of explicit molecule
							 force = connect.SetToLength( m_SforceLength);


#ifdef VMD_OUTPUT
							 if( force.Length() > 0.05)
							 {
							     // vmd visualization of redirected force
//							     LOG << "graphics 1 color green" << std::endl;
							     LOG << "graphics 1 cone {";
							     for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
							     {
								 LOG << *itr << " ";
							     }
							     LOG << "} {";
							     for( std::vector< float>::const_iterator force_itr = force.begin(), itr = atom_ptr->GetPosition().begin(); force_itr != force.end(); ++force_itr, ++itr)
							     {
								 LOG << ( *force_itr + *itr) << " ";
							     }
							     LOG << "} radius 0.1" << std::endl;
							     
							     
//							 // vmd visualization of nsp to which the force points
//							 LOG << "graphics 0 sphere {";
//							 for( std::vector< float>::const_iterator itr = nsp.begin(); itr != nsp.end(); ++itr)
//							 {
//								 LOG << *itr << " ";
//							 }
//							 LOG << "} radius 0.2" << std::endl;
							     
							     
							     // vmd visualization for the geometric center
							     math::Vector3N tmp = MOL_ITR->GetThisMolecule()->GeometricCenter();
//							     LOG << "graphics 2 color magenta" << std::endl;
							     LOG << "graphics 2 sphere {";
							     for( std::vector< float>::const_iterator itr =  tmp.begin(); itr != tmp.end(); ++itr)
							     {
								 LOG << *itr << " ";
							     }
							     LOG << "} radius 0.2" << std::endl;
							 }
#endif
					     }
					     else
					     {
							 ++SFORCE_COUNT;
							 squared_dist = ( m_SurfGrid.CalcPositionFromID( atom_grid_point->GetNSP()) - atom_ptr->GetPosition()).SquaredLength();
							 if( squared_dist > MAX_SQUARED_DIST)
							 {
								 MAX_SQUARED_DIST = squared_dist;
							 }
							 DebugPrint( LOG, "buried non-redirected, dist: " << squared_dist);
#ifdef VMD_OUTPUT
			if( force.Length() > 0.05)
			{
							 // vmd visualization of force before redirection
//							 LOG << "graphics 5 color red" << std::endl;
							 LOG << "graphics 5 cone {";
							 for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
							 {
								 LOG << *itr << "  ";
							 }
							 LOG << "} {";
							 for( std::vector< float>::const_iterator force_itr = force.begin(), itr = atom_ptr->GetPosition().begin(); force_itr != force.end(); ++force_itr, ++itr)
							 {
								 LOG << ( *force_itr + *itr) << " ";
							 }
							 LOG << "} radius 0.1" << std::endl;
			}
			else{ std::cout << "ugly" << std::endl;}
#endif

#ifdef ANALYZE_BURIAL
							 SURF_DISTANCES.InsertValue( squared_dist);
#endif

					     }
					 }
#ifdef VMD_OUTPUT
					 else if( atom_grid_point->GetClassID() == util::e_InteractionGridPoint	&& force.Length() > 0.1)
					 {
//					     LOG << "graphics 3 color green" << std::endl;
					     LOG << "graphics 3 cone {";
					     for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
					     {
						 LOG << *itr << " ";
					     }
					     LOG << "} {";
					     for( std::vector< float>::const_iterator force_itr = force.begin(), itr = atom_ptr->GetPosition().begin(); force_itr != force.end(); ++force_itr, ++itr)
					     {
						 LOG << ( *force_itr + *itr) << " ";
					     }
					     LOG << "} radius 0.1" << std::endl;
					 }
#endif					 
#ifdef TIMERS
					 m_RedirectTimer.Stop();
					 m_CheckTimer.Start();
#endif

//					 LOG << " force: " << atom_ptr->GetPosition()[0] << "  "<< atom_ptr->GetPosition()[1] << "  "<< atom_ptr->GetPosition()[2] << "  " << force[0] << " " << force[1] << " " << force[2] << "moltype: " << mol_type << " gp: " << util::EnumHandler< util::FunctorEnum>().String( gridpoint_type) << std::endl;

					 MOL_ITR->CheckAndWriteForce( WRITE, force, atom_grid_point->GetClassID());
					 
					 ENERGY += atom_grid_point->Energy( *atom_ptr);
					 
					 VIRIAL -= math::MatrixProduct( force, atom_ptr->GetPosition());
					 
#ifdef TIMERS
					 m_CheckTimer.Stop();
#endif

    			} while( MOL_ITR->IterateAtom( LOG));
    		}
    	}
    	else
    	{
			DebugPrint( LOG, "water");
			SimpleCalcAndWriteForcesEnergyVirial( WRITE, MOL_ITR, VIRIAL, ENERGY, SFORCE_COUNT, MAX_SQUARED_DIST, LOG);
    	}
#ifdef TIMERS
    	m_CalcForceVirialEnergyTimer.Stop();
#endif
    } // end CalculateAndWriteForcesEnergyVirial( MOL)
    
    


    void MoleculeForceGrid::SimpleCalcAndWriteForcesEnergyVirial
    (
#ifdef CPPIO
        		std::ostream &WRITE,
#else
        		FILE *WRITE,
#endif
          boost::shared_ptr< mol::MoleculeIterator> &MOL_ITR,
  		  math::Matrix3x3N &VIRIAL,
  		  float &ENERGY,
  		  int &SFORCE_COUNT,
  		  float &MAX_SQUARED_DIST,
  		  std::ostream &LOG
    )
    {
#ifdef TIMERS
    	m_SimpleTimer.Start();
    	m_SimpleGetGPTimer.Start();
#endif

    	DebugPrint( LOG, __FUNCTION__);
    	math::Vector3N
			force,
			zero;
		boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > >
			atom_grid_point;

		std::string
			mol_type = MOL_ITR->GetThisMolecule()->GetAtoms()( 0)->GetResidueType();

		util::FunctorEnum
			gridpoint_type;

		boost::shared_ptr< Atom> 
		    atom_ptr;

#ifdef TIMERS
		m_SimpleGetGPTimer.Stop();
#endif

		do
		{

#ifdef TIMERS
			m_SimpleGetGPTimer.Start();
#endif

			atom_ptr = MOL_ITR->GetThisAtom();

    				// exclude hydrogens before interpolation takes place
    				if
    				(
						m_IgnoreHydrogens
						&& atom_ptr->GetMass() < 2
    				)
    				{
    					util::FunctorEnum
							local_type = m_Grid->GetGridPointType( atom_ptr->GetPosition(), mol_type);

    					if( local_type == util::e_ConstantGridPoint)
    					{
							++SFORCE_COUNT;
							MOL_ITR->CheckAndWriteForce( WRITE, zero, local_type);
							DebugPrint( LOG, "hydrogen excluded in simple: " << atom_ptr->GetMass() << " gridpoint type: " << util::EnumHandler< util::FunctorEnum>().String( local_type));
#ifdef TIMERS
							m_SimpleGetGPTimer.Stop();
#endif
							continue;
					}
    				}


#ifdef FORCE_INDICES
			// contains interpolation
			atom_grid_point = m_Grid->GetGridPoint( atom_ptr->GetPosition(), mol_type, WRITE);
#else
			// contains interpolation
			atom_grid_point = m_Grid->GetGridPoint( atom_ptr->GetPosition(), mol_type);
#endif


			gridpoint_type = atom_grid_point->GetClassID();


#ifdef TIMERS
			m_SimpleGetGPTimer.Stop();
			m_RedirectTimer.Start();
#endif


#ifdef DEBUG
			LOG << "atom GP: " << atom_ptr->GetPosition()( 0) << " " << atom_ptr->GetPosition()( 1) << " " << atom_ptr->GetPosition()( 2) << " "<< util::EnumHandler< util::FunctorEnum>().String( atom_grid_point->GetClassID()) << std::endl;
#endif
#ifdef FORCE_COORDINATES
#ifdef CPPIO
			std::copy( atom_ptr->GetPosition().begin(), atom_ptr->GetPosition().end(), std::ostream_iterator< float>( WRITE, " "));
#else
		math::Vector3N fopo = atom_ptr->GetPosition();
			fprintf( WRITE, "%f  %f  %f  ", fopo[0], fopo[1], fopo[2]);
//			LOG << fopo[0] << "  " << fopo[1] << "  " << fopo[2] << "  ";
#endif
#endif
#ifdef NO_INTERACTION_FORCES
			if(  gridpoint_type == store::InteractionGridPoint().GetClassName())
			{
				MOL_ITR->CheckAndWriteForce( STREAM, zero);
				continue;
			}
#endif


			if( gridpoint_type == util::e_ConstantGridPoint)
			{
				++SFORCE_COUNT;

				float squared_dist = ( m_SurfGrid.CalcPositionFromID( atom_grid_point->GetNSP()) - atom_ptr->GetPosition()).SquaredLength();

				DebugPrint( LOG, "buried in simple, dist: " << sqrt( squared_dist));
				if( squared_dist > MAX_SQUARED_DIST)
				{
					MAX_SQUARED_DIST = squared_dist;
				}
			}

#ifdef TIMERS
			m_RedirectTimer.Stop();
			m_ForceTimer.Start();
#endif


			force = atom_grid_point->operator ()( *atom_ptr);


#ifdef TIMERS
			m_ForceTimer.Stop();
			m_CheckTimer.Start();
#endif
#ifdef VMD_OUTPUT
			if( force.Length() > 0.1)
			{
			    // vmd visualization of force before redirection
//			    LOG << "graphics 4 color red" << std::endl;
			    LOG << "graphics 4 cone {";
			    for( std::vector< float>::const_iterator itr = atom_ptr->GetPosition().begin(); itr != atom_ptr->GetPosition().end(); ++itr)
			    {
			    	LOG << *itr << "  ";
			    }
			    LOG << "} {";
			    for( std::vector< float>::const_iterator force_itr = force.begin(), itr = atom_ptr->GetPosition().begin(); force_itr != force.end(); ++force_itr, ++itr)
			    {
			    	LOG << ( *force_itr + *itr) << " ";
			    }
			    LOG << "} radius 0.1" << std::endl;
			}
#endif
//			LOG << "force: " << force[0] << " " << force[1] << " " << force[2] << " gp: " << util::EnumHandler< util::FunctorEnum>().String( gridpoint_type) << std::endl;
			
//			LOG << " force: " << atom_ptr->GetPosition()[0] << "  "<< atom_ptr->GetPosition()[1] << "  "<< atom_ptr->GetPosition()[2] << "  " << force[0] << " " << force[1] << " " << force[2] << "moltype: " << mol_type << " gp: " << util::EnumHandler< util::FunctorEnum>().String( gridpoint_type) << std::endl;

			MOL_ITR->CheckAndWriteForce( WRITE, force, gridpoint_type);

			ENERGY += atom_grid_point->Energy( *MOL_ITR->GetThisAtom());

			VIRIAL -= math::MatrixProduct( force, atom_ptr->GetPosition());
#ifdef TIMERS
			m_CheckTimer.Stop();
#endif

        } while( MOL_ITR->IterateAtom( LOG));
#ifdef TIMERS
    	m_SimpleTimer.Stop();
#endif

    }  // end SimpleCalcAndWriteForcesEnergyVirial



    /////////////////////////
    //      Read/Write     //
    /////////////////////////


    std::istream &MoleculeForceGrid::Read( std::istream &STREAM)
    {
    	DebugWrite( __PRETTY_FUNCTION__ );
    	Memory();
    	std::string str;
    	STREAM >> str;

    	if( str != mystr::GetClassName( __PRETTY_FUNCTION__))
    	{
    		std::cout << "===> not the correct header <" << str << "> for: " << __PRETTY_FUNCTION__ << std::endl;
    	}

    	m_Grid = boost::shared_ptr< AtomForceGrid>( new AtomForceGrid());

        m_Grid->Read( STREAM);
    	Memory();
    	StandardWrite( "read surf grid" );
    	STREAM >> m_SurfGrid;
    	Memory();
//    	if( m_Grid->GetPositionGrid()->GetClassName() == store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > >().GetClassName())
//    	{
//    		std::cout << "set molecule types in Grid1D" << std::endl;
//    		store::DefaultUniqueMap< std::string, int> map;
//    		int i = 0;
//    		for( std::vector< std::string>::const_iterator itr = m_SurfGrid.GetMoleculeTypes().begin(); itr != m_SurfGrid.GetMoleculeTypes().end(); ++itr, ++i)
//    		{
//    			map.InsertNewKeyAndValue( *itr, i);
//    		}
//		map.SetDefault( "all");
//    		store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > * grid
//				= ( store::Grid1D< boost::shared_ptr< math::Function< mol::Atom, math::Vector3N > > > *) m_Grid->GetSetPositionGrid().get();
//    		grid->SetMolTypes( map);
//    	}
    	if( !STREAM)
    	{
    		std::cout << "===> stream closed after reading surf grid" << std::endl;
    		exit( -1);
    	}
    	StandardWrite( "read angle checker" );
    	STREAM >> m_AngleChecker;
    	Memory();
    	if( !STREAM)
    	{
    		std::cout << "===> stream closed after reading angle checker" << std::endl;
    		exit( -1);
    	}
    	STREAM >> str;
	if( str != "sforce-scale:")
	{
	    std::cout << "(old school grid with the length in the VectorAngle)" << std::endl;
	    STREAM >> str;
	    if( str != "sforce-scale:")
	    {
		std::cout << "expected sforce-scale but got: " << str << std::endl;
		exit( -1);
	    }
	}
    	STREAM >> m_SforceLength;
    	StandardWrite( "read sforce scale: " << m_SforceLength);
    	Memory();
        return STREAM;
    }



    std::ostream &MoleculeForceGrid::Write( std::ostream &STREAM) const
    {
        STREAM << MoleculeForceGrid::GetClassName() << std::endl;
        m_Grid->Write( STREAM);
//        STREAM << m_VdwEpsilonAndRadiusMap;
//        STREAM << m_PartialChargeMap;
//        STREAM << m_MassMap;
//        STREAM << "force-field: " << m_ForceFieldType << std::endl;
        STREAM << m_SurfGrid;
        STREAM << m_AngleChecker;
        STREAM << "sforce-scale: " << m_SforceLength << std::endl;
        return STREAM;
    }

    std::ostream &MoleculeForceGrid::WriteAsPdb( std::ostream &STREAM) const
    {
        return m_Grid->WriteAsPdb( STREAM);
    }


    std::string MoleculeForceGrid::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace mol

