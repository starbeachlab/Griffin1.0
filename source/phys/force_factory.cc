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


#include "../../include/phys/force_factory.h"

namespace phys
{


    namespace factory
    {
        std::istream &ReadForces( std::istream &STREAM, boost::shared_ptr< phys::ForceContainer> &FORCES, const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE)
        {
            std::string str;
            while( STREAM)
            {
                std::getline( STREAM, str);
                if( str.size() > 5)
                {
                    std::vector< std::string> vec( mystr::SplitString( str));
                    std::string type( vec[0]);

                    if( type == "coulomb")
                    {
                        std::cout << " mount coulomb ..." << std::endl;
                        assert( vec.size() == 2);
                        float cutoff( mystr::ConvertStringToNumericalValue< float>( vec[ 1]));
                        FORCES->InsertNewKeyAndValue( "coulomb", boost::shared_ptr< ElectrostaticForce>( new ElectrostaticForce( IMPLICIT_MOLECULE, cutoff)));
                        std::cout << " ... coulomb force with cutoff: " << cutoff << " mounted" << std::endl;
                    }
                    else if( type == "vdw-attractive")
                    {
                        std::cout << " mount vdw ..." << std::endl;
                        assert( vec.size() == 2);
                        float cutoff( mystr::ConvertStringToNumericalValue< float>( vec[ 1]));
                        FORCES->InsertNewKeyAndValue( "vdw-attractive", boost::shared_ptr< AttractiveVanDerWaalsForce>( new AttractiveVanDerWaalsForce( IMPLICIT_MOLECULE, cutoff)));
                        std::cout << " ... attractive vdw force with cutoff: " << cutoff << " is mounted" << std::endl;
                    }
                    else if( type == "vdw-repulsive")
                    {
                        std::cout << " mount vdw ..." << std::endl;
                        assert( vec.size() == 2);
                        float cutoff( mystr::ConvertStringToNumericalValue< float>( vec[ 1]));
                        FORCES->InsertNewKeyAndValue( "vdw-repulsive", boost::shared_ptr<RepulsiveVanDerWaalsForce>( new RepulsiveVanDerWaalsForce( IMPLICIT_MOLECULE, cutoff)));
                        std::cout << " ... repulsive vdw force with cutoff: " << cutoff << " is mounted" << std::endl;
                    }
                    else if( type == "density")
                    {
                        std::cout << " mount density ..." << std::endl;
                        assert( vec.size() == 2);
                        std::string density_file( vec[ 1]);
                //        float threshold( mystr::ConvertStringToNumericalValue< float>( vec[3]));
    //                    FORCES->InsertNewKeyAndValue( "density", boost::shared_ptr< DensityForceHandler>( new DensityForceHandler( weight, density_file, threshold)));
                        std::cout << " ... density force: " << density_file << " is mounted" << std::endl;
                    }
                }
            }
            return STREAM;
        }

        boost::shared_ptr< PotentialForceContainer>
        BuildPotentialForceContainerFromForceContainer( const phys::ForceContainer &FORCES)
        {
            boost::shared_ptr< PotentialForceContainer> potentials( new PotentialForceContainer());
            //store::ShPtrMap< std::string, Force> pots;
            for( std::map< std::string, boost::shared_ptr< Force> >::const_iterator force_itr = FORCES.begin(); force_itr != FORCES.end(); ++force_itr)
            {
                DebugWrite( "insert: " << force_itr->first);
                potentials->InsertNewKeyAndValue( force_itr->first, force_itr->second->GetAssociatedPotentialForceObject());
            }
            return potentials;
        }


//        boost::shared_ptr< PotentialForceContainer>
//        BuildPotentialForces
//        (
//                const boost::shared_ptr< phys::ForceContainer> &FORCES,
//                const math::Vector3N &GRID_POS
//        )
//        {
//            boost::shared_ptr< PotentialForceContainer> potential_forces( BuildPotentialForceContainerFromForceContainer( *FORCES));
//            for( std::map< std::string, boost::shared_ptr< Force> >::const_iterator force_itr = FORCES->begin(); force_itr != FORCES->end(); ++force_itr)
//            {
//                DebugWrite( "build " << force_itr->first);
//                potential_forces->operator()( force_itr->first) = force_itr->second->PotentialForce( GRID_POS);
//            }
//            return potential_forces;
//        }


        PotentialForceContainer
        BuildPotentialForces
        (
                const boost::shared_ptr< phys::ForceContainer> &FORCES,
                const math::Vector3N &GRID_POS,
                const float &MIN_FORCE_MAGNITUDE,
                const float &DISTANCE
        )
        {
            PotentialForceContainer potential_force_container;

            boost::shared_ptr< phys::PotentialForce> potential_force;

            for( std::map< std::string, boost::shared_ptr< Force> >::const_iterator force_itr = FORCES->begin(); force_itr != FORCES->end(); ++force_itr)
            {
                DebugWrite( "insert: " << force_itr->first);

                if( force_itr->second->GetCutoff() < DISTANCE)
                {
                	DebugWrite( "out of limit: dist: " << DISTANCE << " cutoff" << force_itr->second->GetCutoff());
                	continue;
                }

                potential_force = force_itr->second->PotentialForce( GRID_POS);

                DebugWrite( "check: length: " << potential_force->GetSetVector().Length() << " type: " << potential_force->GetClassName());

                if( potential_force->GetSetVector().Length() > MIN_FORCE_MAGNITUDE)
                {
                    DebugWrite( "inserted: " << force_itr->first);
                	potential_force_container.InsertNewKeyAndValue( force_itr->first, potential_force);
                }
            }
            return potential_force_container;
        }

    } // end namespace factory
} // end namespace phys



