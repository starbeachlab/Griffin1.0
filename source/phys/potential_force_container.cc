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


#include "../../include/phys/potential_force_container.h"

namespace phys
{
    /////////////////////////
    //     OPERATORS       //
    /////////////////////////

    math::Vector3N PotentialForceContainer::operator()( const std::string &TYPE, const mol::Atom &ATOM) const
    {
        return store::ShPtrMap< std::string, PotentialForce>::operator()( TYPE)->operator()( ATOM);
    }

    math::Vector3N PotentialForceContainer::operator()( const mol::Atom &ATOM) const
    {
        DebugWrite( __PRETTY_FUNCTION__);
        math::Vector3N result;
        for
        (
                std::map< std::string, boost::shared_ptr< PotentialForce> >::const_iterator itr( store::ShPtrMap< std::string, PotentialForce>::begin());
                itr != store::ShPtrMap< std::string, PotentialForce>::end();
                ++itr
        )
        {
            DebugWrite( itr->first << " " << itr->second->GetClassName() << " " << itr->second->operator()( ATOM));

            result += itr->second->operator()( ATOM);
        }
        DebugWrite( "total_force: " << result);
        return result;
    }

    float PotentialForceContainer::Energy( const mol::Atom &ATOM) const
    {
        DebugWrite( __PRETTY_FUNCTION__);
        float energy = 0.0;
        for
        (
                std::map< std::string, boost::shared_ptr< PotentialForce> >::const_iterator itr( store::ShPtrMap< std::string, PotentialForce>::begin());
                itr != store::ShPtrMap< std::string, PotentialForce>::end();
                ++itr
        )
        {
            DebugWrite( itr->first << " " << itr->second->GetClassName() << " " << itr->second->Energy( ATOM));
            energy += itr->second->Energy( ATOM);
        }
        DebugWrite( "total_energy: " << energy);
        return energy;
    }

//    void PotentialForceContainer::operator += ( const std::pair< std::string, math::Vector3N >  &TYPE_FORCE)
//    {
//        DebugWrite( __PRETTY_FUNCTION__);
//        std::map< std::string, boost::shared_ptr< PotentialForce>>::iterator itr( store::ShPtrMap< std::string, PotentialForce>::find( TYPE_FORCE.first));
//        if( itr == store::ShPtrMap< std::string, PotentialForce>::end())
//        {
//            std::cout << TYPE_FORCE.first << " is not a defined key" << std::endl;
//            exit( -1);
//        }
//        else
//        {
//             itr->second += TYPE_FORCE.second;
//        }
//    }

    /////////////////////////
    //      FUNCTIONS      //
    /////////////////////////

    boost::shared_ptr< PotentialForce> &PotentialForceContainer::PotentialForceFromType( const std::string &TYPE)
    {
        return store::ShPtrMap< std::string, PotentialForce>::operator()( TYPE);
    }

    /////////////////////////
    //      Read/Write     //
    /////////////////////////

    std::istream &PotentialForceContainer::Read( std::istream &STREAM) //, const std::string &EXCLUDE)
    {
       	DebugWrite( __PRETTY_FUNCTION__);
    	std::string str;
    	int size;
    	STREAM >> str;
    	STREAM >> str;
    	DebugWrite(  str);
    	if( str == GetClassName())
    	{
    		DebugWrite(  "class name");
    		STREAM >> size;
    	}
    	else
    	{
    		size = mystr::ConvertStringToNumericalValue< int>( str);
    		DebugWrite(  "convert " << size);
    	}

    	DebugWrite(  "size: " << size);
    	for( int i = 0; i < size; ++i)
    	{
    		str.clear();
    		STREAM >> str;
    		if( str == "vdw-attractive")
    		{
#ifndef NO_VDW
    			DebugWrite(  "vdw-attractive");
    			PotentialAttractiveVanDerWaalsForce vdw;
    			STREAM >> vdw;
    			InsertNewKeyAndValue( str,  boost::shared_ptr< PotentialForce>( new PotentialAttractiveVanDerWaalsForce( vdw)));
#endif
    		}
#ifndef NO_VDW
    		else if( str == "vdw-repulsive")
    		{
    			DebugWrite(  "vdw-repulsive");
    			PotentialRepulsiveVanDerWaalsForce vdw;
    			STREAM >> vdw;
    			InsertNewKeyAndValue( str, boost::shared_ptr< PotentialForce>( new PotentialRepulsiveVanDerWaalsForce( vdw)));
    		}
#endif
#ifndef NO_COULOMB
    		else if( str == "coulomb")
    		{
    			DebugWrite(  "coulomb");
    			PotentialElectrostaticForce coulomb;
    			STREAM >> coulomb;
    			InsertNewKeyAndValue( str, boost::shared_ptr< PotentialForce>( new PotentialElectrostaticForce( coulomb)));
    		}
#endif
    		else
    		{
    			std::cout << "===> undefined interaction type (" << str << ") in: " << __PRETTY_FUNCTION__ << std::endl;
    		}

    	}
    	DebugWrite( __FUNCTION__ << " done");
        return STREAM;
    }

    std::ostream &PotentialForceContainer::Write( std::ostream &STREAM) const
    {
        STREAM << GetClassName() << "  ";
        store::ShPtrMap< std::string, PotentialForce>::Write( STREAM);
        return STREAM;
    }

    std::string PotentialForceContainer::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace phys

