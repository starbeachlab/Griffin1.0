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


#include "../../include/phys/force_container.h"


namespace phys
{


//    boost::shared_ptr< PotentialForceContainer> ForceContainer::CalculatePotentialForceContainer( const math::Vector3N &POSITION, const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL) const
//    {
//        boost::shared_ptr< PotentialForceContainer> container( factory::BuildPotentialForceContainerFromForceContainer( *this));
//        for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator mol_itr( MOL->GetAtoms().begin()); mol_itr != MOL->GetAtoms().end(); ++mol_itr)
//            for
//            (
//                    std::map< std::string, boost::shared_ptr< Force> >::const_iterator itr( std::map< std::string, boost::shared_ptr< Force> >::begin());
//                    itr != std::map< std::string, boost::shared_ptr< Force> >::end();
//                    ++itr
//            )
//            {
//                // for force type (itr->first) calculate force between neutral atom at POSITION and entire molecule
//                container->PotentialForceFromType( itr->first)->AddToPotentialForce( *mol_itr, itr->second->operator()( boost::shared_ptr< mol::Atom>( new mol::Atom( "neutral", POSITION)), *mol_itr));
//            }
//        return container;
//    }

	float ForceContainer::GetMaximumCutoff() const
	{
		float max_cutoff = 0.0;
		for( std::map< std::string, boost::shared_ptr< Force> >::const_iterator itr = this->begin(); itr != this->end(); ++itr)
		{
			if( itr->second->GetCutoff() > max_cutoff)
			{
				max_cutoff = itr->second->GetCutoff();
			}
		}
		return max_cutoff;
	}


    /////////////////////////
    //      Read/Write     //
    /////////////////////////

    std::istream &ForceContainer::Read( std::istream &STREAM)
    {
        return store::ShPtrMap< std::string, Force>::Read( STREAM);
    }

    std::ostream &ForceContainer::Write( std::ostream &STREAM) const
    {
        STREAM << GetClassName() << std::endl;
        return store::ShPtrMap< std::string, Force>::Write( STREAM);
    }

    std::string ForceContainer::GetClassName() const
    {
        return mystr::GetClassName( __PRETTY_FUNCTION__);
    }

} // end namespace phys


