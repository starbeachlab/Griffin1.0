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


#ifndef FORCE_CONTAINER_H_
#define FORCE_CONTAINER_H_

#include "../string/io_string_functions.h"
#include "../storage/shared_pointer_map.t.h"
#include "force.h"

//#include "force_container.hh"
//#include "potential_force_container.hh"
//#include "force_factory.hh"


namespace phys
{

//    /////////////////////////////////
//    ////   FORWARD DECLARATIONS
//    /////////////////////////////////
//
//    namespace factory
//    {
//        boost::shared_ptr< PotentialForceContainer>
//        BuildPotentialForceContainerFromForceContainer( const phys::ForceContainer &FORCES);
//    }

    class ForceContainer
    : public store::ShPtrMap< std::string, Force>
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ForceContainer()
        : store::ShPtrMap< std::string, Force>()
        {}

        //! construct from data
        ForceContainer( const std::vector< std::string> &DEFINED)
        : store::ShPtrMap< std::string, Force>( DEFINED)
        {}

        //! copy constructor
        ForceContainer( const ForceContainer &ORIGINAL)
        : store::ShPtrMap< std::string, Force>( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~ForceContainer(){}

        //! virtual copy constructor
        virtual ForceContainer *Clone() const{ return new ForceContainer( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

        virtual float GetMaximumCutoff() const;

        //    protected:
//        virtual boost::shared_ptr< PotentialForceContainer> EmptyPotentialForceContainer() const
//        {
//            boost::shared_ptr< PotentialForceContainer> container( new PotentialForceContainer( store::ShPtrMap< std::string, Force>::GetKeys()));
//            for
//            (
//                    std::map< std::string, boost::shared_ptr< Force> >::const_iterator itr( store::ShPtrMap< std::string, Force>::begin());
//                    itr != std::map< std::string, boost::shared_ptr< Force> >::end();
//                    ++itr
//            )
//            {
//                container->PotentialForce( itr->first) = itr->second->GetAssociatedPotentialForceObject();
//            }
//            return container;
//        }


    public:
//        virtual boost::shared_ptr< PotentialForceContainer> CalculatePotentialForceContainer( const math::Vector3N &POSITION, const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL) const;

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;


    }; // end class ForceContainer
} // end namespace phys




#endif /* FORCE_CONTAINER_H_ */
