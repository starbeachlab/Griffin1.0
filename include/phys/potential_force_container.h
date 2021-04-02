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


#ifndef Potential_Force_Container_H_
#define Potential_Force_Container_H_

#include "../string/io_string_functions.h"
#include "../storage/shared_pointer_map.t.h"
#include "potential_force.h"
#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"
#include "potential_attractive_vdw_force.h"
#include "potential_repulsive_vdw_force.h"
#include "potential_electrostatic_force.h"


namespace phys
{

//    /////////////////////////////////
//    ////   FORWARD DECLARATIONS
//    /////////////////////////////////
//    class Force;
//    class ForceContainer;
//    class PotentialForceContainer;
//
//    namespace factory
//    {
//        boost::shared_ptr< store::Map< std::string, PotentialForce> >
//        BuildPotentialForceContainerFromForceContainer( const boost::shared_ptr< ForceContainer> &FORCES);
//    }


    ////////////////////////////////////////////////
    //! class PotentialForceContainer
    //!
    //!
    //!  @date Jan 20, 2009
    //!  @author: Rene Staritzbichler
    //!  @example:
    ////////////////////////////////////////////////


    class PotentialForceContainer
    : public store::ShPtrMap< std::string, PotentialForce>
    {

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        PotentialForceContainer()
        : store::ShPtrMap< std::string, PotentialForce>()
        {}

        //! construct from data
        PotentialForceContainer( const std::vector< std::string> &DEFINED)
        : store::ShPtrMap< std::string, PotentialForce>( DEFINED)
        {}

//        //! construct from std::vector< PotentialForce >
//        PotentialForceContainer( const boost::shared_ptr< ForceContainer> &FORCES)
//        : store::ShPtrMap< std::string, PotentialForce>( *( factory::BuildPotentialForceContainerFromForceContainer( FORCES)))
//        {}

        //! copy constructor
        PotentialForceContainer( const PotentialForceContainer &ORIGINAL)
        : store::ShPtrMap< std::string, PotentialForce>( ORIGINAL)
        {}

        //! virtual destructor
        virtual ~PotentialForceContainer(){}

        //! virtual copy constructor
        virtual PotentialForceContainer *Clone() const{ return new PotentialForceContainer( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        virtual math::Vector3N operator()( const std::string &TYPE, const mol::Atom &ATOM) const;

        virtual math::Vector3N operator()( const mol::Atom &ATOM) const;

        virtual const boost::shared_ptr< PotentialForce> &operator()( const std::string &TYPE) const
        { return store::ShPtrMap< std::string, PotentialForce>::operator()( TYPE);}

        virtual boost::shared_ptr< PotentialForce> &operator()( const std::string &TYPE)
        { return store::ShPtrMap< std::string, PotentialForce>::operator()( TYPE);}

//        virtual void operator += ( const std::pair< std::string, math::Vector3N >  &TYPE_FORCE);

//        PotentialForceContainer operator += ( const PotentialForce &FORCE)
//        {
//            m_Container.PushBack( PotentialForce);
//        }
        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////
        virtual float Energy( const mol::Atom &ATOM) const;

        boost::shared_ptr< PotentialForce> &PotentialForceFromType( const std::string &TYPE);

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM);

        virtual std::ostream &Write( std::ostream &STREAM) const;

        virtual std::string GetClassName() const;

    }; // end class PotentialForceContainer
} // end namespace phys




#endif /* PotentialForce_Container_H_ */
