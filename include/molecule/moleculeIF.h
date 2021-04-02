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


#ifndef MOLECULE_IF_H
#define MOLECULE_IF_H

#include "../storage/shared_pointer_vector.t.h"

namespace mol
{

    template< typename t_ATOM>
    class MoleculeIF
    {
    public:
        virtual ~MoleculeIF(){}

        virtual MoleculeIF *Clone() const = 0;

        virtual std::ostream &Write( std::ostream &STREAM) const = 0;

        virtual std::istream &Read( std::istream &STREAM) = 0;

        virtual const store::ShPtrVec< t_ATOM> &GetAtoms() const = 0;

        virtual const std::string &GetType() const = 0;

        virtual const size_t &GetID() const = 0;

        virtual std::string GetClassName() const = 0;


    }; // end class MoleculeIF

} // end namespace mol

#endif
