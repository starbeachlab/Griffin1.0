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


class LipidReader
{
public:
  
  //! Reads a list of atom.
  LipidContainer
  ReadAtomList( std::istream &READ)
  {
    LipidContainer array;
    while( READ)
      {
    // temporary variables
    size_t 
      atom_nr( std::numeric_limits< size_t>::max()), 
      lipid_nr( std::numeric_limits< size_t>::max());
    std::string 
      atom_name, 
      lipid_name;
    
    // read from file into temporaries
    READ >> atom_nr >> atom_name >> lipid_name >> lipid_nr;

    // check for empty lines and exclude 'H'-Atoms
    //    if( atom_name.substr(0,1) != "H" && lipid_name.size() == 4 && lipid_nr < std::numeric_limits< size_t>::max())
    if( lipid_name.size() == 4 && lipid_nr < std::numeric_limits< size_t>::max())
      {
        //        std::cout << "pushed" << std::endl;
        array.push_back( LipidPtr( new LipidAtom( atom_nr, atom_name, lipid_name, lipid_nr)));
      }
      }
    return array;
  }


  void
  ReadNAMDFile( std::istream &STREAM,  LipidContainer &ARRAY)
  {
    // rely on that the number of lines in the NAMD file agree with AtomFile
    // don't worry about that fu..up loop syntax...
    size_t i(1);
    for( std::deque< boost::shared_ptr< LipidAtom> >::iterator itr = ARRAY.begin(); itr != ARRAY.end(); ++itr, ++i)
      {
    float tmp, x, y, z;
    STREAM >> tmp;
    assert( tmp == i);
    STREAM >> tmp >> x >> y >> z;
    (*itr)->SetCharge( tmp);
    (*itr)->SetCoordinates( x, y, z);
      }
  }

  LipidContainer FilterForLipidName( const LipidContainer &LIPIDS, std::string &FILTER)
  {
    LipidContainer tmp;
    for( std::deque< boost::shared_ptr< LipidAtom> >::const_iterator itr = LIPIDS.begin(); itr != LIPIDS.end(); ++itr)
      if( ( *itr)->GetLipidName() == FILTER)
    { tmp.push_back( *itr);}
    return tmp;
  }
  

};  // end class LipidAtomReader
