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
//!	Read/write functions for atoms.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef ATOM_FILE_HANDLER_H_
#define ATOM_FILE_HANDLER_H_


#include "../string/string_functions.h"
#include "../storage/map.t.h"
#include "../storage/map_functions.t.h"

extern float g_PDBFactor;

namespace mol
{
    template< typename t_ATOM>
    class AtomFileHandler
    {
    public:

    	// JUST A REMINDER ABOUT PDB ATOM SECTION FORMAT:

        //            COLUMNS        DATA  TYPE    FIELD        DEFINITION
        //            -------------------------------------------------------------------------------------
        //             1 -  6        Record name   "ATOM  "
        //             7 - 11        Integer       serial       Atom  serial number.
        //            13 - 16        Atom          name         Atom name.
        //            17             Character     altLoc       Alternate location indicator.
        //            18 - 20        Residue name  resName      Residue name.
        //            22             Character     chainID      Chain identifier.
        //            23 - 26        Integer       resSeq       Residue sequence number.
        //            27             AChar         iCode        Code for insertion of residues.
        //            31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        //            39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        //            47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        //            55 - 60        Real(6.2)     occupancy    Occupancy.
        //            61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        //            77 - 78        LString(2)    element      Element symbol, right-justified.
        //            79 - 80        LString(2)    charge       Charge  on the atom.

        boost::shared_ptr< t_ATOM> ReadFromPdbLine( std::string &STR) const
        {
            DebugWrite( __FUNCTION__ << ": " << STR);
            std::string
                atom_name( mystr::TrimString( STR.substr( 12, 4))),
                residue( mystr::TrimString( STR.substr( 17, 4)));
            if( atom_name.length() == 0 || residue.length() == 0)
            {
            	std::cerr << "\ncould not read atom name or residue name!" << std::endl;
            	std::cerr << "residue name: " << residue << std::endl;
            	std::cerr << "atom name: " << atom_name << std::endl;
            	std::cerr << "bailing out" << std::endl;
            	exit(-1);
            }
            size_t
                residue_id = mystr::ConvertStringToNumericalValue< size_t>( STR.substr( 22, 5)),
                atom_id = mystr::ConvertStringToNumericalValue< size_t>( STR.substr( 6, 5));

            DebugWrite( " element: " << atom_name);
            math::Vector3N
                pos
                (
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 30, 8)),
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 38, 8)),
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 46, 8))
                );
            DebugWrite(  __FUNCTION__ << " create atom: ");
            return boost::shared_ptr< t_ATOM>( new t_ATOM( "", atom_name, residue, pos, 0.0, residue_id, atom_id));

        }

        boost::shared_ptr< t_ATOM> ReadInGriffinFormat( std::istream &STREAM) const
        {
        	DebugWrite( __FUNCTION__ );
        	std::string
		    str,
		    atom_name,
			  residue;

        	size_t
			  residue_id,
		    atom_id,
		    mol_id;

			float
			  mass,
			  partial_charge,
			  epsilon,
			  vdw_radius,
			  temperature;

			math::Vector3N
				pos;

			STREAM >> atom_id
				>> atom_name
				>> residue_id
			       >> residue;
//			STREAM >> mol_id;
			STREAM >> str;
			if( mystr::IsNumerical( str))
			{
			    mol_id = mystr::ConvertStringToNumericalValue< size_t>( str);
			}
			else
			{
			    mol_id = -999;
			}
			STREAM >> pos[0] >> pos[1] >> pos[2]
			    >> mass
			    >> partial_charge
			    >> vdw_radius
			    >> epsilon
			    >> temperature;

            DebugWrite(  __FUNCTION__ << " create atom: ");

            return boost::shared_ptr< t_ATOM>( new t_ATOM( "", atom_name, residue, pos, residue_id - 1, atom_id - 1, mass, partial_charge, epsilon, vdw_radius, mol_id - 1));
        }



        boost::shared_ptr< t_ATOM> ReadFromPdbLine
        (
                std::string &STR,
                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR,
                std::ostream &LOG = std::cout
        ) const
        {
            DebugWrite( __FUNCTION__ << ": " << STR);
            std::string
                name( mystr::TrimString( STR.substr( 12, 4))),
                residue( mystr::TrimString( STR.substr( 17, 4)));
            size_t
                residue_id = mystr::ConvertStringToNumericalValue< size_t>( STR.substr( 22, 5)),
                atom_id = mystr::ConvertStringToNumericalValue< size_t>( STR.substr( 6, 5));
            math::Vector3N
                pos
                (
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 30, 8)),
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 38, 8)),
                        mystr::ConvertStringToNumericalValue< float>( STR.substr( 46, 8))
                );

            if( residue == "HSD")
            {
                residue = "HIS";
            }

            std::string
                type( store::ReadFromMap( TRANSLATOR, residue, name));

            if( type != "")
            {
                DebugWrite( "name: " << name << " type: " << type << " residue: " << residue);

                return boost::shared_ptr< t_ATOM>( new t_ATOM( type, name, residue, pos, 0.0, residue_id, atom_id));
            }
//            if( type == "")
//            {
//                type = name.substr( 0, 1);
//                if( type[0] > 47 && type[0] < 58) // number check
//                {
//                    LOG << __FUNCTION__ << " skip number in string " << name << std::endl;
//                    type = name.substr( 1, 1);
//                }
//                LOG << "===> " << __FUNCTION__ << ": simple type recognition fix applied! atom-name " << name;
//                name = type;
//                LOG << " translated to type = " << type << " and name: " << name << " for residue-type: " << residue << std::endl;
//            }
            LOG << "===> " << __FUNCTION__ << ": no matching type found! " /*<< STR*/ << std::endl;
            return boost::shared_ptr< t_ATOM>();
        }


        std::ostream &WriteToPdb
         (
                 std::ostream &STREAM,
                 const boost::shared_ptr< t_ATOM> &ATOM,
                 const size_t ATOM_ID
         )
         {
        	std::string
				line = BasicPDBLine( ATOM, ATOM_ID),
				tmp = mystr::NumericalValueToString( ATOM->GetResidueID() % 9999);

            line.replace( 22 + ( 4 - tmp.size()), tmp.size(), tmp);

        	STREAM << line << std::endl;
        	return STREAM;
         }

        std::ostream &WriteToPdb
        (
        		std::ostream &STREAM,
        		const boost::shared_ptr< t_ATOM> &ATOM,
        		const size_t RESIDUE_ID,
        		const size_t ATOM_ID
        )
        {
        	std::string
				line = BasicPDBLine( ATOM, ATOM_ID),
				tmp = mystr::NumericalValueToString( RESIDUE_ID % 9999);

            line.replace( 22 + ( 4 - tmp.size()), tmp.size(), tmp);

        	STREAM << line << std::endl;
        	return STREAM;
        }


        std::ostream &WriteInGriffinFormat
        (
        		std::ostream &STREAM,
        		const boost::shared_ptr< t_ATOM> &ATOM
        )
        {
        	return WriteInGriffinFormat( STREAM, ATOM, ATOM->GetResidueID() + 1, ATOM->GetAtomID() + 1);
        }


        std::ostream &WriteInGriffinFormat
        (
        		std::ostream &STREAM,
        		const boost::shared_ptr< t_ATOM> &ATOM,
        		const size_t ATOM_ID
        )
        {
        	return WriteInGriffinFormat( STREAM, ATOM, ATOM->GetResidueID() + 1, ATOM_ID);
        }


        std::ostream &WriteInGriffinFormat
        (
        		std::ostream &STREAM,
        		const boost::shared_ptr< t_ATOM> &ATOM,
        		const size_t RESIDUE_ID,
        		const size_t ATOM_ID
        )
        {
//        	char chain = 'x';
        	math::Vector3N pos = ATOM->GetPosition();
        	float temperature = 0.0;

			STREAM.width( 10);
			STREAM << ATOM_ID << " ";

			STREAM.width( 5);
			STREAM << ATOM->GetAtomName() << " ";

			STREAM.width( 8);
			STREAM << RESIDUE_ID << " ";

			STREAM.width( 5);
			STREAM << ATOM->GetResidueType() << " ";

			STREAM.width( 8);
			STREAM << ATOM->GetMolID() << " ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << pos[0] << " ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << pos[1] << " ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << pos[2] << " ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << ATOM->GetMass() << "  ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << ATOM->GetPartialCharge() << "  ";

			STREAM.width( 10);
			STREAM.precision( 4);
			STREAM << ATOM->GetVanDerWaalsRadius() << "  ";

			STREAM.width( 10);
//			STREAM.precision( 5);
			STREAM << ATOM->GetEpsilon() << "  ";

			STREAM.width( 10);
			STREAM << temperature << " ";

			STREAM << std::endl;

        	return STREAM;
         }




	std::string
        BasicPDBLine
        (
                const boost::shared_ptr< t_ATOM> &ATOM,
                const size_t ATOM_ID
        )
        {
            size_t
                first;
            std::string
                tmp,
                line( 67, ' ');

            line.replace( 0, 4, "ATOM");

            if( ATOM_ID > 99999)
            {
            	tmp = "*****";
            }
            else
            {
            	tmp = mystr::NumericalValueToString( ATOM_ID);
            }
            first = 6 + ( 5 - tmp.size());
            line.replace( first, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();

            tmp = ATOM->GetAtomName();
//            first = 12 + ( 4 - tmp.size());
            first = 12;
            if( tmp.size() < 4)
            {
                ++first;
            }
            line.replace( first, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();

            tmp = ATOM->GetResidueType();
//            first = 17 + ( 4 - tmp.size());
            line.replace( 17, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();


////            if( ATOM->GetResidueID() > 9999)
//            if( RESIDUE_ID > 9999)
//            {
//                tmp = "****";
//            }
//            else
//            {
//                tmp = mystr::NumericalValueToString( ATOM->GetResidueID() % 9999);
//            }
//            first = 22 + ( 4 - tmp.size());
//            line.replace( first, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();

            math::Vector3N pos( ATOM->GetPosition());
            tmp = mystr::NumericalValueToString( g_PDBFactor * pos( 0), 8, 3);
            if( tmp.size() > 8)
            {
                tmp = tmp.substr( 0, 8);
            }
            first = 30 + ( 8 - tmp.size());
            line.replace( first, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();

            tmp = mystr::NumericalValueToString( g_PDBFactor * pos( 1), 8, 3);
            if( tmp.size() > 8)
            {
                tmp = tmp.substr( 0, 8);
            }
            first = 38 + ( 8 - tmp.size());
            line.replace( first, tmp.size(), tmp);

//            std::cout << tmp << " ";
//            std::cout.flush();

            tmp = mystr::NumericalValueToString( g_PDBFactor * pos( 2), 8, 3);
            if( tmp.size() > 8)
            {
                tmp = tmp.substr( 0, 8);
            }
            first = 46 + ( 8 - tmp.size());
            line.replace( first, tmp.size(), tmp);

////            line.replace( 56, 4, "1.00");
////
////            line.replace( 62, 4, "0.00");
////

//            std::cout << tmp << " ";
//            std::cout.flush();

//            tmp = ATOM->GetElementName();
//            first = 76 + ( 2 - tmp.size());
//            line.replace( first, tmp.size(), tmp);
//
//            std::cout << tmp << " ";
//            std::cout.flush();

//            std::cout << std::endl;
            return line;
        }


    }; // end class AtomFileHandler
} // end namespace mol




#endif /* ATOM_FILE_HANDLER_H_ */
