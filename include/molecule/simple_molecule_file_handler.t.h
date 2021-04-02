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
//!	Read/write functions for simple molecules.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef SIMPLE_MOLECULE_FILE_HANDLER_H_
#define SIMPLE_MOLECULE_FILE_HANDLER_H_

#include <iterator>

#include "atom_file_handler.t.h"
#include "../external/std_functions.h"
#include "../readwrite/stream_functions.h"
#include "../storage/multiplets.h"
#include "../external/std_functions.h"
#include "../phys/constants.h"
#include "atom.h"
#include "simple_molecule.t.h"
#include "element_data_map.h"

namespace mol
{
    namespace file
    {


        //! a control for reading pdb files, it repairs e.g. problems due to the limited residue numbers
        inline bool LineChecker( const size_t &NR)
        {
            static size_t s_WaterAtomCounter( 0);
            return ++s_WaterAtomCounter % NR == 0;
        }

        //! this routine reads atoms into simple molecules, creating a new molecule for each residue. additionally it tries to repair pdb-free-style-poetry
        inline
        boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
        ReadMoleculesFromPdb( std::istream &STREAM, std::ostream &LOG = std::cout)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            std::string
                str,
                residue_type,
                previous_type;
            int
                residue_id,
                previous_id( std::numeric_limits< int>::max());
            size_t
                molecule_count( 0);

            boost::shared_ptr< SimpleMolecule< Atom> >
                molecule( new SimpleMolecule< Atom>());
            boost::shared_ptr< Atom>
                atom;

            boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
                all_molecules( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());

            // collecting and fitting limits
            while( STREAM)
            {
                std::getline( STREAM, str);
                if( str.length() > 6 && str.substr( 0, 4) == "ATOM")
                {
                    atom = AtomFileHandler< Atom>().ReadFromPdbLine( str);

                    residue_id = atom->GetResidueID();
                    residue_type = atom->GetResidueType();

                    DebugPrint( LOG, "atom: " << residue_id << " " << residue_type << " prev: " << previous_id << " " << previous_type);

                    // if new molecule, add // Linechecker is for the water molecules that have no reliable resid in pdb file
                    if( residue_id != previous_id || residue_type != previous_type || ( atom->GetResidueType() == "TIP3" && molecule->GetAtoms().size() == 3))
//                        if( ( atom->GetResidueType() == "TIP3" && LineChecker( 3)) || residue_type != previous_type || ( atom->GetResidueType() != "TIP3" && residue_id != previous_id))
                    {
                        previous_id = residue_id;
                        previous_type = residue_type;

                        if( molecule_count > 0 && molecule->GetAtoms().size() > 0)
                        {
                            DebugPrint( LOG, "add molecule to all");
                            all_molecules->push_back( molecule);
                            molecule = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>( residue_type, residue_id));
                        }
                        ++molecule_count;
                    }
                    // add atom to molecule
                    molecule->AddAtom( atom);

                    if( atom->GetResidueType() == "TIP3" && molecule->GetAtoms().size() == 3)
                    {
                        all_molecules->push_back( molecule);
                        molecule = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>( residue_type, residue_id));
                    }

                }
            }
            if( molecule->GetAtoms().size() > 0)
            {
                DebugPrint( LOG, "add last");
                all_molecules->push_back( molecule);
            }
            return all_molecules;
        }



        //! this routine reads atoms into simple molecules, creating a new molecule for each residue.
        inline
        boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
        ReadMoleculesInGriffinFormat( std::istream &STREAM, std::ostream &LOG = std::cout)
        {
            DebugWrite( __PRETTY_FUNCTION__);

            if( !STREAM)
            {
            	std::cout << "no input stream opened in " << __FUNCTION__ << std::endl;
            }

            std::string
                str,
                residue_type,
                previous_type;
            int
                residue_id,
                previous_id( std::numeric_limits< int>::max());

            boost::shared_ptr< SimpleMolecule< Atom> >
                molecule( new SimpleMolecule< Atom>());
            boost::shared_ptr< Atom>
                atom;

            boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
                all_molecules( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());

			// read first atom to set type of first molecule (as residue type...) // TODO: make clear protein vs pope/tip3
			atom = AtomFileHandler< Atom>().ReadInGriffinFormat( STREAM);
			previous_type = atom->GetResidueType();
			molecule->SetType( previous_type);
			previous_id = atom->GetResidueID();
			molecule->AddAtom( atom);

            // collecting and fitting limits
            while( STREAM)
            {
                    atom = AtomFileHandler< Atom>().ReadInGriffinFormat( STREAM);

                    residue_id = atom->GetResidueID();
                    residue_type = atom->GetResidueType();

                    DebugPrint( LOG, "atom: " << residue_id << " " << residue_type << " prev: " << previous_id << " " << previous_type);

                    // if new molecule, add
                    if( residue_id != previous_id || residue_type != previous_type || ( atom->GetResidueType() == "TIP3" && molecule->GetAtoms().size() == 3))
                    {
                        previous_id = residue_id;
                        previous_type = residue_type;

                        DebugPrint( LOG, "add molecule to all");
                        all_molecules->push_back( molecule);
                        molecule = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>( residue_type, residue_id));
                    }
                    // add atom to molecule
                    molecule->AddAtom( atom);
            }
            if( molecule->GetAtoms().size() > 1)
            {
                DebugPrint( LOG, "add last");
                all_molecules->push_back( molecule);
            }
            return all_molecules;
        }





        //! this routine reads atoms into simple molecules, creating a new molecule for each residue. additionally it tries to repair pdb-free-style-poetry
        inline
        boost::shared_ptr< mol::SimpleMolecule< mol::Atom> >
        ReadMoleculeInGriffinFormat( std::istream &STREAM, std::ostream &LOG = std::cout)
        {
            DebugWrite( __PRETTY_FUNCTION__);

            boost::shared_ptr< SimpleMolecule< Atom> >
                molecule( new SimpleMolecule< Atom>());
            boost::shared_ptr< Atom>
                atom;

            while( STREAM)
            {
            	atom = AtomFileHandler< Atom>().ReadInGriffinFormat( STREAM);

            	// add atom to molecule
            	molecule->AddAtom( atom);

            }
            if( molecule->GetAtoms().size() > 1)
            {
                DebugPrint( LOG, "remove last");
                molecule->Atoms().pop_back();
            }
            return molecule;
        }


        //! read a simple molecule from pdb file (with atom name/type translator and log)
        template< typename t_ATOM>
        boost::shared_ptr< SimpleMolecule< t_ATOM> >
        ReadFromPdb
        (
                std::istream &STREAM,
                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR,
                std::ostream &LOG = std::cout
        )
        {
            DebugWrite( __PRETTY_FUNCTION__);
            boost::shared_ptr< SimpleMolecule< t_ATOM> > molecule( new SimpleMolecule< t_ATOM>());
            std::string str;
 //           std::cout << __FUNCTION__ << ": NOTE: Griffin only reads PBD ATOM and HETATM lines, others have no effect" << std::endl;
            while( STREAM)
            {
                str.clear();
                std::getline( STREAM, str);
                if( str.size() > 0 && ( str.substr( 0, 4) == "ATOM" || str.substr( 0, 6) == "HETATM"))
                {
                    boost::shared_ptr< t_ATOM> atom = AtomFileHandler< t_ATOM>().ReadFromPdbLine( str, TRANSLATOR);
                    if( atom)
                    {
                        DebugWrite( "name: <" << atom->GetAtomName() <<  "> type: <" << atom->GetType() << ">");
                        molecule->AddAtom( atom);
                    }
                    else
                    {
                        LOG << "===> Atom Line was NOT added!!! in " << __FUNCTION__ << ":\n" << str << std::endl;
                    }
                }
            }
            return molecule;
        }


        //! read a simple molecule from pdb file
        template< typename t_ATOM>
        boost::shared_ptr< SimpleMolecule< t_ATOM> >
        ReadFromPdb( std::istream &STREAM)
        {
            DebugWrite( __PRETTY_FUNCTION__);
            boost::shared_ptr< SimpleMolecule< t_ATOM> > molecule( new SimpleMolecule< t_ATOM>());
            std::string str;
            while( STREAM)
            {
                str.clear();
                std::getline( STREAM, str);
                if( str.size() > 0 && ( str.substr( 0, 4) == "ATOM" || str.substr( 0, 6) == "HETATM"))
                {
                    boost::shared_ptr< t_ATOM> atom( AtomFileHandler< t_ATOM>().ReadFromPdbLine( str));
                    DebugWrite( "name: <" << atom->GetAtomName() <<  "> type: <" << atom->GetType() << ">");
                    molecule->AddAtom( atom);
                }
            }
            molecule->SetType( molecule->GetAtoms()( 0)->GetResidueType());
            molecule->SetID( molecule->GetAtoms()( 0)->GetResidueID());
            return molecule;
        }

        //! update only coordinates of molecules from pdb file
        inline
        void
        ReadNewCoordinatesFromPdb( std::istream &STREAM, store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > &MOLS)
        {
            DebugWrite( __PRETTY_FUNCTION__);

            std::string
                line;

            std::vector< boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > >::iterator
                mol_itr( MOLS.begin());

            std::vector< boost::shared_ptr< mol::Atom> >::iterator
                atom_itr = ( *mol_itr)->Atoms().begin();


            // collecting and fitting limits
            while( STREAM) // TODO: iteration through molecule, not defined by stream!!!
            {
                line.clear();
                std::getline( STREAM, line);
                if( line.length() > 54 && (line.substr( 0, 4) == "ATOM" || line.substr( 0, 6) == "HETATM"))
                {

                    ( *atom_itr)->SetPosition
                    (
                            math::Vector3N
                            (
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 30, 8)),
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 38, 8)),
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 46, 8))
                            )
                    );

                    ++atom_itr;
                    if( atom_itr == ( *mol_itr)->Atoms().end() && ++mol_itr != MOLS.end())
                    {
                        atom_itr = ( *mol_itr)->Atoms().begin();
                    }

                }
            }
        }


        //! update only coordinates of single molecule from pdb file
        inline
        void
        ReadNewCoordinatesFromPdb( std::istream &STREAM, boost::shared_ptr< SimpleMolecule< Atom> > &MOL)
        {
        	DebugWrite( __FUNCTION__);
            std::string line;
            std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin();
            while( STREAM)
            {
                line.clear();
                std::getline( STREAM, line);
                if( line.size() > 0 && ( line.substr( 0, 4) == "ATOM" || line.substr( 0, 6) == "HETATM"))
                {
                    DebugWrite( ( *itr)->GetAtomName()  << " =>  " << line );
                    ( *itr)->SetPosition
                    (
                            math::Vector3N
                            (
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 30, 8)),
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 38, 8)),
                                    mystr::ConvertStringToNumericalValue< float>( line.substr( 46, 8))
                            )
                    );
                    ++itr;
                }
            }
        	DebugWrite( __FUNCTION__ << " done/!-");
        }

        //! update coordinates from gmx dump file
        inline
        void
        ReadNewCoordinatesFromGmxDump( std::istream &STREAM, boost::shared_ptr< SimpleMolecule< Atom> > &MOL, std::ostream &LOG = std::cout)
        {
            // TODO: security checks?

            LOG << __FUNCTION__ << std::endl;

            std::string
                line;
            bool
                found( false);
            size_t
                first,
                second,
                third,
                last;

            while( !found)
            {
                line.clear();
                std::getline( STREAM, line);
                if( line.substr( 0, 3) == "x (")
                {
                    found = true;
                    std::string
                        str( line.substr( 3, line.find( ")") - 5));
                    LOG << "nr atoms in file: "<< str << std::endl;
//                    std::cout << "nr atoms in file: "<< str << std::endl;
                    size_t
                        nr_atoms( mystr::ConvertStringToNumericalValue< size_t>( str));

                    if( nr_atoms != MOL->GetAtoms().size())
                    {
                        LOG << "===> incorrect number of atoms: in file: " << nr_atoms << " in molecule: " << MOL->GetAtoms().size() << std::endl;
                    }
                }
            }
            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                line.clear();
                std::getline( STREAM, line);
                DebugWrite( ( *itr)->GetAtomName()  << "=>  " << line.substr( 13, 12) << " " << line.substr( 27, 12) << " " << line.substr( 41, 12));

                first = line.find( "{") + 1;
                second = line.find( ",", first);
                third = line.find( ",", second + 2);
                last = line.find( "}", third);
                first = line.find_first_not_of( " ", first);


                ( *itr)->SetPosition
                (
                        math::Vector3N
                        (
                                mystr::ConvertStringToNumericalValue< float>( line.substr( first, second-first)),
                                mystr::ConvertStringToNumericalValue< float>( line.substr( second + 1, third - second - 1)),
                                mystr::ConvertStringToNumericalValue< float>( line.substr( third + 1, last - third - 1))
                        )
                );
            }
///            LOG << __FUNCTION__ << " done" << std::endl;
        }

        //! update coordinates from xyz files
        inline
        void
        ReadCoordinatesFromXyz
        (
#ifdef CPPIO
        		std::istream &STREAM,
#else
        		FILE *FILENAME,
#endif
        		boost::shared_ptr< SimpleMolecule< Atom> > &MOL,
        		std::ostream &LOG = std::cout
        )
        {
        	// TODO: security checks?
        	DebugPrint( LOG, __FUNCTION__);

        	float
                x, y, z;

            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
#ifdef CPPIO
                STREAM >> x >> y >> z;
#else
                fscanf( FILENAME, "%f %f %f", &x, &y, &z);
#endif
                ( *itr)->SetPosition(  x, y, z);
            }
        }



        inline
        boost::shared_ptr< SimpleMolecule< Atom> > &
        ReadMoleculeFromNamdFile( std::istream &STREAM, boost::shared_ptr< SimpleMolecule< Atom> > &MOL)
        {
            std::string tmp, atom, residue;
            while( STREAM)
            {
                STREAM >> tmp;
                if( tmp.size() == 0)
                { continue;}

                STREAM >> atom >> residue >> tmp;
//                MOL->AddAtom( boost::shared_ptr< Atom>( new Atom( AtomNameToAtomType()( atom), atom, residue)));
            }
            return MOL;
        }


        inline
        boost::shared_ptr< SimpleMolecule< Atom> > &
        ReadNewPositionsFromNamdFile
        (
#ifdef CPPIO
        		std::istream &STREAM,
#else
        		FILE *FILENAME,
#endif
        		boost::shared_ptr< SimpleMolecule< Atom> > &MOL
        )
        {
        	char tmp[15];
            float x, y, z;
            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {

#ifdef CPPIO
                STREAM >> tmp >> tmp >> x >> y >> z;  // todo: introduce seek
#else
                fscanf( FILENAME, "%s %s %f %f %f", tmp, tmp, &x, &y, &z);   // THIS IS OF LIMITED PRECISION !!!
#endif
                ( *itr)->SetPosition( x, y, z);
            }
            return MOL;
        }

        inline
        boost::shared_ptr< store::Map< std::string, float> > &
        ReadPartialChargeMapFromCharmmTopologyFile( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, float> > &MAP)
        {
            std::string line, residue, atom, block;
            float charge;
            while( STREAM)
            {
                line.clear();
                std::getline( STREAM, line);
                if( line.size() > 5)
                {
                    std::vector< std::string> vec( mystr::SplitString( line));
                    if( vec[0] == "RESI" || vec[0] == "PRES")
                    {
                        block = vec[0];
                        residue = vec[1];
                        if( residue == "HSD") // repair
                        {
                            residue = "HIS";
                        }
                    }
                    else if( vec[0] == "ATOM" && (   block == "RESI"
                            || ( block == "PRES" && residue == "CTER")
                            || ( block == "PRES" && residue == "NTER")))
                    {
                        atom = vec[1];
                        charge = mystr::ConvertStringToNumericalValue< float>( vec[3]);
                        MAP->InsertNewKeyAndValue( residue + ":" + atom, charge);
                    }
//                    if( vec[0] == "PRES")
//                    { break;}
                }
            }
            return MAP;
        }



        inline
        boost::shared_ptr< SimpleMolecule< Atom> > &
        ReadPartialChargesIntoMolecule
        (
                const boost::shared_ptr< store::Map< std::string, float> > &MAP,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOL
        )
        {

            std::string residue, atom;
            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                residue = ( *itr)->GetResidueType();
                atom = ( *itr)->GetAtomName();
                if( residue == "ILE" && atom == "CD1") // repair
                {
                    atom = "CD";
                }
                else if( residue == "LYS" && atom == "OXT") // repair
                {
                    atom = "O";
                }
//                ( *itr)->SetPartialCharge( MAP->operator()( ( *itr)->GetType()));
                ( *itr)->SetPartialCharge( store::ReadFromMap( MAP, residue, atom)); // MAP->operator()( residue + ":" + atom));
                DebugWrite( __FUNCTION__ << " " << residue << " " << atom << "  " << store::ReadFromMap( MAP, residue, atom));
            }
            return MOL;
        }


        inline
        boost::shared_ptr< store::Map< std::string, float> > &
        ReadMassFromCharmm( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, float> > &MAP)
        {
            std::string line, atom;
            float mass;
            while( STREAM)
            {
                std::getline( STREAM, line);
                if( line.size() > 0)
                {
                    std::vector< std::string> vec( mystr::SplitString( line));
                    if( vec[0] == "MASS")
                    {
                        atom = vec[2];
                        mass = mystr::ConvertStringToNumericalValue< float>( vec[3]);
                        MAP->InsertNewKeyAndValue( atom, mass);
                    }
                    if( vec[0] == "DECL")
                    { break;}
                }
            }
            return MAP;
        }



        inline
        boost::shared_ptr< SimpleMolecule< Atom> > &
        ReadMassIntoMolecule
        (
                const boost::shared_ptr< store::Map< std::string, float> > &MAP,
//                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOL
        )
        {
            std::cout  << __FUNCTION__ << std::endl;
            std::string atom, residue;
            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                residue = ( *itr)->GetResidueType();
                atom = ( *itr)->GetAtomName();
//                if( residue == "ILE" && atom == "CD1") // repair
//                {
//                    atom = "CD";
//                }
//                else if( residue == "LYS" && atom == "OXT") // repair
//                {
//                    atom = "O";
//                }
                DebugWrite(  "atom: " << atom << " residue: " << residue);
                DebugWrite(  "type: " << ( *itr)->GetType());
                DebugWrite(  "mass: " << MAP->operator()( ( *itr)->GetType()));
                ( *itr)->SetMass( MAP->operator()( ( *itr)->GetType()));
            }
            std::cout  << __FUNCTION__ <<  " done" <<  std::endl;
            return MOL;
        }

        inline
        boost::shared_ptr< SimpleMolecule< Atom> > &
        ReadMassIntoMolecule
        (
                const boost::shared_ptr< store::Map< std::string, float> > &MAP,
                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOL
        )
        {
//            std::cout << __FUNCTION__ << std::endl;
            std::string atom, residue;
            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                residue = ( *itr)->GetResidueType();
                atom = ( *itr)->GetAtomName();
//                if( residue == "ILE" && atom == "CD1") // repair
//                {
//                    atom = "CD";
//                }
//                else if( residue == "LYS" && atom == "OXT") // repair
//                {
//                    atom = "O";
//                }
//                std::cout << residue << " :: " << atom << " :: ";
//                std::cout.flush();
//                std::cout << TRANSLATOR->operator()( residue + ":" + atom) << " :: ";
//                std::cout.flush();
//                std::cout << MAP->operator()( TRANSLATOR->operator()( residue + ":" + atom)) << std::endl;

                DebugWrite(  "atom: " << atom);
                DebugWrite(  "translated: " << TRANSLATOR->operator()( residue + ":" + atom));
                DebugWrite(  "mass: " << MAP->operator()( TRANSLATOR->operator()( residue + ":" + atom)));
                ( *itr)->SetMass( MAP->operator()( store::ReadFromMap( TRANSLATOR, residue, atom)));
            }
            return MOL;
        }

        inline
        boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &
        ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &MAP)
        {
            std::string line, residue, atom;
            bool not_found( true);
            while( STREAM && not_found)
            {
                std::getline( STREAM, line);
                if( line.find( "epsilon") != std::string::npos)
                {
                    not_found = false;
                }
            }
            while( STREAM)
            {
                JumpOverComments( STREAM, '!');
                std::getline( STREAM, line);
                if( line.size() > 0)
                {
                    std::vector< std::string> vec( mystr::SplitString( line));
                    if( vec.size() > 2 && vec[0].length() <= 4 && mystr::IsCapitolLetter( vec[0][0]) && mystr::IsNumerical( vec[2]))
                    {
                        MAP->InsertNewKeyAndValue
                        (
                                vec[0],
                                std::make_pair
                                (
                                        mystr::ConvertStringToNumericalValue< float>( vec[2]),
                                        mystr::ConvertStringToNumericalValue< float>( vec[3])
                                )
                        );
                    }
                }
            }
            return MAP;
        }


        inline
        void
        ReadDefaultRadii
        (
                boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &MAP
        )
        {
            std::cout << __FUNCTION__ << " subsequent warnings about exotic element types like D from DUM, M from MG, Z from ZN are no problem if the implicit molecule does not contain them! (don't use this flag if it does!)" << std::endl;
            std::string id;
            float radius;

            for( std::map< std::string, std::pair< float, float> >::iterator itr = MAP->begin(); itr != MAP->end(); ++itr)
            {
                id = itr->first.substr( 0, 1);
                if( id[0] > 47 && id[0] < 58)
                {
                    id = itr->first.substr( 1, 1);
                }

                boost::shared_ptr< ElementData> element = mol::ElementDataMap()( id);
                if( element)
                {
                    radius = element->GetVanDerWaalsRadius();

                    itr->second.second = radius;
                }
                else
                {
                    std::cout << "===> " << __FUNCTION__ << " first char of string not found in mol::ElementDataMap - string: " << itr->first << std::endl;
                }
            }
        }


        inline
        void
        ReadVdwEpsilonAndRadiusIntoMolecule
        (
                const boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &MAP,
//                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOL,
                std::ostream &LOG = std::cout
        )
        {
//            std::cout  << __FUNCTION__ << std::endl;
            DebugPrint( LOG, __PRETTY_FUNCTION__);
            std::string
                type;

            for( std::vector< boost::shared_ptr< Atom> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                type = ( *itr)->GetType();

                if( type == "OXT")
                {
                    type = "OC";
                }

#ifdef DEBUG
                    LOG << __FUNCTION__ << std::endl;
                    LOG << "type " << type << std::endl;
                    LOG << "atom: " << *itr << std::endl;
                    LOG << "epsilon: " << MAP->operator()( type).first << std::endl;
                    LOG << "vdw-radius: " << MAP->operator()( type).second << std::endl;
#endif

                ( *itr)->SetEpsilon( MAP->operator()( type).first);
                ( *itr)->SetVanDerWaalsRadius( MAP->operator()( type).second);
            }
        }

        inline
        void
        ReadVdwEpsilonAndRadiusIntoMolecules
        (
                const boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &MAP,
                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS
        )
        {
            std::cout  << __FUNCTION__ << std::endl;
            for( std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator itr = MOLS->begin(); itr != MOLS->end(); ++itr)
            {
                ReadVdwEpsilonAndRadiusIntoMolecule( MAP/*, TRANSLATOR*/, *itr);
            }
        }

        inline
        void
        ReadAtomTypeTranslator( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR)
        {
            std::string line, residue, block;
            std::vector< std::string> vec;
            while( STREAM)
            {
                line.clear();
                std::getline( STREAM, line);
                if( line.size() > 5)
                {
                    vec = mystr::SplitString( line);
                    if( vec[0] == "RESI" || vec[0] == "PRES")
                    {
                        block = vec[0];
                        residue = vec[1];
                        if( residue == "HSD") // repair
                        {
                            residue = "HIS";
                        }
                    }
                    else if( vec[0] == "ATOM" && (   block == "RESI"
                            || ( block == "PRES" && residue == "CTER")
                            || ( block == "PRES" && residue == "NTER")))
                    {
                        DebugWrite( "res: " << residue << " name: " << vec[1] << " type: " << vec[ 2]);
                        TRANSLATOR->InsertNewKeyAndValue( residue + ":" + vec[1], vec[2]);
//                        if( TRANSLATOR->find( vec[1]) == TRANSLATOR->end())
//                        {
//                            DebugWrite( "=> inserted just residue name, but why?");
//                            TRANSLATOR->insert( std::make_pair( vec[1], vec[2]));
//                        }
                    }
                }
            }
            TRANSLATOR->insert( std::make_pair(  "ILE:CD1", "CT3"));
        }

        template< typename t_ATOM>
        void
        ReadPartialChargesFromPsfFile
        (
                std::istream &STREAM,
                boost::shared_ptr< SimpleMolecule< t_ATOM> > &MOL,
                std::ostream &LOG = std::cout
        )
        {
            DebugPrint( LOG,  __FUNCTION__ );
            std::string line;
            std::vector< std::string> vec;

            DebugPrint( LOG, "#atoms: " << MOL->GetAtoms().size());

            for( typename std::vector< boost::shared_ptr< t_ATOM> >::iterator itr = MOL->Atoms().begin(); itr != MOL->Atoms().end(); ++itr)
            {
                std::getline( STREAM, line);
                vec = mystr::SplitString( line);

#ifdef DEBUG
                LOG << __FUNCTION__ << " check: " << vec[4] << " ?= " << ( *itr)->GetAtomName() << std::endl;
                LOG << line << std::endl;
                std::cout << __FUNCTION__ << " check: " << vec[4] << " ?= " << ( *itr)->GetAtomName() << std::endl;
                assert( vec[4] == ( *itr)->GetAtomName());
#endif

                ( *itr)->SetPartialCharge( mystr::ConvertStringToNumericalValue< float>( vec[6]));
            }

        }


        template< typename t_ATOM>
        void
        ReadPartialChargesFromPsfFile
        (
                std::istream &STREAM,
                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > > &MOLS,
                std::ostream &LOG = std::cout
        )
        {
            LOG << __FUNCTION__ << std::endl;
            std::string line;
            bool is_header( true);
            while( STREAM && is_header)
            {
            	line.clear();
                std::getline( STREAM, line);
                DebugPrint( LOG,  "header: " << line);
                if( line.find( "NATOM") != std::string::npos)
                {
                    is_header = false;
                }
            }
            DebugPrint( LOG, "#mols: " << MOLS->size());

            for( typename std::vector< boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > >::iterator itr = MOLS->begin(); itr != MOLS->end(); ++itr)
            {
                ReadPartialChargesFromPsfFile( STREAM, *itr, LOG);
            }
        }

        template< typename t_ATOM>
        void
        ReadMoleculesFromPSF
        (
                std::istream &STREAM,
//                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > > &MOLS,
                store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > &MOLS,
                std::ostream &LOG = std::cout
        )
        {
            LOG << __FUNCTION__ << std::endl;

//            MOLS = boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > >( new store::ShPtrVec< mol::SimpleMolecule< t_ATOM> >());
            MOLS.clear();

            std::string line;
            bool is_header = true;
            while( STREAM && is_header)
            {
            	line.clear();
                std::getline( STREAM, line);
                DebugPrint( LOG,  "header: " << line);
                if( line.find( "NATOM") != std::string::npos)
                {
                    is_header = false;
                }
            }
            DebugPrint( LOG, "#mols: " << MOLS.size());

        	std::vector< std::string>
				vec = mystr::SplitString( line);

			size_t
				atom_id,
				residue_id( 0),
				previous_id = std::numeric_limits< size_t>::max(),
				size = mystr::ConvertStringToNumericalValue< int>( vec[0]);

			if( size == 0 || size > 1e6)
			{
				LOG << "WARNING: check size of molecule: " << size << std::endl;
			}

			float
				partial_charge,
				mass;

			std::string
			        type;

			boost::shared_ptr< mol::SimpleMolecule< t_ATOM> >
				molecule = boost::shared_ptr< mol::SimpleMolecule< t_ATOM> >( new mol::SimpleMolecule< t_ATOM>());

			boost::shared_ptr< t_ATOM>
				atom = boost::shared_ptr< t_ATOM>( new t_ATOM());

			math::Vector3N pos;

			LOG << size << " atoms" << std::endl;

            //while ( STREAM && line.find( "NBOND") != std::string::npos)
			for( size_t i = 0; i < size; ++i)
            {
            	line.clear();
            	std::getline( STREAM, line);
            	vec = mystr::SplitString( line);

            	if( vec.size() < 8)
            	{
            		LOG << __FUNCTION__ << "  WARNING: expects 8 entries per line, got <" << line << ">" << std::endl;
            	}

            	atom_id = mystr::ConvertStringToNumericalValue< int>( vec[0]);

            	residue_id = mystr::ConvertStringToNumericalValue< int>( vec[2]);

		//            	DebugPrint( LOG, "resid: " << residue_id << " atom_id: " << atom_id);

            	if( residue_id != previous_id)
            	{
                    DebugPrint( LOG,  "new mol: " << residue_id);
            		if( i != 0)
            		{
                        molecule->SetType( type);
                        molecule->SetID( previous_id);
                		MOLS.push_back( molecule);
            			molecule = boost::shared_ptr< mol::SimpleMolecule< t_ATOM> >( new mol::SimpleMolecule< t_ATOM>());
            		}
            		previous_id = residue_id;
            	}

            	type = vec[3];

            	partial_charge = mystr::ConvertStringToNumericalValue< float>( vec[6]);

            	mass = mystr::ConvertStringToNumericalValue< float>( vec[7]);

            	atom = boost::shared_ptr< t_ATOM>( new t_ATOM( vec[5], vec[4], type, pos, partial_charge, residue_id, atom_id));
            	atom->SetMass( mass);

            	molecule->AddAtom( atom);

            }
    		MOLS.push_back( molecule);
    		DebugPrint( LOG, __FUNCTION__ << " done");
        }


        template< typename t_ATOM>
        void
        ReadMoleculeFromPSF
        (
                std::istream &STREAM,
                boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > &MOL,
                std::ostream &LOG = std::cout
        )
        {
            LOG << __FUNCTION__ << std::endl;

            MOL = boost::shared_ptr< mol::SimpleMolecule< t_ATOM> >( new mol::SimpleMolecule< t_ATOM>());

            std::string line;
            bool is_header( true);
            while( STREAM && is_header)
            {
            	line.clear();
                std::getline( STREAM, line);
                DebugPrint( LOG,  "header: " << line);
                if( line.find( "NATOM") != std::string::npos)
                {
                    is_header = false;
                }
            }

        	std::vector< std::string>
				vec = mystr::SplitString( line);

			size_t
				atom_id,
				residue_id( 0),
//				previous_id = std::numeric_limits< size_t>::max(),
				size = mystr::ConvertStringToNumericalValue< int>( vec[0]);

			float
				partial_charge,
				mass;

			std::string
			        type;

			boost::shared_ptr< t_ATOM>
				atom = boost::shared_ptr< t_ATOM>( new t_ATOM());

			math::Vector3N pos;

			LOG << size << " atoms" << std::endl;

            //while ( STREAM && line.find( "NBOND") != std::string::npos)
			for( size_t i = 0; i < size; ++i)
            {
            	line.clear();
            	std::getline( STREAM, line);
            	vec = mystr::SplitString( line);

            	atom_id = mystr::ConvertStringToNumericalValue< int>( vec[0]);

            	residue_id = mystr::ConvertStringToNumericalValue< int>( vec[2]);

            	DebugPrint( LOG, "resid: " << residue_id << " atom_id: " << atom_id << " ==> " << line);

            	type = vec[3];

            	partial_charge = mystr::ConvertStringToNumericalValue< float>( vec[6]);

            	mass = mystr::ConvertStringToNumericalValue< float>( vec[7]);

            	atom = boost::shared_ptr< t_ATOM>( new t_ATOM( vec[5], vec[4], type, pos, partial_charge, residue_id, atom_id));

            	DebugPrint( LOG, "atom: " << atom);

            	atom->SetMass( mass);

            	MOL->AddAtom( atom);

            }
    		DebugPrint( LOG, __FUNCTION__ << " done");
        }


        template< typename t_ATOM>
        std::ostream &WriteToPdb
        (
                std::ostream &STREAM,
                const boost::shared_ptr< SimpleMolecule< t_ATOM> > &MOL,
                const size_t MOLECULE_ID,
                size_t &ATOM_COUNT
        )
        {
			DebugWrite( __FUNCTION__);
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = MOL->GetAtoms().begin(); itr != MOL->GetAtoms().end(); ++itr, ++ATOM_COUNT)
            {
                AtomFileHandler< t_ATOM>().WriteToPdb( STREAM, *itr, MOLECULE_ID, ATOM_COUNT);
            }
            return STREAM;
        }


        template< typename t_ATOM>
        std::ostream &WriteToPdb
        (
                std::ostream &STREAM,
                const boost::shared_ptr< SimpleMolecule< t_ATOM> > &MOL
        )
        {
			DebugWrite( __FUNCTION__);
            size_t count( 0);
            for( typename std::vector< boost::shared_ptr< t_ATOM> >::const_iterator itr = MOL->GetAtoms().begin(); itr != MOL->GetAtoms().end(); ++itr, ++count)
            {
                AtomFileHandler< t_ATOM>().WriteToPdb( STREAM, *itr, count);
            }
            return STREAM;
        }


        template< typename t_ATOM>
        std::ostream &WriteToPdb
        (
                std::ostream &STREAM,
                const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > > &MOLS
        )
        {
            size_t mol_count( 1), atom_count = 1;
            for( typename std::vector< boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > >::iterator itr = MOLS->begin(); itr != MOLS->end(); ++itr, ++mol_count)
            {
                WriteToPdb( STREAM, *itr, mol_count, atom_count);
            }
            return STREAM;
        }

	inline
        std::ostream &WriteInGriffinFormat
        (
                std::ostream &STREAM,
                const boost::shared_ptr< SimpleMolecule< mol::Atom> > &MOL,
                const size_t MOLECULE_ID,
                size_t &ATOM_COUNT
        )
        {
			DebugWrite( __FUNCTION__);
            for( std::vector< boost::shared_ptr< Atom> >::const_iterator itr = MOL->GetAtoms().begin(); itr != MOL->GetAtoms().end(); ++itr, ++ATOM_COUNT)
            {
                AtomFileHandler< Atom>().WriteInGriffinFormat( STREAM, *itr, MOLECULE_ID, ATOM_COUNT);
            }
            return STREAM;
        }


	inline
        std::ostream &WriteInGriffinFormat
        (
                std::ostream &STREAM,
                const boost::shared_ptr< SimpleMolecule< Atom> > &MOL
        )
        {
            size_t count( 1);
            for( std::vector< boost::shared_ptr< Atom> >::const_iterator itr = MOL->GetAtoms().begin(); itr != MOL->GetAtoms().end(); ++itr, ++count)
            {
                AtomFileHandler< Atom>().WriteInGriffinFormat( STREAM, *itr, count);
            }
            return STREAM;
        }

	inline
        std::ostream &WriteInGriffinFormat
        (
                std::ostream &STREAM,
                const boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< Atom> > > &MOLS
        )
        {
			DebugWrite( __FUNCTION__);
			size_t mol_count( 1), atom_count = 1;
            for( std::vector< boost::shared_ptr< mol::SimpleMolecule< Atom> > >::iterator itr = MOLS->begin(); itr != MOLS->end(); ++itr, ++mol_count)
            {
                WriteInGriffinFormat( STREAM, *itr, mol_count, atom_count);
            }
            return STREAM;
        }


        template< typename t_ATOM>
        void SortByTypeAndID
        (
                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > > &MOLS
        )
        {
            boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > >
				mols =  boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< t_ATOM> > >( new store::ShPtrVec< mol::SimpleMolecule< t_ATOM> >());

        	std::multimap< std::pair< std::string, int>, boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > > sorted_molecules;

        	for( typename std::vector< boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > >::const_iterator itr = MOLS->begin(); itr != MOLS->end(); ++itr)
        	{
        		DebugWrite( "sorted_molecules.insert( " << ( *itr)->GetAtoms()( 0)->GetResidueType() << ", " << ( *itr)->GetAtoms()( 0)->GetResidueID() << ")");
        		sorted_molecules.insert( std::make_pair( std::make_pair( ( *itr)->GetAtoms()( 0)->GetResidueType(), ( *itr)->GetAtoms()( 0)->GetResidueID()), *itr));
        	}

        	for( typename std::multimap< std::pair< std::string, int>, boost::shared_ptr< mol::SimpleMolecule< t_ATOM> > >::const_iterator itr = sorted_molecules.begin(); itr != sorted_molecules.end(); ++itr)
        	{
        		mols->push_back( itr->second);
        	}
        	DebugWrite( "MOLS: " << MOLS->size() << " sorted: " << sorted_molecules.size() << " mols: " << mols->size());
        	MOLS = mols;
        }




        inline
        void ReadAtomTypeMapFromGromacsCharm( std::istream &STREAM, boost::shared_ptr< store::Map< std::string, std::string> > &ATOM_TYPE_MAP)
        {
            DebugWrite( __FUNCTION__);
            std::string
                line,
                residue;
            std::vector< std::string>
                vec;
            bool
                write_it( false);

            while( !STREAM.eof())
            {
                line.clear();
                std::getline( STREAM, line);
                DebugWrite( "<" << line << ">");

                if( line.size() > 5)
                {
                    vec.clear();
                    vec = mystr::SplitString( line);

                    DebugWrite( "nr-splitted-string: " << vec.size());

                    if( line.substr( 0, 1) == "[")
                    {
                        residue = vec[ 1];
                    }
                    else if( vec[ 0] == "[")
                    {
                        if( vec[ 1] == "atoms")
                        {
                            DebugWrite( "atom-section");
                            write_it = true;
                        }
                        else
                        {
                            write_it = false;
                        }
                    }
                    else if( write_it)
                    {
                        vec = mystr::SplitString( line, "\t");
                        DebugWrite( "nr-splitted-string: " << vec.size());
                        DebugWrite( "insert new key and value: " << residue + ":" + vec[ 0] << " - " << vec[ 1]);
                        DebugWrite( "next");
                        ATOM_TYPE_MAP->InsertNewKeyAndValue( residue + ":" + vec[ 0], vec[ 1]);
                    }
                }
            }
        }

        inline
        void ReadMassAndVdwRadiiAndEpsilonMapsFromGromacsCharmm
        (
                std::ifstream &STREAM,
                boost::shared_ptr< store::Map< std::string, float> > &MASS_MAP,
//                boost::shared_ptr< store::Map< std::string, float> > &PARTIAL_CHARGE_MAP,
                boost::shared_ptr< store::Map< std::string, std::pair< float, float> > > &VDW_EPSILON_AND_RADIUS_MAP,
                const bool &USE_CHARMM_TIP3P = false
        )
        {
            std::cout  << __FUNCTION__ << std::endl;

	    // NOTE: if you want to use heavy atoms the logic here has to change, because the .itp is overriding OT instead of having an #ifdef #else structure

            DebugWrite( __FUNCTION__);
            bool
                read_it( true),
                is_heavy_h_set( false),
                is_charmm_tip3p_set( false);
            float
                mass,
//                charge,
                sigma,
                epsilon,
				radius_conversion_factor = pow( 2.0, 1.0 / 6.0) / 2.0; // ~ 0.56123
            //                radius_conversion_factor = 10.0 * pow( 2.0, 1.0 / 6.0) / 2.0; // ~ 5.6123   // test: leads to same values as reading from namd charmm files
            std::string
                line,
                type;
            std::vector< std::string>
                vec;

            std::getline( STREAM, line);
            std::getline( STREAM, line);

            while( !STREAM.eof())
            {
                line.clear();
                std::getline( STREAM, line);
                DebugWrite( "<" << line << ">");
                vec = mystr::SplitString( line, " ", "\t");
                DebugWrite( "nr-splitted-string: " << vec.size());
                DebugWrite( vec);
                if( vec.size() > 6 && read_it && vec[0] != ";")
                {
                    type    = vec[ 0];
                    mass    = mystr::ConvertStringToNumericalValue< float>( vec[ 2]);

//                    charge  = mystr::ConvertStringToNumericalValue< float>( vec[ 3]);

                    sigma   = radius_conversion_factor * mystr::ConvertStringToNumericalValue< float>( vec[ 5]);

                    //                    epsilon = -1.0 * phys::KjouleToKcal * mystr::ConvertStringToNumericalValue< float>( vec[ 6]);     // test: leads to same values as reading from namd charmm files
                    epsilon = -1.0 * mystr::ConvertStringToNumericalValue< float>( vec[ 6]);

                    MASS_MAP->InsertNewKeyAndValue( type, mass);
//                    PARTIAL_CHARGE_MAP->InsertNewKeyAndValue( type, charge);
                    VDW_EPSILON_AND_RADIUS_MAP->InsertNewKeyAndValue( type, std::make_pair( epsilon, sigma));
                }
                else if( vec.size() == 3 && vec[ 1] == "pairtypes")
                {
                    return;
                }
                else if( line.substr( 0, 1) == "#")
                {
                    if( vec[ 0] == "#ifdef" && vec[ 1] == "CHARMM_TIP3P")
                    {
			if( !USE_CHARMM_TIP3P)
			{
			    read_it = false;
			}
                        is_charmm_tip3p_set = true;
                    }
                    else if( vec[ 0] == "#ifdef" && vec[ 1] == "HEAVY_H")
                    {
                        read_it = false;
                        is_heavy_h_set = true;
                    }
                    else if( vec[ 0] == "#endif")
                    {
                        read_it = true;
                    }
                    else if( vec[ 0] == "#else")
                    {
                        if( is_heavy_h_set)
                        {
                            is_heavy_h_set = false;

			    if( is_charmm_tip3p_set && USE_CHARMM_TIP3P)
			    {
				read_it = true;
			    }
			    else if( !is_charmm_tip3p_set && !USE_CHARMM_TIP3P)
			    {
				read_it = true;
			    }
                        }
                        else if( is_charmm_tip3p_set)
                        {
                            is_charmm_tip3p_set = false;
                        }
 //                      else
 //                       {
 //                           read_it = true;
 //                       }
                    }
                }
            }
//            std::cout  << __FUNCTION__ << " done" << std::endl;
        }

        struct MolBlock
        {
            std::string    m_Type;
            size_t         m_NrMolecules;
            size_t         m_AtomsPerMolecule;
        };


        inline
        std::ostream &operator << ( std::ostream &STREAM, const MolBlock &BLOCK)
        {
            STREAM << BLOCK.m_Type << " nr_mols: " << BLOCK.m_NrMolecules << " atoms_per_mol: " << BLOCK.m_AtomsPerMolecule << std::endl;
            return STREAM;
        }


// why ???????????????????????????????????
//        template< typename T1>
        inline
        std::ostream &operator << ( std::ostream &STREAM, const std::vector< MolBlock> &VECTOR)  // TODO: place structs elsewhere!!
        {
            for( /*typename*/ std::vector< MolBlock>::const_iterator itr = VECTOR.begin(); itr != VECTOR.end(); ++itr)
            {
                STREAM << *itr << std::endl;
            }
            return STREAM;
        }


        inline
        void ReadAtomNameAndTypeAndResidueNameAndPartialChargeFromGmxDump
        (
                std::ifstream &STREAM,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOL,
                std::ostream &LOG = std::cout
        )
        {
            //////  DEFINITIONS  //////

            LOG << __FUNCTION__ << std::endl;

            std::string
                line,
                that_line,
                atom_name,
                value_str,
                atom_type;
            bool
                to_continue( true);
            size_t
                total_nr_atoms = 0,
                count( 0),
                nr,
	      residue_id = std::numeric_limits< size_t>::max() - 10,
                prev_res_id = std::numeric_limits< size_t>::max(),
                mol_id = 1,
                mol_count = 1,
                first,
                last,
                line_id,
                atom_count( 0);
            int
                id( -1);
            float
                charge;

            std::vector< std::string>
                vec,
                residue_types;

            std::vector< MolBlock>
                mol_blocks;
            std::vector< MolBlock>::reverse_iterator
                block_itr;

            std::vector< Quartet< size_t, std::string, std::string, float> >
                atom_list;
            std::vector< Quartet< size_t, std::string, std::string, float> >::iterator
                atom_itr;


            //////  CLEAR MOL  //////

            MOL = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>());

            //////  CHOP INITIAL PIECE  //////

            while( !STREAM.eof() && to_continue)
            {
                std::getline( STREAM, line);
                if( line.substr( 0, 9) == "topology:")
                {
                    to_continue = false;
                }
            }

            std::getline( STREAM, line);
            std::getline( STREAM, line);

            if( STREAM.eof())
            {
            	LOG << __FUNCTION__ << " ERROR: end of file, nothing read!" << std::endl;
            	exit( -1);
            }

            //////  GERERAL BLOCK INFO (moltype, nr mols, atoms per mol)  /////////

            vec = mystr::SplitString( line);

            if( vec[0] == "#atoms")
            {
                total_nr_atoms = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
            }
            else
            {
            	LOG << __FUNCTION__ << " WARNING: nr of atoms not found! line: " << line << std::endl;
            }

            LOG << "nr atoms: " << total_nr_atoms << std::endl;

            to_continue = true;

            while( to_continue)
            {
                line.clear();
                std::getline( STREAM, line);
                vec = mystr::SplitString( line);

                if( vec[0] == "molblock")
                {
                    mol_blocks.push_back( MolBlock());
                    DebugWrite( "new block: " << mol_blocks.size());
                    block_itr = mol_blocks.rbegin();
                }
                else if( vec[0] == "moltype")
                {
                    block_itr->m_Type = vec[3].substr( 1, vec[3].length() - 2);
                    DebugWrite( "mol type: " << block_itr->m_Type);
                }
                else if( vec[0] == "#molecules")
                {
                    block_itr->m_NrMolecules = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
                    DebugWrite( "nr mols: " << block_itr->m_NrMolecules);
                }
                else if( vec[ 0] == "#atoms_mol")
                {
                    block_itr->m_AtomsPerMolecule = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
                    count += block_itr->m_NrMolecules * block_itr->m_AtomsPerMolecule;
                    DebugWrite( "atoms per mol: " << block_itr->m_AtomsPerMolecule);
                }
                else if( vec[0] == "ffparams:")
                {
                    DebugWrite( "end block info");
                    to_continue = false;
                }
            }


            LOG << "mol blocks: " << std::endl;
            LOG << mol_blocks << std::endl;


            if( count != total_nr_atoms)
            {
            	LOG << __FUNCTION__ << " ERROR: atoms counts do not match: " << count << " != " << total_nr_atoms << std::endl;
            	exit( -1);
            }

            count = 0;

            to_continue = true;

	    ///////  READ ATOMS & PARAMETERS  /////////

            while( to_continue)
            {
                line.clear();
                std::getline( STREAM, line);


                if( line.size() > 10 && line.substr( 3, 7) == "moltype")
                {
                    ++id;
                }
                else if( line.size() > 15 && line.substr( 9, 5) == "atom ")
                {
                	// read first block with resid and charge
                	residue_types.clear();

                    DebugPrint( LOG, "new atom block! " << id);
                    that_line = line;
                    vec = mystr::SplitString( mystr::TrimString( line));
                    nr = mystr::ConvertStringToNumericalValue< size_t>( vec[1].substr( 1, vec[1].length() - 3));
                    DebugPrint( LOG,  "nr: " << nr);
                    assert( mol_blocks[ id].m_AtomsPerMolecule == nr);
                    atom_list = std::vector< Quartet< size_t, std::string, std::string, float> >( nr);
                    atom_itr = atom_list.begin();
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                    	line.clear();
                        std::getline( STREAM, line);
                        line_id = line.find( "=", line.find( "res"));
                        value_str = line.substr( line_id + 1, line.find( ',',line_id) - line_id - 1);
                        //			std::cout << ">" << value_str << "<" << std::endl;
                        residue_id = mystr::ConvertStringToNumericalValue< size_t>( value_str);
                        atom_itr->first = residue_id;
                        charge = mystr::ConvertStringToNumericalValue< float>( line.substr( 81, 12));
                        atom_itr->forth = charge;
                    }

                    line.clear();
                    std::getline( STREAM, line);
                    assert( line == that_line);
                    atom_itr = atom_list.begin();
                    // read second block with atom names
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        atom_name = line.substr( first, last-first);
                        // DebugWrite( "atom name:" << atom_name);
                        atom_itr->second = atom_name;
                    }

                    std::getline( STREAM, line);
                    atom_itr = atom_list.begin();
                    // read third block translating atom names to types
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        last = line.find_first_of( "\"", first);
                        atom_type = line.substr( first, last-first);
                        // DebugWrite( "atom type:" << atom_type);
                        atom_itr->third = atom_type;
                    }
                    StandardPrint( LOG, "atom list: " << mol_blocks[id].m_Type << " " << atom_list);

                    std::getline( STREAM, line);
                    for( size_t i = 0; i <= residue_id; ++i)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        last = line.find_first_of( "\"", first);
                        residue_types.push_back( line.substr( first, last-first));
                    }
                    // StandardPrint( LOG, "residue types: " << residue_types);
                   LOG << "residue types: ";
                   std::copy( residue_types.begin(), residue_types.end(), std::ostream_iterator< std::string>( LOG, " "));
                   LOG << std::endl;


                    // information complete, read values into molecule and atom_type_map 
                    StandardPrint( LOG, "read " << mol_blocks[ id].m_NrMolecules << " molecules of type: " << mol_blocks[ id].m_Type);
                    for( size_t i = 0; i < mol_blocks[ id].m_NrMolecules; ++i, ++mol_count, ++mol_id)
                    {
                        for( std::vector< Quartet< size_t, std::string, std::string, float> >::const_iterator itr = atom_list.begin(); itr != atom_list.end(); ++itr, ++atom_count)
                        {
                            DebugWrite( itr->first << " " << residue_types[ itr->first] + ":" + itr->second + " " + itr->third);

							//                        store::Map< std::string, std::string>::iterator
							//                            map_itr( ATOM_TYPE_MAP->find( residue_types[ itr->first] + ":" + itr->second));
							//                        if( map_itr == ATOM_TYPE_MAP->end())
							//                        {
							//                            ATOM_TYPE_MAP->insert( std::make_pair( residue_types[ itr->first] + ":" + itr->second, itr->third));
							//                        }
							//                        else if( map_itr->second != itr->third)
							//                        {
							//                            std::cout << "===> a different atom type was already inserted for this \'residue:atom_name\' key: \'" << map_itr->first << "\': " << map_itr->second << " != " << itr->third << std::endl;
							//                            std::cout << "===> molecule block: " << id << " atom: " << atom_count << std::endl;
							//                            exit( -1);
							//                        }

                            if( itr->first != prev_res_id)
							{
                            	prev_res_id = itr->first;
								if( itr != atom_list.begin())
								{
									++mol_id;
								}
							}

                            MOL->AddAtom( boost::shared_ptr< mol::Atom>( new Atom( itr->third, itr->second, /*mol_blocks[ id].m_Type*/ residue_types[ itr->first], math::Vector3N(), itr->forth, mol_id, atom_count, mol_count))); // TODO: fix residue id!
                        }
                    }

                    if( id + 1 == int( mol_blocks.size()))
                    {
                        to_continue = false;
                    }
                }
            }

            if( total_nr_atoms != MOL->GetAtoms().size())
            {
                LOG << "===> incorrect number of atoms read! molecule has: " << MOL->GetAtoms().size() << " expected: " << total_nr_atoms << std::endl;
            }

            //            LOG << __FUNCTION__ << " done" << std::endl;


        } // end ReadAtomTypeTranslatorAndAtomInformationFromGmxDump  





        inline
        void ReadAtomNameAndTypeAndResidueNameAndPartialChargeFromGmxDump
        (
                std::ifstream &STREAM,
                store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > &MOLS,
//                boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > &MOLS,
                std::ostream &LOG = std::cout
        )
        {
            //////  DEFINITIONS  //////

            LOG << __FUNCTION__ << std::endl;

            std::string
                line,
                that_line,
                atom_name,
                atom_type;
            bool
                to_continue( true);
            size_t
                total_nr_atoms,
                count( 0),
                nr,
                residue_id,
                first,
                last,
                atom_count( 0);
            int
                id( -1);
            float
                charge;

            std::vector< std::string>
                vec,
                residue_types;

            std::vector< MolBlock>
                mol_blocks;
            std::vector< MolBlock>::reverse_iterator
                block_itr;

            std::vector< Quartet< size_t, std::string, std::string, float> >
                atom_list;
            std::vector< Quartet< size_t, std::string, std::string, float> >::iterator
                atom_itr;


            //////  CLEAR MOLS  //////

//            MOLS = boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >());
            MOLS.clear();

            //////  CHOP INITIAL PIECE  //////

            while( !STREAM.eof() && to_continue)
            {
                std::getline( STREAM, line);
                if( line.substr( 0, 9) == "topology:")
                {
                    to_continue = false;
                }
            }

            std::getline( STREAM, line);
            std::getline( STREAM, line);

            //////  HERE WE GO //////

            vec = mystr::SplitString( line);

            if( vec[0] == "#atoms")
            {
                total_nr_atoms = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
            }

            LOG << "nr atoms: " << total_nr_atoms << std::endl;

            to_continue = true;

            while( to_continue)
            {
                line.clear();
                std::getline( STREAM, line);
                vec = mystr::SplitString( line);

                if( vec[0] == "molblock")
                {
                    mol_blocks.push_back( MolBlock());
                    DebugWrite( "new block: " << mol_blocks.size());
                    block_itr = mol_blocks.rbegin();
                }
                else if( vec[0] == "moltype")
                {
                    block_itr->m_Type = vec[3].substr( 1, vec[3].length() - 2);
                    DebugWrite( "mol type: " << block_itr->m_Type);
                }
                else if( vec[0] == "#molecules")
                {
                    block_itr->m_NrMolecules = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
                    DebugWrite( "nr mols: " << block_itr->m_NrMolecules);
                }
                else if( vec[ 0] == "#atoms_mol")
                {
                    block_itr->m_AtomsPerMolecule = mystr::ConvertStringToNumericalValue< size_t>( vec[2]);
                    count += block_itr->m_NrMolecules * block_itr->m_AtomsPerMolecule;
                    DebugWrite( "atoms per mol: " << block_itr->m_AtomsPerMolecule);
                }
                else if( vec[0] == "ffparams:")
                {
                    DebugWrite( "end block info");
                    to_continue = false;
                }
            }


            LOG << "mol blocks: " << std::endl;
            LOG << mol_blocks << std::endl;


            assert( count == total_nr_atoms);

            count = 0;

            to_continue = true;

            // parameter section
            while( !STREAM.eof() && to_continue)
            {
                line.clear();
                std::getline( STREAM, line);


                if( line.substr( 3, 7) == "moltype")
                {
                    ++id;
                }
                else if( line.substr( 9, 5) == "atom ")
                {
                    DebugPrint( LOG, "new atom block! " << id);
                    that_line = line;
                    vec = mystr::SplitString( mystr::TrimString( line));
                    nr = mystr::ConvertStringToNumericalValue< size_t>( vec[1].substr( 1, vec[1].length() - 3));
                    DebugPrint( LOG,  "nr: " << nr);
                    assert( mol_blocks[ id].m_AtomsPerMolecule == nr);
                    atom_list = std::vector< Quartet< size_t, std::string, std::string, float> >( nr);
                    atom_itr = atom_list.begin();
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                        std::getline( STREAM, line);
                        residue_id = mystr::ConvertStringToNumericalValue< size_t>( line.substr( 135, 5));
                        atom_itr->first = residue_id;
                        charge = mystr::ConvertStringToNumericalValue< float>( line.substr( 81, 12));
                        atom_itr->forth = charge;
                    }

                    line.clear();
                    std::getline( STREAM, line);
                    assert( line == that_line);
                    atom_itr = atom_list.begin();
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        atom_name = line.substr( first, last-first);
//                        DebugWrite( "atom name:" << atom_name);
                        atom_itr->second = atom_name;
                    }

                    std::getline( STREAM, line);
                    atom_itr = atom_list.begin();
                    for( size_t i = 0; i < nr; ++i, ++atom_itr)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        last = line.find_first_of( "\"", first);
                        atom_type = line.substr( first, last-first);
//                        DebugWrite( "atom type:" << atom_type);
                        atom_itr->third = atom_type;
                    }
                    DebugPrint( LOG, "atom list: " << atom_list);

                    std::getline( STREAM, line);
                    for( size_t i = 0; i <= residue_id; ++i)
                    {
                        std::getline( STREAM, line);
                        first = 1 + line.find( "\"",12);
                        last = line.find( "\"", first);
                        last = line.find_first_of( "\"", first);
                        residue_types.push_back( line.substr( first, last-first));
                    }
//                    DebugPrint( LOG, "residues: " << residue_types);

                    // information complete, read values into molecule and atom_type_map

                    DebugPrint( LOG, "read " << mol_blocks[ id].m_NrMolecules << " molecules of type: " << mol_blocks[ id].m_Type);
                    for( size_t mol_count = 0; mol_count < mol_blocks[ id].m_NrMolecules; ++mol_count)
                    {
                        for( std::vector< Quartet< size_t, std::string, std::string, float> >::const_iterator itr = atom_list.begin(); itr != atom_list.end(); ++itr, ++atom_count)
                        {
                            DebugWrite( itr->first << " " << residue_types[ itr->first] + ":" + itr->second + " " + itr->third);

    //                        store::Map< std::string, std::string>::iterator
    //                            map_itr( ATOM_TYPE_MAP->find( residue_types[ itr->first] + ":" + itr->second));
    //                        if( map_itr == ATOM_TYPE_MAP->end())
    //                        {
    //                            ATOM_TYPE_MAP->insert( std::make_pair( residue_types[ itr->first] + ":" + itr->second, itr->third));
    //                        }
    //                        else if( map_itr->second != itr->third)
    //                        {
    //                            std::cout << "===> a different atom type was already inserted for this \'residue:atom_name\' key: \'" << map_itr->first << "\': " << map_itr->second << " != " << itr->third << std::endl;
    //                            std::cout << "===> molecule block: " << id << " atom: " << atom_count << std::endl;
    //                            exit( -1);
    //                        }

                            MOLS[ mol_count]->AddAtom( boost::shared_ptr< mol::Atom>( new Atom( itr->third, itr->second, residue_types[ itr->first], math::Vector3N(), itr->forth, itr->first, atom_count))); // TODO: fix residue id!
                        }
                    }

                    if( id + 1 == int( mol_blocks.size()))
                    {
                        to_continue = false;
                    }
                }
            }

//            if( total_nr_atoms != MOL->GetAtoms().size())
//            {
//                LOG << "===> incorrect number of atoms read! molecule has: " << MOL->GetAtoms().size() << " expected: " << total_nr_atoms << std::endl;
//            }

//            LOG << __FUNCTION__ << " done" << std::endl;


        } // end ReadAtomTypeTranslatorAndAtomInformationFromGmxDump MOLS!!





//        inline
//        void
//        ReadGridFromFile(  std::istream &STREAM, boost::shared_ptr< mol::MoleculeForceGrid> &GRID)
//        {
//            if( !GRID)
//            {
////                GRID = boost::shared_ptr< mol::MoleculeForceGrid>( new MoleculeForceGrid());
//            }
//
//            std::string
//                line,
//                another;
//
//            std::getline( STREAM, line);
//            std::getline( STREAM, another);
//
//            if( line != "mol::MoleculeForceGrid" ||  another != "mol::AtomForceGrid")
//            {
//                std::cerr << "wrong header!" << std::endl;
//                exit( -1);
//            }
//
//            std::getline( STREAM, line);
//            if( line == "")
//            {
//
//            }
//            std::getline( STREAM, line);
//
//
//        }

    } // end namespace file
} // end namespace mol


#endif /* SIMPLE_MOLECULE_FILE_HANDLER_H_ */
