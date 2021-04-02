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


#include "../../include/molecule/molecule_factory.h"

namespace mol
{
    namespace factory
    {
        void Read( const store::ShPtrVec< SimpleMolecule< Atom> > &MOL, const std::istream &READ, const std::string &INPUT_FILE_FORMAT)
        {

        }

        void CalculateGridLimits( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOL, math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA, const float &CUTOFF)
        {
            // declarations
            MIN.SetAll( std::numeric_limits< float>::max()),
            MAX.SetAll( std::numeric_limits< float>::min());
            std::vector< float>::iterator min_itr, max_itr;
            std::vector< float>::const_iterator pos_itr, delta_itr;

            // determine min, max from molecule
            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr = MOL->GetAtoms().begin(); atom_itr != MOL->GetAtoms().end(); ++atom_itr)
            {
                pos_itr = ( *atom_itr)->GetPosition().begin();
                min_itr = MIN.begin();
                max_itr = MAX.begin();
                for( ; min_itr != MIN.end(); ++pos_itr, ++min_itr, ++max_itr)
                {
                    if( *pos_itr > *max_itr)
                    { *max_itr = *pos_itr;}
                    else if( *pos_itr < *min_itr)
                    { *min_itr = *pos_itr;}
//                    std::cout  << "min: " << *min_itr << " max: " << *max_itr << " pos: " << *pos_itr << std::endl;
                }
            }

            // add / subtract cutoff to limits
            min_itr = MIN.begin();
            max_itr = MAX.begin();
            for( ; min_itr != MIN.end(); ++min_itr, ++max_itr)
            {
//                std::cout  << "before: min: " << *min_itr << " max: " << *max_itr  << std::endl;
                *max_itr += CUTOFF;
                *min_itr -= CUTOFF;
//                std::cout  << "after:  min: " << *min_itr << " max: " << *max_itr  << std::endl;
            }

            AdjustLimits( 1e-5, MIN, MAX, DELTA);
//            // adjust min and max to match delta
//            min_itr = MIN.begin();
//            max_itr = MAX.begin();
//            delta_itr = DELTA.begin();
//            for( ; min_itr != MIN.end(); ++delta_itr, ++min_itr, ++max_itr)
//            {
//                int max_size( int( *max_itr /= *delta_itr));
//                if( *max_itr - float( max_size) * *delta_itr > 1e-8)
//                { *max_itr = float( max_size + 1) * *delta_itr;}
//
//                int min_size( int( *min_itr /= *delta_itr));
//                if( float( min_size) * *delta_itr - *min_itr > 1e-8)
//                { *min_itr = float( min_size - 1) * *delta_itr;}
//                std::cout  << "min: " << *min_itr << " max: " << *max_itr  << std::endl;
//            }

            std::cout << "AtomForceGrid::SetGridLimits:" << std::endl;
            std::cout << "AtomForceGrid::Minimum:  " << MIN << std::endl;
            std::cout << "AtomForceGrid::Maximum:  " << MAX << std::endl;
            std::cout << "AtomForceGrid::Delta:    " << DELTA << std::endl;
        }

        store::Vector3N< size_t> CalculateNumberOfBins( const math::Vector3N &MIN, const math::Vector3N &MAX, const math::Vector3N &DELTA)
        {
//            std::cout << "... calculate number of bins ... " << std::endl;
            store::Vector3N< size_t> bins;
            std::vector< float>::const_iterator min_itr( MIN.begin()), max_itr( MAX.begin()), delta_itr( DELTA.begin());
            for( std::vector< size_t>::iterator itr( bins.begin()); itr != bins.end() && max_itr != MAX.end(); ++max_itr, ++min_itr, ++delta_itr, ++itr)
            {
//                std::cout << "unrounded: "<< ( *max_itr - *min_itr) / ( *delta_itr);
                *itr = math::Round< size_t>( ( *max_itr - *min_itr) / ( *delta_itr));
//                std::cout << " rounded: " << *itr << std::endl;
            }
            std::cout << "number of bins: (" << bins(0) << " " << bins(1)<< " " << bins(2) << ")" << std::endl;
            return bins;
        } // end CalculateNrOfBins


        void AdjustLimits( const float &PRECISION, const math::Vector3N &MIN, math::Vector3N &MAX, const math::Vector3N &DELTA)
        {
//            std::cout << "...adjust limits ..." << MAX << std::endl;
            for( size_t i = 0; i < 3; ++i)
            {
                float d( ( MAX[i] - MIN[i]) / DELTA[i]);
                float rest( d - size_t( d));
                if( rest > PRECISION && rest < 0.5)
                {
                    MAX[i] = MIN[i] + size_t( d) * DELTA[i];
                    std::cout << "maximum for " << i << " decreased to " << MAX[i] << "!" << std::endl;
                }
                else if( rest != 0)
                {
                    MAX[i] = MIN[i] + ( size_t( d) + 1) * DELTA[i];
                    std::cout << "maximum for " << i << " increased to " << MAX[i] << "!" << std::endl;
                }
            }
//            std::cout << "... limit adjustment done: " << MAX << std::endl;
        }


        void AdjustLimits( math::Vector3N &MIN, math::Vector3N &MAX, math::Vector3N &DELTA, const math::Vector3N &FIXED_MIN, const math::Vector3N &FIXED_MAX, const math::Vector3N &FIXED_DELTA)
        {
            DELTA = FIXED_DELTA;
            for( int i = 0; i < 3; ++i)
            {
                if( MIN[i] < FIXED_MIN[i])
                {
                    MIN[i] = FIXED_MIN[i] - float( int( ( FIXED_MIN[i] - MIN[i]) / DELTA[i]) + 1) * DELTA[i];
                }
                else
                {
                    MIN[i] = FIXED_MIN[i]; // + float( int( ( MIN[i] - FIXED_MIN[i]) / DELTA[i])) * DELTA[i];
                }

                if( MAX[i] > FIXED_MAX[i])
                {
                    MAX[i] = FIXED_MAX[i] + float( int( ( MAX[ i] - FIXED_MAX[ i]) / DELTA[ i]) + 1) * DELTA[i];
                }
                else
                {
                    MAX[i] = FIXED_MAX[i]; // - float( int( ( FIXED_MAX[ i] - MAX[ i]) / DELTA[ i])) * DELTA[i];
                }
            }
        }

        void
        BuildImplicitMolecule
        (
                const CommandLineManager &COMMAND,
                boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE
        )
        {

            // TODO : merge with membrane.cc::Read( read, cmd)

            std::string
                file,
                force_field_type( "charmm");

            boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
                vdw_epsilon_and_radius_map( new store::Map< std::string, std::pair< float, float> >());

            std::ofstream
                write;
            std::ifstream
                read;

            if( COMMAND.IsFlagSet( "force_field"))
            {
                force_field_type = COMMAND.GetArgumentStringForFlag( "force_field");
                StandardWrite( "force field: " << force_field_type);
            }

            if( force_field_type == "charmm")
            {
            	StandardWrite(  "read pdb and psf files");

                std::string file_name = COMMAND.GetArgumentStringForFlag( "implicit_xpsf");

                Open( read, file_name);

                mol::file::ReadMoleculeFromPSF( read, IMPLICIT_MOLECULE);

				Close( read);

				file_name = COMMAND.GetArgumentStringForFlag( "implicit_pdb");

                Open( read, file_name);

                mol::file::ReadNewCoordinatesFromPdb( read, IMPLICIT_MOLECULE);

				Close( read);

				StandardWrite(  "read charm force field parameter file");

                if( COMMAND.IsFlagSet( "parameter_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "parameter_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/charmm/par_all27_prot_lipid.prm");
                    if( !read)
                    {
                        read.open( "par_all27_prot_lipid.prm");
                    }
                    if( !read)
                    {
                        std::cout << "par_all27_prot_lipid.prm not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/charmm/', '.'" << std::endl;
                        std::cout << "try '-parameter_file' flag" << std::endl;
                        exit( -1);
                    }
                }

                mol::file::ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( read, vdw_epsilon_and_radius_map);
                Close( read);

            }
            /*
            else if( force_field_type == "charmm_old")
            {
                StandardWrite(  "read charm force field parameter and topology files");
                StandardWrite(  "read atom type translating map ...");

                   if( COMMAND.IsFlagSet( "topology_file"))
                {
                    file = COMMAND.GetArgumentStringForFlag( "topology_file");
                    Open( read, file);
                }
                else
                {
                    file = "/usr/local/tara/lib/charmm/top_all27_prot_lipid.rtf";
                    read.open( file.c_str());
                    if( !read)
                    {
                        file = "top_all27_prot_lipid.rtf";
                        read.open( file.c_str());
                        if( !read)
                        {
                            std::cout << "top_all27_prot_lipid.rtf not found!" << std::endl;
                            std::cout << "paths searched: '/usr/local/tara/lib/charmm/', '.'" << std::endl;
                            std::cout << "try '-topology_file' flag" << std::endl;
                            exit( -1);
                        }
                    }
                }

                mol::file::ReadAtomTypeTranslator( read, atom_type_map);
                 Close( read);
                DebugWrite( "atom type map: " << atom_type_map);

                StandardWrite(  "read partial charges ...");

                Open( read, file);
                mol::file::ReadPartialChargeMapFromCharmmTopologyFile( read, partial_charge_map);
                DebugWrite( "partial charge map: " <<  partial_charge_map);
                 Close( read);
//                StandardWrite(  "... done");

                StandardWrite(  "read mass values ... ");
                Open( read, file);
                mol::file::ReadMassFromCharmm( read, mass_map);
//                StandardWrite(  "... done");
                 Close( read);

                StandardWrite(  "read epsilon values ... ");
                if( COMMAND.IsFlagSet( "parameter_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "parameter_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/charmm/par_all27_prot_lipid.prm");
                    if( !read)
                    {
                        read.open( "par_all27_prot_lipid.prm");
                    }
                    if( !read)
                    {
                        std::cout << "par_all27_prot_lipid.prm not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/charmm/', '.'" << std::endl;
                        std::cout << "try '-parameter_file' flag" << std::endl;
                        exit( -1);
                    }
                }

                mol::file::ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( read, vdw_epsilon_and_radius_map);
        //            StandardWrite(  partial_charge_map;
                 Close( read);

                if( COMMAND.IsFlagSet( "implicit_molecule"))
                {
                    StandardWrite( "read implicit molecule from pdb ...");
                    Open( read, COMMAND.GetArgumentStringForFlag( "implicit_molecule"));
                    IMPLICIT_MOLECULE = mol::file::ReadFromPdb< mol::Atom>( read, atom_type_map);
                     Close( read);

 //                   StandardWrite( IMPLICIT_MOLECULE);

                    // TODO: getting rid of reading partial charges externally?
                    StandardWrite(  "read partial charges into molecule ...");
                    mol::file::ReadPartialChargesIntoMolecule( partial_charge_map, IMPLICIT_MOLECULE);
//                    StandardWrite(  "... done");
                    DebugWrite( "write partial_charge_map");
                    DebugWrite( partial_charge_map);
                }
//                StandardWrite(  "... done");
            }
            */
            else if( /*COMMAND.IsFlagSet( "build_grid") &&*/ force_field_type == "gromacs_charmm")
            {
                StandardWrite( "use gromacs charm force field");

                boost::shared_ptr< store::Map< std::string, float> >
                    partial_charge_map( new store::Map< std::string, float>());
                boost::shared_ptr< store::Map< std::string, float> >
                    mass_map( new store::Map< std::string, float>());
                boost::shared_ptr< store::Map< std::string, std::string> >
                    atom_type_map( new store::Map< std::string, std::string>());



                ////     READ IMPLICIT MOLECULE   ////
                StandardWrite( "read implicit molecule from gmx dump file...");

//                file = COMMAND.GetArgumentStringForFlag( "implicit_molecule");
                Open( read, COMMAND.GetArgumentStringForFlag( "implicit_molecule"));
                mol::file::ReadAtomNameAndTypeAndResidueNameAndPartialChargeFromGmxDump( read, IMPLICIT_MOLECULE);
                mol::file::ReadNewCoordinatesFromGmxDump( read, IMPLICIT_MOLECULE);
                 Close( read);

                if( COMMAND.IsFlagSet( "parameter_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "parameter_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/gromacs/ffcharmm27nb.itp");
                    if( !read)
                    {
                        read.open( "ffcharmm27nb.itp");
                    }
                    if( !read)
                    {
                        std::cout << "ffcharmm27nb.itp not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/gromacs/', '.'" << std::endl;
                        std::cout << "try '-parameter_file' flag" << std::endl;
                        exit( -1);
                    }
                }

                mol::file::ReadMassAndVdwRadiiAndEpsilonMapsFromGromacsCharmm( read, mass_map, /*partial_charge_map,*/ vdw_epsilon_and_radius_map);
                Close( read);

                StandardWrite(  "read mass into molecule ...");
                mol::file::ReadMassIntoMolecule( mass_map, IMPLICIT_MOLECULE);

                DebugWrite( "write mass_map");
                DebugWrite( mass_map);
            }


            if( COMMAND.IsFlagSet( "use_default_radii"))
            {
                mol::file::ReadDefaultRadii( vdw_epsilon_and_radius_map);
            }
            StandardWrite(  "read vdw radii and epsilon values into molecule ...");
            mol::file::ReadVdwEpsilonAndRadiusIntoMolecule( vdw_epsilon_and_radius_map, IMPLICIT_MOLECULE);


            DebugWrite( "write vdw_epsilon_and_radius_map");
            DebugWrite( vdw_epsilon_and_radius_map);
        //            DebugWrite( "write atom_type_map");
        //            DebugWrite( atom_type_map);


            //TODO: check whether mass_map and partial_charge_map can be read into one object!!

            if( COMMAND.IsFlagSet( "write_implicit_molecule"))
            {
                StandardWrite(  "write implicit molecule to file ...");
                Open( write, COMMAND.GetArgumentStringForFlag( "write_implicit_molecule"));
                write << *IMPLICIT_MOLECULE;
                Close( write);
            }

        }  // end BuildImplicitMolecule


        void
        BuildSimpleImplicitMolecule
        (
                const CommandLineManager &COMMAND,
                boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &IMPLICIT_MOLECULE
        )
        {

            // TODO: get rid of this???

            std::string
                file,
                force_field_type( "charmm");

            boost::shared_ptr< store::Map< std::string, float> >
                partial_charge_map( new store::Map< std::string, float>());
            boost::shared_ptr< store::Map< std::string, float> >
                mass_map( new store::Map< std::string, float>());
            boost::shared_ptr< store::Map< std::string, std::pair< float, float> > >
                vdw_epsilon_and_radius_map( new store::Map< std::string, std::pair< float, float> >());
            boost::shared_ptr< store::Map< std::string, std::string> >
                atom_type_map( new store::Map< std::string, std::string>());


            std::ofstream
                write;
            std::ifstream
                read;

            if( COMMAND.IsFlagSet( "force_field"))
            {
                force_field_type = COMMAND.GetArgumentStringForFlag( "force_field");
                StandardWrite( "force field: " << force_field_type);
            }

            if( force_field_type == "charmm")
            {
                StandardWrite(  "read charm force field parameter and topology files");
                StandardWrite(  "read atom type translating map ...");

                if( COMMAND.IsFlagSet( "topology_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "topology_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/charmm/top_all27_prot_lipid.rtf");
                    if( !read)
                    {
                        read.open( "top_all27_prot_lipid.rtf");
                    }
                    if( !read)
                    {
                        std::cout << "top_all27_prot_lipid.rtf not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/charmm/', '.'" << std::endl;
                        std::cout << "try '-topology_file' flag" << std::endl;
                        exit( -1);
                    }

                }

                mol::file::ReadAtomTypeTranslator( read, atom_type_map);
                Close( read);

                DebugWrite( atom_type_map);

                StandardWrite( "read partial charges ...");

                Open( read, file);
                mol::file::ReadPartialChargeMapFromCharmmTopologyFile( read, partial_charge_map);
        //            StandardWrite(  partial_charge_map;
                   Close( read);
//                StandardWrite(  "... done");

                StandardWrite( "read mass values ... ");
                Open( read, file);
                mol::file::ReadMassFromCharmm( read, mass_map);
//                StandardWrite(  "... done");
                   Close( read);

                StandardWrite( "read epsilon values ... ");
                if( COMMAND.IsFlagSet( "parameter_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "parameter_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/charmm/par_all27_prot_lipid.prm");
                    if( !read)
                    {
                        read.open( "par_all27_prot_lipid.prm");
                    }
                    if( !read)
                    {
                        std::cout << "par_all27_prot_lipid.prm not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/charmm/', '.'" << std::endl;
                        std::cout << "try '-parameter_file' flag" << std::endl;
                        exit( -1);
                    }
                }

                mol::file::ReadVdwEpsilonAndRadiusMapFromCharmmParameterFile( read, vdw_epsilon_and_radius_map);
        //            StandardWrite(  partial_charge_map;
                  Close( read);

                if( COMMAND.IsFlagSet( "implicit_molecule"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "implicit_molecule"));
                    IMPLICIT_MOLECULE = mol::file::ReadFromPdb< mol::Atom>( read, atom_type_map);
                      Close( read);

                // TODO: getting rid of reading partial charges externally?
                StandardWrite( "read partial charges into molecule ...");
                mol::file::ReadPartialChargesIntoMolecule( partial_charge_map, IMPLICIT_MOLECULE);
//                StandardWrite(  "... done");
                DebugWrite( "write partial_charge_map");
                DebugWrite( partial_charge_map);
                }
//                StandardWrite(  "... done");
            }
            else if( force_field_type == "gromacs_charmm")
            {
                StandardWrite( "use gromacs charm force field");

                ////     READ IMPLICIT MOLECULE   ////
                StandardWrite( "read implicit molecule from gmx dump file...");

                file = COMMAND.GetArgumentStringForFlag( "implicit_molecule");
                Open( read, file);
                mol::file::ReadAtomNameAndTypeAndResidueNameAndPartialChargeFromGmxDump( read, IMPLICIT_MOLECULE);
                mol::file::ReadNewCoordinatesFromGmxDump( read, IMPLICIT_MOLECULE);
                  Close( read);


                if( COMMAND.IsFlagSet( "parameter_file"))
                {
                    Open( read, COMMAND.GetArgumentStringForFlag( "parameter_file"));
                }
                else
                {
                    read.open( "/usr/local/tara/lib/gromacs/ffcharmm27nb.itp");
                    if( !read)
                    {
                        read.open( "ffcharmm27nb.itp");
                    }
                    if( !read)
                    {
                        std::cout << "ffcharmm27nb.itp not found!" << std::endl;
                        std::cout << "paths searched: '/usr/local/tara/lib/gromacs/', '.'" << std::endl;
                        std::cout << "try '-parameter_file' flag" << std::endl;
                        exit( -1);
                    }
                }

                mol::file::ReadMassAndVdwRadiiAndEpsilonMapsFromGromacsCharmm( read, mass_map, /*partial_charge_map,*/ vdw_epsilon_and_radius_map);
                  Close( read);
            }

        //        if( COMMAND.IsFlagSet( "build_grid"))
        //        {
            DebugWrite( "write vdw_epsilon_and_radius_map");
            DebugWrite( vdw_epsilon_and_radius_map);
        //            DebugWrite( "write atom_type_map");
        //            DebugWrite( atom_type_map);
            DebugWrite( "write mass_map");
            DebugWrite( mass_map);

            StandardWrite(  "read mass into molecule ...");
            mol::file::ReadMassIntoMolecule( mass_map, IMPLICIT_MOLECULE);
//            StandardWrite(  "... done");

            StandardWrite(  "read vdw epsilon into molecule ...");
            mol::file::ReadVdwEpsilonAndRadiusIntoMolecule( vdw_epsilon_and_radius_map, IMPLICIT_MOLECULE);
 //           StandardWrite(  "... done");
        //        }

            //TODO: check whether mass_map, vdw_epsilon_and_radius_map, vdw_epsilon_and_radius_map can be read into one object!!

            if( COMMAND.IsFlagSet( "write_implicit_molecule"))
            {
                StandardWrite(  "write implicit molecule to file ...");
                Open( write, COMMAND.GetArgumentStringForFlag( "write_implicit_molecule"));
                write << *IMPLICIT_MOLECULE;
                Close( write);
            }

        } // end BuildSimpleImplicitMolecule

		void
		CalcVdwFactors( boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE)
		{
            for( std::vector< boost::shared_ptr< mol::Atom> >::const_iterator atom_itr = MOLECULE->GetAtoms().begin(); atom_itr != MOLECULE->GetAtoms().end(); ++atom_itr)
            {
            	( *atom_itr)->CalcVdwFactors();
            }
		}


    } // end namespace factory
} // end namespace mol
