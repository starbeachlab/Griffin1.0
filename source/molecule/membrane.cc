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


#include "../../include/molecule/membrane.h"


namespace mol
{


        const store::Limits3D &Membrane::GetTotalLimits() const
        {
            return m_TotalLimits;
        }

        std::pair< float, float> Membrane::GetXLimits() const
        {
            return std::make_pair( m_TotalLimits.GetXMin(), m_TotalLimits.GetXMax());
        }

        std::pair< float, float> Membrane::GetYLimits() const
        {
            return std::make_pair( m_TotalLimits.GetYMin(), m_TotalLimits.GetYMax());
        }

        std::vector< float> Membrane::GetZLimits() const
        {
            std::vector< float> limits( 5);

            limits[ 0] = m_TotalLimits.GetZMin();
            limits[ 1] = m_LipidMin;
            limits[ 2] = m_Center;
            limits[ 3] = m_LipidMax;
            limits[ 4] = m_TotalLimits.GetZMax();

            return limits;
        }

        void Membrane::SetXLimits( const float &MIN, const float &MAX)
        {
            m_TotalLimits.SetXMin( MIN);
            m_TotalLimits.SetXMax( MAX);
        }

        void Membrane::SetYLimits( const float &MIN, const float &MAX)
        {
            m_TotalLimits.SetYMin( MIN);
            m_TotalLimits.SetYMax( MAX);
        }

        void Membrane::SetZLimits( const float &LOWER_SOLUTION, const float &LOWER_LIPID, const float &UPPER_LIPID, const float &UPPER_SOLUTION)
        {
            m_TotalLimits.SetZMin( LOWER_SOLUTION);
            m_LipidMin = LOWER_LIPID;
            m_Center = 0.5 * ( LOWER_LIPID + UPPER_LIPID);
            m_LipidMax = UPPER_LIPID;
            m_TotalLimits.SetZMax( UPPER_SOLUTION);
        }

        void Membrane::SetZLimits( const float &LOWER_SOLUTION, const float &LOWER_LIPID, const float &CENTER, const float &UPPER_LIPID, const float &UPPER_SOLUTION)
        {
            m_TotalLimits.SetZMin( LOWER_SOLUTION);
            m_LipidMin = LOWER_LIPID;
            m_Center = CENTER;
            m_LipidMax = UPPER_LIPID;
            m_TotalLimits.SetZMax( UPPER_SOLUTION);
        }

        std::vector< std::string> Membrane::GetMoleculeTypesInLayer( const size_t &LAYER) const
        {
            return m_MoleculeLayer[ LAYER].GetKeys();
        }

        std::vector< std::string> Membrane::GetMoleculeTypes() const
        {
            std::vector< std::string> keys;
            for( size_t i = 0; i < 4; ++i)
            {
                keys = Merge( keys, m_MoleculeLayer[ i].GetKeys());
            }
            return RemoveDuplicates( keys);
        }



        const store::Map< std::string, boost::shared_ptr< store::ShPtrVec< SimpleMolecule< Atom> > > > &
        Membrane::GetMoleculeLayer( const size_t &LAYER) const
        {
            return m_MoleculeLayer[ LAYER];
        }

        const std::vector< std::string> &
        Membrane::GetLipidKeys() const
        {
            return m_LipidKeys;
        }

        const std::vector< std::string> &
        Membrane::GetSolKeys() const
        {
            return m_LipidKeys;
        }


        store::ShPtrVec< SimpleMolecule< Atom> >
        Membrane::GetMolecules() const
        {
            store::ShPtrVec< SimpleMolecule< Atom> > all;
            for( size_t i = 0; i < 4; ++i)
            {
                for( std::map< std::string, boost::shared_ptr< store::ShPtrVec< SimpleMolecule< Atom> > > >::const_iterator itr = m_MoleculeLayer[ i].begin(); itr != m_MoleculeLayer[ i].end(); ++itr)
                {
                    all = Merge( all, *itr->second);
                }
            }
            return all;
        }


        store::Map< std::string, float>
        Membrane::Densities( const size_t &LAYER) const
        {
            store::Map< std::string, float> densities;
            float volume( FourLayerBoxLimits( LAYER).Volume());
            for( std::map< std::string, boost::shared_ptr< store::ShPtrVec< SimpleMolecule< Atom> > > >::const_iterator itr = m_MoleculeLayer[ LAYER].begin(); itr != m_MoleculeLayer[ LAYER].end(); ++itr)
            {
                densities.InsertNewKeyAndValue( itr->first, float( itr->second->size()) / volume);
            }
            return densities;
        }

        store::Limits3D
        Membrane::FourLayerBoxLimits( const size_t &LAYER) const
        {
            store::Limits3D limits( m_TotalLimits);
            switch( LAYER)
            {
            case 0:
                limits.SetZMax( m_LipidMin);
                break;
            case 1:
                limits.SetZMin( m_LipidMin);
                limits.SetZMax( m_Center);
                break;
            case 2:
                limits.SetZMin( m_Center);
                limits.SetZMax( m_LipidMax);
                break;
            case 3:
                limits.SetZMin( m_LipidMax);
                break;
            default:
                exit( -1);
            }
            return limits;
        }


        bool
        Membrane::IsValidKey( const std::string &KEY) const
        {
            for( size_t i = 0; i < m_MoleculeLayer.size(); ++i)
            {
                if( m_MoleculeLayer[ i].IsValidKey( KEY))
                {
                    return true;
                }
            }
            return false;
        }


        size_t
        Membrane::LayerID( const float &Z_COORDINATE) const
        {  // TODO: more efficient by first asking pos or neg?
            if( Z_COORDINATE >= m_TotalLimits.GetZMin() && Z_COORDINATE < m_LipidMin)
            {
                return 0;
            }
            else if( Z_COORDINATE >= m_LipidMin && Z_COORDINATE < m_Center)
            {
                return 1;
            }
            else if( Z_COORDINATE >= m_Center && Z_COORDINATE < m_LipidMax)
            {
                return 2;
            }
            else if( Z_COORDINATE >= m_LipidMax && Z_COORDINATE <= m_TotalLimits.GetZMax())
            {
                return 3;
            }
            return std::numeric_limits< size_t>::max();
        }


        void
        Membrane::AdjustLimits( const std::string &MOLECULE_TYPE, const math::Vector3N &POSITION/*, const float &RADIUS*/)
        {
            DebugWrite( __FUNCTION__);
            if( std::find( m_LipidKeys.begin(), m_LipidKeys.end(), MOLECULE_TYPE) != m_LipidKeys.end())
            {
//                for( size_t i = 0; i < 3; ++i)
//                {
//                    if( POSITION( i) /*+ RADIUS*/ > m_LipidLimits.GetMax( i))
//                    {
//                        DebugWrite( "raise lipid max " <<  i << "  " << POSITION( i)/* + RADIUS */);
//                        m_LipidLimits.SetMax( i, POSITION( i)/* + RADIUS */);
//                    }
//                    else if( POSITION( i) /* - RADIUS */ < m_LipidLimits.GetMin( i))
//                    {
//                        DebugWrite( "lower lipid min " <<  i << " " << POSITION( i) /* - RADIUS */);
//                        m_LipidLimits.SetMin( i, POSITION( i) /* - RADIUS */);
//                    }
//                }
            }
            else
            {
                for( size_t i = 0; i < 3; ++i)
                {
                    if( POSITION( i)/* + RADIUS */ > m_TotalLimits.GetMax( i))
                    {
                        DebugWrite( "raise total max " <<  i << "  " << POSITION( i)/* + RADIUS */);
                        m_TotalLimits.SetMax( i, POSITION( i)/* + RADIUS */);
                    }
                    else if( POSITION( i) /* - RADIUS */ < m_TotalLimits.GetMin( i))
                    {
                        DebugWrite( "lower total min " <<  i << " " << POSITION( i) /* - RADIUS */);
                        m_TotalLimits.SetMin( i, POSITION( i) /* - RADIUS */);
                    }
                }
            }
        }

        void
        Membrane::AdjustLimits( const store::Limits3D &LIMITS)
        {
            for( size_t i = 0; i < 3; ++i)
            {
                m_TotalLimits.SetMin( i, LIMITS.GetMin( i));
            }
            for( size_t i = 0; i < 2; ++i)
            {
//                m_LipidLimits.SetMin( i, LIMITS.GetMin( i));
            }
        }



        std::vector< float>
        Membrane::ZLayerLimits() const
        {
            std::vector< float> layers( 5);
            layers[0] = m_TotalLimits.GetZMin();
            layers[1] = m_LipidMin;
            layers[2] = m_Center;
            layers[3] = m_LipidMax;
            layers[4] = m_TotalLimits.GetZMax();
            return layers;
        }


        //! a control for reading pdb files, it repairs e.g. problems due to the limited residue numbers
        bool Membrane::LineChecker
        (
                const boost::shared_ptr< Atom> &ATOM,
                boost::shared_ptr< SimpleMolecule< Atom> > &MOLECULE
        )
        {
//            DebugWrite( __FUNCTION__);
            DebugWrite( __FUNCTION__ << " previous: " << m_PreviousType << " " << m_PreviousID << " new: " << ATOM->GetResidueType() << " " << ATOM->GetResidueID() << " " << ( ( ATOM->GetResidueID() == m_PreviousID) ? "same": "diff"));

            if( m_PreviousID == std::numeric_limits< int>::max())  
	    {
                DebugWrite( __FUNCTION__ << " first line");
                m_PreviousID =     ATOM->GetResidueID();
                m_PreviousType = ATOM->GetResidueType();
                MOLECULE->SetType( m_PreviousType);
            }

            if( ATOM->GetResidueType() == "TIP3")
            {
                static size_t s_WaterAtomCounter( 0);
                if( ++s_WaterAtomCounter % 3 == 1)
                {
                    DebugWrite( __FUNCTION__ << " new water, #water_atoms: " << s_WaterAtomCounter);
                    m_PreviousType = ATOM->GetResidueType();
                    m_PreviousID =     ATOM->GetResidueID();
                    return true;
                }
            }
            else if( ATOM->GetResidueID() != m_PreviousID || ATOM->GetResidueType() != m_PreviousType)
            {
                m_PreviousID =     ATOM->GetResidueID();
                m_PreviousType = ATOM->GetResidueType();
                DebugWrite( __FUNCTION__ << " new molecule: " << m_PreviousType << " " << m_PreviousID);
                return true;
            }

            return false;
        }

        float Membrane::BilayerThickness() const
        {
            return fabs( m_LipidMax - m_LipidMin);
        }

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        void Membrane::AnalyseAndBuild
        (
                const CommandLineManager &COMMAND
        )
        {
            // TODO: merge with code in molecule_factory.cc::BuildImplicitMolecule( cmd, impl_mol)

            // declarations
            std::ofstream
                write;
            std::ifstream
                read;
            std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >
                lipids;
            std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >
                waters;
            std::string
                line,
                residue_type;
            boost::shared_ptr< SimpleMolecule< Atom> >
                molecule( new SimpleMolecule< Atom>());
            boost::shared_ptr< Atom>
                atom;
            std::string
                mode( "average"),
                file;
            std::pair< float, float>
                z_limits;
            math::Vector3N
                cms;
            float
                z_cms = 0.0,
                total_z_cms = 0.0,
                lipid_min(  std::numeric_limits< float>::max()),
                lipid_max( -std::numeric_limits< float>::max()),
                x_min(  std::numeric_limits< float>::max()),
                x_max( -std::numeric_limits< float>::max()),
                y_min(  std::numeric_limits< float>::max()),
                y_max( -std::numeric_limits< float>::max()),
                z_min(  std::numeric_limits< float>::max()),
                z_max( -std::numeric_limits< float>::max());


            boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >
                all_molecules( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >);

	    Open( read, COMMAND.GetArgumentStringForFlag( "membrane"));

	    all_molecules = file::ReadMoleculesInGriffinFormat( read);

	    Close( read);


	    DebugWrite( "molecule: " << all_molecules);

            math::Vector
                lipid_z_mins,
                lipid_z_maxs;

            store::Vector3N< std::pair< float, float> >
                limits;

            float delta( 1.0), z_pos;

            DebugWrite( "seperate lipids and waters, box limits");

            //separate molecules and calculate z limits
            for( std::vector< boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator mol_itr = all_molecules->begin(); mol_itr != all_molecules->end(); ++mol_itr)
            {
                limits = ( *mol_itr)->Limits();
                z_limits = limits( 2);
                z_cms = ( *mol_itr)->CMS()( 2);
                if( std::find( m_LipidKeys.begin(), m_LipidKeys.end(), ( *mol_itr)->GetAtoms()( 0)->GetResidueType()) != m_LipidKeys.end())
                {
                    lipids.push_back( Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> >( *mol_itr, z_cms, z_limits));
                    total_z_cms += z_cms;
                    if( z_limits.first < lipid_min)
                    {
                        lipid_min = z_limits.first;
                    }
                    if( z_limits.second > lipid_max)
                    {
                        lipid_max = z_limits.second;
                    }
                }
                else
                {
                    waters.push_back( Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> >( *mol_itr, z_cms, z_limits));
                }

                if( limits( 0).first < x_min)
                {
                    x_min = limits( 0).first;
                }
                if( limits( 0).second > x_max)
                {
                    x_max = limits( 0).second;
                }
                if( limits( 1).first < y_min)
                {
                    y_min = limits( 1).first;
                }
                if( limits( 1).second > y_max)
                {
                    y_max = limits( 1).second;
                }
                if( limits( 2).first < z_min)
                {
                    z_min = limits( 2).first;
                }
                if( limits( 2).second > z_max)
                {
                    z_max = limits( 2).second;
                }
            }
            total_z_cms /= lipids.size();
            std::cout << "\nlimits of bilayer box (determined by min/max atom positions): " << std::endl;
            std::cout << "total-box-x-min:     " << x_min << std::endl;
            std::cout << "total-box-x-max:     " << x_max << std::endl;
            std::cout << "total-box-y-min:     " << y_min << std::endl;
            std::cout << "total-box-y-max:     " << y_max << std::endl;
            std::cout << "total-box-z-min:     " << z_min << std::endl;
            std::cout << "total-box-z-max:     " << z_max << std::endl;
            std::cout << "lipid-box-z-min:     " << lipid_min << std::endl;
            std::cout << "lipid-box-z-center:  " << 0.5 * ( lipid_max + lipid_min) << std::endl;
            std::cout << "lipid-box-z-cms:     " << total_z_cms << std::endl;
            std::cout << "lipid-box-z-max:     " << lipid_max << std::endl << std::endl;

            // insert lipid keys into the lipid layers
            for( std::vector< std::string>::const_iterator itr = m_LipidKeys.begin(); itr != m_LipidKeys.end(); ++itr)
            {
                m_MoleculeLayer[ 1].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
                m_MoleculeLayer[ 2].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
            }

            // insert solution keys into all layers
            for( std::vector< std::string>::const_iterator itr = m_SolKeys.begin(); itr != m_SolKeys.end(); ++itr)
            {
                m_MoleculeLayer[ 0].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
                m_MoleculeLayer[ 1].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
                m_MoleculeLayer[ 2].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
                m_MoleculeLayer[ 3].insert( std::make_pair( *itr, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< mol::SimpleMolecule< mol::Atom> >())));
            }


            DebugWrite( "insert molecules into layers");

            // sort lipids into layers and calculate max values
            for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = lipids.begin(); itr != lipids.end(); ++itr)
            {
                // lower layer
                if( itr->second < total_z_cms)
                {
                    m_MoleculeLayer[ 1][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                    lipid_z_mins.push_back( itr->third.first);
                }
                // upper layer
                else
                {
                    m_MoleculeLayer[ 2][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                    lipid_z_maxs.push_back( itr->third.second);
                }
            }

            DebugWrite( "lipid layer statistics");


            std::pair< float, float>
                min_statistic = lipid_z_mins.RmsdOfElements(),
                max_statistic = lipid_z_maxs.RmsdOfElements();

            std::cout << "\nbilayer z dimensions as determined by statistics of the maximum (minimum) coordinates of every lipid: " << std::endl;
            std::cout << "mean-min: " << min_statistic.first << " - deviation: " << min_statistic.second << " = " << min_statistic.first - min_statistic.second << std::endl;
            std::cout << "mean-max: " << max_statistic.first << " + deviation: " << max_statistic.second << " = " << max_statistic.first + max_statistic.second << std::endl << std::endl;

            // assign limits

            m_TotalLimits.SetXMin( x_min);
            m_TotalLimits.SetXMax( x_max);
            m_TotalLimits.SetYMin( y_min);
            m_TotalLimits.SetYMax( y_max);
            m_TotalLimits.SetZMin( z_min);
            m_TotalLimits.SetZMax( z_max);



            if( COMMAND.IsFlagSet( "auto_dimensions"))
            {
                mode = COMMAND.GetArgumentStringForFlag( "auto_dimensions");
            }

            if( mode == "max")
            {
                m_LipidMin = lipid_min;
                m_LipidMax = lipid_max;
            }
            else if( mode == "average")
            {
                m_LipidMin = min_statistic.first;
                m_LipidMax = max_statistic.first;
            }
            else if( mode == "standard-deviation")
            {
                m_LipidMin = min_statistic.first - min_statistic.second;
                m_LipidMax = max_statistic.first + max_statistic.second;
            }
            m_Center = 0.5 * ( m_LipidMin + m_LipidMax);

            if( COMMAND.IsFlagSet( "dimensions"))
            {
            	StandardWrite( "set bilayer dimensions");

			    std::vector< float> values = COMMAND.GetArgumentsForFlag< float>( "dimensions");

			    float
				    xmin = values[0],
				    xmax = values[1],
				    ymin = values[2],
				    ymax = values[3],
				    zmin = values[4],
				    zmax = values[5],
				    z_thickness = 0.0;

			    if( values.size() == 8)
			    {
			    	m_Center = values[6];
			    	z_thickness = values[7];
			    }
			    else if( values.size() == 7)
			    {
			    	m_Center = 0.5 * ( z_min + z_max);
			    	z_thickness = values[6];
			    }

			    m_TotalLimits = store::Limits3D( xmin, xmax, ymin, ymax, zmin, zmax);
			    m_LipidMin = m_Center - 0.5 * z_thickness;
			    m_LipidMax = m_Center + 0.5 * z_thickness;

			    std::pair< float, float> lim = GetXLimits();
			    std::cout << "\nused limits of bilayer box (after adjustment): " << std::endl;
			    std::cout << "total-box-x-min:     " << lim.first  << std::endl;
			    std::cout << "total-box-x-max:     " << lim.second << std::endl;

			    lim = GetYLimits();
			    std::cout << "total-box-y-min:     " << lim.first  << std::endl;
			    std::cout << "total-box-y-max:     " << lim.second << std::endl;

			    std::vector< float> lims = GetZLimits();
			    std::cout << "total-box-z-min:     " << lims[0] << std::endl;
			    std::cout << "lipid-box-z-min:     " << lims[1] << std::endl;
			    std::cout << "lipid-box-z-center:  " << lims[2] << std::endl;
			    std::cout << "lipid-box-z-max:     " << lims[3] << std::endl;
			    std::cout << "total-box-z-max:     " << lims[4] << std::endl << std::endl;

            }



            for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = waters.begin(); itr != waters.end(); ++itr)
            {
            	DebugWrite( itr->first->GetAtoms()( 0)->GetResidueType());
                // lower layer
                if( itr->second < m_LipidMin)
                {
                    m_MoleculeLayer[ 0][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                }
                // upper layer
                else if( itr->second < m_Center)
                {
                    m_MoleculeLayer[ 1][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                }
                else if( itr->second < m_LipidMax)
                {
                    m_MoleculeLayer[ 2][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                }
                else
                {
                    m_MoleculeLayer[ 3][ itr->first->GetAtoms()( 0)->GetResidueType()]->push_back( itr->first);
                }
            }


            if( COMMAND.IsFlagSet( "write_histograms"))
            {

                StandardWrite( "write histograms");


                std::string prefix = COMMAND.GetArgumentStringForFlag( "write_histograms");

                math::Histogram< float>
                    histogram( m_TotalLimits.GetZMin(), delta, size_t( m_TotalLimits.ZDelta() / delta));

                for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = lipids.begin(); itr != lipids.end(); ++itr)
                {
                    for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = itr->first->GetAtoms().begin(); atom_itr != itr->first->GetAtoms().end(); ++atom_itr)
                    {
                        z_pos = ( *atom_itr)->GetPosition()( 2);
                        histogram.InsertValue( z_pos);
                    }
                }

                WriteObject( write, histogram, prefix + "_lipid_z_histogram.txt");

                StandardWrite( "lipid histogram in z written");

                math::Distribution< float>
    //                atom_nr_histogram( m_TotalLimits.GetZMin(), delta, size_t( m_TotalLimits.GetZDelta() / delta)),
                    mass_histogram( m_TotalLimits.GetZMin(), delta, size_t( m_TotalLimits.ZDelta() / delta));

                for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = lipids.begin(); itr != lipids.end(); ++itr)
                {
                    for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = itr->first->GetAtoms().begin(); atom_itr != itr->first->GetAtoms().end(); ++atom_itr)
                    {
                        z_pos = ( *atom_itr)->GetPosition()( 2);
                        mass_histogram.InsertValue( z_pos, ( *atom_itr)->GetMass());
                    }
                }

                WriteObject( write, mass_histogram, prefix + "_lipid_mass_z_histogram.txt");

                StandardWrite( "lipid mass histogram in z written");


                math::Distribution< float>
                    mass_histogram_x( m_TotalLimits.GetXMin(), delta, size_t( m_TotalLimits.XDelta() / delta)),
                    mass_histogram_y( m_TotalLimits.GetYMin(), delta, size_t( m_TotalLimits.YDelta() / delta)),
                    mass_histogram_z( m_TotalLimits.GetZMin(), delta, size_t( m_TotalLimits.ZDelta() / delta));

                for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = lipids.begin(); itr != lipids.end(); ++itr)
                {
                    for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = itr->first->GetAtoms().begin(); atom_itr != itr->first->GetAtoms().end(); ++atom_itr)
                    {
                        mass_histogram_x.InsertValue( ( *atom_itr)->GetPosition()( 0), ( *atom_itr)->GetMass());
                        mass_histogram_y.InsertValue( ( *atom_itr)->GetPosition()( 1), ( *atom_itr)->GetMass());
                        mass_histogram_z.InsertValue( ( *atom_itr)->GetPosition()( 2), ( *atom_itr)->GetMass());
                    }
                }

                for( std::vector< Triplet< boost::shared_ptr< SimpleMolecule< Atom> >, float, std::pair< float, float> > >::iterator itr = waters.begin(); itr != waters.end(); ++itr)
                {
                    for( std::vector< boost::shared_ptr< Atom> >::const_iterator atom_itr = itr->first->GetAtoms().begin(); atom_itr != itr->first->GetAtoms().end(); ++atom_itr)
                    {
                        mass_histogram_x.InsertValue( ( *atom_itr)->GetPosition()( 0), ( *atom_itr)->GetMass());
                        mass_histogram_y.InsertValue( ( *atom_itr)->GetPosition()( 1), ( *atom_itr)->GetMass());
                        mass_histogram_z.InsertValue( ( *atom_itr)->GetPosition()( 2), ( *atom_itr)->GetMass());
                    }
                }

                std::string name = prefix + "_mass_histogram";

                WriteObject( write, mass_histogram_x, name + "_x.txt");
                WriteObject( write, mass_histogram_y, name + "_y.txt");
                WriteObject( write, mass_histogram_z, name + "_z.txt");


                Open( write, name + "_gnuplot.txt");
                write << "plot \"" << name <<  "_x.txt\" us 1:2 w l, \"" << name << "_y.txt\" us 1:2 w l, \"" << name << "_z.txt\" us 1:2 w l" << std::endl;
                write << "pause(-1)" << std::endl;
                write.close();
                write.clear();

                StandardWrite( "mass histograms in x,y,z written");
            }
        }


        std::istream &Membrane::Read
        (
                std::istream &STREAM,
                const boost::shared_ptr< store::Map< std::string, float> > &MASS_MAP,
                const boost::shared_ptr< store::Map< std::string, std::string> > &TRANSLATOR
        )
        {
            std::string
                str,
                residue_type;
            float
                z_cms,
                vdw;
            boost::shared_ptr< SimpleMolecule< Atom> >
                molecule( new SimpleMolecule< Atom>());
            boost::shared_ptr< Atom>
                atom;
            std::multimap< float, boost::shared_ptr< SimpleMolecule< Atom> > >
                collector;

//            bool
//                carry_on( true);
//
//            // find first atom line
//            while( carry_on)
//            {
//                std::getline( STREAM, str);
//                if( str.length() > 6 && str.substr( 0, 4) == "ATOM")
//                {
//                    atom = AtomFileHandler< Atom>().ReadFromPdbLine( str);
//                    atom->SetMass( MASS_MAP->operator()( TRANSLATOR->operator()( atom->GetAtomName())));
//
//                    atom->SetResidueID( molecule_count);
//                    atom->SetAtomID( ++atom_count);
//
//                    DebugWrite( atom);
//
//                    LineChecker( atom);  // for correct internal numbering within LineChecker
//
//
////                    previous_id = mystr::ConvertStringToNumericalValue< size_t>( str.substr( 22, 4));
////                    previous_type = atom->GetResidueType();
//
//                    vdw = atom->GetVanDerWaalsRadius();
//
//                    molecule->SetType( previous_type);
//
//                    DebugWrite( "mol-type altered: " << molecule->GetType());
//
//                    // if residue type is one of the expected lipid types, adjust limits
//                    AdjustLimits( previous_type, atom->GetPosition(), vdw);
//                    molecule->AddAtom( atom);
//                    carry_on = false;
//                }
//            }

            // collecting and fitting limits
            while( STREAM)
            {
                std::getline( STREAM, str);
                if( str.length() > 6 && str.substr( 0, 4) == "ATOM")
                {
                    atom = AtomFileHandler< Atom>().ReadFromPdbLine( str);
                    atom->SetMass( MASS_MAP->operator()( TRANSLATOR->operator()( atom->GetAtomName())));

//                    residue_id = mystr::ConvertStringToNumericalValue< size_t>( str.substr( 22, 4));
                    residue_type = atom->GetResidueType();
                    vdw = atom->GetVanDerWaalsRadius();

//                    atom->SetResidueID( molecule_count);
//                    atom->SetAtomID( ++atom_count);

                    DebugWrite( atom);

                    // if residue type is one of the expected lipid types, adjust limits
                    AdjustLimits( residue_type, atom->GetPosition()/*, vdw*/);
                    // build collector
                    // if new molecule, add old to collector and reset molecule
                    if( LineChecker( atom, molecule)) // note: LineChecker has to be first!
                    {
//                        previous_id = residue_id;
//                        previous_type = residue_type;

//                        ++molecule_count;

                        z_cms = molecule->CMS()( 2);
//                        DebugWrite( "cms-z: " << z_cms << " mol-type: " << residue_type);
                        if( z_cms != std::numeric_limits< float>::max())
                        {
                            collector.insert( std::make_pair( z_cms, molecule));
                        }
                        molecule = boost::shared_ptr< SimpleMolecule< Atom> >( new SimpleMolecule< Atom>( residue_type));
                    }
                    // add atom to molecule
                    molecule->AddAtom( atom);
                }
            }

            DebugWrite( "atoms read, now reading into layers");

            if( m_MoleculeLayer.size() == 0)
            {
                DebugWrite( "add molecule layers");
                for( size_t i = 0; i < 4; ++i)
                {
                    m_MoleculeLayer.push_back( store::Map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > > >());
                }
            }

            // read collector into layers, according to limits
            for( std::multimap< float, boost::shared_ptr< SimpleMolecule< Atom> > >::const_iterator itr = collector.begin(); itr != collector.end(); ++itr)
            {
                size_t
                    layer_id( LayerID( itr->first));
                if( layer_id > 3)
                {
                    std::cout << "===> LAYER NOT IDENTIFIED" << std::endl;
                    exit( -1);
                }
                std::string
                    molecule_type( itr->second->GetType());

                DebugWrite( "layer: " << layer_id << " molecule type: <"<< molecule_type << "> z: " << itr->first);


                if( !m_MoleculeLayer[ layer_id].IsValidKey( molecule_type))
                {
                    DebugWrite( "new molecule type" << molecule_type);
                    m_MoleculeLayer[ layer_id].InsertNewKeyAndValue( molecule_type, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< mol::Atom> > >( new store::ShPtrVec< SimpleMolecule< Atom> >()));
                }
//                std::cout << "push back: " << itr->second;
//                std::cout << "cms: " << itr->second->CMS() << std::endl;
                m_MoleculeLayer[ layer_id][ molecule_type]->push_back( itr->second);
            }

            DebugWrite( "atoms read into layers");
            return STREAM;
        }



        //                        //            COLUMNS        DATA  TYPE    FIELD        DEFINITION
        //                        //            -------------------------------------------------------------------------------------
        //                        //             1 -  6        Record name   "ATOM  "
        //                        //             7 - 11        Integer       serial       Atom  serial number.
        //                        //            13 - 16        Atom          name         Atom name.
        //                        //            17             Character     altLoc       Alternate location indicator.
        //                        //            18 - 20        Residue name  resName      Residue name.
        //                        //            22             Character     chainID      Chain identifier.
        //                        //            23 - 26        Integer       resSeq       Residue sequence number.
        //                        //            27             AChar         iCode        Code for insertion of residues.
        //                        //            31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        //                        //            39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        //                        //            47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        //                        //            55 - 60        Real(6.2)     occupancy    Occupancy.
        //                        //            61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        //                        //            77 - 78        LString(2)    element      Element symbol, right-justified.
        //                        //            79 - 80        LString(2)    charge       Charge  on the atom.



        std::ostream &Membrane::Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << std::endl;
            STREAM << m_TotalLimits,
            STREAM << m_MoleculeLayer;
            return STREAM;
        }

        std::ostream &Membrane::WriteLimits( std::ostream &STREAM) const
        {
            STREAM << "box-limits: " << m_TotalLimits;
            STREAM << "lipid-limits: " << m_LipidMin << " " << m_LipidMax << std::endl;
            STREAM << "center: " << m_Center << std::endl;
            return STREAM;
        }

        std::ostream &Membrane::Write( std::ostream &STREAM, const CommandLineManager &COMMAND) const
        {
            WriteLimits( STREAM);
            return STREAM;
        }

        std::ostream &Membrane::WriteLayerOverview( std::ostream &STREAM) const
        {
            STREAM << "total box sizes: " << std::endl;
            STREAM << m_TotalLimits;
            store::Limits3D limits;
            for( size_t i = 0; i < 4; ++i)
            {
                limits = FourLayerBoxLimits( i);
                for( std::map< std::string, boost::shared_ptr< store::ShPtrVec< mol::SimpleMolecule< Atom> > > >::const_iterator itr = m_MoleculeLayer[ i].begin(); itr != m_MoleculeLayer[ i].end(); ++itr)
                {
                    STREAM << "layer: " << i << " z_min: " << limits.GetZMin() << " z_max: " << limits.GetZMax() << " molecule_type: " << itr->first << " nr_atoms: " << itr->second->size() << std::endl;
                }
            }
            return STREAM;
        }


        std::string Membrane::GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }

} // end namespace mol


