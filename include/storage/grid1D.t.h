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


#ifndef GRID1D_T_H
#define GRID1D_T_H

#include <vector>
#include <iterator>

#include "position_grid.t.h"
#include "atom_grid_factory.t.h"
#include "typemapped_gridpoint.t.h"
#include "recursive_typemapped_gridpoint.t.h"
#include "default_unique_map.t.h"
#include "../utilities/energy_enum.h"

namespace store
{


    template< typename t_RETURN>
    class Grid1D
    : public PositionGrid< t_RETURN>
    {
    protected:
    	std::vector< std::pair< util::FunctorEnum, int> >  		              	   	    m_GridPoint;
    	std::vector< std::vector< std::pair< util::FunctorEnum, int> > >                m_RecursiveTypeMapped;
    	std::vector< std::vector< int> >	  							        	    m_TypeMapped;
    	std::vector< std::vector< std::pair< util::EnergyEnum, int> > >  	            m_Interaction;
    	std::vector< math::Vector3N> 										            m_Force;
    	std::vector< float> 												            m_Energy;
    	std::vector< int>													            m_NSP;
    	store::DefaultUniqueMap< std::string, int>	                                                m_MolTypes;
    	store::Map< std::string, int>	                                                m_EnergyTypes;

    public:
    	Grid1D()
    	: PositionGrid< t_RETURN>(),
		m_GridPoint(),
		m_RecursiveTypeMapped(),
		m_TypeMapped(),
		m_Interaction(),
		m_Force(),
		m_Energy(),
		m_NSP(),
		m_MolTypes( "all", 0),
		m_EnergyTypes()
		{}


    	Grid1D
    	(
    			const math::Vector3N &MIN,
    			math::Vector3N &MAX,
    			const math::Vector3N &DELTA,
    			const store::Vector3N< size_t> &NR_BINS,
    			const store::DefaultUniqueMap< std::string, int> &MOL_TYPES,// = std::map< std::string, int>(),
    			const store::Map< std::string, int> &ENERGY_TYPES// = std::map< std::string, int>()
    	)
    	: PositionGrid< t_RETURN>( MIN, MAX, DELTA, NR_BINS),
		m_GridPoint(),
		m_RecursiveTypeMapped(),
		m_TypeMapped(),
		m_Interaction(),
		m_Force(),
		m_Energy(),
		m_NSP(),
		m_MolTypes( MOL_TYPES),
		m_EnergyTypes( ENERGY_TYPES)
		{
		    Allocate();
		}


    	Grid1D( const Grid1D &ORIG)
    	: PositionGrid< t_RETURN>( ORIG),
		m_GridPoint( ORIG.m_GridPoint),
		m_RecursiveTypeMapped( ORIG.m_RecursiveTypeMapped),
		m_TypeMapped( ORIG.m_TypeMapped),
		m_Interaction( ORIG.m_Interaction),
		m_Force( ORIG.m_Force),
		m_Energy( ORIG.m_Energy),
		m_NSP( ORIG.m_NSP),
		m_MolTypes( ORIG.m_MolTypes),
		m_EnergyTypes( ORIG.m_EnergyTypes)
		{}

    	virtual ~Grid1D(){}

    	virtual Grid1D * Clone() const{ return new Grid1D( *this);}


		void Allocate()
		{
			int size = PositionGrid< t_RETURN>::TotalNrBins();
				m_GridPoint.resize( size);
				m_Force.reserve( size); // try to allocate expected memory, estimate memory for each element and assign
				m_Energy.reserve( size);
				m_NSP.reserve( size);
		}

		const PositionGrid< t_RETURN> &
        GetPositionGrid() const
        {
        	return *this;
        }

        PositionGrid< t_RETURN> &
        GetSetPositionGrid()
        {
        	return *this;
        }


       void SetPositionGrid( const PositionGrid< t_RETURN> &POSGRID)
       {
    	   PositionGrid< t_RETURN>::operator = (POSGRID);
       }


        void SetMolTypes( const store::DefaultUniqueMap< std::string, int> &MOLTYPES)
    	{
    		m_MolTypes = MOLTYPES;
    	}

    	const store::DefaultUniqueMap< std::string, int> &GetMolTypes() const
    	{
    		return m_MolTypes;
    	}

    	void SetEnergyTypes( const store::Map< std::string, int> &ENERGIES)
    	{
    		m_EnergyTypes = ENERGIES;
    	}

    	const store::Map< std::string, int> &GetEnergyTypes() const
    	{
    		return m_EnergyTypes;
    	}

#ifdef FORCE_INDICES
        virtual
        t_RETURN
        GetGridPoint
	    ( 
		const math::Vector3N &POSITION, 
		const std::string &MOL_TYPE, 
#ifdef CPPIO
		std::ostream &STREAM
#else
		FILE *STREAM
#endif
	    ) const 
//        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE, std::ostream &STREAM) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION), MOL_TYPE, STREAM);
        }
#endif


        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION, const std::string &MOL_TYPE) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << " moltype: " << MOL_TYPE);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION), MOL_TYPE);
        }




        virtual t_RETURN 
	GetGridPoint
	( 
	    const store::Vector3N< int> &IDS, 
	    const std::string &MOL_TYPE
	) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")");
            int id = PositionGrid< t_RETURN>::IDFromIDs( IDS);

            if( id == g_IntMax)
            {
            	return t_RETURN(  new VoidGridPoint< mol::Atom, math::Vector3N>());
            }
            return ConstructGridPoint( m_GridPoint[ id], MOL_TYPE);
        }

#ifdef FORCE_INDICES
        virtual t_RETURN 
	GetGridPoint
	( 
	    const store::Vector3N< int> &IDS, 
	    const std::string &MOL_TYPE,
#ifdef CPPIO
		std::ostream &STREAM
#else
		FILE *STREAM
#endif
	) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")");
            int id = PositionGrid< t_RETURN>::IDFromIDs( IDS);

#ifdef CPPIO
	    STREAM << << "  " << IDS[0] << "  " << IDS[1] << "  " << IDS[2] << "    " << id << "    ";
#else
	    fprintf( STREAM, "  %d  %d %d  %d  ", IDS[0], IDS[1], IDS[2], id);
#endif

            if( id == g_IntMax)
            {
            	return t_RETURN(  new VoidGridPoint< mol::Atom, math::Vector3N>());
            }
            return ConstructGridPoint( m_GridPoint[ id], MOL_TYPE);
        }
#endif


        virtual t_RETURN ConstructGridPoint( const std::pair< util::FunctorEnum, int> &INFO, const std::string &MOL_TYPE) const
        {
            int locator = INFO.second;
            switch( INFO.first)
            {
            case util::e_VoidGridPoint:
				return t_RETURN( new VoidGridPoint< mol::Atom, math::Vector3N>());
            case util::e_ConstantGridPoint:
				return t_RETURN( new ConstantGridPoint< mol::Atom, math::Vector3N>( m_Force[ locator], m_NSP[ locator], m_Energy[ locator]));
            case util::e_InteractionGridPoint:
            {
				InteractionGridPoint inter;
				for( std::vector< std::pair< util::EnergyEnum, int> >::const_iterator itr = m_Interaction[locator].begin(); itr != m_Interaction[locator].end(); ++itr)
				{
					switch( itr->first)
					{
					case util::e_Coulomb:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "coulomb", boost::shared_ptr< phys::PotentialElectrostaticForce>( new phys::PotentialElectrostaticForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					case util::e_AttractiveVdw:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-attractive", boost::shared_ptr< phys::PotentialAttractiveVanDerWaalsForce>( new phys::PotentialAttractiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					case util::e_RepulsiveVdw:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-repulsive", boost::shared_ptr< phys::PotentialRepulsiveVanDerWaalsForce>( new phys::PotentialRepulsiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					default: std::cout << "====> " << __FUNCTION__ << " undefined type here " << itr->first << " " << itr-> second << std::endl;
					return t_RETURN();
					}
				}
				return t_RETURN( new InteractionGridPoint( inter));
            }
            case util::e_TypeMappedGridPoint:
            {
             	int const_locator = m_TypeMapped[ locator][ m_MolTypes( MOL_TYPE)];
				return t_RETURN( new ConstantGridPoint< mol::Atom, math::Vector3N>( m_Force[ const_locator], m_NSP[ const_locator], m_Energy[ const_locator]));
            }
            case util::e_RecursiveTypeMappedGridPoint:
            {
				return ConstructGridPoint( m_RecursiveTypeMapped[ locator][ m_MolTypes( MOL_TYPE)], MOL_TYPE);
            }
            default: std::cerr << "not possible to construct grid point" << std::endl; exit( -1);
            }
        }

        virtual
        util::FunctorEnum
        GetGridPointType( const math::Vector3N &POSITION, const std::string &MOL_TYPE)
        {
        	int id = PositionGrid< t_RETURN>::IDFromPosition( POSITION);
        	if( id == g_IntMax)
        	{
        		return util::e_VoidGridPoint;
        	}
        	// check for gp type and access moltype if necessary
        	std::pair< util::FunctorEnum, int> pair = m_GridPoint[ id];
        	switch( pair.first)
        	{
        	case util::e_TypeMappedGridPoint:
        	{
        		return util::e_ConstantGridPoint;
        	}
        	break;
        	case util::e_RecursiveTypeMappedGridPoint:
        	{
        		return m_RecursiveTypeMapped[ pair.second][ m_MolTypes( MOL_TYPE)].first;
        	}
        	break;
        	default: return pair.first;
        	}
        }





        virtual t_RETURN GetGridPoint( const math::Vector3N &POSITION) const
        {
            DebugWrite( __PRETTY_FUNCTION__);
            return GetGridPoint( PositionGrid< t_RETURN>::IDsFromPosition( POSITION));
        }


        virtual t_RETURN GetGridPoint( const store::Vector3N< int> &IDS) const
        {
            DebugWrite( __PRETTY_FUNCTION__ << "   (" << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << ")");
	

            int id = PositionGrid< t_RETURN>::IDFromIDs( IDS);

//	    std::cout << IDS(0) << " " << IDS( 1) << " " << IDS( 2) << " => " << id <<  std::endl;
	
            if( id == g_IntMax)
            {
//            	return t_RETURN( factory::BuildNeutralObject< t_RETURN>());
            	return t_RETURN(  new VoidGridPoint< mol::Atom, math::Vector3N>());
            }
            return ConstructGridPoint( m_GridPoint[ id]);
        }


        virtual t_RETURN ConstructGridPoint( const std::pair< util::FunctorEnum, int> &INFO) const
        {
            int locator = INFO.second;
            switch( INFO.first)
            {
            case util::e_VoidGridPoint:
				return t_RETURN( new VoidGridPoint< mol::Atom, math::Vector3N>());
            case util::e_ConstantGridPoint:
				return t_RETURN( new ConstantGridPoint< mol::Atom, math::Vector3N>( m_Force[ locator], m_NSP[ locator], m_Energy[ locator]));
            case util::e_InteractionGridPoint:
            {
				std::vector< std::pair< util::EnergyEnum, int> >::const_iterator 
				    itr = m_Interaction[locator].begin();
				InteractionGridPoint 
				    inter( m_NSP[ itr->second]);

				for( ; itr != m_Interaction[locator].end(); ++itr)
				{
					switch( itr->first)
					{
					case util::e_Coulomb:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "coulomb", boost::shared_ptr< phys::PotentialElectrostaticForce>( new phys::PotentialElectrostaticForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					case util::e_AttractiveVdw:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-attractive", boost::shared_ptr< phys::PotentialAttractiveVanDerWaalsForce>( new phys::PotentialAttractiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					case util::e_RepulsiveVdw:
						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-repulsive", boost::shared_ptr< phys::PotentialRepulsiveVanDerWaalsForce>( new phys::PotentialRepulsiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
						break;
					default:
					{ 
					    std::cout << "====> " << __FUNCTION__ << " undefined type at position: " << locator << ": " << itr->first << " " << itr-> second << std::endl;
					    std::cout << "====> " << __FUNCTION__ << " entire map: ";
					    WriteTable< util::EnergyEnum>( std::cout, m_Interaction[locator]);
					    return t_RETURN();
					}
					}
				}
				return t_RETURN( new InteractionGridPoint( inter));
            }
            case util::e_TypeMappedGridPoint:
            {
            	TypeMappedGridPoint< mol::Atom, math::Vector3N>
					gp;
            	std::vector< int>
					ids = m_TypeMapped[ locator];
            	int
					id;

            	if( ids.size() != m_MolTypes.size())
            	{
            		std::cerr << "====> mol types size " << m_MolTypes.size() << " != " << ids.size() << " typemapped size" << std::endl;
            		std::cerr << "typemapped grid point size moltypes: " << gp.GetMolTypes().size() << std::endl << "ids: ";
            		std::copy( ids.begin(), ids.end(), std::ostream_iterator< int>( std::cout, " "));
            		std::cout << std::endl;
            		exit( -1);
            	}

            	for( std::map< std::string, int>::const_iterator itr = m_MolTypes.begin(); itr != m_MolTypes.end(); ++itr)
            	{
            		id = ids[ itr->second];
            		gp.InsertNewKeyAndValues( itr->first, m_Force[ id], m_NSP[ id], m_Energy[ id]);
            	}

            	return t_RETURN( new TypeMappedGridPoint< mol::Atom, math::Vector3N>( gp));
            }
            case util::e_RecursiveTypeMappedGridPoint:
            {
            	RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N>
					gp;

             	std::vector< std::pair< util::FunctorEnum, int> >
					type_id = m_RecursiveTypeMapped[ locator];
            	boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> >
					function;

               	if( type_id.size() != m_MolTypes.size())
               	{
               		std::cerr << "recursive typemapped grid point size moltypes: " << m_MolTypes.size() << " != " << type_id.size() << " data size" << std::endl;
					std::cerr << "grid1D size moltypes: " << gp.GetMolTypes().size() << std::endl;
					exit( -1);
               	}

            	for( std::map< std::string, int>::const_iterator itr = m_MolTypes.begin(); itr != m_MolTypes.end(); ++itr)
            	{
            		switch( type_id[ itr->second].first)
            		{
                    case util::e_VoidGridPoint:
         				function = t_RETURN( new VoidGridPoint< mol::Atom, math::Vector3N>());
         				break;
            		case util::e_ConstantGridPoint:
            		{
            			int id = type_id[ itr->second].second;
            			function = boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> >
            			(
            					new ConstantGridPoint< mol::Atom, math::Vector3N>
								(
										m_Force[ id],
										m_NSP[ id],
										m_Energy[ id]
								)
            			);
            		}
            		break;
            		case util::e_InteractionGridPoint:
            		{
        				int
							id = type_id[ itr->second].second;
        				std::vector< std::pair< util::EnergyEnum, int> >::const_iterator
        				    itr = m_Interaction[ id].begin();

					InteractionGridPoint 
					    inter( m_NSP[ itr->second]);
					
//					std::cout << "=> NSP of interaction gp inserted in recursive: " << m_NSP[ itr->second] << " id: " << itr->second << std::endl;

        				for( ; itr != m_Interaction[ id].end(); ++itr)
        				{
        					switch( itr->first)
        					{
        					case util::e_Coulomb:
        						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "coulomb", boost::shared_ptr< phys::PotentialElectrostaticForce>( new phys::PotentialElectrostaticForce( m_Force[ itr->second], m_Energy[ itr->second])));
        						break;
        					case util::e_AttractiveVdw:
        						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-attractive", boost::shared_ptr< phys::PotentialAttractiveVanDerWaalsForce>( new phys::PotentialAttractiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
        						break;
        					case util::e_RepulsiveVdw:
        						inter.GetSetPotentialForcesContainer().InsertNewKeyAndValue( "vdw-repulsive", boost::shared_ptr< phys::PotentialRepulsiveVanDerWaalsForce>( new phys::PotentialRepulsiveVanDerWaalsForce( m_Force[ itr->second], m_Energy[ itr->second])));
        						break;
        					default:
        					{
        					    std::cout << "====> " << __FUNCTION__ << " undefined type at position: " << locator << ": " << type_id[ itr->second].first << " " << id << std::endl;
        					    std::cout << "====> " << __FUNCTION__ << " entire map: ";
        					    WriteTable< util::EnergyEnum>( std::cout, m_Interaction[ id]);
        					    return t_RETURN();
        					}
        					}
        				}
            			function = boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> >( new InteractionGridPoint( inter));
            		}
            		break;
            		default:
            		{
            			std::cout << "===> undefined grid point mounted on the recursive type mapped one: " << util::EnumHandler< util::FunctorEnum>().String( type_id[ itr->second].first) << " : " << type_id[ itr->second].second << std::endl;
            			exit( 1);
            		}
            		break;
            		}
            		gp.InsertNewKeyAndValue( itr->first, function);
            	}
            	return t_RETURN( new RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N>( gp));
            }
            default: std::cerr << "not possible to construct grid point" << std::endl; exit( -1);
            }
        }


	virtual void UpdateGridPoint( const store::Vector3N< int> &INDICES, const t_RETURN &DATA)
	{
//	    std::cout << __FUNCTION__ << std::endl;
	        std::pair< util::FunctorEnum, int>
		    pair = m_GridPoint[ PositionGrid< t_RETURN>::IDFromIDs( INDICES)];


		if( DATA->GetClassID() != pair.first)
		{
		    std::cout << "====> types do not math in " << __FUNCTION__ << ": " << util::EnumHandler< util::FunctorEnum>().String( DATA->GetClassID()) << " != " << util::EnumHandler< util::FunctorEnum>().String( pair.first) << std::endl;
		    exit( -1);
		}

        	switch( DATA->GetClassID())
        	{
		    case util::e_ConstantGridPoint:
			UpdateConstantGridPoint( pair.second, DATA);
			break;
		    case util::e_TypeMappedGridPoint:
			UpdateTypeMappedGridPoint( pair.second, DATA);
			break;
		    case util::e_RecursiveTypeMappedGridPoint:
			UpdateRecursiveTypeMappedGridPoint( pair.second, DATA);
			break;
		    default:
//			std::cout << "====> " << __FUNCTION__ << " is only defined for Constant, TypeMapped and Recursive GridPoints so far" << std::endl;
			break;
		}

	}

		
	virtual void UpdateConstantGridPoint( const int &ID, const t_RETURN &DATA)
	{
//	    std::cout << __FUNCTION__ << std::endl;
	    store::ConstantGridPoint< mol::Atom, math::Vector3N > *
		ptr = ( store::ConstantGridPoint< mol::Atom, math::Vector3N > *) DATA.get();
	    m_Force[ ID] = ptr->ReturnValue();
	    m_Energy[ ID] = ptr->Energy();
	}


	virtual void UpdateTypeMappedGridPoint( const int &ID, const t_RETURN &DATA)
	{
//	    std::cout << __FUNCTION__ << std::endl;

	    int
		id;

	    std::vector< int>
		ids = m_TypeMapped[ ID];

	    store::TypeMappedGridPoint< mol::Atom, math::Vector3N > * 
		tptr = ( store::TypeMappedGridPoint< mol::Atom, math::Vector3N > *) DATA.get();

	    if( ids.size() != tptr->TypeMap().size())
	    {
		std::cerr << "====> " <<  __FUNCTION__ << ": number of contained ConstantGridPoints doesn't match: " << ids.size() << " != " << tptr->TypeMap().size() << std::endl;
		exit( -1);
	    }

	    for( store::Map< std::string, store::ConstantGridPoint< mol::Atom, math::Vector3N> >::iterator itr = tptr->TypeMap().begin(); itr != tptr->TypeMap().end(); ++itr)
	    {
		id = ids[ m_MolTypes( itr->first)];
		m_Force[ id] = itr->second.ReturnValue();
		m_Energy[ id] = itr->second.Energy();							
	    }
	}



	virtual void UpdateRecursiveTypeMappedGridPoint( const int &ID, const t_RETURN &DATA)
	{
//	    std::cout << __FUNCTION__ << std::endl;

	    std::pair< util::FunctorEnum, int>
		id_pair;

	    std::vector< std::pair< util::FunctorEnum, int> > 
		ids = m_RecursiveTypeMapped[ ID];

	    store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > * 
		rptr = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N > *) DATA.get();

	    for( store::Map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> > >::iterator itr = rptr->TypeMap().begin(); itr != rptr->TypeMap().end(); ++itr)
	    {
		if( itr->second->GetClassID() == util::e_ConstantGridPoint)
		{
		    id_pair = ids[ m_MolTypes( itr->first)];

		    if( id_pair.first != util::e_ConstantGridPoint)
		    {
			std::cerr << "====> " << __FUNCTION__ << ": grid point types do not match for mol-type: " << itr->first << ": " << id_pair.first << " != " << util::e_ConstantGridPoint << std::endl;
			exit( -1);
		    }

		    store::ConstantGridPoint< mol::Atom, math::Vector3N > *
			cptr = ( store::ConstantGridPoint< mol::Atom, math::Vector3N > *) itr->second.get();
		    m_Force[ id_pair.second] = cptr->ReturnValue();
		    m_Energy[ id_pair.second] = cptr->Energy();
		}
	    }
	}



        virtual void
        SetGridPoint( const math::Vector3N &POSITION, const t_RETURN &DATA)
        {
        	store::Vector3N< int> ids = PositionGrid< t_RETURN>::IDsFromPosition( POSITION);
        	SetGridPoint( ids, DATA);
        }


        virtual void SetGridPoint( const store::Vector3N< int> &INDICES, const t_RETURN &DATA)
        {
        	int
				id = PositionGrid< t_RETURN>::IDFromIDs( INDICES);

        	switch( DATA->GetClassID())
        	{
            case util::e_VoidGridPoint:
            	InsertVoidGridPoint( id);
            	break;
            case util::e_ConstantGridPoint:
            	m_GridPoint[ id] = std::make_pair( util::e_ConstantGridPoint, m_Force.size());
            	InsertConstantGridPoint( /*id,*/ DATA);
            	break;
            case util::e_InteractionGridPoint:
            	m_GridPoint[ id] = std::make_pair( util::e_InteractionGridPoint,  m_Interaction.size());
            	InsertInteractionGridPoint( /*id,*/ DATA);
            	break;
            case util::e_TypeMappedGridPoint:
            	InsertTypeMappedGridPoint( id, DATA);
            	break;
            case util::e_RecursiveTypeMappedGridPoint:
            	InsertRecursiveTypeMappedGridPoint( id, DATA);
            	break;
            default:
            	std::cerr << "====> undefined grid point type in " << __PRETTY_FUNCTION__ << std::endl;
//            	exit( -1);
            	break;
        	}
        }

	

        void InsertVoidGridPoint( const int &ID)
        {
        	m_GridPoint[ ID] = std::make_pair( util::e_VoidGridPoint, 0);
        }


        void InsertConstantGridPoint( /*const int &ID,*/ const t_RETURN &DATA)
        {
//	    std::cout <<  __FUNCTION__ << std::endl;
//        	m_GridPoint[ ID] = std::make_pair( util::e_ConstantGridPoint, m_Force.size());
        	ConstantGridPoint< mol::Atom, math::Vector3N>
				*gp = ( ConstantGridPoint< mol::Atom, math::Vector3N> *) DATA.get();
        	m_Force.push_back( gp->GetReturnValue());
        	m_Energy.push_back( gp->GetEnergy());
        	m_NSP.push_back( DATA->GetNSP());
        }


        void InsertInteractionGridPoint( /*const int &ID,*/ const t_RETURN &DATA)
        {
//	 	   std::cout <<  __FUNCTION__ << std::endl;

        	std::vector< std::pair< util::EnergyEnum, int> >
				list;

        	// avoiding copy businness when pushing back
        	list.reserve( m_EnergyTypes.size());

//        	m_GridPoint[ ID] = std::make_pair( util::e_InteractionGridPoint,  m_Interaction.size());

        	store::InteractionGridPoint
				*gp = ( store::InteractionGridPoint *) DATA.get();

        	phys::PotentialForceContainer
				forces = gp->GetSetPotentialForcesContainer();

        	int
				force_id = m_Force.size();

        	for( std::map< std::string, boost::shared_ptr< phys::PotentialForce> >::const_iterator itr = forces.begin(); itr != forces.end(); ++itr, ++force_id)
        	{
        		if( m_EnergyTypes.find( itr->first) == m_EnergyTypes.end())
        		{
        			int i = -1;
        			for( std::map< std::string, int>::const_iterator ftr = m_EnergyTypes.begin(); ftr != m_EnergyTypes.end(); ++ftr)
        			{
        				if( ftr->second > i)
        				{
        					i = ftr->second;
        				}
        			}
					m_EnergyTypes.InsertNewKeyAndValue( itr->first, ++i);
				std::cout << "\n==> add " << itr->first << " " << i << std::endl;
				}

				list.push_back
				(
					std::make_pair
					(
						util::EnumHandler< util::EnergyEnum>().ID( itr->first),
						force_id
					)
				);

				m_Force.push_back( itr->second->GetSetVector());
				m_Energy.push_back( itr->second->GetSetEnergy());
				m_NSP.push_back( gp->GetNSP());
        		

        	}

        	m_Interaction.push_back( list);
        }



    	void InsertTypeMappedGridPoint( const int &ID, const t_RETURN &DATA)
        {
//	    std::cout <<  __FUNCTION__ << std::endl;
        	std::vector< int>
				list( m_MolTypes.size(), g_IntMax);

        	m_GridPoint[ ID] = std::make_pair( util::e_TypeMappedGridPoint,  m_TypeMapped.size());

        	store::TypeMappedGridPoint< mol::Atom, math::Vector3N>
				*gp = ( store::TypeMappedGridPoint< mol::Atom, math::Vector3N> *) DATA.get();

        	store::Map< std::string, ConstantGridPoint< mol::Atom, math::Vector3N> >
				type_mapped = gp->TypeMap();

        	int
				force_id = m_Force.size();

        	for( std::map< std::string, ConstantGridPoint< mol::Atom, math::Vector3N> >::const_iterator itr = type_mapped.begin(); itr != type_mapped.end(); ++itr, ++force_id)
        	{
				list[ m_MolTypes( itr->first)] = force_id;
				m_Force.push_back( itr->second.GetReturnValue());
				m_Energy.push_back( itr->second.GetEnergy());
				m_NSP.push_back( itr->second.GetNSP());
        	}

        	m_TypeMapped.push_back( list);
        }


    	void InsertRecursiveTypeMappedGridPoint( const int &ID, const t_RETURN &DATA)
        {
//	    std::cout <<  __FUNCTION__ << std::endl;
    		std::vector< std::pair< util::FunctorEnum, int> >
				list( m_MolTypes.size(), std::make_pair( util::e_UNDEFINED_Functor, g_IntMax));

    		//list.reserve( m_MolTypes.size());

    		m_GridPoint[ ID] = std::make_pair( util::e_RecursiveTypeMappedGridPoint,  m_RecursiveTypeMapped.size());

        	store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N>
				*gp = ( store::RecursiveTypeMappedGridPoint< mol::Atom, math::Vector3N> *) DATA.get();

        	int
		    id;


        	for( std::map< std::string, boost::shared_ptr< math::Function< mol::Atom, math::Vector3N> > >::const_iterator itr = gp->TypeMap().begin(); itr != gp->TypeMap().end(); ++itr)
        	{
		        id = m_MolTypes( itr->first);

        		switch( itr->second->GetClassID())
        		{
        		case util::e_ConstantGridPoint:
            		list[id] = std::make_pair( util::e_ConstantGridPoint, m_Force.size());
                	InsertConstantGridPoint( itr->second);
                	break;
                case util::e_InteractionGridPoint:
            		list[ id] = std::make_pair( util::e_InteractionGridPoint, m_Interaction.size());
                	InsertInteractionGridPoint( itr->second);
                	break;
                case util::e_VoidGridPoint: break;
                default:
                	std::cerr << "====> only constant or interaction grid point should be mounted to the recursive type mapped (got:" << itr->first << ")" << std::endl;
                	std::cerr << "====> undefined grid point type in " << __PRETTY_FUNCTION__ << std::endl;
    //            	exit( -1);
                	break;
        		}
//        		m_RecursiveTypeMapped.push_back( list);
        	}


    		m_RecursiveTypeMapped.push_back( list);
        }



    	std::ostream &Write( std::ostream& STREAM) const
    	{
	   
    		STREAM << "Grid1D" << std::endl;
    		std::cout<< "Grid1D" << std::endl;
    		PositionGrid< t_RETURN>::Write( STREAM);
     		STREAM << "grid_points:" << std::endl;
     		std::cout<< "grid_points:" << std::endl;
    		WriteTable< util::FunctorEnum>( STREAM, m_GridPoint);
    		STREAM << "recursive_type_mapped:" << std::endl;
    		std::cout<< "recursive_type_mapped:" << std::endl;
    		WriteTableOfTables< util::FunctorEnum>( STREAM, m_RecursiveTypeMapped);
    		STREAM << "type_mapped:" << std::endl;
    		std::cout<< "type_mapped:" << std::endl;
    		STREAM << m_TypeMapped << std::endl;
    		STREAM << "interaction:" << std::endl;
    		std::cout<< "interaction:" << std::endl;
    		WriteTableOfTables< util::EnergyEnum>( STREAM, m_Interaction);
    		STREAM << "force:" << std::endl;
    		std::cout<< "force:" << std::endl;
    		STREAM << m_Force << std::endl;
    		STREAM << "energy:" << std::endl;
    		std::cout<< "energy:" << std::endl;
    		STREAM << m_Energy << std::endl;
    		STREAM << "nsp:" << std::endl;
    		std::cout<< "nsp:" << std::endl;
    		STREAM << m_NSP << std::endl;
    		STREAM << "mol_types:" << std::endl;
    		std::cout<< "mol_types:" << std::endl;
    		STREAM << m_MolTypes << std::endl;
    		STREAM << "energy_types:" << std::endl;
    		std::cout<< "energy_types:" << std::endl;
    		STREAM << m_EnergyTypes << std::endl;
    		return STREAM;
    	}


    	template< typename t_TYPE>
	    std::ostream &WriteTableOfTables( std::ostream &STREAM,  const std::vector< std::vector< std::pair< t_TYPE, int> > > &TABLE) const
	    {
			STREAM << "TableOfTables" << "  ";
			STREAM << TABLE.size() << std::endl;

			for( typename std::vector< std::vector< std::pair< t_TYPE, int> > >::const_iterator itr = TABLE.begin(); itr != TABLE.end(); ++itr)
			{
				WriteTable< t_TYPE>( STREAM, *itr);
			}
			STREAM << std::endl;
			return STREAM;
	    }

    	template< typename t_TYPE>
	    std::ostream &WriteTable( std::ostream &STREAM,  const std::vector< std::pair< t_TYPE, int> > &TABLE) const
	    {
			STREAM << "Table" << "  ";
			STREAM << TABLE.size() << std::endl;

			for( typename std::vector< std::pair< t_TYPE, int> >::const_iterator line_itr = TABLE.begin(); line_itr != TABLE.end(); ++line_itr)
			{
//				    STREAM /* << "std::pair "*/ <<  line_itr->first << "  " << line_itr->second << std::endl;
				STREAM /* << "std::pair "*/ <<  util::EnumHandler< t_TYPE>().String( line_itr->first) << "  " << line_itr->second << std::endl;
			}
			return STREAM;
	    }

    	std::istream &Read( std::istream &STREAM)
    	{

    		std::string
				str;
//   		STREAM >> str;
//
//		std::cout << "read grid point table: " << std::endl;
//
//   		if( str != "Grid1D")
//   		{
//   			std::cout << "====> not the correct id, expected: Grid1D, got: " << str << std::endl;
//   		}

        	PositionGrid< t_RETURN>::Read( STREAM);


    		STREAM >> str;
    		if( str != "grid_points:")
    		{
    			std::cout << "====> not the correct id, expected: grid_points, got: " << str << std::endl;
    		}
    		m_GridPoint = ReadTable< util::FunctorEnum>( STREAM);

    		std::cout << "read recursive type mapped table" << std::endl;
		
    		STREAM >> str;
    		if( str != "recursive_type_mapped:")
    		{
    			std::cout << "====> not the correct id, expected: recursive_type_mapped, got: " << str << std::endl;
    		}
    		m_RecursiveTypeMapped = ReadTableOfTables< util::FunctorEnum>( STREAM);

    		std::cout << "read type mapped table" << std::endl;

    		STREAM >> str;
    		if( str != "type_mapped:")
    		{
    			std::cout << "====> not the correct id, expected: type_mapped, got: " << str << std::endl;
    		}
    		STREAM >> m_TypeMapped;

    		std::cout << "read interaction table" << std::endl;

    		STREAM >> str;
    		if( str != "interaction:")
    		{
    			std::cout << "====> not the correct id, expected: interaction, got: " << str << std::endl;
    		}
    		m_Interaction = ReadTableOfTables< util::EnergyEnum>( STREAM);

    		std::cout << "read force table" << std::endl;

     		STREAM >> str;
    		if( str != "force:")
    		{
    			std::cout << "====> not the correct id, expected: force, got: " << str << std::endl;
    		}
    		STREAM >> m_Force;

    		std::cout << "read energy table" << std::endl;

    		STREAM >> str;
    		if( str != "energy:")
    		{
    			std::cout << "====> not the correct id, expected: energy, got: " << str << std::endl;
    		}
    		STREAM >> m_Energy;

    		std::cout << "read nsp table" << std::endl;

    		STREAM >> str;
    		if( str != "nsp:")
    		{
    			std::cout << "====> not the correct id, expected: nsp, got: " << str << std::endl;
    		}
    		STREAM >> m_NSP;

    		std::cout << "read mol types" << std::endl;

    		STREAM >> str;
    		if( str != "mol_types:")
    		{
    			std::cout << "====> not the correct id, expected: mol_types, got: " << str << std::endl;
    		}
    		STREAM >> m_MolTypes;

    		std::cout << "moltypes: " << m_MolTypes;

    		std::cout << "read energy types" << std::endl;

   		STREAM >> str;
    		if( str != "energy_types:")
    		{
    			std::cout << "====> not the correct id, expected: energy_types, got: " << str << std::endl;
    		}
    		STREAM >> m_EnergyTypes;

    		std::cout << "read grid1D done" << std::endl;

    		return STREAM;
    	}


        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }


        ////////////////////////////////////////////
        //                                        //
        //   BLOCK OF PRIVATE READING FUNCTIONS   //
        //                                        //
        ////////////////////////////////////////////

    private:
        template< typename t_TYPE>
        std::vector< std::vector< std::pair< t_TYPE, int> > >
        ReadTableOfTables( std::istream &STREAM)
    	{

    		std::string
				str;
    		int
				size;
    		std::vector< std::vector< std::pair< t_TYPE, int> > >
				table;
       		std::vector< std::pair< t_TYPE, int> >
				sub_table;

    		m_Interaction.clear();

    		// do output directly here to avoid to call operator << pair< enum, X>
    		STREAM >> str;
    		if( str != "TableOfTables")
    		{
    			std::cout << "====> not the correct id, expected std::vector, got: " << str << std::endl;
    		}
    		STREAM >> size;
    		if( size > 0)
    		{
    			sub_table = ReadTable< t_TYPE>( STREAM);
    			table.resize( size, std::vector< std::pair< t_TYPE, int> >( sub_table.size()));
    			table[0] = sub_table;

    			for( typename std::vector< std::vector< std::pair< t_TYPE, int> > >::iterator itr = table.begin() + 1; itr != table.end(); ++itr)
    			{
    				*itr = ReadTable< t_TYPE>( STREAM);
    			}
    		}
    		return table;
    	}



        template< typename t_TYPE>
    	std::vector< std::pair< t_TYPE, int> >
    	ReadTable( std::istream &STREAM)
    	{

    		std::string
				str;
    		int
				value,
				size;
    		std::vector< std::pair< t_TYPE, int> >
				table;

    		STREAM >> str;
    		if( str != "Table")
    		{
    			std::cout << "====> not the correct id, expected std::vector, got: " << str << std::endl;
    		}
    		STREAM >> size;
    		if( size > 0)
    		{
    			table.resize( size);

    			typename std::vector< std::pair< t_TYPE, int> >::iterator itr = table.begin();
				for( int i = 0; i < size; ++i, ++itr)
				{
//					STREAM >> str;
//					if( str != "std::pair")
//					{
//						std::cout << "====> not the correct id, expected std::pair, got: " << str << std::endl;
//					}
					STREAM >> str >> value;
					*itr = std::make_pair( util::EnumHandler< t_TYPE>().ID( str), value);
				}
    		}

    		return table;
    	}



    }; // end class

} // end namespace

#endif /* GRID1D_T_H */
