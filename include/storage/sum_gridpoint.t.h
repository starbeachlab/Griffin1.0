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


#ifndef SUM_GRIDPOINT_T_H_
#define SUM_GRIDPOINT_T_H_


namespace store
{

    template< typename t_INPUT, typename t_RETURN>
    class SumGridPoint
    : public math::Function< t_INPUT, t_RETURN>
    {
    protected:
        std::vector< std::pair< float, boost::shared_ptr< math::Function< t_INPUT, t_RETURN > > > >   m_WeightsAndGridPoints;
        util::FunctorEnum   							      m_Type;
        int   									      m_NSP;       //!< nearest surface point id

    public:
        SumGridPoint
        (
        		const util::FunctorEnum &TYPE,
				const int &NSP
        )
        : m_WeightsAndGridPoints(),
          m_Type( TYPE),
          m_NSP( NSP)
        {}

        SumGridPoint( const SumGridPoint &SGP)
        : m_WeightsAndGridPoints( SGP.m_WeightAndGridPoints),
          m_Type( SGP.m_Type),
          m_NSP( SGP.m_NSP)
        {}

        virtual ~SumGridPoint(){}


        void Insert( const float &WEIGHT, const boost::shared_ptr< math::Function< t_INPUT, t_RETURN > > &GP)
        {
        	m_WeightsAndGridPoints.push_back( std::make_pair( WEIGHT, GP));
        }

        virtual const int &GetNSP() const
        {
            return m_NSP;
        }


        virtual float Energy( const mol::Atom &ATOM) const
        {
        	DebugWrite( __PRETTY_FUNCTION__);
        	float energy = 0;
        	for( typename std::vector< std::pair< float, boost::shared_ptr< math::Function< t_INPUT, t_RETURN > > > >::const_iterator itr = m_WeightsAndGridPoints.begin(); itr != m_WeightsAndGridPoints.end(); ++itr)
        	{ 
        		energy += itr->first * itr->second->Energy( ATOM);
        	}
        	return energy;
        }

        virtual t_RETURN  & operator()( const t_INPUT &INPUT) const
	{
        	DebugWrite( __PRETTY_FUNCTION__);
        	t_RETURN result;
		int i = 0;
        	for( typename std::vector< std::pair< float, boost::shared_ptr< math::Function< t_INPUT, t_RETURN > > > >::const_iterator itr = m_WeightsAndGridPoints.begin(); itr != m_WeightsAndGridPoints.end(); ++itr, ++i)
         	{
//		    std::cout << i << ". gp" << std::endl;
//		    std::cout <<  i << ":  " << itr->first << std::endl;
//		    std::cout <<  i << ":  " << itr->first << " *  " <<  *itr->second << std::endl;
//		    std::cout <<  i << ":  " << itr->first << " *  " <<  itr->second->GetClassName() << std::endl;
//		    std::cout << i << " now: " << std::endl;
		    result += itr->first * itr->second->operator()( INPUT);
		}
//		std::cout << "sumGP: " <<  result[0] << " " << result[1] << "  " << result[2] << std::endl;
        	return math::Function< t_INPUT, t_RETURN>::s_Tmp = result;
         }
 

         virtual std::ostream & Write( std::ostream &STREAM) const
         {
             STREAM << "SumGridPoint: core-type: " << util::EnumHandler< util::FunctorEnum>().String( m_Type) << " nsp: " << m_NSP << std::endl;
	     int i = 0;
	     for( typename std::vector< std::pair< float, boost::shared_ptr< math::Function< t_INPUT, t_RETURN > > > >::const_iterator itr = m_WeightsAndGridPoints.begin(); itr != m_WeightsAndGridPoints.end(); ++itr, ++i)
	     {
//		 STREAM <<  i << " <=> " << itr->first << "  " << itr->second->GetClassName() << std::endl;
		 STREAM << itr->first << "  " << *itr->second << std::endl;
	     }
             return STREAM;
         }


         virtual util::FunctorEnum GetClassID() const
         {
         	return m_Type;
         }


    }; // end class
} // end namespace

#endif /* SUM_GRIDPOINT_T_H_ */
