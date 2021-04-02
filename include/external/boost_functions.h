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


#ifndef BOOST_FUNCTIONS_H_
#define BOOST_FUNCTIONS_H_

#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "../macro/griffin_definitions.h"
#include "../macro/macro_functions_read_write.h"


template< typename t_DATA>
inline
std::ostream &operator << ( std::ostream &STREAM, const boost::shared_ptr< t_DATA> &STR)
{
    STREAM << *STR << std::endl;
    return STREAM;
}

template< typename t_DATA>
inline
std::istream &operator >> ( std::istream &STREAM, boost::shared_ptr< t_DATA> &STR)
{
    t_DATA data;
    STREAM >> data;
    STR = boost::shared_ptr< t_DATA>( new t_DATA( data));
    return STREAM;
}

template< typename t_DATA>
inline
std::ostream &operator << ( std::ostream &STREAM, const boost::scoped_ptr< t_DATA> &STR)
{
    STREAM << *STR << std::endl;
    return STREAM;
}

template< typename t_DATA>
inline
std::istream &operator >> ( std::istream &STREAM, boost::scoped_ptr< t_DATA> &STR)
{
    t_DATA tmp;
    STREAM >> tmp;
    STR( new t_DATA( tmp));
    return STREAM;
}

//template< typename t_DATA>
//inline
//boost::shared_ptr< t_DATA> &operator += ( boost::shared_ptr< t_DATA> &PTR, const boost::shared_ptr< t_DATA> &B_PTR)
//{
//
//#ifdef DEBUG
//    if( PTR.use_count() == 0)
//    {
//        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
//    }
//    else
//    {
//        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
//    }
//    if( B_PTR.use_count() == 0)
//    {
//        std::cout << "===> no object behind second pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
//    }
//    else
//    {
//        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
//    }
//#endif
//
//    if( PTR.use_count() == 0)
//    {
//        PTR = B_PTR;
//    }
//    else
//    {
//        *PTR += *B_PTR;
//    }
//    return PTR;
//}

template< typename t_DATA>
inline
boost::shared_ptr< t_DATA> operator + ( const boost::shared_ptr< t_DATA> &A_PTR, const boost::shared_ptr< t_DATA> &B_PTR)
{
#ifdef DEBUG
    if( A_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  A_PTR->GetClassName() << std::endl;
    }
    if( B_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind second pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
    }
#endif
    if( A_PTR.use_count() == 0 && B_PTR.use_count() == 0)
    {
        std::cout << "===> no objects behind both pointers! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
        return boost::shared_ptr< t_DATA>();
    }
    else if( B_PTR.use_count() == 0)
    {
        return boost::shared_ptr< t_DATA>( A_PTR);
    }
    else if( A_PTR.use_count() == 0)
    {
        return boost::shared_ptr< t_DATA>( B_PTR);
    }
    return boost::shared_ptr< t_DATA>( new t_DATA( *A_PTR + *B_PTR));
}

template< typename t_DATA>
inline
boost::shared_ptr< t_DATA> &operator -= ( boost::shared_ptr< t_DATA> &PTR, const boost::shared_ptr< t_DATA> &B_PTR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
    if( B_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
    }
#endif
    *PTR -= *B_PTR;
    return PTR;
}

template< typename t_DATA>
inline
boost::shared_ptr< t_DATA> operator - ( const boost::shared_ptr< t_DATA> &A_PTR, const boost::shared_ptr< t_DATA> &B_PTR)
{
#ifdef DEBUG
    if( A_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  A_PTR->GetClassName() << std::endl;
    }
    if( B_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind second pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
    }
#endif
    boost::shared_ptr< t_DATA> tmp( new t_DATA());
    *tmp = *A_PTR - *B_PTR;
    return tmp;
}

template< typename t_DATA>
inline
boost::shared_ptr< t_DATA> &operator *= ( boost::shared_ptr< t_DATA> &PTR, const boost::shared_ptr< t_DATA> &B_PTR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
    if( B_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind second pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
    }
#endif
    *PTR *= *B_PTR;
    return PTR;
}

template< typename t_DATA>
inline
boost::shared_ptr< t_DATA> operator * ( const boost::shared_ptr< t_DATA> &A_PTR, const boost::shared_ptr< t_DATA> &B_PTR)
{
#ifdef DEBUG
    if( A_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind first pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind first boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  A_PTR->GetClassName() << std::endl;
    }
    if( B_PTR.use_count() == 0)
    {
        std::cout << "===> no object behind second pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind second boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  B_PTR->GetClassName() << std::endl;
    }
#endif
    boost::shared_ptr< t_DATA> tmp( new t_DATA());
    *tmp = *A_PTR * *B_PTR;
    return tmp;
}

template< typename t_DATA, typename t_SCALAR>
inline
boost::shared_ptr< t_DATA> &operator *= ( boost::shared_ptr< t_DATA> &PTR, const t_SCALAR &SCALAR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
#endif
    *PTR *= SCALAR;
    return PTR;
}

template< typename t_DATA, typename t_SCALAR>
inline
boost::shared_ptr< t_DATA> operator * ( const boost::shared_ptr< t_DATA> &PTR, const t_SCALAR &SCALAR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
#endif
    boost::shared_ptr< t_DATA> tmp( new t_DATA());
    *tmp = *PTR * SCALAR;
    return tmp;
}

template< typename t_SCALAR, typename t_DATA>
inline
boost::shared_ptr< t_DATA> &operator *= ( const t_SCALAR &SCALAR, boost::shared_ptr< t_DATA> &PTR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
#endif
    *PTR *= SCALAR;
    return PTR;
}

template< typename t_SCALAR, typename t_DATA>
inline
boost::shared_ptr< t_DATA> operator * ( const t_SCALAR &SCALAR, const boost::shared_ptr< t_DATA> &PTR)
{
#ifdef DEBUG
    if( PTR.use_count() == 0)
    {
        std::cout << "===> no object behind pointer! in " << __PRETTY_FUNCTION__ << "  "  << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    else
    {
        std::cout << "object behind boost pointer: " << __PRETTY_FUNCTION__ << ":  " <<  PTR->GetClassName() << std::endl;
    }
#endif
    boost::shared_ptr< t_DATA> tmp( new t_DATA( *PTR * SCALAR));
//    tmp *= SCALAR;
    return tmp;
}


//template< typename t_INPUT, typename t_RETURN>
//inline
//t_RETURN operator()( const boost::shared_ptr< math::Function< t_INPUT, t_RETURN> > &FUNCTION_POINTER, const t_INPUT &INPUT)
//{
//    return FUNCTION_POINTER->operator()( INPUT);
//}

#endif /* BOOST_FUNCTIONS_H_ */
