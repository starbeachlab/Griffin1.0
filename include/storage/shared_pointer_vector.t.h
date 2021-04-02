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


#ifndef SHARED_POINTER_VECTOR_T_H_
#define SHARED_POINTER_VECTOR_T_H_

#include <vector>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "../readwrite/class_name.h"
#include "../readwrite/stream_operator.h"

namespace store
{
    template< typename t_DATA>
    class ShPtrVec
    : public std::vector< boost::shared_ptr< t_DATA> >,
    public StreamOperator, public readwrite::ClassName
    {
    protected:
        ///////////////
        //    DATA     //
        ///////////////

    public:
        //////////////////////////////////
        //  CONSTRUCTION & DESTRUCTION  //
        //////////////////////////////////

        //! default constructor
        ShPtrVec()
        : std::vector< boost::shared_ptr< t_DATA> >(), StreamOperator(), readwrite::ClassName()
        {}

        //! default constructor
        ShPtrVec( const size_t &SIZE)
        : std::vector< boost::shared_ptr< t_DATA> >( SIZE), StreamOperator(), readwrite::ClassName()
        {}

        //! construct from parent data
        ShPtrVec( const std::vector< boost::shared_ptr< t_DATA> > &DATA)
        : std::vector< boost::shared_ptr< t_DATA> >( DATA), StreamOperator(), readwrite::ClassName()
        {}

        //! copy constructor - soft copy!!
        ShPtrVec( const ShPtrVec< t_DATA> &ORIGINAL)
        : std::vector< boost::shared_ptr< t_DATA> >( ORIGINAL.size()), StreamOperator(), readwrite::ClassName()
        {
            SoftCopy( ORIGINAL);
        }

        //! virtual destructor
        virtual ~ShPtrVec(){}

        //! virtual copy constructor
        //!(needed for copying pointers to derived classes)
        virtual ShPtrVec *Clone() const{ return new ShPtrVec( *this);}

        /////////////////////////
        //     OPERATORS       //
        /////////////////////////

        ShPtrVec< t_DATA> &operator = ( const ShPtrVec< t_DATA> &VEC)
        { /*std::cout << "operator =" << std::endl;*/ return SoftCopy( VEC);}

        ShPtrVec< t_DATA> &operator += ( const ShPtrVec< t_DATA> &VEC)
        {
            std::vector< boost::shared_ptr< t_DATA> >::insert( std::vector< boost::shared_ptr< t_DATA> >::end(), VEC.begin(), VEC.end());
            return *this;
        }

        virtual boost::shared_ptr< t_DATA> & operator()( const size_t &ID)
        {
            assert( ID < std::vector< boost::shared_ptr< t_DATA> >::size());
            return std::vector< boost::shared_ptr< t_DATA> >::operator[]( ID);
        }

        virtual const boost::shared_ptr< t_DATA> & operator()( const size_t &ID) const
        {
            assert( ID < std::vector< boost::shared_ptr< t_DATA> >::size());
            return std::vector< boost::shared_ptr< t_DATA> >::operator[]( ID);
        }

        virtual const boost::shared_ptr< t_DATA> & operator[]( const size_t &ID) const
        {
            assert( ID < std::vector< boost::shared_ptr< t_DATA> >::size());
            return std::vector< boost::shared_ptr< t_DATA> >::operator[]( ID);
        }

        /////////////////////////
        //      FUNCTIONS      //
        /////////////////////////

//        //! hard copy original to this (new pointers to new objects)
//        ShPtrVec< t_DATA> &HardCopy( const ShPtrVec< t_DATA> &ORIGINAL)
//        {
//            std::cout << "ShPtrVec Hard copy" << std::endl;
//            resize( ORIGINAL.size());
//            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator orig_itr( ORIGINAL.begin());
//            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator itr( this->begin()); itr != this->end(); ++itr, ++orig_itr)
//            {
//                *itr = boost::shared_ptr< t_DATA>( ( *orig_itr)->Clone());
//            }
//            return *this;
//        }

        //! hard copy this (new pointers to new objects)
        ShPtrVec< t_DATA> HardCopy() const
        {
            std::cout << "ShPtrVec Hard copy" << std::endl;
            ShPtrVec< t_DATA> copy( this->size());
            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator orig_itr( this->begin());
            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator itr( copy.begin()); itr != copy.end(); ++itr, ++orig_itr)
            {
                *itr = boost::shared_ptr< t_DATA>( ( *orig_itr)->Clone());
            }
            return copy;
        }

        //! soft copy this (new pointers to new objects)
        ShPtrVec< t_DATA> SoftCopy() const
        {
            std::cout << "ShPtrVec Hard copy" << std::endl;
            ShPtrVec< t_DATA> copy( this->size());
            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator orig_itr( this->begin());
            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator itr( copy.begin()); itr != copy.end(); ++itr, ++orig_itr)
            {
                *itr = *orig_itr;
            }
            return copy;
        }

        //! soft copy original to this (new pointers to old objects
        ShPtrVec< t_DATA> &SoftCopy( const ShPtrVec< t_DATA> &ORIGINAL)
        {
//            std::cout << __PRETTY_FUNCTION__ << std::endl;
//            std::cout << mystr::GetClassName( __PRETTY_FUNCTION__) << " soft copy" << std::endl;
            resize( ORIGINAL.size());
            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator orig_itr( ORIGINAL.begin());
            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator itr( this->begin()); itr != this->end(); ++itr, ++orig_itr)
            {
                *itr = *orig_itr;
            }
            return *this;
        }


        //! hard copy this to new (new pointers to new objects)
        ShPtrVec< t_DATA> GetHardCopy() const
        {
            ShPtrVec< t_DATA> copy( this->size());
            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator this_itr( this->begin());
            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator copy_itr( copy.begin()); this_itr != this->end(); ++copy_itr, ++this_itr)
            {
                *copy_itr = boost::shared_ptr< t_DATA>( ( *this_itr)->Clone());
            }
            return copy;
        }

        //! soft copy this to new (new pointers to old objects
        ShPtrVec< t_DATA> GetSoftCopy() const
        {
            ShPtrVec< t_DATA> copy( this->size());
            typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator this_itr( this->begin());
            for( typename std::vector< boost::shared_ptr< t_DATA> >::iterator copy_itr( copy.begin()); this_itr != this->end(); ++copy_itr, ++this_itr)
            {
                *copy_itr = *this_itr;
            }
            return copy;
        }



        //! data access
        virtual const std::vector< boost::shared_ptr< t_DATA> > &GetData() const{ return *this;}

        virtual std::vector< boost::shared_ptr< t_DATA> > &Data(){ return *this;}

        /////////////////////////
        //      Read/Write     //
        /////////////////////////

        virtual std::istream &Read( std::istream &STREAM)
        {
            std::string str;
            STREAM >> str;
            assert( str == GetClassName());
            size_t size;
            STREAM >> size;
            this->clear(); // empty vector
            for( size_t i = 0; i < size; ++i)
            {
                exit( -1);
//                t_DATA tmp;
//                STREAM >> tmp;
//                this->push_back( boost::shared_ptr< t_DATA>( new t_DATA( tmp)));
            }
            return STREAM;
        }

        virtual std::ostream &Write( std::ostream &STREAM) const
        {
            STREAM << GetClassName() << std::endl;
            STREAM << this->size() << std::endl;
            for( typename std::vector< boost::shared_ptr< t_DATA> >::const_iterator itr (this->begin()); itr != this->end(); ++itr)
            { STREAM << **itr << std::endl;}
            return STREAM;
        }

        virtual std::string GetClassName() const
        {
            return mystr::GetClassName( __PRETTY_FUNCTION__);
        }


    }; // end class ShPtrVec
} // end namespace store




#endif /* SHARED_POINTER_VECTOR_T_H_ */
