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


#ifndef SURF_GRID_H
#define SURF_GRID_H

#include <iostream>
#include <cmath>
#include <iterator>
#include <set>

#include <boost/shared_ptr.hpp>

#include "../math/double_functions.h"
#include "../molecule/simple_molecule.t.h"
#include "../molecule/atom.h"
#include "../math/tupleXD.t.h"
#include "../external/std_functions.h"
#include "../readwrite/stream_functions.h"
#include "../storage/multiplets.h"
#include "../math/tupleXD.t.h"
#include "../storage/map.t.h"
#include "../storage/limits3D.h"

#include "../geometry/object.h"

namespace geom
{

    class SurfGrid
    : public math::TupleXD< short, 3>
   {
   protected:
        mutable std::vector< std::multimap< short, std::pair< short, short> > >    m_SurfIndices;  // outer for moltype,
        mutable std::vector< std::string> m_MolTypes;


   public:

        SurfGrid()
        : math::TupleXD< short, 3>(),
          m_SurfIndices( 1),
          m_MolTypes( 1, "all")
          {}

          SurfGrid
          (
                  const float &VOXELSIZE,
                  const store::Limits3D & LIMITS = store::Limits3D()
          )
          : math::TupleXD< short, 3>( VOXELSIZE, LIMITS.GetMin(), LIMITS.GetMax(), short( 1)),
            m_SurfIndices( 1),
            m_MolTypes( 1, "all")
            {
              //            BuildFromMolecule( MOLECULE, SMOOTH, PROBE_RADIUS);
            }

          virtual ~SurfGrid(){}



        virtual
        const std::vector< std::string> &
        GetMoleculeTypes() const
        {
            return m_MolTypes;
        }

        virtual
        const std::vector< std::multimap< short, std::pair< short, short> > > &
        GetSurfaceIndices() const
        {
            return m_SurfIndices;
        }

        virtual
        std::vector< int>
        GetIndicesOfSurfacePoint( const std::multimap< short, std::pair< short, short> >::const_iterator &ITR) const
        {
            std::vector< int> indices( 3);
            indices[0] = ITR->first;
            indices[1] = ITR->second.first;
            indices[2] = ITR->second.second;
            return indices;
        }

        virtual void
        ClearData()
        {
            math::TupleXD< short, 3>::m_Data.clear();
            m_SurfIndices.clear();
        }

        virtual
        void
        AnalyzeObjectMap(  const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS) const;

        virtual
        void
        AddSurfaceObjects( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS);

        virtual
        void
        AddSurfaceObject( const boost::shared_ptr< geom::Object> &OBJECT, const short &MOL_TYPE_ID);

        virtual
        void 
        AdjustLimitsForSurfaceObjects( const boost::shared_ptr< store::Map< std::string, store::ShPtrVec< geom::Object> > > &OBJECTS);

        virtual
        void 
        AdjustLimitsForSurfaceObject( const boost::shared_ptr< geom::Object> &OBJECT);
        

        virtual
        void
        BuildFromMolecule
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
                const bool &SMOOTH = true,
                const float &PROBE_RADIUS = 1.4
        );

        virtual
        void
        Prepare
        (
                const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE,
                float *X_COO,
                float *Y_COO,
                float *Z_COO,
                float *VDW_RADII,
                std::vector< float> &CENTER,
                const float &FACTOR,
                const float &OFFSET,
		const float &MARGIN
        );

        virtual
        void
        Build
        (
                const bool &SMOOTH,
                const float &PROBE_RADIUS,
                const size_t &SIZE,
                const float *X_COO,
                const float *Y_COO,
                const float *Z_COO,
                const float *VDW_RADII,
                const std::vector< float> &CENTER
        );

//        virtual
//        void
//        BuildVanDerWaalsGrid( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE);

//        virtual
//        void
//        SmoothGrid( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE);

        //! to distinguish different molecule types indices are translated to keys: 'all': -1, NEXT: -2, ..
        virtual
        int
        IndexToKey( const int &ID) const;

        //! returns the vector indices of a grid value below zero
        virtual
        std::vector< int>
        KeyToIndices( const int &KEY) const;

        virtual
//        std::vector< store::Vector3N< short> >
        void
        IndicesOfSurfaceGridPoints() const;

//        virtual
//        store::Vector3N< short>
//        FindNearestGridPoint( const math::Vector3N &POSITION) const;
//
//        virtual
//        float
//        DistanceBetweenGridPoints( const store::Vector3N< short> &IDS_A, const store::Vector3N< short>  &IDS_B) const;

        virtual
        float
        CalcVoxelVolume() const;

        virtual
        store::Vector3N< int>
        AdjustMaxIndices( const store::Vector3N< int> &MAX) const;

        virtual
        store::Vector3N< int>
        AdjustMinIndices( const store::Vector3N< int> &MIN) const;

        virtual
        std::vector< float>
        Volume( const store::Vector3N< int> &MIN, const store::Vector3N< int> &MAX) const;

        virtual
        store::Map< std::string, float>
        Volume( const math::Vector3N &MIN, const math::Vector3N &MAX) const;

        virtual
        size_t
        CountAtomsInside( const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > &MOLECULE) const;

        //! return distance and vector id of closest point on the surface
        virtual
        std::pair< float, int>
        ClosestSurfaceGridPoint( const math::Vector3N &POSITION, const std::string &TYPE) const;

        //! return distance and vector id of closest point on the surface
        virtual
        std::pair< float, int>
        BasicClosestSurfPoint( const math::Vector3N &POSITION, const std::string &TYPE) const;

        virtual
        bool
        IsInside( const short &VALUE, const std::string &NAME) const;

        virtual
        std::ofstream &
        WriteSurfAsPdb( std::ofstream &STREAM) const;

        virtual
        std::ostream &
        Write( std::ostream &STREAM) const;

        virtual
        std::istream &
        Read( std::istream &STREAM);

   };
    
    std::istream & operator >> ( std::istream &STREAM, SurfGrid &GRID);

    std::ostream & operator << ( std::ostream &STREAM, const SurfGrid &GRID);


    //! Generate random point on the sphere surface (used in MAYER)
    void GenerateRandomPointOnSphereSurface
    (
            float &XP,
            float &YP,
            float &ZP
    );

    //  Segment part formula (used in MAYER)
    float
    SegmentPartFormula
    (
            const float &FIRST_VDW_PLUS_PROBE_RADIUS,
            const float &SECOND_VDW_PLUS_PROBE_RADIUS,
            const float &ATOM_DISTANCE
    );
    //TRANX = HALF*(NCLX-1)*DCEL;


    void
    MayerContactAndReentrance
    (
         const  int      &NTPRP,  /// nr atoms
         const  float   *X, // atom pos of protein
         const  float   *Y,
         const  float   *Z,
         const  float   *PBRAD, // vdw radii
         const  float   &RADW, // probe radius
         const  int      &NCLX,  // number of cells of grid
         const  int      &NCLY,
         const  int      &NCLZ,
         const  float   &DCEL,   // cell size
         const  float   &TRANX,  // translation to middle of grid from min
         const  float   &TRANY,
         const  float   &TRANZ,
         const  float   &XBCEN,  // center of the grid
         const  float   &YBCEN,
         const  float   &ZBCEN,
         std::vector< short>  &RMAY,  // the actual grid
                float   *POX,     // sphere surface positions
                float   *POY,
                float   *POZ,
         const  int      &MAPT   // number of spheres points
    );



    //  when the dielectric boundary is defined by the van der Waals surface
    //  M = 0 : inside solute
    //      1 : outside solute
    void
    MayerMStep
    (
         const  int    &NTPRP,
         const  float *X,
         const  float *Y,
         const  float *Z,
         const  float *PBRAD,
         const  float &RADW,
         const  int    &NCLX,
         const  int    &NCLY,
         const  int    &NCLZ,
         const  float &DCEL,
         const  float &TRANX,
         const  float &TRANY,
         const  float &TRANZ,
         const  float &XBCEN,
         const  float &YBCEN,
         const  float &ZBCEN,
         std::vector< short> &RMAY
     );


    void
    PrepareSurfGrid
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            //            std::vector< short> &RMAY,
            std::vector< float> &MIN,
            std::vector< float> &MAX,
            std::vector< int> &NR_BINS,
            const float &DELTA,
            //            const float &SMOOTH, 
            //            const float &PROBE_RADIUS,
            float *X_COO,
            float *Y_COO,
            float *Z_COO,
            float *VDW_RADII,
            std::vector< float> &CENTER,
            const float &FACTOR,
            const float &OFFSET,
	    const float &MARGIN
    );


    void
    BuildSurfGrid
    (
            std::vector< short> &RMAY,
            std::vector< float> &MIN,
            std::vector< int> &NR_BINS,
            const float &DELTA,
            const float &SMOOTH, 
            const float &PROBE_RADIUS,
            const size_t &SIZE,
            const float *X_COO,
            const float *Y_COO,
            const float *Z_COO,
            const float *VDW_RADII,
            const std::vector< float> &CENTER
    );


    void
    BuildSurfGrid
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            std::vector< short> &RMAY,
            std::vector< float> &MIN,
            std::vector< int> &NR_BINS,
            const float &DELTA,
            const float &SMOOTH = true, 
            const float &PROBE_RADIUS = 1.4
    );

    float
    GetGridPoint
    (
            const math::Vector3N &POSITION,
            const std::vector< float> &RMAY,
            const math::Vector3N &MIN,
            const store::Vector3N< size_t> &NR_BINS,
            const float &DELTA
    );

    void
    CalcSurfGrid
    (
            const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
            const std::string &NAME,
            const size_t &NR_POINTS = 1e6
    );



} // end namespace charmm


#endif //SURF_GRID_H
