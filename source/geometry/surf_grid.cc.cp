#include "../../include/geometry/surf_grid.h"

namespace geom
{
	
	//! Generate random point on the sphere surface (used in MAYER)
	void GenerateRandomPointOnSphereSurface
	(
			float &XP,
			float &YP,
			float &ZP
	)
	{
//		StandardWrite( __FUNCTION__);
		float rp;
		do
		{
			XP = math::BasicRandom( -1.0, 1.0); // equal behaviour of random generators???
			YP = math::BasicRandom( -1.0, 1.0);
			ZP = math::BasicRandom( -1.0, 1.0);
			rp = XP*XP + YP*YP + ZP*ZP;
//			StandardWrite( "rp: " << rp);
		}
		while( rp == 0 || rp > 1.0);
		XP=XP/rp;
		YP=YP/rp;
		ZP=ZP/rp;
	}



	//  Quick sorting of arr (used in MAYER)
	void QuickSort
	(
			int &n,
			std::vector< float> &arr,
			std::vector< int> &ibrr
	)
	{

		StandardWrite( __FUNCTION__);
		int
			m=7,
			nstack=50,
			jstack,
			l,
			ir,
			ib,
			j,
			i,
			iq;
		float
			fm=7875.,
			fa=211.,
			fc=1663.,
			fmi = 1. / fm,
			a,
			fx;

		int istack[ nstack];


		jstack=0;
		l=1;
		ir=n;
		fx=0.;
		while( true)
		{
			if(ir - l < m)
			{
				for( j = l + 1; j <= ir; ++j)
				{
					a = arr[j];
					ib = ibrr[j];
					for( i = j - 1; i >= 1; --i) /// CHECK!
					{
						if(arr[i] <= a)
						{
							break;
							//goto 12
						}
						arr[i+1]=arr[i];
						ibrr[i+1]=ibrr[i];
					}
					if( i > 0 && arr[i] > a)
					{
						i=0;
					}
					arr[i+1]=a;  // 12
					ibrr[i+1]=ib;
				}

				if(jstack == 0)
				{
					return;
				}
				ir=istack[jstack];
				l=istack[jstack-1];
				jstack=jstack-2;
			}
			else
			{
				i=l;
				j=ir;
				fx = math::Modulo( fx*fa+fc, fm);  // mod(a,b) : a%b ??
				iq= int( l+(ir-l+1)*(fx*fmi)); // conversion to int??
				a=arr[iq];
				ib=ibrr[iq];
				arr[iq]=arr[l];
				ibrr[iq]=ibrr[l];
				while( true)
				{
					while(j > 0 && a < arr[j])
					{
						j=j-1;
					}

					if(j <= i)
					{
						arr[i]=a;
						ibrr[i]=ib;
						break;
					}

					arr[i]=arr[j];
					ibrr[i]=ibrr[j];
					i=i+1;

					while(i <= n && a > arr[i])
					{
						i=i+1;
					}

					if(j <= i)
					{
						arr[j]=a;
						ibrr[j]=ib;
						i=j;
						break;
					}

					arr[j]=arr[i];
					ibrr[j]=ibrr[i];;
					j=j-1;

				}

				jstack=jstack+2;
				if(jstack > nstack)
				{
					std::cerr << __FUNCTION__ << ": nstack must be larger!" << std::endl;
					exit( -1);
				}
				if(ir-i >= i-l)
				{
					istack[jstack]=ir;
					istack[jstack-1]=i+1;
					ir=i-1;
				}
				else
				{
					istack[jstack]=i-1;
					istack[jstack-1]=l;
					l=i+1;
				}
			}
		}
	}





	//  Segment part formula (used in MAYER)
	float
	SegPar
	(
			float &raneti,
			float &ranetj,
			float &aij
	)
	{
		float v1, v2, segpar;

		v1=raneti+ranetj;
		v2=fabs(raneti-ranetj);
		if(aij > v1)
		{
			segpar=0.;
		}
		else if(aij < v2)
		{
			if(raneti > ranetj)
			{
				segpar=0.;
			}
			else
			{
				segpar=1.;
			}
		}
		else
		{
			segpar = ( ranetj * ranetj - pow( ( raneti - aij), 2)) / (4 * aij * raneti);
		}
		return segpar;
	}

	//TRANX = HALF*(NCLX-1)*DCEL;


	void
	MayerContactAndReentrance
	(
	 const  int     &NTPRP,  /// nr atoms
	 const  float  *X, // atom pos of protein
	 const  float  *Y,
	 const  float  *Z,
	 const  float  *PBRAD, // vdw radii
	 const  float  &RADW, // probe radius
	 const  int     &NCLX,  // number of cells of grid
	 const  int     &NCLY,
	 const  int     &NCLZ,
	 const  float  &DCEL,   // cell size
	 const  float  &TRANX,  //
	 const  float  &TRANY,
	 const  float  &TRANZ,
 	 const  float  &XBCEN,  // center of the grid
	 const  float  &YBCEN,
	 const  float  &ZBCEN,
	 std::vector< float>  &RMAY,
			float  *POX,     // sphere surface positions
			float  *POY,
			float  *POZ,
	 const  int     &MAPT,   // number of spheres points
	 std::vector< int>     &LISTR,  // index list of neighbors
	 std::vector< float>  &FISTR   // neighbor list with score of ability to influence
	)
	{

		StandardWrite( __FUNCTION__);

		int  iseed, jl, i, j, l, m, nl, ml, kl, ll, nfil, ix, iy, iz;
		int  jx1, jx2, jy1, jy2, jz1, jz2, k, ipx, ipy, ipz, ncyz;
		float  wrsf, raneti, ranetj, xi, yi, zi, xj, yj, zj, aij, fistr0;
		float  xp, yp, zp, bisq, xc, yc, zc, xsq, ysq, zsq, dsq;
		//  float  segpar;

		float rsmall = 1e-6; // as of squantm/sqnt_qm2_mopac.src90

		iseed = 314159265;
		ncyz = NCLY * NCLZ;

//#ifdef NEWRNG
		srand( iseed);
		//iseed = 1;
		//Cmh050712      ISEED=MYNODP              !##PARALLEL
//#endif

		// generate points on the sphere surface and adjust parameter
		wrsf = 0.0;
		for( i = 0; i < NTPRP; ++i)
		{

			raneti = PBRAD[i] + RADW;
			if( wrsf < raneti) // get max atom radius plus probe radius
			{
				wrsf = raneti;
			}
		}
		wrsf = MAPT / ( wrsf * wrsf) + 0.000001; // ???
		StandardWrite( "sphere surface point creation ...");
		for( i = 0; i < MAPT; ++i)
		{
			GenerateRandomPointOnSphereSurface( POX[ i], POY[ i], POZ[ i]);
		}

		StandardWrite( "main loop ...");

		// main loop by atoms
		for( i = 0; i < NTPRP; ++i)
		{
			nl = 0;

//			FISTR.clear();
//			LISTR.clear();
			std::multimap< float, int> neighbor_map;

			xi = X[ i];
			yi = Y[ i];
			zi = Z[ i];
			raneti = PBRAD[ i] + RADW;
			StandardWrite( "atom " << i << " sort closest neighbors");

			// create the list of closest neighbors for the atom from the main loop
			// according to their screening ability of the neighbors
			// 0<Segpar<1
			// Segpar=0  if the neighbor does not reduce the accessible surface of the atom
			//           and it means that it is not a neighbor at all.
			// Segpar=1  if the neighbor buries the atom totally.
			for( j = 0; j < NTPRP; ++j)
			{

				if( i == j)
				{
					continue;
				}
				xj = X[ j];
				yj = Y[ j];
				zj = Z[ j];
				ranetj = PBRAD[ j] + RADW;
				aij = sqrt( ( xi - xj) * ( xi - xj) + ( yi - yj) * ( yi - yj) + ( zi - zj) * ( zi - zj));
				fistr0 = SegPar( raneti, ranetj, aij);
//				std::cout << i << ": " << raneti << ", " << j << ": " << ranetj << " dist: " << aij << " segpar: " << fistr0 << std::endl;
				if( fistr0 == 0.)
				{
					continue;
				}
				if( fistr0 == 1.)
				{
					std::cout << "get out of loop ..." << std::endl;
					break;
				}

				++nl;
//				FISTR.push_back( fistr0);
//				LISTR.push_back( j);
				neighbor_map.insert( std::make_pair( fistr0, j));
			}

			if( fistr0 == 1.)
			{
				std::cout << "... and to the end of this one" << std::endl;
				continue;
			}
			// sort neighbors according to their ability to
			// reduce an accessible surface of the atom
//			{
//					std::vector< float> fist_cp( FISTR);
//					std::vector< int> list_cp( LISTR);
//
//			if( nl > 0)
//			{
//				QuickSort( nl, FISTR, LISTR);
//			}
//
//			std::vector< int>::const_iterator  l_itr = list_cp.begin(),  L_itr = LISTR.begin();
			std::cout << "neigbours of atom " << i << ": " << std::endl;

			for( std::multimap< float, int>::const_iterator itr = neighbor_map.begin(); itr != neighbor_map.end(); ++itr)
			{
				std::cout <<  itr->first << " : " << itr->second << std::endl;
			}
//			}
			ml = int( wrsf * raneti * raneti);
			if( ml > MAPT)
			{
				ml = MAPT;
			}


			StandardWrite( "loop over sphere surface");

			// loop over points on the sphere surface
			for( jl = 0; jl < ml; ++jl)
			{
				xp = xi + raneti * POX[ jl];
				yp = yi + raneti * POY[ jl];
				zp = zi + raneti * POZ[ jl];

				// reject the points which are not on the accessible surface.
				if( nl > 0)
				{
					std::multimap< float, int>::const_iterator m_itr = neighbor_map.begin();
					for( kl = 0; kl < nl; ++kl, ++m_itr)
					{
						ll = m_itr->second; // LISTR[ kl];
						aij = pow( ( xp - X[ ll]), 2) + pow( ( yp - Y[ ll]), 2) + pow( ( zp - Z[ ll]), 2);  // distance of the surface point to a neighbor atom
						if( aij < pow( ( PBRAD[ll] + RADW), 2))
						{
							std::cout << "surf point closer to neighbor than probe radius: " << jl << " neigb: " << kl << " id: " << ll << " dist: " << aij << " < (" << PBRAD[ll] << " + " <<  RADW << ")^2 =  " << pow( ( PBRAD[ll] + RADW), 2) << std::endl;
							break;
						}
					}
					if( aij < pow( ( PBRAD[ll] + RADW), 2))
					{
						std::cout << "... sneak out of this iteration" << std::endl;
						continue;
					}
					else
					{
						std::cout << "this one seems just fine! " << std::endl;
					}
				}

				// reset grid points inside a probe sphere around any
				// of random accessible points as a dielectric media.
				bisq = RADW;
				nfil = int( bisq / DCEL) + 2;
				xp = xp + TRANX - XBCEN;
				yp = yp + TRANY - YBCEN;
				zp = zp + TRANZ - ZBCEN;
				bisq = bisq * bisq + rsmall;

				ix = int( xp / DCEL) + 1;
				iy = int( yp / DCEL) + 1;
				iz = int( zp / DCEL) + 1;

				jx1 = ix - nfil + 1;
				if( jx1 < 1){ jx1 = 1;}
				jx2 = ix + nfil;
				if( jx2 > NCLX){ jx2 = NCLX;}
				jy1 = iy - nfil + 1;
				if( jy1 < 1){ jy1 = 1;}
				jy2 = iy + nfil;
				if( jy2 > NCLY){ jy2 = NCLY;}
				jz1 = iz - nfil + 1;
				if( jz1 < 1){ jz1 = 1;}
				jz2 = iz + nfil;
				if( jz2 > NCLZ){ jz2 = NCLZ;}

				for( k = jx1; k <= jx2; ++k)
				{
					ipx = ( k - 1) * ncyz;
					xc  = ( k - 1) * DCEL;
					xsq = ( xc - xp) * ( xc - xp);
					for( l = jy1; l <= jy2; ++l)
					{
						ipy = ( l - 1) * NCLZ;
						yc  = ( l - 1) * DCEL;
						ysq = ( yc - yp) * ( yc - yp);
						for( m = jz1; m <= jz2; ++m)
						{
							ipz = (m - 1) + ipy + ipx;
							zc = ( m - 1) * DCEL;
							zsq = ( zc - zp) * ( zc - zp);
							dsq = xsq + ysq + zsq;


							if( RMAY[ ipz] < 0.)
							{
								if( dsq <= bisq)
								{
									std::cout << "invert rmay value " << k << " " << l << " " << m << " " << ipz << std::endl;
									assert( ipz < RMAY.size());
									RMAY[ ipz] = -RMAY[ ipz]; // this defines points in the probe shell to be outside of smoothed surface // ! bulk kappa restored
								}
							}

						}
					}
				}
			}
		}
		StandardWrite( "set remaining negative values in rmay to zero = inside smoothed surf");
		int cc =0;
		for( k = 0; k < NCLX; k++)
		{
			ipx = k * ncyz;
			for( l = 0; l < NCLY; l++)
			{
				ipy = l * NCLZ;
				for( m = 0; m < NCLZ; ++m, ++cc)
				{
					ipz = m + ipy + ipx;
//					std::cout << k << " " << l << " " << m << " " << ipz << " " << cc << std::endl;
					if(  RMAY[ipz] < 0.)
					{
						std::cout << "generously set to zero" << std::endl;
						RMAY[ipz] = 0.;
					}
				}
			}
		}
		std::cout << RMAY.size() << std::endl;
	} // end MayerContactAndReentrance



	//  when the dielectric boundary is defined by the van der Waals surface
	//  M = 0 : inside solute
	//      1 : outside solute
	void
	MayerMStep
	(
   	   const  int    &NTPRP,  // number of atoms
   	   const  float *X,
   	   const  float *Y,
   	   const  float *Z,
   	   const  float *PBRAD,  // radii
   	   const  float &RADW,  // water radius
   	   const  int    &NCLX,  // number of cells in x
   	   const  int    &NCLY,
   	   const  int    &NCLZ,
   	   const  float &DCEL,  // cell size
   	   const  float &TRANX,
   	   const  float &TRANY,
   	   const  float &TRANZ,
   	   const  float &XBCEN, // box center x coordin
   	   const  float &YBCEN,
   	   const  float &ZBCEN,
  std::vector< float> &RMAY
	)
	{
		StandardWrite( __FUNCTION__);

		int      i, k, l, m, ix, iy, iz, jx1, jx2, jy1, jy2, jz1, jz2;
		int      ipx, ipy, ipz;
		int      nfil, ncyz;
		float   xi, yi, zi, xc, yc, zc, dsq, xsq, ysq, zsq;
		float   sqr, sqrw;


		//cwi      ncyz = NCLX * NCLY
		ncyz = NCLY * NCLZ;
		for( i = 0; i < NTPRP; ++i)
		{
			sqr  = PBRAD[i];
			if(sqr <= 0.0)
			{
				break;
			}
		  xi = X[i] + TRANX - XBCEN;
		  yi = Y[i] + TRANY - YBCEN;
		  zi = Z[i] + TRANZ - ZBCEN;
		  sqrw = sqr + RADW;
		  nfil = int( sqrw / DCEL) + 2;
		  sqr  = sqr  * sqr;
		  sqrw = sqrw * sqrw;
		  ix = int(xi / DCEL) + 1;
		  iy = int(yi / DCEL) + 1;
		  iz = int(zi / DCEL) + 1;

		  jx1 = ix - nfil + 1;
		  if(jx1 < 1){ jx1 = 1;}
		  jx2 = ix + nfil - 1;
		  if(jx2 > NCLX){ jx2 = NCLX;}
		  jy1 = iy - nfil + 1;
		  if(jy1 < 1){ jy1 = 1;}
		  jy2 = iy + nfil - 1;
		  if(jy2 > NCLY){ jy2 = NCLY;}
		  jz1 = iz - nfil + 1;
		  if(jz1 < 1){ jz1 = 1;}
		  jz2 = iz + nfil - 1;
		  if(jz2 > NCLZ){ jz2 = NCLZ;}

		  for( k = jx1; k <= jx2; ++k)
		  {
			  ipx = (k - 1) * ncyz;
			  xc = (k - 1) * DCEL;
			  for( l = jy1; l <= jy2; ++l)
			  {
				  ipy = (l - 1) * NCLZ;
				  yc = (l - 1) * DCEL;
				  for( m = jz1; m <= jz2; ++m)
				  {
//					  std::cout << k << "  " << l << "  " << m << "  " << ipz << std::endl;
					  ipz = m + ipy + ipx - 1;  // the -1 solves the z-shift, arrays start at zero!!
					  zc = (m - 1) * DCEL;
					  xsq = (xc - xi) * (xc - xi);
					  ysq = (yc - yi) * (yc - yi);
					  zsq = (zc - zi) * (zc - zi);

					  if(RMAY[ipz] != 0.)
					  {
						  dsq = xsq + ysq + zsq;  // distance squared to neares atom center
						  if(dsq <= sqr)   // if square of distance is smaller than square of the radius => zero = inside of the atom
						  {
							  // Zero the Debye-Huckel factor inside the solute
//										  std::cout << "// Zero the Debye-Huckel factor inside the solute" << std::endl;
							  RMAY[ipz] = 0.;
						  }
						  else if(dsq > sqr && dsq <= sqrw)
						  {
							  if(RMAY[ipz] > 0.)
							  {
//								  std::cout << "invert rmay" << std::endl;
								  RMAY[ipz] = -RMAY[ipz];
							  }
						  }
					  }
				  }
			  }
		  }
		}
	}// end MayerMStep


	void
	BuildSurfGrid
	(
			const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
			std::vector< float> &RMAY,
			math::Vector3N &MIN,
			store::Vector3N< size_t> &NR_BINS,
			float &DELTA
	)
	{
		StandardWrite( __FUNCTION__);


		int size = MOLECULE->GetAtoms().size();

		float x_coo[ size];
		float y_coo[ size];
		float z_coo[ size];
		float vdw_radii[ size];

		std::vector< int> nr_cells( 3);
		std::vector< float> center( 3);

		float
			probe_radius = 1.4;

		math::Vector3N pos, min, max, total_delta;

		min.SetAll( std::numeric_limits< float>::max());
		max.SetAll( -std::numeric_limits<float>::max());

		for( int i = 0; i < size; ++i)
		{
			pos = MOLECULE->GetAtoms()( i)->GetPosition();
			x_coo[ i] = pos( 0);
			y_coo[ i] = pos( 1);
			z_coo[ i] = pos( 2);

			for( int j = 0; j < 3; ++j)
			{
				if( pos( j) < min( j))
				{
					min( j) = pos( j);
				}
				if( pos( j) > max( j))
				{
					max( j) = pos( j);
				}
			}

			vdw_radii[ i] = MOLECULE->GetAtoms()( i)->GetVanDerWaalsRadius();
//			std::cout <<  x_coo[i] << "  " << y_coo[ i] << "  " << z_coo[ i] << "  " << vdw_radii[i] << std::endl;
		}

		for( int i = 0; i < 3; ++i)
		{
			max( i) += 3.0;
			min( i) -= 3.0;
			total_delta(i) = max( i) - min( i);
			nr_cells[ i] = int( total_delta(i) / DELTA) + 1;
			total_delta(i) = nr_cells[i] * DELTA;
			center[i] = 0.5 * ( max( i) + min( i));
			NR_BINS(i) = nr_cells[i]; // rubbish
			MIN( i) = center[i] - 0.5 * total_delta( i);
		}

		const int total_array_length = nr_cells[0] * nr_cells[1] * nr_cells[2];
		RMAY. resize( total_array_length, 1);

		int number_surface_points = 10000;
		float surf_point_x[ number_surface_points];
		float surf_point_y[ number_surface_points];
		float surf_point_z[ number_surface_points];

		std::vector< float> neighbor_score_list;
		std::vector< int>    neighbor_index_list;

		MayerMStep
		(
				size,
				x_coo,
				y_coo,
				z_coo,
				vdw_radii,
				probe_radius,
				nr_cells[0],
				nr_cells[1],
				nr_cells[2],
				DELTA,
				0.5 * ( nr_cells[ 0] - 1) * DELTA,
				0.5 * ( nr_cells[ 1] - 1) * DELTA,
				0.5 * ( nr_cells[ 2] - 1) * DELTA,
				center[0],
				center[1],
				center[2],
				RMAY
		);

		std::vector< float> copy( RMAY);
		std::ofstream write;
		Open( write, "rmay_grid_mstep.txt");
		std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( write , " "));
		write << std::endl;
		Close( write);

		MayerContactAndReentrance
		(
				size,
				x_coo,
				y_coo,
				z_coo,
				vdw_radii,
				probe_radius,
				nr_cells[0],
				nr_cells[1],
				nr_cells[2],
				DELTA,
				0.5 * ( nr_cells[ 0] - 1) * DELTA,
				0.5 * ( nr_cells[ 1] - 1) * DELTA,
				0.5 * ( nr_cells[ 2] - 1) * DELTA,
				center[0],
				center[1],
				center[2],
				RMAY,
				surf_point_x,     // sphere surface positions
				surf_point_y,
				surf_point_z,
				number_surface_points,   // number of spheres points
				neighbor_index_list,
				neighbor_score_list
		);

		Open( write, "rmay_grid_nreen.txt");
		std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( write , " "));
		write << std::endl;
		Close( write);

		int changed = 0;
		for( std::vector< float>::const_iterator r_itr = RMAY.begin(), c_itr = copy.begin(); r_itr != RMAY.end() && c_itr != copy.end(); ++r_itr, ++c_itr)
		{
			if( *r_itr != *c_itr)
			{
				++changed;
			}
		}
		std::cout << changed << " positions in array have changed due to smoothing" << std::endl;

		std::cout << "some info:\n#cells:" << std::endl;
		std::copy( nr_cells.begin(), nr_cells.end(), std::ostream_iterator< size_t>( std::cout , " "));
		std::cout << "delta: " << DELTA << std::endl;
		std::cout << "product #cells: " << total_array_length << std::endl;

//		std::cout << "rmay: " << std::endl;
//		std::copy( RMAY.begin(), RMAY.end(), std::ostream_iterator< float>( std::cout , " "));
//		std::cout << std::endl;

	}


	float
	GetGridPoint
	(
			const math::Vector3N &POSITION,
			const std::vector< float> &RMAY,
			const math::Vector3N &MIN,
			const store::Vector3N< size_t> &NR_BINS,
			const float &DELTA
	)
	{
		std::vector< size_t> ids( 3);
		math::Vector3N tmp = (POSITION - MIN) / DELTA;
		for( int i = 0; i < 3; ++i)
		{
			ids[i] = size_t( tmp( i));
		}
		size_t id = math::TupleXD< float>( NR_BINS).CalcID( ids);

		return RMAY[ id];
	}


	void
	SurfGrid
	(
			const boost::shared_ptr< mol::SimpleMolecule< mol::Atom> > & MOLECULE,
			const std::string &NAME,
			const size_t &NR_POINTS
	)
	{
		float value;
		std::vector< float>  rmay;
		math::Vector3N pos, min, max;
		store::Vector3N< size_t> nr_cells;
		float delta = 0.5;
		std::cout << __FUNCTION__ << ": build surf grid" << std::endl;

		std::string name = NAME.substr( 0, NAME.size() - 4);

		std::ofstream
			write;

		Open( write, name + "_grid.pdb");

		BuildSurfGrid( MOLECULE, rmay, min, nr_cells, delta);

		max = min + delta * nr_cells;
		int cc = 0;
		for( int i = 0; i < nr_cells[0]; ++i)
			for( int j = 0; j < nr_cells[1]; ++j)
				for( int k = 0; k < nr_cells[2]; ++k, ++cc)
				{
					int l = k + j * nr_cells[2] + i * nr_cells[2] * nr_cells[1];
//					std::cout.width(4);
//					std::cout << i << "  ";
//					std::cout.width(4);
//					std::cout << j << "  " ;
//					std::cout.width(4);
//					std::cout << k << "  ";
//					std::cout.width(8);
//					std::cout << l << "   ";
//					std::cout.width(8);
//					std::cout << cc << std::endl;
					value = rmay[ l];
					if( value == 0)
					{
						pos = min;
						pos( 0) += (i + 0.5) * delta;
						pos( 1) += (j + 0.5) * delta;
						pos( 2) += (k + 0.5) * delta;
						geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(cc) / 1e5), cc);
					}
					if( value == -1)
					{
						pos = min;
						pos( 0) += (i + 0.5) * delta;
						pos( 1) += (j + 0.5) * delta;
						pos( 2) += (k + 0.5) * delta;
						geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(cc) / 1e5), cc, "O");
					}
				}


//		for( int i = 0; i < NR_POINTS; ++i)
//		{
//			pos.Randomize( min, max);
//			value = GetGridPoint( pos, rmay, min, nr_cells, delta);
//			if( value == 0)
//			{
//				geom::DirectedSurfacePoint( pos, math::Vector3N()).WriteAsPdb( write, int( float(i) / 1e5), i);
//			}
//		}

		Close( write);
	}



} // end namespace geom




