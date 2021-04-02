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


#ifndef POLYNOMIAL_FITTING_H_
#define POLYNOMIAL_FITTING_H_


namespace NR
{
    inline
    void PolynomialInterpolation( const math::Vector &XA, const math::Vector &YA, const float X, float &y, float &dy)
    {
            int i,m,ns=0;
            float den,dif,dift,ho,hp,w;

            int n=XA.size();
            math::Vector c(n),d(n);
            dif=fabs(X-XA[0]);
            for (i=0;i<n;i++) {
                    if ((dift=fabs(X-XA[i])) < dif) {
                            ns=i;
                            dif=dift;
                    }
                    c[i]=YA[i];
                    d[i]=YA[i];
            }
            y=YA[ns--];
            for (m=1;m<n;m++) {
                    for (i=0;i<n-m;i++) {
                            ho=XA[i]-X;
                            hp=XA[i+m]-X;
                            w=c[i+1]-d[i];
                            if ((den=ho-hp) == 0.0) nrerror("Error in routine polint");
                            den=w/den;
                            d[i]=hp*den;
                            c[i]=ho*den;
                    }
                    y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
            }
    }

    inline // TODO : cpp file!
    void PolynomialInterpolation2Dim( const math::Vector &X1A, const math::Vector &X2A, const math::Matrix &YA, const float X1,
            const float x2, float &y, float &dy)
    {
            int j,k;

            int m=X1A.size();
            int n=X2A.size();
            math::Vector ymtmp(m),ya_t(n);
            for (j=0;j<m;j++) {
                    for (k=0;k<n;k++) ya_t[k]=YA[j][k];
                    PolynomialInterpolation(X2A,ya_t,x2,ymtmp[j],dy);
            }
            PolynomialInterpolation(X1A,ymtmp,X1,y,dy);
    }

    void FitPolynomialCoefficients(Vec_I_DP &x, Vec_I_DP &y, Vec_O_DP &cof)
    {
            int k,j,i;
            float phi,ff,b;

            int n=x.size();
            Vec_DP s(n);
            for (i=0;i<n;i++) s[i]=cof[i]=0.0;
            s[n-1]= -x[0];
            for (i=1;i<n;i++) {
                    for (j=n-1-i;j<n-1;j++)
                            s[j] -= x[i]*s[j+1];
                    s[n-1] -= x[i];
            }
            for (j=0;j<n;j++) {
                    phi=n;
                    for (k=n-1;k>0;k--)
                            phi=k*s[k]+x[j]*phi;
                    ff=y[j]/phi;
                    b=1.0;
                    for (k=n-1;k>=0;k--) {
                            cof[k] += b*ff;
                            b=s[k]+x[j]*b;
                    }
            }
    }


    void InterpolatePolynomialCoefficients(Vec_I_DP &xa, Vec_I_DP &YA, Vec_O_DP &cof)
    {
            int k,j,i;
            float xmin,dy;

            int n=xa.size();
            Vec_DP x(n),y(n);
            for (j=0;j<n;j++) {
                    x[j]=xa[j];
                    y[j]=YA[j];
            }
            for (j=0;j<n;j++) {
                    Vec_DP x_t(n-j),y_t(n-j);
                    for (k=0;k<n-j;k++) {
                            x_t[k]=x[k];
                            y_t[k]=y[k];
                    }
                    PolynomialInterpolation(x_t,y_t,0.0,cof[j],dy);
                    xmin=1.0e38;
                    k = -1;
                    for (i=0;i<n-j;i++) {
                            if (fabs(x[i]) < xmin) {
                                    xmin=fabs(x[i]);
                                    k=i;
                            }
                            if (x[i] != 0.0)
                                    y[i]=(y[i]-cof[j])/x[i];
                    }
                    for (i=k+1;i<n-j;i++) {
                            y[i-1]=y[i];
                            x[i-1]=x[i];
                    }
            }
    }


    void BicubicInterpolation(Vec_I_DP &y, Vec_I_DP &y1, Vec_I_DP &y2, Vec_I_DP &y12,
            const float x1l, const float x1u, const float x2l, const float x2u,
            const float X1, const float x2, float &ansy, float &ansy1, float &ansy2)
    {
            int i;
            float t,u,d1,d2;
            Mat_DP c(4,4);

            d1=x1u-x1l;
            d2=x2u-x2l;
            bcucof(y,y1,y2,y12,d1,d2,c);
            if (x1u == x1l || x2u == x2l)
                    nrerror("Bad input in routine BicubicInterpolation");
            t=(X1-x1l)/d1;
            u=(x2-x2l)/d2;
            ansy=ansy2=ansy1=0.0;
            for (i=3;i>=0;i--) {
                    ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
                    ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
                    ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
            }
            ansy1 /= d1;
            ansy2 /= d2;
    }




} // end namespace NR




#endif /* POLYNOMIAL_FITTING_H_ */
