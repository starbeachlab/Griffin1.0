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


#include "../../include/utilities/time.h"

double
TimeDiff
(
		struct timeval *DIFFERENCE,
		const struct timeval *END,
		const struct timeval *START
)
{
  struct timeval tmp_diff;

  if( DIFFERENCE == NULL)
  {
	  DIFFERENCE=&tmp_diff;
  }

  DIFFERENCE->tv_sec = END->tv_sec - START->tv_sec ;
  DIFFERENCE->tv_usec= END->tv_usec - START->tv_usec;

  /* Using while instead of if below makes the code slightly more robust. */

  while( DIFFERENCE->tv_usec < 0)
  {
	  DIFFERENCE->tv_usec += 1000000;
	  DIFFERENCE->tv_sec -= 1;
  }

  return double( 1000000LL * DIFFERENCE->tv_sec + DIFFERENCE->tv_usec) / 1.0e6;

} /* timeval_diff() */

double
TimeSum
(
		struct timeval *SUM,
		const struct timeval *END,
		const struct timeval *START
)
{
  struct timeval tmp_sum;

  if( SUM == NULL)
  {
	  SUM=&tmp_sum;
  }

  SUM->tv_sec = END->tv_sec + START->tv_sec ;
  SUM->tv_usec= END->tv_usec + START->tv_usec;

  return double( 1000000LL * SUM->tv_sec + SUM->tv_usec) / 1.0e6;

} /* timeval_sum() */


double
Seconds
(
		const struct timeval *TIME
)
{
  return double( 1000000LL * TIME->tv_sec + TIME->tv_usec) / 1.0e6;
} 




void SumTimer::Reset()
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	m_Sum.tv_sec = 0;
	m_Sum.tv_usec = 0;
//	m_Own.tv_sec = 0;
//	m_Own.tv_usec = 0;
//	m_Started = false;
//	m_Calls = 0;
}


void SumTimer::Start()
{
	gettimeofday( &m_Start, NULL);
//	if( m_Started)
//	{
//		std::cout << __PRETTY_FUNCTION__ << " called on already started timer! (" << m_ID << ")" << std::endl;
//	}
//	m_Started = true;
}


void SumTimer::Stop()
{

	gettimeofday( &m_End, NULL);
	TimeDiff( &m_Diff, &m_End, &m_Start);
	TimeSum( &m_Sum, &m_Diff, &m_Sum);
/*
	gettimeofday( &m_Diff, NULL);
	TimeDiff( &m_Diff, &m_Diff, &m_End);
	TimeSum( &m_Own, &m_Diff, &m_Own);
	++m_Calls;

	if( !m_Started)
	{
		std::cout << __PRETTY_FUNCTION__ << " called on non-started timer! (" << m_ID << ")" << std::endl;
	}
	m_Started = false;
*/
		//	m_Sum = m_Diff;
}


double SumTimer::Status() const
{
	return Seconds( &m_Sum);
}


std::ostream & SumTimer::WriteDiff( std::ostream &STREAM) const
{
    STREAM << m_ID << ": diff: " << Seconds( &m_Diff) << "s" << std::endl;
    return STREAM;
}


std::ostream & SumTimer::WriteStatus( std::ostream &STREAM) const
{
    STREAM << m_ID << ": " << Seconds( &m_Sum) << " s" << std::endl;
//    STREAM << "own: " << Seconds( &m_Own) << " s, #calls: " << m_Calls << std::endl;
    return STREAM;
}




void ContinuousTimer::Reset()
{
//	std::cout << __PRETTY_FUNCTION__ << std::endl;
	m_Sum.tv_sec = 0;
	m_Sum.tv_usec = 0;
}


void ContinuousTimer::Start()
{
	gettimeofday( &m_Start, NULL);
	m_Last = m_Start;
}


void ContinuousTimer::Total()
{
	gettimeofday( &m_End, NULL);
	TimeDiff( &m_Diff, &m_End, &m_Start);
//	TimeSum( &m_Sum, &m_Diff, &m_Sum);
	m_Sum = m_Diff;
	m_Last = m_End;
}

void ContinuousTimer::Intermediate()
{
	gettimeofday( &m_End, NULL);
	TimeDiff( &m_Diff, &m_End, &m_Last);
	TimeSum( &m_Sum, &m_Diff, &m_Sum);
	m_Last = m_End;
}


double ContinuousTimer::Status() const
{
	return Seconds( &m_Sum);
}


std::ostream & ContinuousTimer::WriteDiff( std::ostream &STREAM) const
{
    STREAM << m_ID << ": diff: " << Seconds( &m_Diff) << "s" << std::endl;
    return STREAM;
}


std::ostream & ContinuousTimer::WriteStatus( std::ostream &STREAM) const
{
    STREAM << m_ID << ": " << Seconds( &m_Sum) << " s" << std::endl;
    return STREAM;
}




