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


#ifndef TIME_H_
#define TIME_H_

#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>

double
TimeDiff
(
		struct timeval *DIFFERENCE,
		const struct timeval *END,
		const struct timeval *START
);

double
TimeSum
(
		struct timeval *SUM,
		const struct timeval *END,
		const struct timeval *START
);

double
Seconds
(
		const struct timeval *TIME
);


class Timer
{
private:
	struct timeval m_Start;
	struct timeval m_End;
	std::string    m_ID;
	std::ostream   &m_Stream;
//	struct timeval m_Own;

public:
	Timer( const std::string &STR = "", std::ostream &STREAM = std::cout)
	: m_Start(),
	  m_End(),
	  m_ID( STR),
	  m_Stream( STREAM)
//	, m_Own()
	{
		gettimeofday( &m_Start, NULL);
	}

	~Timer()
	{
		gettimeofday( &m_End, NULL);
		m_Stream << m_ID << ": " << TimeDiff( NULL, &m_End, &m_Start) << " s" << std::endl;
//		gettimeofday( &m_Own, NULL);
//		m_Stream << "own: " << TimeDiff( NULL, &m_Own, &m_End) << " s" << std::endl;
	}

};

class SumTimer
{
private:
	struct timeval m_Start;
	struct timeval m_End;
	struct timeval m_Diff;
	struct timeval m_Sum;
	std::string    m_ID;
//	struct timeval m_Own;
//	int            m_Calls;
//	bool           m_Started;

public:
	SumTimer( const std::string &STR)
	: m_Start(),
	  m_End(),
	  m_Diff(),
	  m_Sum(),
	  m_ID( STR) 
//	  ,m_Own(), m_Calls(0),
//	  m_Started( false)
	{}

	SumTimer( const SumTimer &TIM)
        : m_Start( TIM.m_Start),
	  m_End( TIM.m_End),
	  m_Diff( TIM.m_Diff),
	  m_Sum( TIM.m_Sum),
	  m_ID( TIM.m_ID)
//	,m_Own( TIM.m_Own), m_Calls( TIM.m_Calls), m_Started( TIM.m_Started)
	{}

	~SumTimer()
	{}

	void Reset();

	void Start();

	void Stop();

	double Diff() const;

	double Status() const;

	std::ostream & WriteDiff( std::ostream &STREAM = std::cout) const;

	std::ostream & WriteStatus( std::ostream &STREAM = std::cout) const;

};


class ContinuousTimer
{
private:
	struct timeval m_Start;
	struct timeval m_End;
	struct timeval m_Diff;
	struct timeval m_Sum;
	struct timeval m_Last;
	std::string    m_ID;

public:
	ContinuousTimer( const std::string &STR)
	: m_Start(),
	  m_End(),
	  m_Diff(),
	  m_Sum(),
	    m_Last(),
	  m_ID( STR)
	{}

	ContinuousTimer( const ContinuousTimer &TIM)
        : m_Start( TIM.m_Start),
	  m_End( TIM.m_End),
	  m_Diff( TIM.m_Diff),
	  m_Sum( TIM.m_Sum),
	    m_Last( TIM.m_Last),
	  m_ID( TIM.m_ID)
	{}

	~ContinuousTimer()
	{}

	void Reset();

	void Start();

	void Intermediate();

	void Total();

	double Diff() const;

	double Status() const;

	std::ostream & WriteDiff( std::ostream &STREAM = std::cout) const;

	std::ostream & WriteStatus( std::ostream &STREAM = std::cout) const;

};


#endif /* TIME_H_ */
