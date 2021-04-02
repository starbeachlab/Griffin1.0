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


#ifndef COMMAND_OPTION_H
#define COMMAND_OPTION_H

//#include <sstream>
//#include <map>
//#include <vector>
#include <string>
//#include <cassert>

enum e_Level { e_Required, e_Basic, e_Advanced};

const std::string e_LevelString[3] = { "required", "basic", "advanced"};

const int e_NrLevels = 3;



class CommandOption
{
 private:
  std::string      m_Identifier;
  std::string      m_Argument;
  std::string      m_Help;
  size_t           m_MinimumNumberOfArguments;
  size_t           m_MaximumNumberOfArguments;
  e_Level    	   m_Level;

 public:

  CommandOption
    (
     const std::string &IDENTIFIER,
     const std::string &ARGUMENT,
     const std::string &HELP,
     const size_t MINIMUM = 0,
     const size_t MAXIMUM = 0,
     const e_Level &REQUIRED = e_Basic
    )
    : m_Identifier( IDENTIFIER),
      m_Argument( ARGUMENT),
      m_Help( HELP),
      m_MinimumNumberOfArguments( MINIMUM),
      m_MaximumNumberOfArguments( MAXIMUM),
      m_Level( REQUIRED)
    {}

    virtual ~CommandOption(){}

    std::string const &GetIdentifier() const
    { return m_Identifier;}

    std::string const& GetHelp() const
    { return m_Help;}

    std::string const& GetArgument() const
	{ return m_Argument;}

    bool IsRequired() const
    { return m_Level == e_Required;}

    e_Level GetLevel() const
    {
    	return m_Level;
    }

    size_t GetMinimumNumberOfParameters() const
    { return m_MinimumNumberOfArguments;}

    size_t GetMaximumNumberOfParameters() const
    { return m_MaximumNumberOfArguments;}

};

#endif // COMMAND_OPTION_H

