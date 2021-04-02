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
//!	CommandLineManager allows to define flags, either as required, standard or debug.
//!	User command line input is tested for consistency.
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef COMMAND_LINE_MANAGER_H
#define COMMAND_LINE_MANAGER_H


#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <iterator>
#include "command_option.h"
#include "../string/string_functions.h"

//extern const std::string e_LevelString[];

class CommandLineManager
{
private:

    std::vector< CommandOption>                        m_Flags;
    std::map< std::string, std::vector< std::string> > m_FlagsWithOptions;
    std::string                                        m_GeneralIntroduction;

public:

    CommandLineManager()
    : m_Flags(),
    m_FlagsWithOptions(),
    m_GeneralIntroduction()
    { DefineFlag( "help", "TYPE", "TYPE can be \'standard\' or \'debug\' (default \'standard\'), writes this help", 0, 1);}


    CommandLineManager( const CommandLineManager& ORIGINAL)
    : m_Flags( ORIGINAL.m_Flags),
    m_FlagsWithOptions( ORIGINAL.m_FlagsWithOptions),
    m_GeneralIntroduction( ORIGINAL.m_GeneralIntroduction)
    {}


    virtual ~CommandLineManager(){}


    void DefineFlag( const std::string &STRING, const std::string &HELP, const e_Level &LEVEL)
    { m_Flags.push_back( CommandOption( STRING, "", HELP, 0, 0, LEVEL));}

    void DefineFlag( const std::string &STRING, const std::string &HELP, const size_t MINIMUM = 0, const e_Level &LEVEL = e_Basic)
    { m_Flags.push_back( CommandOption( STRING, "", HELP, MINIMUM, MINIMUM, LEVEL));}


    void DefineFlag( const std::string &STRING, const std::string &HELP, const size_t MINIMUM, const size_t MAXIMUM, const e_Level &LEVEL = e_Basic)
    { m_Flags.push_back( CommandOption( STRING, "", HELP, MINIMUM, MAXIMUM, LEVEL));}

    void DefineFlag( const std::string &STRING, const std::string &ARGUMENTS, const std::string &HELP, const e_Level &LEVEL)
    { m_Flags.push_back( CommandOption( STRING, ARGUMENTS, HELP, 0, 0, LEVEL));}

    void DefineFlag( const std::string &STRING, const std::string &ARGUMENTS, const std::string &HELP, const size_t MINIMUM = 0, const e_Level &LEVEL = e_Basic)
    { m_Flags.push_back( CommandOption( STRING, ARGUMENTS, HELP, MINIMUM, MINIMUM, LEVEL));}


    void DefineFlag( const std::string &STRING, const std::string &ARGUMENTS, const std::string &HELP, const size_t MINIMUM, const size_t MAXIMUM, const e_Level &LEVEL = e_Basic)
    { m_Flags.push_back( CommandOption( STRING, ARGUMENTS, HELP, MINIMUM, MAXIMUM, LEVEL));}



//    void DefineRequiredFlag( const std::string &STRING, const std::string &HELP, const size_t MINIMUM, const size_t MAXIMUM)
//    { m_Flags.push_back( CommandOption( STRING, HELP, MINIMUM, MAXIMUM, true));}
//
//    void DefineRequiredFlag( const std::string &STRING, const std::string &HELP, const size_t MINIMUM)
//    { m_Flags.push_back( CommandOption( STRING, HELP, MINIMUM, MINIMUM, true));}

    void SetGeneralIntroduction( const std::string &INTRO)
    { m_GeneralIntroduction = "About this program:\n" + INTRO;}

    void AppendToGeneralIntroduction( const std::string &INTRO)
    { m_GeneralIntroduction += INTRO;}

    std::string const &GetGeneralIntroduction() const
    { return m_GeneralIntroduction;}

    bool IsValidFlag( const std::string &STRING) const
    {
        for( std::vector< CommandOption>::const_iterator itr( m_Flags.begin()); itr != m_Flags.end(); ++itr)
            if( STRING == itr->GetIdentifier())
            { return true;}
        std::cout << STRING << " is invalid flag" << std::endl;
        return false;
    }

    bool IsFlagSet( const std::string &STRING) const
    {
        return(  IsValidFlag( STRING) &&  m_FlagsWithOptions.find( STRING) != m_FlagsWithOptions.end());
    }

    bool ReadAndCheckFlags( const size_t ARGC, const char * ARGV[])
    {
//        std::cout << m_GeneralIntroduction << std::endl;
        if( ARGC > 1)
        {
            std::string flag, prev_flag( ARGV[1]);
            std::vector< std::string> options;
            //    if( ARGC == 1)
            //      { WriteHelp();}
            // iterate through arguments

//            assert( prev_flag.substr( 0, 1) == "-");

            prev_flag = prev_flag.substr( 1,  prev_flag.size() - 1);
            //    std::cout << "reading flag: " << prev_flag << std::endl;

            for( size_t i = 2; i < ARGC; ++i)
            {
                std::string local( ARGV[i]);
                // if it is a flag:
                if( local.substr( 0,1) == "-" && !mystr::IsNumerical( local))
                {
                    // if it is not the first flag:
                    flag = local.substr( 1,  local.size() - 1);

                    if( Check( prev_flag, options))
                    { /*std::cout << "flag " << prev_flag << " inserted" << std::endl;*/ m_FlagsWithOptions.insert( std::make_pair( prev_flag, options));}
                    else{ std::cout << "incorrect flag/options" << std::endl; return false;}
                    /*std::cout << "reading flag: " << flag << std::endl;*/
                    prev_flag = flag;
                    options = std::vector< std::string>();
                }
                else
                { options.push_back( local);}
            }
            if( Check( prev_flag, options))
            { /*std::cout << "last flag " << prev_flag << " inserted" << std::endl;*/ m_FlagsWithOptions.insert( std::make_pair( prev_flag, options));}
            else{ std::cout << "incorrect flag/options" << std::endl; return false;}

            if( IsFlagSet("help"))
            {
            	std::vector< std::string> str = this->GetAllOptions( "help");
            	if( str.size() == 0 || ( str.size() == 1 && str[0] == "standard"))
            	{
            		WriteHelp( std::cout, e_Basic);
            	}
            	else if( ( str.size() == 1 && str[0] == "debug"))
            	{
            		WriteHelp( std::cout, e_Advanced);
            	}
                exit( -1);
            }

            WriteFlagsWithOptions( std::cout);

            return CheckForRequired();
        }
        return false;
    }



    bool Check( const std::string &FLAG, const std::vector< std::string> &OPTIONS) const
    {

        for( std::vector< CommandOption>::const_iterator itr = m_Flags.begin(); itr != m_Flags.end(); ++itr)
        {
            //    std::cout << "check flag <" << FLAG << "> with allowed option: " << itr->GetIdentifier() << std::endl;
            if( FLAG == itr->GetIdentifier())
            {
                //        std::cout << "identified" << std::endl;
                if( OPTIONS.size() < itr->GetMinimumNumberOfParameters() || OPTIONS.size() > itr->GetMaximumNumberOfParameters())
                {
                	std::cout << "recognized flag: <" << FLAG << ">, wrong number of arguments: " << OPTIONS.size() << " should be between (" << itr->GetMinimumNumberOfParameters() << " " << itr->GetMaximumNumberOfParameters() << ")" << std::endl;
                	std::cout << "given arguments: < ";
                	std::copy( OPTIONS.begin(), OPTIONS.end(), std::ostream_iterator<std::string>( std::cout, " "));
                	std::cout << ">" << std::endl << std::endl;
                	return false;
                }
                else{ return true;}
            }
        }
        std::cout << "unrecognized flag: " << FLAG << std::endl;
        return false;
    }


    bool CheckForRequired() const
    {
        for( std::vector< CommandOption>::const_iterator itr = m_Flags.begin(); itr != m_Flags.end(); ++itr)
        {
            if( itr->IsRequired() && !IsFlagSet( itr->GetIdentifier()))
            {
                std::cout << "required flag " << itr->GetIdentifier() << " is not given!" << std::endl;
                return false;
            }
        }
        return true;
    }


    template< typename T>
    T GetArgumentForFlag( const std::string &STRING) const
    {
        assert( IsFlagSet( STRING));
        std::vector< std::string> options( m_FlagsWithOptions.find( STRING)->second);
        assert( options.size() == 1);
        std::istringstream iss( options[0]);
        T value;
        iss >> value;
        return value;
    }

    template< typename T>
    std::vector< T> GetArgumentsForFlag( const std::string &STRING) const
    {
        assert( IsFlagSet( STRING));
        std::vector< std::string> options( m_FlagsWithOptions.find( STRING)->second);
        std::vector< T> values( options.size());
        typename std::vector< T>::iterator itr( values.begin());
        for( std::vector< std::string>::const_iterator str_itr = options.begin(); str_itr != options.end(); ++itr, ++str_itr)
        {
            std::istringstream iss( *str_itr);
            iss >> *itr;
        }
        return values;
    }


    std::string GetArgumentStringForFlag( const std::string &STRING)  const
    {
        if( !IsFlagSet( STRING))
        {
        	std::cout << __FUNCTION__ << " flag is not defined: " << STRING << std::endl;
        }
        std::vector< std::string> options( m_FlagsWithOptions.find( STRING)->second);
        assert( options.size() == 1);
        return options[0];
    }


    std::vector< std::string> GetAllOptions( const std::string &STRING) const
    {
        assert( IsFlagSet( STRING));
        return m_FlagsWithOptions.find( STRING)->second;
    }

    std::ostream& WriteFlagsWithOptions( std::ostream &STREAM) const
    {
        STREAM << "Following arguments were given:" << std::endl;

        //      std::pair< std::vector< std::string>, std::vector< std::string> > required( GetRequiredAndOptionalFlags());

        for( std::map< std::string, std::vector< std::string> >::const_iterator itr = m_FlagsWithOptions.begin(); itr != m_FlagsWithOptions.end(); ++itr)
        {
            STREAM << "-";
            STREAM.width( 30);
            STREAM.setf( std::ios::left);
            STREAM << itr->first;
	    if( itr->second.size() != 0)
	    {
		STREAM << ":      <";
	    }
            for( std::vector< std::string>::const_iterator opt_itr = itr->second.begin(); opt_itr != itr->second.end(); ++opt_itr)
            {
                STREAM << *opt_itr << ">  ";
                if( opt_itr + 1 != itr->second.end()){ STREAM << "<";}
            }
            STREAM << std::endl;
        }
        return STREAM;
    }

    std::ostream& WriteHelp( std::ostream &STREAM, const e_Level &LEVEL = e_Advanced) const
    {
//        STREAM << std::endl << m_GeneralIntroduction << std::endl << std::endl;

        STREAM << "========================  HELP  ===========================\n" << std::endl;

        int min = std::min( e_NrLevels - 1, (int) LEVEL);
        for( int i = 0; i <= min; ++i)
        {
        	std::string str = e_LevelString[ i];

        	if( i == 0)
        	{
        		STREAM << "\n============  required flags  ============\n" << std::endl;
        	}
           	else if( i == 1)
           	{
               		STREAM << "\n============  optional flags: standard  ============\n" << std::endl;
           	}
           	else if( i == 2)
           	{
               		STREAM << "\n============  optional flags: debug/development  ============\n" << std::endl;
           	}

			for( std::vector< CommandOption>::const_iterator itr = m_Flags.begin(); itr != m_Flags.end(); ++itr)
			{
				if( itr->GetLevel() == i)
				{
					STREAM << "-" << itr->GetIdentifier() << "    " << itr->GetArgument() << std::endl;
					STREAM << "     " << itr->GetHelp() << "" << std::endl;
					if( itr->GetMinimumNumberOfParameters() ==  itr->GetMaximumNumberOfParameters())
					{
					    STREAM << "     " << itr->GetMinimumNumberOfParameters();
					    if( itr->GetMinimumNumberOfParameters() == 1)
					    {
					    	STREAM << " parameter" << std::endl;
					    }
					    else
					    {
					    	STREAM << " parameters" << std::endl;
					    }
					}
					else
					{
					    STREAM << "     from " << itr->GetMinimumNumberOfParameters() << " to " << itr->GetMaximumNumberOfParameters() << " parameters" << std::endl;
					}
					STREAM << std::endl;
				}
			}
        }
        STREAM << "\n=====================  HELP done  ========================\n" << std::endl;
        return STREAM;
    }

    size_t GetNumberOfSetFlags() const
    { return m_FlagsWithOptions.size();}

    size_t GetNumberOfRequiredFlags() const
    {
        size_t nr;
        for( std::vector< CommandOption>::const_iterator itr = m_Flags.begin(); itr != m_Flags.end(); ++itr)
        {
            if( itr->IsRequired())
            { ++nr;}
        }
        return nr;
    }

    std::ostream &IfNoFlagIsSetWriteHelpAndExit( const size_t &ARGC, std::ostream &STREAM) const
    {
        if( ARGC == 1)
        {
            WriteHelp( std::cout, e_Basic);
            std::cout << std::endl << "\n=====================  No flags given!  =====================\n" << std::endl << std::endl;
            exit( -1);
        }
        return STREAM;
    }

}; // end class CommandLineManager

#endif // COMMAND_LINE_MANAGER_H
