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
//!	Writing headers, logos, system checks and quotes.
//!									 
//!									 
//! @author: Rene Staritzbichler					 
//! @date:   1.2.2011						 
/////////////////////////////////////////////////////////////////////////


#ifndef QUOTES_H_
#define QUOTES_H_

#include <cstdlib>

inline
std::string
ThrowQuote()
{
    std::vector< std::string>
        quote;

    quote.push_back( "Arthur Schopenhauer:\n\nEvery nation ridicules other nations, and all are right.");
    quote.push_back( "Arthur Schopenhauer:\n\nThe difficulty is to try and teach the multitude that something can be true and untrue at the same time.");
    quote.push_back( "Arthur Schopenhauer:\n\nMartyrdom is the only way a man can become famous without ability.");
    quote.push_back( "Arthur Schopenhauer:\n\nAll truth passes through three stages. First, it is ridiculed. Second, it is violently opposed. Third, it is accepted as being self-evident.");
    quote.push_back( "Arthur Schopenhauer:\n\nChange alone is eternal, perpetual, immortal.");
    quote.push_back( "Arthur Schopenhauer:\n\nEvery man takes the limits of his own field of vision for the limits of the world.");
    quote.push_back( "Arthur Schopenhauer:\n\nThe doctor sees all the weakness of mankind; the lawyer all the wickedness, the theologian all the stupidity.");
    quote.push_back( "Friedrich Nietzsche:\n\nA subject for a great poet would be God's boredom after the seventh day of creation.");
//    quote.push_back( "Friedrich Nietzsche:\n\nAh, women. They make the highs higher and the lows more frequent.");
    quote.push_back( "Friedrich Nietzsche:\n\nAnd we should consider every day lost on which we have not danced at least once. And we should call every truth false which was not accompanied by at least one laugh.");
    quote.push_back( "Michel de Montaigne:\n\nA man who fears suffering is already suffering from what he fears.");
    quote.push_back( "Michel de Montaigne:\n\nConfidence in the goodness of another is good proof of one's own goodness.");
    quote.push_back( "Dwight D. Eisenhower:\n\nAn intellectual is a man who takes more words than necessary to tell more than he knows.");
    quote.push_back( "Friedrich Nietzsche:\n\nA casual stroll through the lunatic asylum shows that faith does not prove anything.");
    quote.push_back( "Moliere:\n\nDon't appear so scholarly (..). Humanize your talk, and speak to be understood.");
    quote.push_back( "Ernest Hemingway:\n\nHappiness in intelligent people is the rarest thing I know.");
    quote.push_back( "Albert Schweitzer:\n\nHappiness is nothing more than good health and a bad memory.");
    quote.push_back( "Friedrich Nietzsche:\n\nThe advantage of a bad memory is that one enjoys several times the same good things for the first time.");
    quote.push_back( "Janis Joplin:\n\nFreedom's just another word for nothing left to lose.");
//    quote.push_back( "Dr. Seuss:\n\nWhen a fox is in the bottle where the tweetle beetles battle with their paddles in a puddle on a noodle-eating poodle, THIS is what they call ... a tweetle beetle noodle poodle bottled paddled muddled fuddled wuddled fox in socks, sir!");
    quote.push_back( "Dr. Seuss:\n\nThank goodness for all of the things you are not! Thank goodness you're not something someone forgot, and left all alone in some punkerish place like a rusty tin coat hanger in space.\nThat's why I say, \"Duckie! Don't grumble! Don't stew! Some critters are much-much, oh, ever so much-much, so muchly much-much more unlucky than you!\"");
    quote.push_back( "Immanuel Kant: \n\nIt is beyond a doubt that all our knowledge begins with experience.");
    quote.push_back( "Immanuel Kant: \n\nI can because I want what I have to.");
    quote.push_back( "Malcom Muggeridge: \n\nOnly dead fish swim with the stream.");

    std::srand( time( NULL));

    size_t
        randy( std::rand() % quote.size());
    size_t
        first( 3 + quote[ randy].find( ":"));
    std::string
        top( quote[ randy].size() - first, '=');

    return "\n\n\n" + top + "\n\n" + quote[ randy] + "\n\n" +  top + "\n\n\n\n";
}

inline
std::string
ThrowTaraLogo()
{
    std::string intro;
    intro  = "\n==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "================                          ================\n";
    intro +=   "================         Tara 1.0         ================\n";
    intro +=   "================                          ================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "=========================        =========================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    return intro;
}

inline
std::string
ThrowGriffinLogo()
{
    std::string intro;
    intro  = "\n==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "==================                       =================\n";
    intro +=   "=================       Griffin 1.0       ================\n";
    intro +=   "================                           ===============\n";
    intro +=   "================        ===========        ===============\n";
    intro +=   "================        ===========        ===============\n";
    intro +=   "================        ==================================\n";
    intro +=   "================        ==================================\n";
    intro +=   "================        ==================================\n";
    intro +=   "================        ======             ===============\n";
    intro +=   "================        ======             ===============\n";
    intro +=   "================        ===========        ===============\n";
    intro +=   "================        ===========        ===============\n";
    intro +=   "================         =========         ===============\n";
    intro +=   "=================                         ================\n";
    intro +=   "==================                       =================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    intro +=   "==========================================================\n";
    return intro;
}



inline
std::string
ThrowSystemTest()
{
    std::string intro( "\n  system test:\n");
    intro += "  USER = " + SystemVariable( "USER") + "\n";
    intro += "  HOSTNAME = " + SystemVariable( "HOSTNAME") + "\n";
    intro += "  PWD = " + SystemVariable( "PWD") + "\n";
    return intro;
}

inline
std::string
ThrowTaraGridAuthors()
{
    std::string intro( "\n");
    intro += std::string("  GRId-based Force-Field INput for molecular dynamics simulations\n");
    intro += std::string("  Authors: Rene Staritzbichler, Jose D. Faraldo-Gomez, Lucy R. Forrest\n");
    intro += std::string("  Designed and programmed by Rene Staritzbichler in C++\n");
    intro += std::string("  Please cite: R. Staritzbichler, C.Anselmi, L. Forrest, J. Faraldo-Gomez; J.Chem.Theor.Comp., 2011\n");
//    intro += std::string("  Please cite: R. Staritzbichler, C.Anselmi, L. Forrest, J. Faraldo-Gomez; \"GRIFFIN: A versatile methodology for optimization of protein-lipid interfaces for membrane protein simulations\", J.Chem.Theor.Comp., 2011\n");
    intro += std::string("  Max Planck Institute of Biophysics, Frankfurt, Germany\n");
    intro += std::string("  Source code is distributed under GNU general public license agreement\n");
    intro += std::string("  Compiled on: ") + std::string(  __DATE__) + std::string( "\n");
    return intro;
}

inline
std::string
ThrowFullHeader()
{
    return ThrowGriffinLogo() + ThrowTaraGridAuthors() + ThrowSystemTest(); // + ThrowQuote();
}


#endif /* QUOTES_H_ */
