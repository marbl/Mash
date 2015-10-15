// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include <iostream>

#include "CommandList.h"
#include "version.h"
#include <string.h>

using namespace::std;

CommandList::CommandList(string nameNew)
{
    name = nameNew;
}

CommandList::~CommandList()
{
    for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
    {
        delete i->second;
    }
}

void CommandList::addCommand(Command * command)
{
    commands[command->name] = command;
}

void CommandList::print()
{
    vector<vector<string>> columns(1);
    
    cout << endl << "Mash version " << version << endl << endl;
    
    cout << "Copyright © 2015, Battelle National Biodefense Institute (BNBI); all rights" << endl;
	cout << "reserved. Authored by: Brian Ondov, Todd Treangen, and Adam Phillippy." << endl;
	
    cout << "This program is free software and comes with ABSOLUTELY NO WARRANTY, though you" << endl;
    cout << "are welcome to redistribute it under certain conditions. Source code is" << endl;
    cout << "available at github.com/marbl/mash. For more details, type 'mash --license'." << endl;
    
    cout << endl << "Usage:" << endl << endl;
    
    columns[0].push_back(name + " <command> [options] [arguments ...]");
    printColumns(columns);
    
    cout << "Commands:" << endl << endl;
    
    int lengthMax = 0;
    
    for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
    {
        if ( i->first.length() > lengthMax )
        {
            lengthMax = i->first.length();
        }
    }
    
    columns.clear();
    columns.resize(2);
    
    for ( map<string, Command *>::iterator i = commands.begin(); i != commands.end(); i++ )
    {
        columns[0].push_back(i->first);
        columns[1].push_back(i->second->summary);
    }
    
    printColumns(columns);
}

int CommandList::run(int argc, const char ** argv)
{
	if ( argc > 1 && strcmp(argv[1], "--version") == 0 )
	{
		cout << version << endl;
		return 0;
	}
	
	if ( argc > 1 && strcmp(argv[1], "--license") == 0 )
	{
		showLicense();
		return 0;
	}
	
    if ( argc < 2 || commands.count(argv[1]) == 0 )
    {
        print();
        return 0;
    }
    
    return commands.at(argv[1])->run(argc - 2, argv + 2);
}

void CommandList::showLicense()
{
	cout << "\n\
PURPOSE\n\
\n\
Mash is a fast sequence distance estimator that uses the MinHash\n\
algorithm and is designed to work with genomes and metagenomes in the\n\
form of assemblies or reads. It is implemented in C++ and is\n\
distributed with KSeq (lh3lh3.users.sourceforge.net/kseq.shtml) and\n\
MurmurHash3 (code.google.com/p/smhasher/wiki/MurmurHash3).\n\
\n";
#ifdef DIST_LICENSE
	cout << "\n\
This binary includes binary redistributions of software covered by\n\
the Boost and MIT licenses, included below.\n\
\n";

#endif

cout << "\n\
COPYRIGHT LICENSE\n\
\n\
Copyright © 2015, Battelle National Biodefense Institute (BNBI);\n\
all rights reserved. Authored by: Brian Ondov, Todd Treangen, and\n\
Adam Phillippy\n\
\n\
This Software was prepared for the Department of Homeland Security\n\
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as\n\
part of contract HSHQDC-07-C-00020 to manage and operate the National\n\
Biodefense Analysis and Countermeasures Center (NBACC), a Federally\n\
Funded Research and Development Center.\n\
\n\
Redistribution and use in source and binary forms, with or without\n\
modification, are permitted provided that the following conditions are\n\
met:\n\
\n\
1. Redistributions of source code must retain the above copyright\n\
notice, this list of conditions and the following disclaimer.\n\
\n\
2. Redistributions in binary form must reproduce the above copyright\n\
notice, this list of conditions and the following disclaimer in the\n\
documentation and/or other materials provided with the distribution.\n\
\n\
3. Neither the name of the copyright holder nor the names of its\n\
contributors may be used to endorse or promote products derived from\n\
this software without specific prior written permission.\n\
\n\
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n\
\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n\
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n\
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n\
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n\
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n\
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n\
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n\
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n\
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n\
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\
\n";

#ifdef DIST_LICENSE
	cout << "\n\
Boost Software License - Version 1.0 - August 17th, 2003\n\
\n\
Permission is hereby granted, free of charge, to any person or organization\n\
obtaining a copy of the software and accompanying documentation covered by\n\
this license (the \"Software\") to use, reproduce, display, distribute,\n\
execute, and transmit the Software, and to prepare derivative works of the\n\
Software, and to permit third-parties to whom the Software is furnished to\n\
do so, all subject to the following:\n\
\n\
The copyright notices in the Software and this entire statement, including\n\
the above license grant, this restriction and the following disclaimer,\n\
must be included in all copies of the Software, in whole or in part, and\n\
all derivative works of the Software, unless such copies or derivative\n\
works are solely in the form of machine-executable object code generated by\n\
a source language processor.\n\
\n\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n\
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT\n\
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE\n\
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,\n\
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER\n\
DEALINGS IN THE SOFTWARE.\n\
\n\
\n\
MIT License\n\
\n\
Permission is hereby granted, free of charge, to any person obtaining a copy\n\
of this software and associated documentation files (the \"Software\"), to deal\n\
in the Software without restriction, including without limitation the rights\n\
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n\
copies of the Software, and to permit persons to whom the Software is\n\
furnished to do so, subject to the following conditions:\n\
\n\
The above copyright notice and this permission notice shall be included in\n\
all copies or substantial portions of the Software.\n\
\n\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n\
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n\
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n\
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n\
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n\
THE SOFTWARE.\n\
\n";
#endif
}
