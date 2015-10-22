// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "HashPriorityQueue.h"

void HashPriorityQueue::clear()
{
    if ( use64 )
    {
        while ( queue64.size() )
        {
            queue64.pop();
        }
    }
    else
    {
        while ( queue32.size() )
        {
            queue32.pop();
        }
    }
}

hash_u HashPriorityQueue::top() const
{
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = queue64.top();
    }
    else
    {
        hash.hash32 = queue32.top();
    }
    
    return hash;
}
