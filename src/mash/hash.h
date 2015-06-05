#ifndef hash_h
#define hash_h

#include <inttypes.h>

static const int seed = 42; // TODO: better seed???

typedef uint32_t hash32_t;
typedef uint64_t hash64_t;

union hash_u
{
    hash32_t hash32;
    hash64_t hash64;
};

hash_u getHash(const char * seq, int length);
bool hashLessThan(hash_u hash1, hash_u hash2, bool use64);

#endif
