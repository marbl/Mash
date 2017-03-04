#ifndef MASH_HASH_LIST_C_H
#define MASH_HASH_LIST_C_H

#include "../mash/HashList.h"

extern "C" {

MASH_EXTERN_API int mash_create_hash_list(HashList* hlptr, int kmerSize);

MASH_EXTERN_API int mash_hash_list_size(HashList* hlptr, int *size);

} // extern "C"

#endif /* ifndef MASH_HASH_LIST_C_H */
