
#include "base.h"
#include "HashList.h"

extern "C" {

MASH_EXTERN_API int mash_create_hash_list(HashList* hlptr, int kmerSize) {
    
    hlptr = new HashList(kmerSize);
    return 0;

}

MASH_EXTERN_API int mash_hash_list_size(HashList* hlptr, int *size) {
    
    *size = hlptr->size();
    return 0;
    
}

} // extern "C"

