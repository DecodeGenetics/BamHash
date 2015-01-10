#ifndef BAMHASH_CHECKSUM_COMMON_H
#define BAMHASH_CHECKSUM_COMMON_H

#define BAMHASH_VERSION "1.0"

#include <string>
#include <stdint.h>

union hash_t {
  unsigned char c[16];
  struct {
    uint64_t low;
    uint64_t high;
  } p;
};

hash_t str2md5(const char *str, int length);
void hexSum(hash_t out, uint64_t & sum);


#endif // BAMHASH_CHECKSUM_COMMON_H
