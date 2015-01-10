#include <string>
#include <sstream>
#include <openssl/md5.h>
#include <stdint.h>

#include "bamhash_checksum_common.h"



hash_t str2md5(const char *str, int length)
{
	hash_t out;
	MD5((unsigned char*)str, length, (unsigned char*)(out.c));
	return out;
}

void hexSum(hash_t out, uint64_t & sum)
{
	sum += out.p.low;
}

