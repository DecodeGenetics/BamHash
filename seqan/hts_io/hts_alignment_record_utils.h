#ifndef INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_
#define INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_

#include <seqan/bam_io/bam_alignment_record.h>
#include <seqan/hts_io/hts_file.h>

namespace seqan {


CharString inline
getContigName(int32_t const rID, HtsFile const & file)
{
    SEQAN_ASSERT(file.hdr);

    if (rID < 0 || rID >= file.hdr->n_targets)
        return CharString("");

    return CharString(file.hdr->target_name[rID]);
}


} // namespace seqan

#endif // INCLUDE_SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_UTILS_H_
