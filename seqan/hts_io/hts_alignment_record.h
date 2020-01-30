#ifndef SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_
#define SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_

#include <seqan/basic.h>

#include <cstdio>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

namespace seqan {

class HtsSequenceRecord
{
  public:
    String<char> qName;
    String<Iupac> seq;
    String<char> qual;

    static __int32 const INVALID_POS = -1;
    static __int32 const INVALID_REF_ID = -1;
    static __int32 const INVALID_LEN = 0;
    static __uint32 const INVALID_QID = 4294967295u;

    HtsSequenceRecord()
      : qName(), seq(), qual() {}

    HtsSequenceRecord(bam1_t * hts_record)
    {
        HtsSequenceRecord::parse(hts_record);
    }

    void parse(bam1_t * hts_record)
    {
        qName = bam_get_qname(hts_record);
        int32_t lqseq = hts_record->core.l_qseq;
        resize(seq, lqseq, Exact());
        uint8_t* seqptr = bam_get_seq(hts_record);

        for (int i = 0; i < lqseq; ++i)
        {
          seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
        }

        uint8_t* qualptr = bam_get_qual(hts_record);
        resize(qual, lqseq, Exact());

        for (int i = 0; i < lqseq; ++i, ++qualptr)
        {
            qual[i] = static_cast<char>(*qualptr + 33);
        }
    }
};

} // namespace seqan

#endif // SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_
