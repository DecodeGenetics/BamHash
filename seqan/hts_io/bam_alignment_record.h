#ifndef INCLUDE_SEQAN_HTS_IO_BAM_ALIGNMENT_RECORD_H_
#define INCLUDE_SEQAN_HTS_IO_BAM_ALIGNMENT_RECORD_H_

#include <seqan/bam_io/bam_alignment_record.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

namespace seqan
{

inline void
syncCores(bam1_t * hts_record, BamAlignmentRecord const & record)
{
  hts_record->core.tid     = record.rID;
  hts_record->core.pos     = record.beginPos;
  hts_record->core.l_qname = record._l_qname;
  hts_record->core.qual    = record.mapQ;
  hts_record->core.bin     = record.bin;
  hts_record->core.n_cigar = record._n_cigar;
  hts_record->core.flag    = record.flag;
  hts_record->core.l_qseq  = record._l_qseq;
  hts_record->core.mtid    = record.rNextId;
  hts_record->core.mpos    = record.pNext;
  hts_record->core.isize   = record.tLen;
}

inline void
syncCores(BamAlignmentRecord & record, bam1_t * hts_record)
{
  record.rID      = hts_record->core.tid;
  record.beginPos = hts_record->core.pos;
  record._l_qname = hts_record->core.l_qname;
  record.mapQ     = hts_record->core.qual;
  record.bin      = hts_record->core.bin;
  record._n_cigar = hts_record->core.n_cigar;
  record.flag     = hts_record->core.flag;
  record._l_qseq  = hts_record->core.l_qseq;
  record.rNextId  = hts_record->core.mtid;
  record.pNext    = hts_record->core.mpos;
  record.tLen     = hts_record->core.isize;
}

inline std::vector<std::string>
split(std::string const & s, char const & delimiter)
{
  std::vector<std::string> splitted_str(0);
  std::string item;
  std::size_t pos = 0;

  while (true)
  {
    std::size_t new_pos = s.find(delimiter, pos);
    std::string item(s, pos, new_pos - pos);
    // std::cerr << "pos = " << pos << " " << " new_pos = " << new_pos << std::endl;
    // std::cerr << "item = '" << item << "'" << std::endl;

    splitted_str.push_back(item);

    if (new_pos == std::string::npos)
    {
      break;
    }
    else
    {
      pos = new_pos + 1;
    }
  }

  return splitted_str;
}

inline std::string
join(std::vector<std::string> const & splitted_str, char const & delimiter)
{
  // std::string s;
  std::stringstream ss;

  ss << splitted_str.front();

  for (auto it = splitted_str.begin()+1; it != splitted_str.end(); ++it)
  {
    ss << delimiter << *it;
  }

  return ss.str();
}


template<typename TTagType>
inline void
printTag(TTagType & read_int, CharString const & tags, unsigned & i, std::stringstream & ss)
{
  memcpy(&read_int, &tags[i], sizeof(TTagType));
  ss << static_cast<int64_t>(read_int);
  i += sizeof(TTagType);
}


inline std::string
toString(BamAlignmentRecord const & record, bam_hdr_t * hdr)
{
  std::stringstream ss;

  if (!hdr)
  {
    std::cerr << "ERROR: No header. Did you forget to read header?" << std::endl;
    std::exit(1);
  }
  else if (record.rID >= hdr->n_targets)
  {
    std::cerr << "ERROR: Invalid contig in BamAlignmentRecord with rID " << record.rID << std::endl;
    std::exit(1);
  }

  std::string const chrom = record.rID == -1 ? "*" : hdr->target_name[record.rID];
  ss << record.qName << '\t' << record.flag << '\t' << chrom << '\t' << record.beginPos+1 << '\t' << record.mapQ << '\t';

  if (length(record.cigar) > 0)
  {
    for (auto print_it = begin(record.cigar); print_it != end(record.cigar); ++ print_it)
      ss << print_it->count << print_it->operation;
  }
  else
  {
    ss << "*";
  }

  if (record.rNextId == -1)
    ss << "\t*\t";
  else if (record.rNextId == record.rID)
    ss << "\t=\t";
  else
    ss << '\t' << hdr->target_name[record.rNextId] << '\t';

  ss << record.pNext+1 << '\t' << record.tLen << '\t' << record.seq << '\t' << record.qual;

  for (unsigned i = 0; i < length(record.tags);)
  {
    ss << '\t' << record.tags[i] << record.tags[i+1] << ':';
    i += 3;

    switch(record.tags[i-1])
    {
      case 'A':
      {
        /* A printable character */
        ss << "A:" << record.tags[i];
        ++i;
        break;
      }

      case 'Z':
      {
        ss << "Z:";

        /* A string! Let's loop it until qNULL */
        while (record.tags[i] != '\0' && record.tags[i] != '\t' && record.tags[i] != '\n')
        {
          ss << record.tags[i];
          ++i;
        }

        ++i;
        break;
      }

      case 'c':
      {
        ss << "i:";
        int8_t read_int8 = 0;
        printTag(read_int8, record.tags, i, ss);
        break;
      }

      case 'C':
      {
        ss << "i:";
        uint8_t read_uint8 = 0;
        printTag(read_uint8, record.tags, i, ss);
        break;
      }

      case 's':
      {
        ss << "i:";
        int16_t read_int16 = 0;
        printTag(read_int16, record.tags, i, ss);
        break;
      }

      case 'S':
      {
        ss << "i:";
        uint16_t read_uint16 = 0;
        printTag(read_uint16, record.tags, i, ss);
        break;
      }

      case 'i':
      {
        ss << "i:";
        int32_t read_int32 = 0;
        printTag(read_int32, record.tags, i, ss);
        break;
      }

      case 'I':
      {
        ss << "I:";
        uint32_t read_uint32 = 0;
        printTag(read_uint32, record.tags, i, ss);
        break;
      }

      case 'f':
      {
        ss << "f:";
        float read_float = 0.0;
        memcpy(&read_float, &record.tags[i], sizeof(float));
        ss << read_float;
        i += sizeof(float);
        break;
      }

      default:
      {
        i = length(record.tags); // Unkown tag, stop
        break;
      }
    }
  }

  return ss.str();
}


inline bool
parse(BamAlignmentRecord & record, bam1_t * hts_record)
{
  syncCores(record, hts_record);

  // Parse qName
  record.qName = hts_record->data;
  auto it = hts_record->data + record._l_qname;

  // Parse CIGAR
  resize(record.cigar, record._n_cigar, Exact());
  static char const * CIGAR_MAPPING = "MIDNSHP=X*******";

  for (auto cig = begin(record.cigar, Standard()); cig != end(record.cigar, Standard()); ++cig)
  {
    uint32_t opAndCnt;
    memcpy(&opAndCnt, it, sizeof(uint32_t));
    it += sizeof(uint32_t);
    SEQAN_ASSERT_LEQ(opAndCnt & 15, 8u);
    cig->operation = CIGAR_MAPPING[opAndCnt & 15];
    cig->count = opAndCnt >> 4;
  }

  // Parse sequence
  resize(record.seq, record._l_qseq, Exact());

  for (int i = 0; i < record._l_qseq; ++i)
  {
    record.seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(it, i)];
  }

  it += (hts_record->core.l_qseq + 1) >> 1;

  // Parse qualities
  resize(record.qual, record._l_qseq, Exact());

  for (int i = 0; i < record._l_qseq; ++i, ++it)
  {
    record.qual[i] = static_cast<char>(*it + 33);
  }

  it = bam_get_aux(hts_record);
  resize(record.tags, bam_get_l_aux(hts_record), Exact());
  arrayCopyForward(it, it + bam_get_l_aux(hts_record), begin(record.tags, Standard()));

  return true;
}


inline bool
parse(bam1_t * hts_record, bam_hdr_t * hdr, BamAlignmentRecord const & record)
{
  kstring_t * s = static_cast<kstring_t*>(calloc(1, sizeof(kstring_t)));
  std::string str = toString(record, hdr);
  ksprintf(s, "%s", str.data());
  int ret = sam_parse1(s, hdr, hts_record);

  free(s->s);
  free(s);

  if (ret != 0)
  {
    std::cerr << "[seqan::hts_io.bam_alignment_record] ERROR parsing record:\n";
    std::cerr << str << std::endl;
    return false;
  }

  return true;
}

} // namespace seqan


#endif // INCLUDE_SEQAN_HTS_IO_BAM_ALIGNMENT_RECORD_H_
