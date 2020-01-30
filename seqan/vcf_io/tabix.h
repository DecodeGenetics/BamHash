#ifndef SEQAN_VCF_IO_TABIX_H_
#define SEQAN_VCF_IO_TABIX_H_

#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <cstring>
#include <vector>

#include <seqan/vcf_io/vcf_record.h>


namespace seqan {


class Tabix
{
 public:
  htsFile* fp = nullptr;
  tbx_t* tbx = nullptr;
  hts_itr_t* hts_iter = nullptr;
  const tbx_conf_t *idxconf;
  String<String<char> > chroms;
  unsigned rID = 0;
  StringSet<CharString> samples;
};


inline void
clear(Tabix & index)
{
    if (index.fp)
    {
        hts_close(index.fp);
        index.fp = NULL;
    }

    if (index.hts_iter)
    {
        tbx_itr_destroy(index.hts_iter);
        index.hts_iter = NULL;
    }

    if (index.tbx)
    {
        tbx_destroy(index.tbx);
        index.tbx = NULL;
    }

    clear(index.chroms);
    index.rID = 0;
    clear(index.samples);
}

inline void
getHeader(seqan::CharString& header_string, Tabix & index)
{
  // clear(header);
  kstring_t str = {0,0,0};

  while ( hts_getline(index.fp, KS_SEP_LINE, &str) >= 0 )
  {
    if ( !str.l || str.s[0] != index.tbx->conf.meta_char )
    {
      break;
    }
    else
    {
      if (length(str.s) > 5 && strncmp(str.s, "#CHROM", 6) == 0)
      {
        seqan::CharString h_line(str.s);
        seqan::strSplit(index.samples, h_line, seqan::EqualsChar<'\t'>());
        size_t end = 8;
        if (length(index.samples) > 8 && index.samples[8] == "FORMAT") end = 9;
        erase(index.samples, 0, end);
      }
      else
      {
        append(header_string, str.s);
        append(header_string, "\n");
      }
    }
  }

  free(str.s);

  //set back to start
  // index.rID = 0;
  //
  // if (index.hts_iter)
  // {
  //   tbx_itr_destroy(index.hts_iter);
  //   index.hts_iter = NULL;
  // }
  //
  // index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[rID]));
}

inline bool
_onLastRId(Tabix & index)
{
  return index.rID == length(index.chroms) - 1;
}

inline void
_nextRId(Tabix & index)
{
    ++index.rID;

    if (index.hts_iter)
    {
        tbx_itr_destroy(index.hts_iter);
        index.hts_iter = NULL;
    }

    index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[index.rID]));
}

inline bool
_extractLineOfThisRId(String<char> & line, Tabix & index)
{
  kstring_t str = {0,0,0};

  if (index.hts_iter && tbx_itr_next(index.fp, index.tbx, index.hts_iter, &str) >= 0)
  {
    line = str.s;
    free(str.s);
    return true;
  }

  free(str.s);
  return false;
}

inline bool
_extractLineOfThisRId(std::string & line, Tabix & index)
{
  kstring_t str = {0,0,0};

  if (index.hts_iter && tbx_itr_next(index.fp, index.tbx, index.hts_iter, &str) >= 0)
  {
    line = str.s;
    free(str.s);
    return true;
  }

  free(str.s);
  return false;
}

inline void
_insertDataToVcfRecord(VcfRecord & record, String<char> const & line, unsigned const & rID)
{
    seqan::StringSet< seqan::String<char> > splitted_line;
    seqan::strSplit(splitted_line, line, seqan::EqualsChar<'\t'>());
    record.rID = rID;
    lexicalCast(record.beginPos, splitted_line[1]);
    --record.beginPos; // Change from 1-based to 0-based indexing
    record.id = splitted_line[2];
    record.ref = splitted_line[3];
    record.alt = splitted_line[4];
    lexicalCast(record.qual, splitted_line[5]);
    record.filter = splitted_line[6];
    record.info = splitted_line[7];

    if (length(splitted_line) > 8)
    {
        record.format = splitted_line[8];
    }

    if (length(splitted_line) > 9)
    {
        erase(splitted_line, 0, 9);
        record.genotypeInfos = std::move(splitted_line);
    }
}

/**
 * @brief Reads a VCF record to a single string.
 * @details [long description]
 *
 * @param[in,out] line [description]
 * @param[in] index [description]
 */
inline bool
readRawRecord(String<char> & line, Tabix & index)
{
  while (!_onLastRId(index))
  {
    if(_extractLineOfThisRId(line, index))
    {
      return true;
    }

    _nextRId(index);
  }

  // We are on the last rID.
  if(_extractLineOfThisRId(line, index))
  {
    return true;
  }

  return false;
}


inline bool
readRawRecord(std::string & line, Tabix & index)
{
  while (!_onLastRId(index))
  {
    if(_extractLineOfThisRId(line, index))
    {
      return true;
    }

    _nextRId(index);
  }

  // We are on the last rID.
  if(_extractLineOfThisRId(line, index))
  {
    return true;
  }

  return false;
}

/**
 * @brief Read a VCF record from a tabix file.
 * @details The tabix file needs to have previously been opened. If no record can be read false is returned,
 * this can happen either because we've reached the end of the specified region or end of the BCF file.
 *
 * @param record A VCF record.
 * @param index Tabix index.
 *
 * @return True means a new record was read, false is returned otherwise.
 */
inline bool
readRecord(VcfRecord & record, Tabix & index)
{
    seqan::String<char> line;
    if (!seqan::readRawRecord(line, index))
        return false;

    _insertDataToVcfRecord(record, line, index.rID);
    return true;
}

/**
 * @brief Changes the region of the index.
 * @details The
 *
 * @param[in,out] index Tabix index.
 * @param[in] region The region to change to. It should be on the format chrX or chrX:Y-Z.
 */
inline void
setRegion(Tabix & index, const char * region)
{
    if (index.hts_iter)
    {
        tbx_itr_destroy(index.hts_iter);
        index.hts_iter = NULL;
    }

    index.hts_iter = tbx_itr_querys(index.tbx, region);
}

/**
 * @brief Reads a region of a VCF file.
 * @details Note: The records variable is not clear before adding VCF records, to allow extractions of multiple regions.
 *
 * @param[in,out] records A list of VCF records.
 * @param index [description]
 * @param region [description]
 */
inline void
readRegion(seqan::String<VcfRecord> & records, Tabix & index, const char * region)
{
    setRegion(index, region);
    seqan::String<char> line;

    while(_extractLineOfThisRId(line, index))
    {
        VcfRecord record;
        _insertDataToVcfRecord(record, line, index.rID);
        append(records, record);
    }
}


inline bool
readRegion(VcfRecord & record, Tabix & index)
{
    seqan::String<char> line;

    if (!_extractLineOfThisRId(line, index))
        return false;

    _insertDataToVcfRecord(record, line, index.rID);
    return true;
}


inline void
open(Tabix & index, char const * vcfFilename, const char * fileMode = "r")
{
  clear(index);
  struct stat stat_tbi,stat_vcf;
  char *fnidx = (char*) calloc(strlen(vcfFilename) + 5, 1);
  strcat(strcpy(fnidx, vcfFilename), ".tbi");

  if (bgzf_is_bgzf(vcfFilename)!=1 )
  {
    SEQAN_FAIL("File '%s' was not identified as bgzipped. Please use bgzip to compress the file.", vcfFilename);
    std::free(fnidx);
  }

  // Common source of errors: new VCF is used with an old index
  stat(fnidx, &stat_tbi);
  stat(vcfFilename, &stat_vcf);

  if (stat_vcf.st_mtime > stat_tbi.st_mtime )
  {
    SEQAN_FAIL("The index file is older than the bcf file. Please reindex the bcf file.");
  }

  std::free(fnidx);

  if ((index.fp = hts_open(vcfFilename, fileMode)) == 0)
    SEQAN_FAIL("Fail to open the VCF file.");

  if ((index.tbx = tbx_index_load(vcfFilename)) == NULL)
    SEQAN_FAIL("Failed to load the VCF index file.");

  int nseq;
  const char** seq = tbx_seqnames(index.tbx, &nseq);

  for (int i = 0; i < nseq; ++i)
  {
    appendValue(index.chroms, seq[i]);
  }

  std::free(seq);
  index.idxconf = &tbx_conf_vcf;

  // set up the iterator, defaults to the beginning
  index.rID = 0;

  if (seqan::length(index.chroms) > 0)
    index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[0]));
}

}  // namespace seqan

#endif  // SEQAN_VCF_IO_TABIX_H_
