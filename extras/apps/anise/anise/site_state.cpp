#include "site_state.h"

#include <seqan/basic.h>

#include "file_name_tokens.h"
#include "streaming_exception.h"

// ----------------------------------------------------------------------------
// Class AssemblySiteState
// ----------------------------------------------------------------------------

void AssemblySiteState::load(TemporaryFileManager & tmpMgr, int _siteID)
{
    std::fstream f;
    tmpMgr.open(f, std::ios::binary | std::ios::in, SITE_STATE_TOKEN, SITE_STATE_EXT, _siteID);
    if (!f.good())
        throw AniseIOException() << "Could not open site state file for site " << _siteID << " for reading.";

    std::string buffer;
    buffer.resize(1024);
    while (f.good() && f.peek() == '#')  // Skip comments.
        f.getline(&buffer[0], buffer.size());

    f >> buffer;  // "SITE_ID"
    f >> siteID;

    f >> buffer;  // "RID"
    f >> rID;

    f >> buffer;  // "REF_NAME"
    f >> refName;

    f >> buffer;  // "POS"
    f >> pos;
    pos -= 1;

    f >> buffer;  // "STEP_NO"
    f >> stepNo;

    f >> buffer;  // "ACTIVE"
    f >> active;

    f >> buffer;  // "NUM_NEW_ALIGNMENTS"
    f >> numNewAlignments;

    f >> buffer;  // "COMMENT"
    getline(f, comment);
    // trim leading/trailing whitespace
    while (!comment.empty() && isspace(comment.front()))
        comment.erase(0, 1);
    while (!comment.empty() && isspace(comment.back()))
        comment.pop_back();
}
void AssemblySiteState::save(TemporaryFileManager & tmpMgr)
{
    SEQAN_CHECK(siteID >= 0, "Invalid siteID.");

    std::fstream f;
    tmpMgr.open(f, std::ios::binary | std::ios::out, SITE_STATE_TOKEN, SITE_STATE_EXT, siteID);
    if (!f.good())
        throw AniseIOException() << "Could not open site state file for site " << siteID << " for writing.";

    f << "#ANISE SITE STATE\n"
      << "SITE_ID\t" << siteID << "\n"
      << "RID\t" << rID << "\n"
      << "REF_NAME\t" << refName << "\n"
      << "POS\t" << (pos + 1) << "\n"
      << "STEP_NO\t" << stepNo << "\n"
      << "ACTIVE\t" << active << "\n"
      << "NEW_ALIGNMENTS\t" << numNewAlignments << "\n"
      << "COMMENT\t" << comment << "\n";
}
