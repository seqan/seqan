#ifndef SEQAN_APPS_NGS_ROI_VERSION_H_
#define SEQAN_APPS_NGS_ROI_VERSION_H_

#ifdef SEQAN_APP_VERSION
    #ifdef SEQAN_REVISION
        #define VERSION SEQAN_APP_VERSION " [" SEQAN_REVISION "]"
    #else
        #define VERSION SEQAN_APP_VERSION
    #endif
#else
    #define VERSION "0.2.2"
#endif
#ifdef SEQAN_DATE
    #define DATE SEQAN_DATE
#else
    #define DATE    "October 2013"
#endif

#endif  // #ifndef SEQAN_APPS_NGS_ROI_VERSION_H_
