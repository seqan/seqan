#!/usr/bin/awk -f

BEGIN {

    INHEADER = 0
    num_iorev = 0
    num_doc = 0
    num_bug = 0
    num_nodoc = 0
    num_duplicate = 0
    num_delcandidate = 0
    num_hascref = 0
    num_hascomment = 0
    num_notio = 0
    num_noop = 0
    num_stub = 0
    num_tested = 0
    num_nottested = 0
    num_notused = 0
    

    # due to some awkward awk-ness FILENAME is not known in BEGIN
    needtoprinthead = 1 

}


needtoprinthead == 1 {
    print ""
    print ""
    print "============ " FILENAME " ============"
    print ""
    needtoprinthead = 0
}


$0 == "/* IOREV" {
    INHEADER = 1

    print "-------- General File Information --------"
    print ""
}

INHEADER == 1 {

    if ($0 == " */")
        INHEADER = 0
    else {
        s = substr($0, 4, length($0))
        if ((substr(s, 1,1) != "_") && (s != "IOREV"))
            print " " s
            
        if (s == "_windows_")
            print " contains windows specific code"
    }
}


$0 ~ "//IOREV" {
    num_iorev = num_iorev +1
    if ($0 ~ /[[:blank:]][[:alnum:]]+[[:blank:]]/)
        num_hascomment = num_hascomment + 1

    if ($0 ~ "_doc_")
        num_doc = num_doc +1

    if ($0 ~ "_nodoc_" )
        num_nodoc = num_nodoc +1

    if ($0 ~ "_bug_" )
        num_bug = num_bug +1

    if ($0 ~ "_duplicate_") 
        num_duplicate = num_duplicate +1

    if ($0 ~ "_delcandidate_") 
        num_delcandidate = num_delcandidate +1

    if ($0 ~ "_hascref_" )
        num_hascref = num_hascref +1

    if ($0 ~ "_notio_")
        num_notio = num_notio +1

    if ($0 ~ "_noop_")
        num_noop = num_noop +1

    if ($0 ~ "_stub_")
        num_stub = num_stub +1

    if ($0 ~ "_tested_")
        num_stub = num_stub +1

}

END {
    print ""
    print ""
    print "--------      Statistics         --------"
    print ""
    print " # of tagged functions, structs etc:                " num_iorev
    print " # falsely tagged (regex false positive):           " num_notio
    print ""
    print " # well documented:                                 " num_doc
    print " # not documented but should be:                    " num_nodoc
#    print " # well tested / used by other stuff:               " num_tested
    print ""
    print " # that don't do anything, but maybe should:        " num_noop
    print " # that seem to be incomplete/not done:             " num_stub
    print ""
    print " # that contain an obvious bug:                     " num_bug
    print " # that duplicate existing functionality:           " num_duplicate
    print " # that can be removed for other reasons:           " num_delcandidate
    print ""
    print " # of parsing functions that get \"TChar c&\":        " num_hascref
    print ""
    print " # with free-text comment:                          " num_hascomment
    print ""
    print "--------         END              --------"
    print ""
    print ""

}