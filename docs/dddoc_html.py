# coding: latin-1
import os
import datetime
import dddoc
import dddoc_html_trans
import xml.sax
import operator
import sys
from os import F_OK


# Number of warnings.
WARNING_COUNT = 0
OUT_PATH = 'html'

################################################################################

def createDocs(path, buildfull, indexonly, include_dirs):
    global globalDocsPath
    globalDocsPath = path

    global globalBuildFull
    globalBuildFull = buildfull

    global includeDirs
    includeDirs = include_dirs
    
    if not os.access(path, os.F_OK): 
        os.mkdir(path)
    copyFile(path, "dddoc_html.css")
    copyFile(path, "seqan_logo.gif")
    copyFile(path, "dddoc_empty.gif")
    copyFile(path, "dddoc_plus.gif")
    copyFile(path, "dddoc_minus.gif")
    copyFile(path, "dddoc.js")

    gatherGlossary()
        
    if buildfull or indexonly:
        createIndexes(path)
        createSearchfile(path)

    if not indexonly:
        createPages(path)

################################################################################

def createIndexes(path):
    cats = dddoc.DATA["globals.indexes"].keys()
    for cat in cats:
        print 'Indexes for ' + cat,
        
        entries = collectIndexEntries(cat)
        
        subcats = entries.keys()
        subcats.sort()
        for subcat in subcats:
            subcat_entries = entries[subcat]
            filename = os.path.join(path, getIndexname(cat, subcat))
            fl = file(filename, "w")
            print '.',
            pageIndex(fl, path, cat, subcat, entries, subcats)
            fl.close()
 
            
        filename = os.path.join(path, getIndexpage(cat))
        fl = file(filename, "w")
        print '.',
        pageIndexpage(fl, cat)
        fl.close()
        print


################################################################################

def createPages(path):
    cats = dddoc.DATA["globals.categories"].keys()
    for cat in cats:
        print 'Pages for ' + cat,
        
        entries = dddoc.DATA[cat]
        for key in entries.keys():
            data = entries[key]
            filename = os.path.join(path, getFilename(data.name(0), data.name(1)))
            fl = file(filename, "w")
            print '.',
            pageContent(fl, data)
            warningPage(cat, key, data)
            fl.close()
            
        print



################################################################################

def copyFile(path, filename):
    out_path = os.path.join(path, filename)
    if not os.access(out_path, os.F_OK):
        in_fl = file(filename, "rb")
        out_fl = file(out_path, "wb")
        out_fl.write(in_fl.read())
        in_fl.close()
        out_fl.close()
            
#######################################################################

def escapeFiles(text):
    text = text.replace("_", "__")

    ret = ""
    for i in range(len(text)):
    	if (text[i] >= 'A') and (text[i] <= 'Z'):
    		ret += "_"
    	ret += text[i]

    ret = ret.replace("\t", "_09")
    ret = ret.replace("\n", "_0a")
    ret = ret.replace("!", "_21")
    ret = ret.replace("\"", "_22")
    ret = ret.replace("#", "_23")
    ret = ret.replace("$", "_24")
    ret = ret.replace("%", "_25")
    ret = ret.replace("&", "_26")
    ret = ret.replace("'", "_27")
    ret = ret.replace("(", "_28")
    ret = ret.replace(")", "_29")
    ret = ret.replace("*", "_2a")
    ret = ret.replace("+", "_2b")
    ret = ret.replace("/", "_2f")
    ret = ret.replace(":", "_3a")
    ret = ret.replace(",", "_2c")
    ret = ret.replace("<", "_3c")
    ret = ret.replace(">", "_3e")
    ret = ret.replace("?", "_3f")
    ret = ret.replace("\\", "_5c")
    ret = ret.replace("|", "_7c")
    ret = ret.replace(" ", "+")
    
    if (len(ret) == 0) or (ret[0] == '_'): return ret
    else: return '.'+ret

################################################################################

def escapeHTML(text):
    text = text.replace("&", "&amp;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")
    
    if (text.find("\\") >= 0):
        text = text.replace("\\\\", "&backslash;")
    	text = dddoc_html_trans.translate(text);
        text = text.replace("&backslash;", "\\")

    return text

################################################################################

def escapeJavaScript(text):
    text = text.replace("\\", "\\\\")
    text = text.replace("'", "\\'")
    text = text.replace("\n", " ")
    text = text.replace("\r", "")

    return text

################################################################################

def getFilename(cat, item):
    return cat.upper() + escapeFiles(item) + ".html"

def getIndexpage(cat):
    s = dddoc.DATA["globals.project.indexcategory"].text()
    
    if (s[0:(len(cat))] == cat): 
    	return "index.html"
    else:
    	return "INDEXPAGE" + escapeFiles(cat) + ".html"

def getIndexname(cat, subcat = ""):
    if len(subcat) == 0: return "INDEX" + escapeFiles(cat) + ".html"
    else: return "INDEX" + escapeFiles(cat) + escapeFiles(subcat) + ".html"

def getIndexnameLink(cat, item, subcat = ""):
    if len(subcat) == 0: return "INDEX" + escapeFiles(cat) + ".html#" + item
    else: return "INDEX" + escapeFiles(cat) + escapeFiles(subcat) + ".html#" + item

def getDemoFilename(sourcefile):
    return "DEMO" + escapeFiles(sourcefile) + ".html"

    
################################################################################
    
def translateText(text, line=None):
    ret = ''
    str = ''
    in_code = False
    in_link = False
    in_escaped = False
    pos = 0
    while (pos < len(text)):
        c = text[pos]
        if in_escaped:
            str += c;
            in_escaped = False;
        elif in_code:
            if c == '$':
                if str != '': 
                    ret += translateCode(str)
                    str = ''
                else:
                     str = '$'
                in_code = False
            else: str += c
        elif in_link:
            if c == '@': 
                ret += translateLink(str, line=line)
                str = ''
                in_link = False
            else: str += c
        else:
            if c == '$': 
                if str != '' or pos == 0: 
                    ret += escapeHTML(str)
                    str = ''
                else:
                    str = '$'
                in_code = True
            elif c == '@':
                ret += escapeHTML(str)
                str = ''
                in_link = True
#            elif c == '\\':
#                in_escaped = True
            else:
                str += c;
                
        pos += 1
                
    if str != '':
        ret += escapeHTML(str)

    return ret
    
################################################################################
    
def translateCode(text):
    text = text.strip(" \n\r")
    text = escapeHTML(text)
    text = text.replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
    text = text.replace(" ", "&nbsp;")
    text = text.replace("\n", "<br >")
    return '<span class=code><nobr>' + text + '</nobr></span>'
    
################################################################################
    
def translateTooltip(text):
    text = text.replace("\t", " ")
    text = text.replace("\n", " ")
    text = text.replace("$", "")
    text = text.replace("@", "")
    text = text.replace("\"", "&quot;")
    text = text.replace("<", "&lt;")
    text = text.replace("<", "&gt;")
    return text

################################################################################

def brokenLinkText(text):
    """Returns the HTML for a broken link.

    Args:
      text  String to format.

    Returns:
      String, formatted as broken link in HTML.
    """
    return '<span class="broken_link">' + text + '</span>'


def brokenLink(text, line=None):
    """Format the given text as a broken link and print warning.

    If line is not None then the warning includes the source of the
    line to help debugging.

    Args:
      text  String, text to format as broken link.
      line  dddoc.Line object to use as the source.

    Returns:
      String with the HTML for the broken link.
    """
    global WARNING_COUNT
    WARNING_COUNT += 1
    print
    print '!!  WARNING: broken link "' + text + '"'
    if line:
        print '        Location: %s:%d' % (line.file_name, line.line_no)
    # The following is here to debug the source of broken link tags.
    # If you get a broken link message without a source, uncomment the
    # following two lines and add enough "line" parameters to the code
    # so these errors are printed with source the next time.
    #import traceback
    #traceback.print_stack()
    return brokenLinkText(text)


################################################################################

def findGlossary(text):
    global globalGlossary
    
    if globalGlossary.has_key(text):
        return globalGlossary[text][0]
       
    text2 = text.lower() 
    for key in globalGlossary.keys():
        key2 = key.lower()
        if text2.find(key2) == 0:
            return globalGlossary[key][0]
            
    return False
    
################################################################################

def translateLinkDisplaytext(text):
    pos = text.find(':')
    if pos >= 0: 
        protocol = text[:pos]
        rest = text[pos+1:]
    else: 
        protocol = ""
        rest = text

    if (protocol == 'http') or (protocol == 'ftp'):
        arr = dddoc.splitUrl(text)
        return arr[len(arr) - 1]
        
    if (protocol == 'glos'):
        arr = dddoc.splitName(rest)
        if len(arr) == 0: return brokenLinkText(text)
        glos = findGlossary(arr[0])
        if not glos: return brokenLinkText(text)
        if len(arr) == 1: return arr[0]
        return arr[1]
        
    if (protocol == 'nolink'):
        return translateID(rest)

    arr = dddoc.splitName(text)
    if len(arr) == 0: return brokenLinkText(text)

    if (dddoc.DATA["globals.categories"][arr[0]].empty()): return brokenLinkText(text)

    if text.find('.') < 0: #Link to indexpage
        if (dddoc.DATA["globals.indexes"][arr[0]].empty()): return brokenLinkText(text)
        if len(arr) == 1: return getCategoryTitle(arr[0])
        else: return arr[1] 

    if (len(arr) < 2): return brokenLinkText(text)
    
    return translateID(arr[len(arr) - 1])

################################################################################


def translateLink(text, attribs = "", line=None, cat=False):
    global globalDocsPath
    
    pos = text.find(':')
    if pos >= 0: 
        protocol = text[:pos]
        rest = text[pos+1:]
    else: 
        protocol = ""
        rest = text

    #external Link
    if (protocol == 'http') or (protocol == 'ftp'): 
        arr = dddoc.splitUrl(text)
        return '<a href="' + arr[0] + '" ' + attribs + '>' + arr[len(arr) - 1] + '</a>'

    #Glossary Link
    if protocol == 'glos':
        arr = dddoc.splitName(rest)
        if len(arr) == 0: return brokenLink(text, line=line)
        glos = findGlossary(arr[0])
        if not glos: return brokenLink(text, line=line)
        if len(arr) == 1: t = arr[0]
        else: t = arr[1]
        return '<a class=glossary_link ' + glos[2] + t + '</a>'

    #no link    
    if (protocol == 'nolink'):
        return translateID(rest, cat=cat)

    arr = dddoc.splitName(text)
    if len(arr) == 0: return brokenLink(text, line=line)
    
    if (dddoc.DATA["globals.categories"][arr[0]].empty()): return brokenLink(text, line=line)
    
    #Link to Indexpage
    if text.find('.') < 0:
        if (dddoc.DATA["globals.indexes"][arr[0]].empty()): return brokenLink(text, line=line)
        if len(arr) == 1: t = getCategoryTitle(arr[0])
        else: t = arr[1]
        return '<a href="' + getIndexpage(arr[0]) + '" ' + attribs + '>' + t + '</a>'


    #Link to Page
    if (len(arr) < 2): return brokenLink(text, line=line)

    href = getFilename(arr[0], arr[1])

    obj = dddoc.DATA[arr[0]][arr[1]];
    if obj.empty():
        #test existing file
        doc_path = os.path.join(globalDocsPath, href)
        if not os.access(doc_path, os.F_OK): 
            return brokenLink(text, line=line)
        summary = '';
    else:
        summary = translateTooltip(obj["summary"].text());
        
    if len(summary) > 0:
        summary = 'title="' + summary + '"'
    
    return '<a href="' + href + '" ' + summary + ' ' + attribs + '>' + translateID(arr[len(arr) - 1], cat=cat) + '</a>'


################################################################################
    
def translateImage(text):
    in_path = 'img/' + text + '.png'
    global OUT_PATH
    out_path = OUT_PATH + '/' + text + '.png'
    
    if os.access(in_path, os.F_OK):
        in_fl = file(in_path, "rb")
        out_fl = file(out_path, "wb")
        out_fl.write(in_fl.read())
        in_fl.close()
        out_fl.close()
    else:
        global WARNING_COUNT
        WARNING_COUNT += 1
        print
        print "\n!!  WARNING: image not found: \"" + text + ".png\""

    text = escapeHTML(text)
    text = text.replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
    text = text.replace(" ", "&nbsp;")
    text = text.replace("\n", "<br >")
    return '<img class=image src="./' + text + '.png" border=0 />'


################################################################################
    
def translateID(text, line=None, cat=False):
    i = text.find('#');
    if (i >= 0):
        if cat:
            text = '%s (%s)' % (text[i + 1:], text[:i])
        else:
            text = text[i + 1:]
    i = text.find('|');
    if (i >= 0):
        text = text[i + 1:]
    return translateText(text, line)

################################################################################
         
def sortingKey(data, key):
    if data["order"].empty():
        return translateID(key)
    
    return data["order"].text();

 
################################################################################

def getBeforeColon(text):
    pos = text.find(':')
    if pos < 0: return text
    else: return text[0: pos]

################################################################################

def getAfterColon(text):
	return text[text.find(':') + 1:len(text)]

################################################################################

def getPageTitle(data):
	s = data["title"].text()
	if not s:
		return translateID(data.name(1))
	else:
		return s
	
################################################################################
	
def getCategoryTitle(cat, show_cat=False):
	s = dddoc.DATA["Indexpage"][cat]["title"].text()
	if (s == ''): 
		return translateID(dddoc.DATA["globals.indexes"][cat].text(), cat=show_cat)
	else:
		return s


################################################################################

def addCollectIndexEntry(data, cat, subcat, key, entries):
    subcat2 = subcat
    while subcat2 != '':
        if not entries.has_key(subcat2): entries[subcat2] = {}
        pos = subcat2.rfind('.')
        if pos < 0: break
        subcat2 = subcat2[:pos]
        
    key4sorting = sortingKey(data, key)
    if not entries[subcat].has_key(key4sorting): entries[subcat][key4sorting] = []
    entries[subcat][key4sorting].append(cat + "." +key)

################################################################################

def collectIndexEntries(cat):
    entries = {}
    entries[''] = {}
    data = dddoc.DATA[cat]
    for key in data.keys():
        entry = data[key]
        if not entry["hidefromindex"].empty(): continue
        if entry["cat"].empty(): addCollectIndexEntry(entry, cat, '', key, entries)
        else:
            for subcat_line in entry["cat"].lines:
                subcat = subcat_line.text()
                addCollectIndexEntry(entry, cat, subcat, key, entries)
                
    return entries

################################################################################

def pageIndexPrintMembers(fl, data):
    # Note that data is not a dddoc.Data object!
    entrs = data.keys()
    entrs.sort()
    for entr in entrs:
        links = data[entr]
        for link in links:
            fl.write('<div class=index_item>' + translateLink(link, "target=_top", cat=True) + '</div>')

################################################################################

def pageIndex(fl, path, cat, subcat, entries, subcats):
    fl.write('<html>')
    fl.write('<head>')
    fl.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8">');
    fl.write('<link rel="stylesheet" href="dddoc_html.css" type="text/css" />')
    fl.write('</head>\n')
    fl.write('<script src="searchfile.js"></script>\n') 
    fl.write('<script src="dddoc.js"></script>\n') 
    fl.write('<body id=index_body>')
    
    lines = dddoc.DATA["globals.indexes"]
    for cat2 in lines.keys_by_occ():
        if cat2 == cat:
            fl.write('<div class=index_section_high>')
        else:
            fl.write('<div class=index_section>')
        fl.write('<div class=index_cat><a class=index_link target=_top href="' + getIndexpage(cat2) + '">' + lines[cat2].text() + '</a></div>')
        if cat2 == cat:
            
            # print subfolder
            this_reached = False
            members_printed = False
            for subcat2 in subcats:
                is_this = (subcat2 == subcat)
                is_sub = ((subcat2.find(subcat) == 0) and (((len(subcat2) > len(subcat)) and (subcat2[len(subcat)] == '.')) or (len(subcat) == 0)))
                is_child = (is_sub and (subcat2.find('.', len(subcat)+1) < 0))
                is_super = (subcat.find(subcat2) == 0)
                has_depth_greater_2 = (subcat2.find('.') >= 0)
                is_sister = not has_depth_greater_2 or (subcat.find(subcat2[0:subcat2.rfind('.')]) == 0)
                
                if not members_printed and this_reached and not is_sub:
                    pageIndexPrintMembers(fl, entries[subcat])
                    members_printed = True
                    
                if is_this:
                    this_reached = True
                
                if (subcat2 != "") and (is_child or is_super or is_sister or is_this):
                    #print out folder
                    indent = ''
                    display_text = subcat2
                    while True:
                        pos = display_text.find('.')
                        if pos < 0: break
                        indent += '&nbsp;&nbsp;'
                        display_text = display_text[pos+1:]
                        
                    if is_super: image = 'dddoc_minus.gif'
                    else: image = 'dddoc_plus.gif'
                        
                    fl.write('<div class=index_subcat>' + indent + '<img src="' + image + '" border=0><a class=index_link href="' + getIndexname(cat, subcat2) + '">' + display_text + '</a></div>')
                    
            if not members_printed:
                pageIndexPrintMembers(fl, entries[subcat]) 
                    
        fl.write('</div>') #index_section or index_section_high
        
    printSearchmask(fl)
    fl.write('</body>')
    fl.write('</html>')

################################################################################

def printSearchmask(fl):
    fl.write('<div id=searchmask>')
    fl.write('<div id=searchtitle>Searching</div>')
    fl.write('<input id=search name=search onKeyUp="updateSearch(this.value);" onBlur="updateSearch(this.value);" autocomplete="off">')
    fl.write('<div id=result></div>')
    fl.write('</div>')

################################################################################

def pageIndexpage(fl, cat):
    fl.write('<html>')
    fl.write('<head>')
    fl.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8">');
    fl.write('<link rel="stylesheet" href="dddoc_html.css" type="text/css" />')
    fl.write('<title>' + getCategoryTitle(cat) + '</title>')
    fl.write('</head>')
    fl.write('<body>')
    fl.write('<table id=main_table cellspacing=0 cellpadding=0>')
    fl.write('<tr><td valign=top>')
    fl.write('<iframe frameborder=0 id=navigation src="' + getIndexname(cat) + '"></iframe>')
    fl.write('</td><td valign=top>')
    fl.write('<div id=content>')

    #s = dddoc.DATA["globals.indexes"][cat].text()
    fl.write('<div class=indexpage_title>' + getCategoryTitle(cat) + '</div>')

    data = dddoc.DATA["Indexpage"][cat]

    printSummary(fl, data, "summary")

    for line in data.at_level(0).lines:
        s = translateText(line.text())
        if (s): fl.write('<div class=text>' + s + '</div>')

    printTextblock(fl, data, "description")
    printTextblock(fl, data, "remarks")
    printIndexpageMembers(fl, dddoc.DATA[cat])
    printTextblock(fl, data, "example")
    printLink(fl, data, "demo")
    printLink(fl, data, "see")

    pageEnd(fl, data)
    fl.write('</div>')
    fl.write('</td></tr>')
    fl.write('</table>')
    fl.write('<p style="font-size:50%%; color: #909090">Page built @%s</p>' %
             datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
    fl.write('</body>')
    fl.write('</html>')

################################################################################

def addIndexPageMembers(data, key, entries, subcat):
    if not entries.has_key(subcat): 
        entries[subcat] = {}
        
    s = '<tr><td class=value_key valign=top><nobr>'
    s += '<a href="' + getFilename(data.name(0), key) + '">' + translateID(key, cat=True) + '</a>'
    s += '</nobr></td><td class=value_text valign=top>'

    summary = translateText(data[key]["summary"].text())
    if summary: s += summary
    else: s += '&nbsp;'
    
    s += '</td></tr>'
    
    key4sorting = sortingKey(data[key], key)
    if not entries[subcat].has_key(key4sorting): entries[subcat][key4sorting] = ''
    entries[subcat][key4sorting] += s

################################################################################

def printIndexpageMembers(fl, data):
    keys = data.keys()
    if len(keys) > 0:

        entries = {}
        for key in keys:
            if not data[key]["hidefromindex"].empty():
                continue

            linelist = data[key]["cat"].lines
            if (linelist == []):
                addIndexPageMembers(data, key, entries, "")
            else:
                for line in linelist:
                    addIndexPageMembers(data, key, entries, line.text())
                           
        keys2 = entries.keys();
        keys2.sort();

        for key in keys2:
            fl.write('<div class=section>')
            
            if key != "":
                s = key.replace('.', ': ')
                fl.write('<div class=section_headline>' + s + '</div>')
            else:
                s = dddoc.DATA["globals.indexes"][data.name(0)].text()
                if s: fl.write('<div class=section_headline>' + s + '</div>')
                
            fl.write('<table class=indexpage_members_tab cellspacing=0 cellpadding=0>')
            
            entry_keys = entries[key].keys();
            entry_keys.sort();
            for entry in entry_keys:
                fl.write(entries[key][entry])
            fl.write('</table>')
        
            fl.write('</div>')

################################################################################

def pageContent(fl, data):
    pageBegin(fl, data)
    fl.write('<table id=main_table cellspacing=0 cellpadding=0>')
    fl.write('<tr><td valign=top>')
    
    cat = "Class"
    item = ""
    if ((data.name(0) == 'Memfunc') or (data.name(0) == 'Memvar') or (data.name(0) == 'Typedef')):
        arr = dddoc.splitName(data["class"].text())
        if (len(arr) > 2): item = arr[1]
#    elif data.name(0) == 'Spec':
#        arr = dddoc.splitName(data["general"].text())
#        if (len(arr) > 2): item = arr[1]
    else:
        cat = data.name(0)
        item = data.name(1)
    
    subcats = dddoc.DATA[cat][item]["cat"].lines
    subcat = ""
    if (len(subcats) > 0): subcat = subcats[0].text()
    
    fl.write('<iframe frameborder=0 id=navigation src="' + getIndexnameLink(cat, item, subcat) + '"></iframe>')
    fl.write('</td><td valign=top>')
    fl.write('<div id=content>')
    writePage(fl, data)
    pageEnd(fl, data)
    fl.write('</div>')
    fl.write('</td></tr>')
    fl.write('</table>')
    fl.write('<p style="font-size:50%%; color: #909090">Page built @%s</p>' %
             datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
    fl.write('</body>')
    fl.write('</html>')

################################################################################

def pageBegin(fl, data):
    fl.write('<html>')
    fl.write('<head>')
    fl.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8">');
    fl.write('<link rel="stylesheet" href="dddoc_html.css" type="text/css" />')
    fl.write('<title>' + getPageTitle(data) + '</title>')
    fl.write('</head>')
    fl.write('<body>')

################################################################################

def pageEnd(fl, data):
    s = dddoc.DATA["globals.project.footline"].text()
    fl.write('<div class=page_footline>' + s + '</div>')
    fl.write('<div id=page_widthblock>&nbsp;</div>')
    
################################################################################

def writePage(fl, data):
    printTitle(fl, data)
    printSummary(fl, data, "summary")
    
    for line in data.at_level(0).lines:
        s = translateText(line.text())
        if (s): fl.write('<div class=text>' + s + '</div>')
           
    printTextblock(fl, data, "description")
    printGlossary(fl, data, "glossary")
    
    printTree(fl, data)

    printConcept(fl, data)

    printSignature(fl, data, "signature")
    printList(fl, data, "include")
    printTable(fl, data, "param")
    printTextblock(fl, data, "remarks")
    printTextblock(fl, data, "status")
    printTextblock(fl, data, "returns")
    printLink(fl, data, "class")
    printLink(fl, data, "general")
    printShortcutfor(fl, data, "shortcutfor")
    printLinkRek(fl, data, "implements")

    printMemberRek(fl, data, "conceptimplements")
    printLinkRek(fl, data, "conceptusedbymeta")
    printLinkRek(fl, data, "conceptusedbyfunc")

    printMember(fl, data, "spec")
    printMemberRek(fl, data, "type")
    printMemberRek(fl, data, "typedef")
    printMemberRek(fl, data, "memvar")
    printMemberRek(fl, data, "memfunc")
    printMemberRek(fl, data, "function")
    
    printTable(fl, data, "value")
    printTable(fl, data, "tag")

    printMember(fl, data, "shortcut")

    printTextblock(fl, data, "example")
    printLinkRek(fl, data, "demo")
    printFile(fl, data, "file")
    printSnippet(fl, data, "snippet")
    printTextblock(fl, data, "output")
    
    printLink(fl, data, "concept")
    printLink(fl, data, "demofor")
    printLink(fl, data, "see")

################################################################################

def printConcept(fl, data):
    
    if data["baseconcept"].empty() and data["childconcept"].empty() and data["conceptmetafunc"].empty() and data["conceptmemvar"].empty() and data["conceptmemfunc"].empty() and data["conceptfunc"].empty() and data["concepttypedef"].empty(): return
    
    fl.write('<div id=define_concept>')
    
    s = dddoc.DATA["globals.sections.conceptdefinition"].text()
    if s: fl.write('<div class=define_concept_headline>' + s + '</div>')
    
    printLink(fl, data, "baseconcept")
    printMemberRek(fl, data, "conceptmetafunc")
    printMemberRek(fl, data, "concepttypedef")
    printMemberRek(fl, data, "conceptmemvar")
    printMemberRek(fl, data, "conceptmemfunc")
    printMemberRek(fl, data, "conceptfunc")
    printMember(fl, data, "childconcept")
    
    fl.write('</div>')
    
################################################################################

def printTitle(fl, data):
    s = dddoc.DATA["globals.categories"][data.name(0)].text()
    if s: fl.write('<div class=page_category>' + s + '</div>')
    fl.write('<div class=page_title>' + getPageTitle(data) + '</div>')


################################################################################

def printText(fl, data, category):
    lines = data[category]
    if not lines.empty():
        for line in lines.lines:
            s = translateText(line.text())
            if s: fl.write('<div class=text>' + s + '</div>')
   

################################################################################

def printSummary(fl, data, category):
    fl.write('<div id=summary>')
    printText(fl, data, category)
    fl.write('</div>')


################################################################################

def printSignature(fl, data, category):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')
        
        for line in lines.lines:
            s = translateCode(line.text())
            if s: fl.write('<div class=signature_block><nobr>' + s + '</nobr></div>')

        fl.write('</div>');


################################################################################

def printTableContent(fl, data, category):
    lines = data[category]
    if not lines.empty():
        
        keys = lines.keys_by_occ()
        if len(keys) > 0:         
            fl.write('<table class=value_tab cellspacing=0 cellpadding=0>')
    
            for key in keys:
                fl.write('<tr><td class=value_key valign=top><nobr>' + key + '</nobr></td><td class=value_text valign=top>')
                subprintText(fl, lines[key])
                fl.write('</td></tr>')
            fl.write('</table>')
            
################################################################################


def printTable(fl, data, category, showheadline = True):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')
        
        if showheadline:
            s = dddoc.DATA["globals.sections"][category].text()
            if s: fl.write('<div class=section_headline>' + s + '</div>')
        
        s = translateText(lines.text())
        if s: fl.write('<div class=text_block>' + s + '</div>')
    
        printTableContent(fl, data, category)

        fl.write('</div>');
  
################################################################################

def findDataRek(data, field, lines, map):
    highlight = False
    
    lines.extend(data[field].lines)

    for line in lines:
        map[line.text()] = ''
        
    follows = dddoc.DATA["globals.inherit"][field]
    for follow in follows.keys():
        follow_field = follows[follow].text()
        if follow_field.strip() == "": follow_field = field
        
        followups = data[follow].lines
        for followup in followups:
            dataup = dddoc.DATA[followup.text()]
            
            submap = {}
            findDataRek(dataup, follow_field, lines, submap)
            
            for key in submap.keys():
                if not map.has_key(key): 
                    origin = submap[key]
                    if origin == '': origin = followup.text()
                    map[key] = origin
                    highlight = True
        
    return highlight
      
################################################################################

def printMemberOut(fl, data, category, lines, derivedfrom, highlight, showheadline = True):
    if len(lines) > 0:
        fl.write('<div class="section" id="' + category + '">')

        if showheadline:
            s = dddoc.DATA["globals.sections"][category].text()
            if s: fl.write('<div class="section_headline">' + s + '</div>')
        
        fl.write('<table class="value_tab" cellspacing="0" cellpadding="0">')

        # Sort lines by their display link text and iterate over them.
        keyed_lines = dict([(translateLinkDisplaytext(line.text()).lower(), line) for line in lines])
        for key in sorted(keyed_lines.keys()):
            line = keyed_lines[key]
            text = line.text()
            origin = ''
            do_highlight = highlight
            if derivedfrom.has_key(text) and len(derivedfrom[text]) > 0: 
                origin = ' (' + translateLink(derivedfrom[text], line=line) + ')'
                do_highlight = False
                
            if do_highlight: tag_class = 'value_key_high'
            else: tag_class = 'value_key'

            link = translateLink(text, line=line)
            fl.write('<tr><td class="' + tag_class + '" valign="top"><nobr>' + link + '</nobr></td><td class="value_text" valign="top">')
                
            summary = translateText(dddoc.DATA[text]["summary"].text()) + origin
            if summary: fl.write(summary)
            else: fl.write('&nbsp;')
            
            fl.write('</td></tr>')
                
        fl.write('</table>')
            
        fl.write('</div>');

################################################################################

def printMember(fl, data, category, showheadline = True):
    lines = data[category].lines
    printMemberOut(fl, data, category, lines, {}, False, showheadline)

################################################################################

def printMemberRek(fl, data, category, showheadline = True):
    lines = []
    map = {}
    highlight = findDataRek(data, category, lines, map)
    
    printMemberOut(fl, data, category, lines, map, highlight, showheadline)
    

################################################################################

def printLinkOut(fl, data, category, lines, derivedfrom, highlight, showheadline = True):
    if not lines:
        return  # Nothing to do if no lines are given.
    fl.write('<div class="section" id="%s">' % category)

    if showheadline:
        s = dddoc.DATA["globals.sections"][category].text()
        if s: fl.write('<div class="section_headline">' + s + '</div>')

    my_dict = dict([(translateLinkDisplaytext(line.text()).lower(), line) for line in lines])
    str = ''
    for key, line in my_dict.iteritems():
        link = line.text()
        if (link == ''): continue
        origin = ''
        do_highlight = highlight
        if derivedfrom.has_key(link) and len(derivedfrom[link]) > 0: 
            origin = ' (' + translateLink(derivedfrom[link], line=line) + ')'
            do_highlight = False

        if do_highlight: tag_class = 'link_text_high'
        else: tag_class = 'link_text'

        s = translateLink(link, line=line)
        if s:
            if (str != ''): str += ', '
            str += '<span class="' + tag_class + '">' + s + '</span>'

    if (str != ''):
        fl.write('<div class="text_block">' + str + '</div>')

    fl.write('</div>')

################################################################################

def printLink(fl, data, category, showheadline = True):
    lines = data[category].at_level().lines
    printLinkOut(fl, data, category, lines, {}, False, showheadline)

################################################################################

def printLinkRek(fl, data, category, showheadline = True):
    lines = []
    map = {}
    highlight = findDataRek(data, category, lines, map)
    
    printLinkOut(fl, data, category, lines, map, highlight, showheadline)


################################################################################

def printList(fl, data, category, showheadline = True):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')

        if showheadline:
            s = dddoc.DATA["globals.sections"][category].text()
            if s: fl.write('<div class=section_headline>' + s + '</div>')

        map = {}
        for line in lines.lines:
            map[line.text()] = 1
        
        str = ''
        texts = map.keys()
        texts.sort()
        for text in texts:
            if (str != ''): str += ', '
            str += text
                
        if (str != ''):
            fl.write('<div class=text_block>' + str + '</div>');

        fl.write('</div>');

################################################################################

def printShortcutfor(fl, data, category, showheadline = True):
    printLink(fl, data, category, showheadline)
    printSignature(fl, data[category], "signature")

################################################################################

def printTextblock(fl, data, category, showheadline = True):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')
        
        if showheadline:
            s = dddoc.DATA["globals.sections"][category].text()
            if s: fl.write('<div class=section_headline>' + s + '</div>')
        
        fl.write('<div class=text_block>')
        subprintText(fl, lines)
        fl.write('</div>')
            
        fl.write('</div>');


################################################################################

def printGlossary(fl, data, category):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')
        
        keys = lines.keys()
        if len(keys) > 0:         
            for key in keys:
                fl.write('<div class=glossary_entry>')
                fl.write('<a name="GLOSSARY' + escapeFiles(key) + '"></a>')
                fl.write('<div class=glossary_title>' + translateText(key) + '</div>')
                fl.write('<div class=glossary_content>')
                subprintText(fl, lines[key])
                fl.write('</div>')
                fl.write('</div>')

        fl.write('</div>');


################################################################################

def printFile(fl, data, category):
    # Note: This somehow works on the demos page.
    global globalDocsPath
    global includeDirs
    
    filename = data[category].text()
    
    filename = filename.replace("\n", "")
    filename = filename.replace("\\", "/")

    # Return if the file name is empty.
    if not filename:
        return

    # Try to build the file name from the include dirs.
    for prefix in ['.'] + includeDirs:
        filenameCandidate = os.path.join(prefix, filename)
        if os.access(filenameCandidate, F_OK):
            filename = filenameCandidate
            break

    # Return if we cannot open the file.
    if (not os.access(filename, F_OK)):
        global WARNING_COUNT
        WARNING_COUNT += 1
        print
        print '!!  WARNING: unknown file "' + filename + '"'
        return

    # Read in file...
    f = open(filename)
    lines = f.readlines()
    f.close()

    linenumber = 0  # Of code, non-comment.
    codemode = False

    pos = filename.rfind("/")
    if (pos >= 0): 
        s = filename[pos+1:]
    else: 
        s = filename

    #copy file
    f_out = open(os.path.join(globalDocsPath, s), "w")

    fl.write('<div class=section_headline>File "<a href="' + s + '">' + s + '</a>"</div>')

    fl.write('<div class=codefile >')
    line_no = 0  # Absolute in file.
    for line in lines:
        line_no += 1
        is_comment = (line[0:3] == '///')

        if is_comment:
            if codemode:
                fl.write('</table><div class=comment>')    
                codemode = False
            line_obj = dddoc.Line([], line, filename, line_no)
            fl.write(translateText(line[3:], line_obj))

        else:
            if not codemode:
                if (len(line) <= 1): continue
                if linenumber: fl.write('</div>')
                fl.write('<table cellspacing=0 cellpadding=0 class=codefiletab>')
                codemode = True

            linenumber += 1

            fl.write('<tr>')
            fl.write('<td align=right class=linenumber>')
            fl.write(str(linenumber))
            fl.write('</td>')

            text = escapeHTML(line)
            text = text.replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
            text = text.replace(" ", "&nbsp;")
            text = text.replace("\n", "<br >")
            fl.write('<td class=content><nobr>')
            fl.write(text)
            fl.write('</nobr></td>')        

            fl.write('</tr>')

            f_out.write(line)

    if codemode:
        fl.write('</table>')
    else:
        fl.write('</div>')

    fl.write('</div>')    

    f_out.close

################################################################################

def _loadSnippet(path, snippet_key):
    result = []
    current_key = None
    current_lines = []
    with open(path, 'rb') as f:
        fcontents = f.read()
    for line in fcontents.splitlines():
        line = line.rstrip()  # Strip line ending and trailing whitespace.
        if line.strip().startswith('//![') and line.strip().endswith(']'):
            key = line.strip()[4:-1].strip()
            if key == current_key:
                if key == snippet_key:
                    result = current_lines
                current_lines = []
                current_key = None
            else:
                current_key = key
        elif current_key:
            current_lines.append(line)
    return result

def printSnippet(fl, data, category):
    # Note: This somehow works on the demos page.
    global globalDocsPath
    global includeDirs
    
    filename = data[category].text()
    
    filename = filename.replace("\n", "")
    filename = filename.replace("\\", "/")

    # Return if the file name is empty.
    if not filename:
        return

    snippet_id = '<none>'
    if '|' in filename:
        filename, snippet_id = filename.split('|', 1)

    # Try to build the file name from the include dirs.
    for prefix in ['.'] + includeDirs:
        filenameCandidate = os.path.join(prefix, filename)
        if os.access(filenameCandidate, F_OK):
            filename = filenameCandidate
            break

    # Return if we cannot open the file.
    if (not os.access(filename, F_OK)):
        global WARNING_COUNT
        WARNING_COUNT += 1
        print
        print '!!  WARNING: unknown file "' + filename + '"'
        return

    # Read in file...
    lines = _loadSnippet(filename, snippet_id)

    if not lines:
        print
        print '!!  WARNING: unknown snippet "' + snippet_id + '" in "' + filename + '"'
        return

    linenumber = 0  # Of code, non-comment.
    codemode = False

    pos = filename.rfind("/")
    if (pos >= 0): 
        s = filename[pos+1:]
    else: 
        s = filename

    #copy file
    with open(os.path.join(globalDocsPath, s), "w") as f_out:
        with open(filename, 'rb') as f2:
            f_out.write(f2.read())

    fl.write('<div class=codefile >')
    line_no = 0  # Absolute in file.
    for line in lines:
        line_no += 1
        is_comment = (line[0:3] == '///')

        if is_comment:
            if codemode:
                fl.write('</table><div class=comment>')    
                codemode = False
            line_obj = dddoc.Line([], line, filename, line_no)
            fl.write(translateText(line[3:], line_obj))

        else:
            if not codemode:
                if (len(line) <= 1): continue
                if linenumber: fl.write('</div>')
                fl.write('<table cellspacing=0 cellpadding=0 class=codefiletab>')
                codemode = True

            linenumber += 1

            fl.write('<tr>')
            fl.write('<td align=right class=linenumber>')
            fl.write(str(linenumber))
            fl.write('</td>')

            text = escapeHTML(line)
            text = text.replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
            text = text.replace(" ", "&nbsp;")
            text = text.replace("\n", "<br >")
            fl.write('<td class=content><nobr>')
            fl.write(text)
            fl.write('</nobr></td>')        

            fl.write('</tr>')

    if codemode:
        fl.write('</table>')
    else:
        fl.write('</div>')

    fl.write('</div>')    
    fl.write('<div class=section_headline>Snippet from "<a href="' + s + '">' + s + '</a>"</div>')


################################################################################

def subprintText(fl, data, subcategory = False):
    if data.empty():
        return
        
    headline = ''
    if subcategory:
        s = dddoc.DATA["globals.subsections"][subcategory].text()
        if s: headline = '<span class=section_sub_headline>' + s + '</span>'
    
    subprintText(fl, data["summary"])    

    for line in data.at_level(0).lines:
        s = translateText(line.text())
        if s: 
            fl.write('<div class=text_sub_block>' + headline + ' ' + s + '</div>')
            headline = ''
            
    in_table = False
#    in_ol = False
            
    section_cout = 0
    subsection_cout = 0
    
    for line in data.at_level(1).by_occ().lines:
        name = line.name(data.level)
        
        if (name == 'table') or (name == 'tableheader'):
            if not in_table:
                fl.write('<table class=table_explicite cellspacing=0 cellpadding=0>')
            subprintTableLine(fl, line.text(), (name == 'tableheader'), line)
            in_table = True
        else:
            if in_table:
                fl.write('</table>')
            in_table = False

#        if name == 'enumerate':
#            if not in_ol:
#                fl.write('<ol>')
#            fl.write('<li>' + translateText(line.text()))
#            in_ol = True
#        else:
#            if in_ol:
#                fl.write('</ol>')
#            in_ol = False
            
        if name == 'contents':
            subprintContents(fl, data)
            
        elif name == 'section': 
            section_cout += 1
            subsection_cout = 0
            t = translateSection(line.text(), section_cout)
            fl.write('<div class=section_headline_explicite><a name="' + t + '"></a>' + t + '</div>')
            headline = ''
            
        elif name == 'subsection': 
            subsection_cout += 1
            t = translateSubsection(line.text(), section_cout, subsection_cout)
            fl.write('<div class=section_sub_headline_explicite><a name="' + t + '"></a>' + t + '</div>')
            headline = ''
            
        elif name == 'text': 
            s = translateText(line.text(), line=line)
            fl.write('<div class=text_sub_block>' + headline + ' ' + s + '</div>')
            headline = ''
            
        elif name == 'code': 
            s = translateCode(line.text())
            fl.write('<div class=code_sub_block>' + s + '</div>')
            
        elif name == 'snippet':
            printSnippet(fl, data, 'snippet')
            #s = translateCode(line.text())
            #fl.write('<div class=code_sub_block>' + s + '</div>')
            
        elif name == 'output':
            s = translateCode(line.text())
            fl.write('<div class=output_sub_block>' + s + '</div>')
            
        elif name == 'image': 
            subprintImage(fl, line.text())
            
        elif name == 'note':
            s = dddoc.DATA["globals.subsections.note"].text()
            s = '<span class=section_sub_headline>' + s + '</span> '
            s += translateText(line.text())
            fl.write('<div class=note_sub_block>' + s + '</div>')
            
        elif name == 'field':
            subprintField(fl, line.text())

            
    if in_table:
        fl.write('</table>')
       
    subprintLink(fl, data["metafunction"], "metafunction")
    subprintConceptAndType(fl, data)
    subprintText(fl, data["value"], "value")
    subprintText(fl, data["default"], "default")
    printTableContent(fl, data, "param")
    subprintText(fl, data["remarks"], "remarks")    
    subprintLink(fl, data["see"], "see")


################################################################################

def translateSection(text, section_count):
    i = text.find('#')
    if (i >= 0):
        text = text[:i] + str(section_count) + text[i+1:]
    return translateText(text);

################################################################################

def translateSubsection(text, section_count, subsection_count):
    text = translateText(text)
    i = text.find('#')
    if (i >= 0):
        j = text.find('#', i+1)
        if (j >= 0):
            text = text[:i] + str(section_count) + text[i+1:j] + str(subsection_count) + text[j+1:]
        else:
            text = text[:i] + str(subsection_count) + text[i+1:]
    return translateText(text);

################################################################################

def getLinkList(data_lines, not_name_types = {}):
    pairs = [(translateLinkDisplaytext(line.text()).lower(), line) for line in data_lines]
    
    str = ''
    
    for key, line in sorted(pairs):
        link = line.text()
        if not_name_types.has_key(link): continue
        s = translateLink(link, line=line)
        if s:
            if (str != ''): str += ', '
            str += s
            
    return str

################################################################################

def subprintLink(fl, data, subcategory):                         
    if data.empty():
        return
        
    headline = ''
    if subcategory:
        s = dddoc.DATA["globals.subsections"][subcategory].text()
        if s: headline = '<span class=section_sub_headline>' + s + '</span>'
        
    str = getLinkList(data.lines)
            
    if (str != ''):
        fl.write('<div class=text_sub_block>' + headline + ' ' + str + '</div>')

################################################################################


def subprintConceptAndType(fl, data):    
    not_name_types = {}
    type_title = ''

    #display concepts
    c_data = data["concept"]
    t_data = data["type"]
    
    if not c_data.empty():
        if len(c_data.lines) + len(t_data.lines) == 1: tag = 'span'
        else: tag = 'div'
        
        fl.write('<div class=text_sub_block id=concept_block>');

        title = dddoc.DATA["globals.subsections.concept"].text()
        if title: fl.write('<' + tag + ' class=section_sub_headline>' + title + '</' + tag + '> ')
        
        pairs = [(translateLinkDisplaytext(line.text()).lower(), line) for line in c_data.lines]
        for key, line in sorted(pairs):
            link = line.text()
            str_concept = translateLink(link, line=line)
            if str_concept:
                lines = []
                findDataRek(dddoc.DATA[link], "conceptimplements", lines, {})
                str_types = getLinkList(lines)
                if (str_types != ''): str_types = ': ' + str_types
                fl.write('<' + tag + ' class=text_sub_block>' + str_concept + str_types + '</' + tag + '>')
                
                for lin in lines:
                    not_name_types[lin.text()] = 1
                
        subprintType(fl, t_data, '', not_name_types)
        
        fl.write('</div>')
        
    else:
        s = dddoc.DATA["globals.subsections.type"].text()
        if s: type_title = '<span class=section_sub_headline>' + s + '</span> '
        
        subprintType(fl, t_data, type_title, not_name_types)


################################################################################

def subprintType(fl, data, headline, not_name_types):                         
    if data.empty():
        return
        
    str = getLinkList(data.lines, not_name_types)
            
    if (str != ''):
        fl.write('<div class=text_sub_block>' + headline + ' ' + str + '</div>')


################################################################################

def subprintTableLine(fl, text, is_header, line=None):
    fl.write('<tr>')

    t = ''

    while len(text) > 0:
        i = text.find('|');
        j = text.find('@');

        if (j >= 0) and (j < i):
            j2 = text.find('@', j+1)
            if (j2 >= 0):
                t += text[:j2+1]
                text = text[j2+1:]
                continue;

        if (i >= 0):
            t += text[:i]
            text = text[i+1:]
        else:
            t += text
            text = ''
        
        if len(t) > 0: 
            t = translateText(t, line=line)
        else: 
            t = '&nbsp;'

        if is_header:
            fl.write('<td class=table_header_explicite valign=top><center>' + t + '</center></td>')
        else:
            fl.write('<td class=table_cell_explicite valign=top>' + t + '</td>')
        
        t = ''
        
    fl.write('</tr>')

################################################################################

def subprintImage(fl, text):
    t = ''
    s1 = ''
    s2 = ''
    is_image = True
    
    while len(text) > 0:
        i = text.find('|');
        j = text.find('@');
        
        if (j >= 0) and (j < i):
            j2 = text.find('@', j+1)
            if (j2 >= 0):
                t += text[:j2+1]
                text = text[j2+1:]
                continue;
        
        if i >= 0:
            t += text[:i]
            text = text[i+1:]
        else:
            t += text
            text = ''
        
        if len(t) > 0:
            if is_image:
                s1 += '<td>' + translateImage(t) + '</td>'
            else:
                s2 += '<td valign=top align=left class=image_sub_block_caption>' + translateText(t) + '</td>'
        else:
            if is_image:
                s1 += '<td>&nbsp;</td>'
            else:
                s2 += '<td>&nbsp;</td>'
        
        t = ''
        is_image = not is_image
        
    s = '<center><table cellspacing=0 cellpadding=0 border=0><tr>' + s1 + '</tr><tr>' + s2 + '</tr></table></center>'
    fl.write('<div class=image_sub_block>' + s + '</div>')

################################################################################

def subprintContents(fl, data):
    s = ''
    
    section_cout = 0
    subsection_cout = 0
    
    for line in data.at_level(1).by_occ().lines:
        name = line.name(data.level)
        
        if name == 'section':
            section_cout += 1
            subsection_cout = 0
            t = translateSection(line.text(), section_cout)
            s += '<div class=contents_section><a href="#' + t + '">' + t + '</a></div>'
            headline = ''
            
        if name == 'subsection': 
            subsection_cout += 1
            t = translateSubsection(line.text(), section_cout, subsection_cout)
            s += '<div class=contents_subsection><a href="#' + t + '">' + t + '</a></div>'
            headline = ''
    
    
    t = dddoc.DATA["globals.subsections.contents"].text()
    if t: s = '<div class=contents_headline><center>' + t + '</center></div>' + s
    
    fl.write('<div class=contents>' + s + '</div>');

################################################################################

def subprintField(fl, text):
    i = text.rfind('.')
    entry = text[:i]
    field = text[i+1:]

    data = dddoc.DATA[entry]

    if   (field == "description"): printTextblock(fl, data, "description", False)
    elif (field == "signature"): printSignature(fl, data, "signature")
    elif (field == "param"): printTable(fl, data, "param", False)
    elif (field == "returns"): printTextblock(fl, data, "returns", False)
    elif (field == "class"): printLink(fl, data, "class", False)
    elif (field == "general"): printLink(fl, data, "general", False)
    elif (field == "shortcutfor"): printShortcutfor(fl, data, "shortcutfor", False)
    elif (field == "implements"): printLinkRek(fl, data, "implements", False)
    elif (field == "baseconcept"): printLink(fl, data, "baseconcept", False)

    elif (field == "spec"): printMember(fl, data, "spec", False)
    elif (field == "shortcut"): printMember(fl, data, "shortcut", False)
    elif (field == "type"): printMemberRek(fl, data, "type", False)
    elif (field == "typedef"): printMemberRek(fl, data, "typedef", False)
    elif (field == "memvar"): printMemberRek(fl, data, "memvar", False)
    elif (field == "memfunc"): printMemberRek(fl, data, "memfunc", False)
    elif (field == "function"): printMemberRek(fl, data, "function", False)
    
    elif (field == "childconcept"): printMember(fl, data, "childconcept", False)
    elif (field == "conceptimplements"): printMemberRek(fl, data, "conceptimplements", False)
    elif (field == "conceptmetafunc"): printMemberRek(fl, data, "conceptmetafunc", False)
    elif (field == "concepttypedef"): printMemberRek(fl, data, "concepttypedef", False)
    elif (field == "conceptmemvar"): printMemberRek(fl, data, "conceptmemvar", False)
    elif (field == "conceptmemfunc"): printMemberRek(fl, data, "conceptmemfunc", False)
    elif (field == "conceptfunc"): printMemberRek(fl, data, "conceptfunc", False)
    
    elif (field == "value"): printTable(fl, data, "value", False)

    elif (field == "conceptusedbymeta"): printLinkRek(fl, data, "conceptusedbymeta", False)
    elif (field == "conceptusedbyfunc"): printLinkRek(fl, data, "conceptusedbyfunc", False)

    elif (field == "remarks"): printTextblock(fl, data, "remarks", False)
    elif (field == "example"): printTextblock(fl, data, "example", False)
    elif (field == "demo"): printLinkRek(fl, data, "demo", False)
    elif (field == "file"): printFile(fl, data, "file")
    elif (field == "snippet"): printSnippet(fl, data, "snippet")
    
    elif (field == "concept"): printLink(fl, data, "concept", False)
    elif (field == "status"): printTextblock(fl, data, "status", False)
    elif (field == "include"): printList(fl, data, "include", False)
    elif (field == "demofor"): printLink(fl, data, "demofor", False)
    elif (field == "see"): printLink(fl, data, "see", False)


        

################################################################################

def gatherGlossary():
    global globalGlossary 
    
    globalGlossary = {}
    got_it = {}
    
    print "Gather Glossary:",
    
    lines = dddoc.DATA.lines;
    for line in lines:
        key = line.name(3)
        if (line.name(2) == 'glossary') and (key != '(unknown)'):
            fname = line.name(0) + '.' + line.name(1) + '.glossary.' + key
            if not got_it.has_key(fname):
                got_it[fname] = 1
                print ".",
               
                href = getFilename(line.name(0), line.name(1)) + '#GLOSSARY' + escapeFiles(key)
                link = 'href="' + href + '" title="' + escapeJavaScript(dddoc.DATA[fname].text()) + '">'
                
                if not globalGlossary.has_key(key):  globalGlossary[key] = []
                globalGlossary[key].append([key, '(Glossary)', link])
                
    print

################################################################################

def createSearchfile(path):
    print 'Create Searchfile'
    
    global globalGlossary
    
    db = globalGlossary
    
    cats = dddoc.DATA["globals.categories"].keys()
    for cat in cats:
        entries = dddoc.DATA[cat]
        title = getCategoryTitle(cat).lower()
        pushSearchResult(db, title, cat, "")
       
        for key in entries.keys():
            data = entries[key]
            title = getPageTitle(data).lower()
            pushSearchResult(db, title, cat, key)

    searchfile = os.path.join(path, "searchfile.js")
    fl = file(searchfile, "w")
    
    fl.write('var DB = [\n')

    keys = db.keys()
    keys.sort()
    for key in keys:
        if len(key) > 0:
            for entry in db[key]:
                fl.write('[\'' + escapeJavaScript(entry[0]) + '\', \'' + escapeJavaScript(entry[1]) + '\', \'' + escapeJavaScript(entry[2]) + '\'],\n')
            
    fl.write('false];\n')
    fl.close()

################################################################################

def pushSearchResult(db, title, cat, name):
    if (len(name) == 0):
        link = '<a target=_parent href="' + getIndexpage(cat) + '">'
        key = getCategoryTitle(cat)
        text = ''
        
    else:
        obj = dddoc.DATA[cat][name]
        href = getFilename(cat, name)
        if obj.empty(): 
            summary = '';
        else:
            summary = translateTooltip(obj["summary"].text())
        if len(summary) > 0:
            summary = 'title="' + summary + '"'
        
        link = 'href="' + href + '" ' + summary + '>'
        key = translateID(name, cat=True)
        text =  '(' + cat + ')'

    if not db.has_key(title): db[title] = []
    db[title].append([key, text, link])

################################################################################

def warningPage(cat, key, data):
    # Warning for multiple or no summary field.
    global globalBuildFull
    global WARNING_COUNT
    if globalBuildFull: # only if full documentation is built
        desc = data["summary"]
        if desc.empty():
            WARNING_COUNT += 1
            print
            print "\n!!  WARNING: no summary field for \"" + cat + "." + key + "\""
            print '        Location: %s:%d' % (data.lines[0].file_name, data.lines[0].line_no)
        elif len(desc.lines) > 1:
            WARNING_COUNT += 1
            print
            print "\n!!  WARNING: multiple summary fields for \"" + cat + "." + key + "\""
            for line in desc.lines:
                print '        Location: %s:%d' % (line.file_name, line.line_no)

    # Warning for break of naming convention.
    convention = dddoc.DATA["globals.namingconventions"][cat].text()
    title = getPageTitle(data)
    
    if (convention == 'bigsmall' and ((title[0] < 'A') or (title[0] > 'Z'))):
        # Ugly hack: Allowing to violate bigsmall convention for tag groups that
        # have the substring ' Tags'
        if ' Tags' not in title:
            WARNING_COUNT += 1
            print
            print "\n!!  WARNING: \"" + title + "\" breaks naming convention: " + cat + " must start with capital letter."
            print '        Location: %s:%d' % (data.lines[0].file_name, data.lines[0].line_no)
    elif (convention == 'smallbig' and ((title[0] < 'a') or (title[0] > 'z'))):
        WARNING_COUNT += 1
        print
        print "\n!!  WARNING: \"" + title + "\" breaks naming convention: " + cat + " must start with lower case."
        print '        Location: %s:%d' % (data.lines[0].file_name, data.lines[0].line_no)
    elif (title[len(title)-1] == '_'):
        WARNING_COUNT += 1
        print
        print "\n!!  WARNING: \"" + title + "\" breaks naming convention: public identifiers must not end with \"_\"." # and private identifiers should not be documented
        print '        Location: %s:%d' % (data.lines[0].file_name, data.lines[0].line_no)


    #warning for unknown params
    sigs = data["signature"].text()
    sigs = sigs.strip(" \n")
    if len(sigs) > 0:
        params = data["param"]
        if params.empty():
            return
        keys = params.keys_by_occ()
        for k in keys:
            if sigs.find(k) < 0:
                # Maybe ignore the "unknown param" warning for this parameter
                # if the "nowarn" child includes "unknown param".
                if "unknown param" in data["param.%s.nowarn" % k].text():
                    continue
                WARNING_COUNT += 1
                print
                print "\n!!  WARNING: unknown param \"" + k + "\" in \"" + cat + "." + key + "\""
                line = params.find(k).lines[0]
                print '        Location: %s:%d' % (line.file_name, line.line_no)
                        
                        
################################################################################

def printTree(fl, data):
    MAX_NUMBER_OF_NODES = 6

    cat = data.name(0)
    treedown_fields = dddoc.DATA["globals.treedown"][cat]
    treeup_fields= dddoc.DATA["globals.treeup"][cat]
    
    if treedown_fields.empty() and treeup_fields.empty(): return
    
    treedown = followTree(data, treedown_fields)
    treeup = followTree(data, treeup_fields)

    if len(treedown) + len(treeup) == 0: return;
    fl.write("<div class=tree><center><table cellspacing=0 cellpadding=0>")
    if len(treeup) > 0:
        fl.write("<tr><td class=tree_td_subtree colspan=2 align=middle>" + translateTree(treeup, False, MAX_NUMBER_OF_NODES) + "</td></tr>")
        fl.write("<tr><td class=tree_td_line><img src=dddoc_empty.gif border=0 /></td><td class=tree_td_none><img src=dddoc_empty.gif border=0 /></td></tr>")

    fl.write("<tr><td class=tree_td_central colspan=2><center><div class=tree_central id=\"tree_node_" + cat + "\">")
    fl.write(getPageTitle(data))
    fl.write("</div></center></td></tr>")
    
    if len(treedown) > 0:
        fl.write("<tr><td class=tree_td_line><img src=dddoc_empty.gif border=0 /></td><td class=tree_td_none><img src=dddoc_empty.gif border=0 /></td></tr>")
        fl.write("<tr><td class=tree_td_subtree colspan=2 align=middle>" + translateTree(treedown, True, MAX_NUMBER_OF_NODES) + "</td></tr>")
        
        
    fl.write("</table></center></div>")

################################################################################

def followTree(data, follow_fields):
    tree = {}
    if data.empty(): return tree
    
    for field_line in follow_fields.lines:
        field = field_line.text()
        children = data[field]
        if not children.empty():
            for child_line in children.lines:
                child = child_line.text()
                if not tree.has_key(child):
                    tree[child] = followTree(dddoc.DATA[child], follow_fields)
                    
    return tree

################################################################################

def translateTree(tree, print_down, max_num_of_nodes):
    
    tr1 = ''
    tr2 = ''
    tr3 = ''
    tr4 = ''
    multiblock = ''
    i = 0

    keys = tree.keys()
    keys.sort()
    
    new_max_num_of_nodes = round((max_num_of_nodes + len(keys) - 1)/ len(keys))
    
    if print_down: valign = "valign=top"
    else: valign = "valign=bottom"

    for key in keys:
        if (i > 0) and (i % max_num_of_nodes == 0):
            separator_colspan = (2 *max_num_of_nodes).__str__()
            if print_down: 
                multiblock = multiblock + '<tr><td class=tree_td_connector rowspan = 5><img src=dddoc_empty.gif border=0 /></td>' + tr1 + '</tr><tr>' + tr2 + '</tr><tr>' + tr3 + '</tr><tr>' + tr4 + '</tr><tr><td colspan=' + separator_colspan + '><img id=tree_block_separator src=dddoc_empty.gif border=0/></td></tr>'
            else: 
                multiblock = '<tr><td class=tree_td_connector rowspan = 5><img src=dddoc_empty.gif border=0 /></td><td colspan=' + separator_colspan + '><img id=tree_block_separator src=dddoc_empty.gif border=0/></td></tr><tr>' + tr4 + '</tr><tr>' + tr3 + '</tr><tr>' + tr2 + '</tr><tr>' + tr1 + '</tr>' + multiblock
            tr1 = ''
            tr2 = ''
            tr3 = ''
            tr4 = ''
            
        i += 1

        if len(tree) > 1:
            if (i % max_num_of_nodes == 1): left_class = "tree_td_first"
            else: left_class = "tree_td_left"
            if (i % max_num_of_nodes == 0) or (i == len(tree)): right_class = "tree_td_last"
            else: right_class = "tree_td_right"
            tr1 += "<td class=" + left_class + "><img src=dddoc_empty.gif border=0 /></td><td class=" + right_class + "><img src=dddoc_empty.gif border=0 /></td>"

        
        pos = key.find('.')
        if pos > 0: node_cat = key[:pos]
        else: node_cat = key
        tr2 += "<td class=tree_td_node colspan=2><nobr><center><span class=tree_node id=\"tree_node_" + node_cat + "\">"
        tr2 += translateLink(key)
        tr2 += "</span></center></nobr></td>"
        
        subtree = tree[key]
        if len(subtree) == 0: 
            tr3 += "<td class=tree_td_none></td><td class=tree_td_none></td>"
            tr4 += "<td class=tree_td_subtree colspan=2 align=middle></td>"
        else: 
            tr3 += "<td class=tree_td_line><img src=dddoc_empty.gif border=0 /></td><td class=tree_td_none><img src=dddoc_empty.gif border=0 /></td>"
            tr4 += "<td class=tree_td_subtree colspan=2 align=middle " + valign + ">" + translateTree(subtree, print_down, new_max_num_of_nodes) + "</td>"
            

    if multiblock: connector = '<td class=tree_td_connector_end rowspan=4><img src=dddoc_empty.gif border=0 /></td>'
    else: connector = ''
               
    if print_down: 
        str = '<table class=tree_table_down cellspacing=0 cellpadding=0>'
        str += multiblock + '<tr>' + connector + tr1 + '</tr><tr>' + tr2 + '</tr><tr>' + tr3 + '</tr><tr>' + tr4 + '</tr>'
        str += '</table>'
    else: 
        str = '<table class=tree_table_up cellspacing=0 cellpadding=0>'
        str += '<tr>' + connector + tr4 + '</tr><tr>' + tr3 + '</tr><tr>' + tr2 + '</tr><tr>' + tr1 + '</tr>' + multiblock
        str += '</table>'
        
    if multiblock: return "<div id=tree_multiblock>" + str + "</div>"
    else:return str
                

    
