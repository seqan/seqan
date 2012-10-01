import os
import dddoc
import dddoc_html_trans
import xml.sax
import operator
from os import F_OK


################################################################################

def createDocs(path):
    fl = file("temp.tex", "w")
    createPages(fl)
    fl.close()

################################################################################

def createPages(fl):
    cats = dddoc.DATA["globals.categories"].keys_by_occ()
    for cat in cats:
        print 'Pages for ' + cat,
        
        pageIndex(fl, cat)
        
        entries = dddoc.DATA[cat]
        keys = entries.keys().sort()
        for key in keys:
            data = entries[key]
            print '.',
            pageContent(fl, data)
            
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
            
################################################################################

def escapeTeX(text):

    return text
    
################################################################################
    
def translate(text):
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
                ret += translateCode(str)
                str = ''
                in_code = False
            else: str += c
        elif in_link:
            if c == '@': 
                ret += translateLink(str)
                str = ''
                in_link = False
            else: str += c
        else:
            if c == '$': 
                ret += translateText(str)
                str = ''
                in_code = True
            elif c == '@':
                ret += translateText(str)
                str = ''
                in_link = True
#            elif c == '\\':
#                in_escaped = True
            else:
                str += c;
                
        pos += 1
                
    if str != '':
        ret += translateText(str)

    return ret
    
################################################################################
    
def translateText(text):

    return escapeTeX(text)
 
################################################################################
   
def translateCode(text):
    return '\\begin{lstlisting}' + text + '\\end{lstlisting}'
    
################################################################################

def brokenLink(text):
    print
    print '   WARNING: broken link "' + text + '"'
    return text
    
################################################################################

def translateLink(text):
    arr = dddoc.splitName(text)
    if len(arr) < 2: return brokenLink(text)

    obj = dddoc.DATA[arr[0]][arr[1]];
    if obj.empty(): return brokenLink(text);
    
    summary = translateTooltip(obj["summary"].text());
    if len(summary) > 0:
        summary = 'title="' + summary + '"'
    
    href = getFilename(arr[0], arr[1])

    return '\\url(<a href="' + href + '" ' + summary + '>' + translateID(arr[len(arr) - 1]) + '</a>'

################################################################################
    
def translateImage(text):
    text = escapeHTML(text)
    text = text.replace("\t", "&nbsp;&nbsp;&nbsp;&nbsp;")
    text = text.replace(" ", "&nbsp;")
    text = text.replace("\n", "<br >")
    return '<img class=image src="../img/' + text + '.png" border=0 />'


################################################################################
    
def translateID(text):
    i = text.find('#');
    if (i >= 0):
        text = text[i + 1: len(text)]
    return translateText(text)
        
################################################################################

def getBeforeColon(text):
	return text[0: text.find(':')]

################################################################################

def getAfterColon(text):
	return text[text.find(':') + 1:len(text)]

################################################################################

def getPageTitle(data):
	s = data["title"].text()
	if (s == ''): 
		return translateID(data.name(1))
	else:
		return s
	
################################################################################
	
def getCategoryTitle(cat):
	s = dddoc.DATA["Indexpage"][cat]["title"].text()
	if (s == ''): 
		return translateID(dddoc.DATA["globals.indexes"][cat].text())
	else:
		return s

################################################################################

def pageIndex(fl, cat, subcat = ""):
    fl.write('<html>')
    fl.write('<head>')
    fl.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8">');
    fl.write('<link rel="stylesheet" href="dddoc_html.css" type="text/css" />')
    fl.write('</head>')
    fl.write('<body id=index_body>')
    
    lines = dddoc.DATA["globals.indexes"]
    for cat2 in lines.keys_by_occ():
        if cat2 == cat:
            fl.write('<div class=index_section_high>')
        else:
            fl.write('<div class=index_section>')
        fl.write('<div class=index_cat><a class=index_link target=_top href="' + getIndexpage(cat2) + '">' + lines[cat2].text() + '</a></div>')
        if cat2 == cat:
            keys = dddoc.DATA[cat2].keys_by_occ()
            entries = {}
            for key in keys:
                data = dddoc.DATA[cat2][key]
                if not data["hidefromindex"].empty():
                    continue
                    
                subcats2 = data["cat"].lines
                
                map = {}
                if (len(subcats2) == 0): map[""] = 1
                else:
                    for line in subcats2:
                        map[line.text()] = 1
        
                for subcat2 in map.keys():
                    if not entries.has_key(subcat2): entries[subcat2] = {}
                    if subcat == subcat2:
                        s = '<div class=index_item>'
                        s += '<a class=index_link id="' + key + '" target=_top title="' + translateID(key) + '" href="' + getFilename(cat, key) + '">' + translateID(key) + '</a>'
                        s += '</div>'
                        if not entries[subcat2].has_key(translateID(key)): entries[subcat2][translateID(key)] = ''
                        entries[subcat2][translateID(key)] += s
                        
            keys = entries.keys()
            keys.sort()
            for key in keys:
            
                text = ''
                subkeys = entries[key].keys()
                subkeys.sort()
                for subkey in subkeys:
                    text += entries[key][subkey]

                if key != "":
                    fl.write('<div class=index_subsection>')
                    fl.write('<div class=index_subcat><a class=index_link href="' + getIndexname(cat2, key) + '"></div>')
                    if key == subcat:
                        fl.write('<img src="dddoc_minus.gif" border=0>')
                    else:
                        fl.write('<img src="dddoc_plus.gif" border=0>')
                    fl.write(key + '</a>')
                    fl.write(text)
                    fl.write('</div>')
                else:
                    fl.write(text)
        fl.write('</div>')
    fl.write('</body>')
    fl.write('</html>')

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
    fl.write('</body>')
    fl.write('</html>')

################################################################################

def addIndexPageMembers(data, key, entries, subcat):
    if not entries.has_key(subcat): 
        entries[subcat] = {}
        
    s = '<tr><td class=value_key valign=top><nobr>'
    s += '<a href="' + getFilename(data.name(0), key) + '">' + translateID(key) + '</a>'
    s += '</nobr></td><td class=value_text valign=top>'

    summary = translateText(data[key]["summary"].text())
    if summary: s += summary
    else: s += '&nbsp;'
    
    s += '</td></tr>'
    
    if not entries[subcat].has_key(translateID(key)): entries[subcat][translateID(key)] = ''
    entries[subcat][translateID(key)] += s

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
                fl.write('<div class=section_headline>' + key + '</div>')
            else:
                s = dddoc.DATA["globals.indexes"][data.name(0)].text()
                if s: fl.write('<div class=section_headline>' + s + '</div>')
                
            fl.write('<table class=value_tab cellspacing=0 cellpadding=0>')
            
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
    if ((data.name(0) == 'Memfunc') or (data.name(0) == 'Memvar')):
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
    printTree(fl, data)
    printSignature(fl, data, "signature")
    printTable(fl, data, "param")
    printTextblock(fl, data, "returns")
    printLink(fl, data, "class")
    printLink(fl, data, "general")
    printShortcutfor(fl, data, "shortcutfor")
    printLinkRek(fl, data, "implements", "general")
    printLink(fl, data, "baseconcept")

    printMember(fl, data, "spec")
    printMember(fl, data, "shortcut")
    printMemberRek(fl, data, "type", "general")
    printMemberRek(fl, data, "memvar", "base")
    printMemberRek(fl, data, "memfunc", "base")
    printMemberRek(fl, data, "function", "general")
    
    printMember(fl, data, "childconcept")
    printMemberRek(fl, data, "conceptimplements", "childconcept")
    printMemberRek(fl, data, "conceptmetafunc", "baseconcept")
    printMemberRek(fl, data, "conceptmemvar", "baseconcept")
    printMemberRek(fl, data, "conceptmemfunc", "baseconcept")
    printMemberRek(fl, data, "conceptfunc", "baseconcept")
    
    printTable(fl, data, "value")

    printLinkRek(fl, data, "conceptusedby", "childconcept")

    printTextblock(fl, data, "remarks")
    printTextblock(fl, data, "example")
    printLinkRek(fl, data, "demo", "general")
    printFile(fl, data, "file")
    
    printLink(fl, data, "concept")
    printList(fl, data, "include")
    printLink(fl, data, "demofor")
    printLink(fl, data, "see")

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


def printTable(fl, data, category):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')
        
        s = dddoc.DATA["globals.sections"][category].text()
        if s: fl.write('<div class=section_headline>' + s + '</div>')
        
        s = translateText(lines.text())
        if s: fl.write('<div class=text_block>' + s + '</div>')
    
        printTableContent(fl, data, category)

        fl.write('</div>');
  
################################################################################

def findDataRek(fl, data, category, follow, lines, map):
    highlight = False
    
    lines.extend(data[category].lines)

    for line in lines:
        map[line.text()] = ''
    
    followups = data[follow].lines
    for followup in followups:
        dataup = dddoc.DATA[followup.text()]
        
        submap = {}
        findDataRek(fl, dataup, category, follow, lines, submap)
        
        for key in submap.keys():
            if not map.has_key(key): 
                origin = submap[key]
                if origin == '': origin = followup.text()
                map[key] = origin
                highlight = True
        
    return highlight
      
################################################################################

def printMemberOut(fl, data, category, lines, derivedfrom, highlight):
    if len(lines) > 0:
        fl.write('<div class=section id=' + category + '>')

        s = dddoc.DATA["globals.sections"][category].text()
        if s: fl.write('<div class=section_headline>' + s + '</div>')
        
        fl.write('<table class=value_tab cellspacing=0 cellpadding=0>')
        
        map = {}
        for line in lines:
            map[line.text()] = 1
        
        texts = map.keys()
        texts.sort()
        for text in texts:
            origin = ''
            do_highlight = highlight
            if derivedfrom.has_key(text) and len(derivedfrom[text]) > 0: 
                origin = ' (' + translateLink(derivedfrom[text]) + ')'
                do_highlight = False
                
            if do_highlight: tag_class = 'value_key_high'
            else: tag_class = 'value_key'

            link = translateLink(text)        
            fl.write('<tr><td class=' + tag_class + ' valign=top><nobr>' + link + '</nobr></td><td class=value_text valign=top>')
                
            summary = translateText(dddoc.DATA[text]["summary"].text()) + origin
            if summary: fl.write(summary)
            else: fl.write('&nbsp;')
            
            fl.write('</td></tr>')
                
        fl.write('</table>')
            
        fl.write('</div>');

################################################################################

def printMember(fl, data, category):
    lines = data[category].lines
    printMemberOut(fl, data, category, lines, {}, False)

################################################################################

def printMemberRek(fl, data, category, follow):
    lines = []
    map = {}
    highlight = findDataRek(fl, data, category, follow, lines, map)
    
    printMemberOut(fl, data, category, lines, map, highlight)
    

################################################################################

def printLinkOut(fl, data, category, lines, derivedfrom, highlight):
    if len(lines) > 0:
        fl.write('<div class=section id=' + category + '>')

        s = dddoc.DATA["globals.sections"][category].text()
        if s: fl.write('<div class=section_headline>' + s + '</div>')

        map = {}
        for line in lines:
            map[line.text()] = 1
        
        str = ''
        links = map.keys()
        links.sort()
        for link in links:
            if (link == ''): continue
            origin = ''
            do_highlight = highlight
            if derivedfrom.has_key(link) and len(derivedfrom[link]) > 0: 
                origin = ' (' + translateLink(derivedfrom[link]) + ')'
                do_highlight = False
                
            if do_highlight: tag_class = 'link_text_high'
            else: tag_class = 'link_text'

            s = translateLink(link)
            if s:
                if (str != ''): str += ', '
                str += '<span class=' + tag_class + '>' + s + '</span>'
                
        if (str != ''):
            fl.write('<div class=text_block>' + str + '</div>');

################################################################################

def printLink(fl, data, category):
    lines = data[category].at_level().lines
    printLinkOut(fl, data, category, lines, {}, False)

################################################################################

def printLinkRek(fl, data, category, follow):
    lines = []
    map = {}
    highlight = findDataRek(fl, data, category, follow, lines, map)
    
    printLinkOut(fl, data, category, lines, map, highlight)


################################################################################

def printList(fl, data, category):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')

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

################################################################################

def printShortcutfor(fl, data, category):
    printLink(fl, data, category)
    printSignature(fl, data[category], "signature")

################################################################################

def printTextblock(fl, data, category):
    lines = data[category]
    if not lines.empty():
        fl.write('<div class=section id=' + category + '>')

        s = dddoc.DATA["globals.sections"][category].text()
        if s: fl.write('<div class=section_headline>' + s + '</div>')
        
        fl.write('<div class=text_block>')
        subprintText(fl, lines)
        fl.write('</div>')
            
        fl.write('</div>');

################################################################################

def printFile(fl, data, category):
    filename = data[category].text()
    
    filename = filename.replace("\n", "")
    
    if (filename != ''):
        if (not os.access(filename, F_OK)):
            print '  WARNING: unknown file "' + filename + '"'
        else:
            f = open(filename)
            lines = f.readlines()
            f.close()
            
            linenumber = 0
            codemode = False

            pos = filename.rfind("/")
            if (pos >= 0): s = filename[pos+1:]
            else: s = filename
            fl.write('<div class=section_headline>File "' + s + '"</div>')
        
            fl.write('<div class=codefile >')
            for line in lines:
                is_comment = (line[0:3] == '///')
                
                if is_comment:
                    if codemode:
                        fl.write('</table><div class=comment>')    
                        codemode = False
                    fl.write(translateText(line[3:]))
                    
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
            
    for line in data.at_level(1).by_occ().lines:
        name = line.name(data.level)
        if name == 'text': 
            s = translateText(line.text())
            fl.write('<div class=text_sub_block>' + headline + ' ' + s + '</div>')
            headline = ''
            
        elif name == 'code': 
            s = translateCode(line.text())
            fl.write('<div class=code_sub_block>' + s + '</div>')
            
        elif name == 'image': 
            img = translateImage(getBeforeColon(line.text()))
            s = '<center><table width=0%><tr><td>' + img + '</td></tr>'
            
            caption = translateText(getAfterColon(line.text()))
            if (len(caption) >= 2):
            	s += '<tr><td lign=left class=image_sub_block_captionn>' + caption + '</td></tr>'
            s += '</table></center>'
            
            fl.write('<div class=image_sub_block>' + s + '</div>')
            
        elif name == 'note':
            s = dddoc.DATA["globals.subsections.note"].text()
            s = '<span class=section_sub_headline>' + s + '</span> '
            s += translateText(line.text())
            fl.write('<div class=note_sub_block>' + s + '</div>')
            

    subprintLink(fl, data["metafunction"], "metafunction")
    subprintLink(fl, data["type"], "type")
    subprintLink(fl, data["concept"], "concept")
    subprintText(fl, data["value"], "value")
    subprintText(fl, data["default"], "default")
    printTableContent(fl, data, "param")
    subprintText(fl, data["remarks"], "remarks")    
    subprintLink(fl, data["see"], "see")

################################################################################

def subprintLink(fl, data, subcategory):                         
    if data.empty():
        return
        
    headline = ''
    if subcategory:
        s = dddoc.DATA["globals.subsections"][subcategory].text()
        if s: headline = '<span class=section_sub_headline>' + s + '</span>'
        
    map = {}
    for line in data.lines:
        map[line.text()] = 1
    
    str = ''
    
    links = map.keys()
    links.sort()
    for link in links:
        s = translateLink(link)
        if s:
            if (str != ''): str += ', '
            str += s
            
    if (str != ''):
        fl.write('<div class=text_sub_block>' + headline + ' ' + str + '</div>');

