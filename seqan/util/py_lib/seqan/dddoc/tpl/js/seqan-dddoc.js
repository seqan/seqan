// Namespace seqan.
var seqan = {};
// Namespace seqan.doc.
seqan.doc = {};

/*
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
    
    if not ret or ret[0] == '_':
        return ret
    else:
        return '.' + ret*/

seqan.doc.FILENAME_KEEP = ['-']

seqan.doc.FILENAME_SUBS = {
    "_": "__",
    "\t": "_09",
    "\n": "_0a",
    "!": "_21",
    "\"": "_22",
    "#": "_23",
    "$": "_24",
    "%": "_25",
    "&": "_26",
    "'": "_27",
    "(": "_28",
    ",": "_29",
    "*": "_2a",
    "+": "_2b",
    "/": "_2f",
    ":": "_3a",
    ":": "_2c",
    "<": "_3c",
    ">": "_3e",
    "?": "_3f",
    "\\": "_5c",
    "|": "_7c",
    " ": "+"
}

seqan.doc.escapeFilename = function(str)
{
    ret = [];
    for (var i = 0; i < str.length; ++i)
    {
        if ((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= '0' && str[i] <= '9') || seqan.doc.FILENAME_KEEP.indexOf(str[i]) != -1)
        {
            ret.push(str[i]);
        }
        else if (str[i] >= 'A' && str[i] <= 'Z')
        {
            ret.push('_');
            ret.push(str[i]);
        }
        else if (seqan.doc.FILENAME_SUBS[str[i]] !== undefined)
        {
            ret.push(seqan.doc.FILENAME_SUBS[str[i]]);
        }
    }

    ret = ret.join('');

    if (ret.length == 0 || ret[0] == '_')
        return ret;
    else
        return '.' + ret;
}

seqan.doc.getPagePath = function(cat, subcat, prefix)
{
    result = cat.toUpperCase() + seqan.doc.escapeFilename(subcat) + ".html"
    if (prefix)
        result = prefix + '/' + result;
    return result
}

// hack in startsWith()
if (typeof String.prototype.startsWith != 'function') {
  String.prototype.startsWith = function (str){
    return this.indexOf(str) == 0;
  };
}