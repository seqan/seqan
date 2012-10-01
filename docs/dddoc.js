var MAX_RESULT = 50;

function updateSearch(text)
{
    if (DB == null) return;
    
    s = '';
    count = 1;
        
    if (text.length >= 2)
    {
        if (text.length < 3)
            reg = new RegExp('^(' + text.toLowerCase() + ')', "gi");
        else
            reg = new RegExp('(' + text.toLowerCase() + ')', "gi");
            
        for (i = 0; i < DB.length-1; ++i)
        {
            entry = DB[i];
            key = entry[0];
            
            if (key.match(reg))
            {
                displaytext = entry[0].replace(reg, '<b>$1</b>');
                s += '<div><nobr><a target=_parent ' + entry[2] + displaytext + ' ' + entry[1] + '</a></nobr></div>';
                ++count;
                if (count >= MAX_RESULT) break;
            }
        }
    }
        
    result.innerHTML = s ;
}
