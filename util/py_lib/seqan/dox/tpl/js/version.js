/* Jongkyu Kim (j.kim@fu-berlin.de), 2016.01.12 */
DOCUMENT_URL = "http://docs.seqan.de/seqan/";
VERSION_JSON_URL = "http://docs.seqan.de/version.php";

function changeVersion(formid)
{
    var frm = document.getElementById(formid);
    window.top.location.href = DOCUMENT_URL + frm.options[frm.selectedIndex].value;
}    

function addVersionSelection(arr)
{
    // add HTMLs 
    var version_select = document.createElement("select");
    var version_div = document.createElement("div");

    version_select.setAttribute("id","version_select");
    version_div.setAttribute("align","center");
    version_div.appendChild(document.createElement("br"));
    version_div.appendChild(document.createTextNode("Version: "));
    version_div.appendChild(version_select);
    document.body.appendChild(version_div);
            
    version_select.addEventListener("change", function(){changeVersion(this.id);}, false);

    // current selection is.. 
    cur_sel = window.location.pathname.split("/")[2];
    cur_sel = "master"; // should be deleted

    for(i=0; i < arr.length; ++i)
    {
        var op = document.createElement("option");
        op.value = arr[i];
        op.text = arr[i];
        op.selected = ( arr[i] == cur_sel ) ? true : false;
        version_select.add(op);
    }
}

// get JSON data & add selection form
var req = new XMLHttpRequest();
req.open("GET", VERSION_JSON_URL, true);
req.setRequestHeader("Content-type", "application/json");
req.onreadystatechange = function()
{
    if( req.readyState == 4 && req.status == 200 )
    {
       var response = JSON.parse(req.responseText);
        addVersionSelection(response); // add selection form
    }
}
req.send();
