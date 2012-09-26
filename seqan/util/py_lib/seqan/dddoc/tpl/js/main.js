function toggleSource( id )
{
  var $src = $('#' + id).toggle();
  $('#l_' + id).html($src.css('display') == 'none'  ? 'show' : 'hide');
}

function openCode( url )
{
  window.open( url, "SOURCE_CODE", "resizable=yes,scrollbars=yes,toolbar=no,status=no,height=480,width=750" ).focus();
}


window.highlight = function(url) {
    var hash = url.match(/#([^#]+)$/)
    if(hash) {
        $('a[name=' + hash[1] + ']').parent().effect('highlight', {}, 'slow')
    }
}

$(function() {
    highlight('#' + location.hash);
});

function setPermalink(path)
{
    var loc = window.parent.location.origin + window.parent.location.pathname;
    $('#permalink-field')[0].value = loc + '?i=' + path;
}

function togglePermalink()
{
    var linkPos = $('.permalink').offset();
    var linkWidth = $('.permalink').width();
    var width = $('#permalink-div').width();
    var height = $('#permalink-div').height();
    $('#permalink-div').toggle();
    $('#permalink-div').css({'left': (linkPos.left - width + linkWidth + 10) + 'px', 'top': (linkPos.top - height - 15) + 'px'})
    if ($('#permalink-div').is(':visible'))
    {
        $('#permalink-field')[0].focus();
        $('#permalink-field')[0].select();
    }
}