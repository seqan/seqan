// Add "more / less" toggle link to each element with folded class.
$( document ).ready(function() {
  $( ".foldable" ).each(function() {
    $( this ).wrap( '<div class="fold-out-wrapper folded-in"></div>' );
    var toggleFunc = function() {
        var parent = $( this ).parent();
        if (parent.hasClass('folded-in')) {
            parent.removeClass('folded-in');
            parent.addClass('folded-out');
        } else {
            parent.removeClass('folded-out');
            parent.addClass('folded-in');
        }
    };
    var foldOutMore = jQuery('<a />', {class: 'fold-out-more', html: '&#9654; more...'});
    foldOutMore.click(toggleFunc);
    $( this ).parent().prepend(foldOutMore);
    var foldOutLess = jQuery('<a />', {class: 'fold-out-less', html: '&#9660; less...'});
    foldOutLess.click(toggleFunc);
    $( this ).parent().prepend(foldOutLess);
  });  
});

