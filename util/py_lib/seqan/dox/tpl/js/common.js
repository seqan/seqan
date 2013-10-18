/**
 * Various
 */
(function ($) {

	$(document).ready(function () {
	
	    // only keep top-level nav items
	    $('#toc ol ol').filter(function() { return $(this).find('a[href=#Examples]').length == 0; }).remove();

	    // hightlight nav items based on scroll area
	    $('body').scrollspy({ target: '#toc', offset: 50 });
	    
	    var id;
	    $(window).resize(function() {
            clearTimeout(id);
            id = setTimeout(function() { $('body').scrollspy('refresh'); }, 500);
        });
        
        // tooltips
        $('[title]:not([href])').tooltip({ container: 'body' });
        console.log($('[title]:not([href])'));

        // smooth scrolling
        $('a[href*=#]:not([href=#])').click(function() {
            if (location.pathname.replace(/^\//,'') == this.pathname.replace(/^\//,'') || location.hostname == this.hostname) {
                var target = $(this.hash);
                target = target.length ? target : $('[name=' + this.hash.slice(1) +']');
                if (target.length) {
                    $('html,body').animate({ scrollTop: target.offset().top-20 }, 350);
                    return false;
                }
            }
        });

    });

})(jQuery);

/**
 * Language Entity Labels
 */
(function ($) {

	$.fn.extend({
		createLangEntityLabel: function(langEntity) {
            if(!window.langEntities) return;
            var entry = window.langEntities[langEntity];
			if (!entry) entry = window.langEntities['unknown'];
					
			return $('<a href="language_entities.html#' + langEntity + '">' + entry.ideogram + '</a>');
		},
		
		pimpLangEntityLabels: function() {
			return this.each(function() {
				$(this).find('[data-lang-entity]').each(function () {
					var $this = $(this);
					if($this.attr('data-pimped')) return true;
            
					var langEntity = $this.attr('data-lang-entity');
					
					// if dealing with list items
					// wrap the inner nodes to not destroy the li
					if($this.prop('tagName') == 'LI') {
						$this.wrapInner('<span data-lang-entity="' + langEntity + '"/>')
							.removeAttr('data-lang-entity');
						$this = $this.children();
					}
					
					if($this.hasClass("signature")) {
    					$this.prepend($().createLangEntityLabel(langEntity));
					} else {					    
    					$this.wrap('<span data-lang-entity="' + langEntity + '" data-pimped="true"/>')
    						.removeAttr('data-lang-entity')
    						.before($().createLangEntityLabel(langEntity));
                    }
				});
			});
		}
	});

    $(document).ready(function () {
        $('body').on('mouseover', '[data-lang-entity] > a:first-child', function() {
            var langEntity = $(this).attr('data-lang-entity') || $(this).parent().attr('data-lang-entity');
            showPopOver(this, langEntity, true);
        });
        
        $('#search').on('mouseover', '[data-lang-entity-container] label', function() {
            var langEntity = $(this).parents('[data-lang-entity-container]').attr('data-lang-entity-container');
            showPopOver(this, langEntity, false);
        });
        
        function showPopOver(el, langEntity, showMore) {
    		var $this = $(el);
    		if($this.attr('data-original-title')) return; // already set up

    		var langEntityData = window.langEntities[langEntity];
    		
    		$this.popover({
				html: true,
				trigger: 'hover',
				template: '<div class="popover" data-lang-entity-container="' + langEntity + '"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>',
				title: langEntityData.name,
				content: function() {
					return '<div class="description">' + langEntityData.description + '</div>' + (showMore ? '<p class="more">Click now for more information</p>' : '');
				},
				container: 'body',
				placement: function(tip, element) {
                    var $element, above, actualHeight, actualWidth, below, boundBottom, boundLeft, boundRight, boundTop, elementAbove, elementBelow, elementLeft, elementRight, isWithinBounds, left, pos, right;
                    isWithinBounds = function(elementPosition) {
                    return boundTop < elementPosition.top && boundLeft < elementPosition.left && boundRight > (elementPosition.left + actualWidth) && boundBottom > (elementPosition.top + actualHeight);
                    };
                    $element = $(element);
                    pos = $.extend({}, $element.offset(), {
                        width: element.offsetWidth,
                        height: element.offsetHeight
                    });
                    actualWidth = 283;
                    actualHeight = 117;
                    boundTop = $(document).scrollTop() + 40; // DIRTY: takes the small data-lang-entity window "hat" into account
                    boundLeft = $(document).scrollLeft();
                    boundRight = boundLeft + $(window).width();
                    boundBottom = boundTop + $(window).height();
                    
                    elementAbove = {
                        top: pos.top - actualHeight,
                        left: pos.left + pos.width / 2 - actualWidth / 2
                    };
                    elementBelow = {
                        top: pos.top + pos.height,
                        left: pos.left + pos.width / 2 - actualWidth / 2
                    };
                    elementLeft = {
                        top: pos.top + pos.height / 2 - actualHeight / 2,
                        left: pos.left - actualWidth
                    };
                    elementRight = {
                        top: pos.top + pos.height / 2 - actualHeight / 2,
                        left: pos.left + pos.width
                    };
                    above = isWithinBounds(elementAbove);
                    below = isWithinBounds(elementBelow);
                    left = isWithinBounds(elementLeft);
                    right = isWithinBounds(elementRight);
                    if (above) {
                        return "top";
                    } else {
                        if (below) {
                            return "bottom";
                        } else {
                            if (left) {
                                return "left";
                            } else {
                                if (right) {
                                    return "right";
                                } else {
                                    return "bottom";
                                }
                            }
                        }
                    }
                }
			}).popover('show');
    	}
    	
    	if(!$('html').hasClass('list')) {
    		$('html').pimpLangEntityLabels();
    	}
    });

})(jQuery);





/**
 * Search Bar
 */
(function ($) {

    function createFilterableSearch($el) {
        // make filter box look fancy
        $el.find('select').multiselect({
            buttonClass: 'btn btn-primary',
            includeSelectAllOption: false,
            selectAllText: "Select all",
            selectAllValue: 'all',
            dropRight: false,
            buttonText: function (checkedOptions) {
                var options = $(arguments[1][0]).find('option').map(function () {
                    return $(this).val();
                });

                //return (options.length - checkedOptions.length) + "/" + options.length + " excluded";

                if (options.length == checkedOptions.length) return "all visible";
                else if (options.length - checkedOptions.length == 1) return "1 excluded";
                else if (checkedOptions.length == 0) return "all excluded!";
                else return (options.length - checkedOptions.length) + " excluded";
            },
            onChange: function (element, checked) {}
        });

        // copies the options value to <li>'s data-lang-entity-container attribute classes of the parent li element (for easier styling)
        $el.find('.multiselect-container input[value]').each(function () {
            $this = $(this);
            $this.parents('li').attr('data-lang-entity-container', $this.val());
            $this.parents('a').click(function (e) {
                if (e.target == this) {
                    // link and not the label or the input was clicked
                    $(this).find('input').click();
                }
            });
        });

        if($el.jsonsearch) $el.jsonsearch({
            numElementsPerPage: -1,
            target: 'main',
            raw: true,
            showUrl: false,
            minimumLength: 1,
            descriptiveWords: 25,
            highlightTerms: true,
            highlightEveryTerm: true,
            output: $("#results").css('display', 'none'),
            searchOnKeyPress: true,
            data: window.searchData,
            stopWords: ["and", "be", "by", "do", "for", "he", "how", "if", "is", "it", "my", "not", "of", "or", "the", "to", "up", "what", "when"], // filtered out of query
            replaceWords: [ // words replaced in the query
                ["test_A", "test_B"]
            ],
            stemWords: [ // silently adds the stem if the corresponding word was found in the query
                ["javascript", "script"]
            ],
            langEntityGroups: [
                ["grouped_typedef", "typedef"],
                ["global_typedef", "typedef"],
                ["member_typedef", "typedef"],
                ["interface_meta_function", "meta_function"],
                ["global_function", "function"],
                ["interface_function", "function"],
                ["member_function", "function"],
                ["grouped_tag", "tag"],
                ["global_variable", "variable"],
                ["local_variable", "variable"],
                ["member_variable", "variable"]
            ],
            callback: function($form, $results) {
            	if($form.find('input[type=search]').val().length == 0) {
            	    $("#results").slideUp();
            	    $el.find('.pre-action').slideDown();
                } else {
                    $("#results").show();
                    $el.find('.pre-action').slideUp();
                }
            }
        });
    }

    $(document).ready(function () {
        createFilterableSearch($('#search'));
    });

})(jQuery);