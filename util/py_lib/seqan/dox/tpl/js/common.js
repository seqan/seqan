/**
 * Various
 */
(function ($) {

	$(document).ready(function () {
	
	    // only keep top-level nav items
	    // exceptions: examples and mainpage
	    $('html:not(.page_mainpage) #toc ol ol').filter(function() { return $(this).find('a[href=#Examples], a[href=#Example]').length == 0; }).remove();
		$('html.page_languageentities #toc ol ol').remove();

	    // hightlight nav items based on scroll area
	    $('body').scrollspy({ target: '#toc', offset: 50 });
	    
	    var id;
	    $(window).resize(function() {
            clearTimeout(id);
            id = setTimeout(function() { $('body').scrollspy('refresh'); }, 500);
        });
        
        // shows 'open in frameset' link if opened separately
		if(window == window.parent && window.name != 'list') {
        	$('#content').prepend('<div class="open-in-frame alert alert-info"><a href="index.html#' + window.location.href.replace(/^.*\/|\.[^.]*$/g, '') + '"><strong>Looking for a different entry?</strong> Unhide the navigation bar and start your search.</a></div>'); 
        }
        
        // if loaded in a frame, checks the addresses hash and uses it to load the
        // specified page in the right frame (e.g. docs.seqan.de/index.html#abc will open abc.html in the main frame) 
        if(window != window.parent && window.name == 'list') {
        	try {
        		var redirectTo = null;
        		if($.urlParam('p', window.parent.location)) {
        			redirectTo = $.urlParam('p', window.parent.location) + '.html';
        		} else {
        			var hash = window.parent.location.hash;
        			if(typeof hash === 'string' && hash.length > 1) {
        				redirectTo = hash.substr(1) + '.html';
        			}
        		}
        		
        		if(redirectTo) {
        			window.parent['main'].location = redirectTo;
        		}
    		} catch(e) {
    		    // some browsers like Chrome don't allow this cross-frame access if using file://
    		}
        }
 
        // tooltips
        $('[title]:not([href])').tooltip({ container: 'body' });

        // smooth scrolling
        $('a[href*=#]:not([href=#])').smoothScroll({ offset: -20 });
        
        // autofocus search field
        if($('html').hasClass('list')) {
        	window.setTimeout(function() {
        		$('input[type=search]').focus();
        	}, 50);
        }

    });

})(jQuery);

/**
 * Get URL parameter functionality
 */
(function ($) {
	$.extend({
		urlParam: function(name, location) {
			if(!location) location = window.location;
            return decodeURIComponent((new RegExp('[?|&]' + name + '=' + '([^&;]+?)(&|#|;|$)').exec(location.search) || [, ""])[1].replace(/\+/g, '%20')) || null;
		}
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
			
			if(langEntity == 'tutorial') return $('<span>' + entry.ideogram + '</span>');
			return $('<a href="page_LanguageEntities.html#' + langEntity + '">' + entry.ideogram + '</a>');
		},
		
		pimpLangEntityLabel: function(langEntity) {
    		$this = $(this);
    		if(jQuery.inArray($this.prop('tagName'), ['A']) != -1) {
    		    // a tags may and should not be nested
			    $this.wrap('<span data-lang-entity="' + langEntity + '" data-pimped="true"/>')
					.removeAttr('data-lang-entity')
					.before($().createLangEntityLabel(langEntity));
			} else {
    			$this.wrapInner('<span/>').prepend($().createLangEntityLabel(langEntity));
			}
		},
		
		/**
		 * Annotates all tags with a data-lang-entity(-container) attribute the following way:
		 * 1) if the tag has the data-lang-entity attribute,
		 *      it will be prefixed with a colored lang entity label
		 * 2) if the tag has the data-lang-entity-container attribute,
         *      the highest headline level will be prefixed,
         *      whereas all other headlines won't be prefixed (by removing their eventually set data-lang-entity attribute)
		 */
		pimpLangEntityLabels: function() {
			return this.each(function() {
				$(this).find('[data-lang-entity]').each(function () {
					var $this = $(this);
					if($this.attr('data-pimped')) return true;
					//if($this.parents('[data-lang-entity-container]').length > 0) return true;
            
					var langEntity = $this.attr('data-lang-entity');
					$this.pimpLangEntityLabel(langEntity);				
				});
			}).each(function() {
			    $(this).find('[data-lang-entity-container]').each(function () {
    			    // only use one headline level
    			    var headlineHandled = false;
    			    for(var i=1; i<=6; i++) {
        			    $el = $(this).find('h' + i);
        			    if(!headlineHandled) {
            			    if($el.length > 0) {
            					if($el.attr('data-pimped')) return true;
            					
            					var langEntity = $el.parents('[data-lang-entity-container]').attr('data-lang-entity-container');
            					$el.attr('data-lang-entity', langEntity);
            					$el.pimpLangEntityLabel(langEntity);
            					
                			    headlineHandled = true;
            			    }
        			    } else {
            			    //$el.removeAttr('data-lang-entity');
        			    }
    			    }
    			 });
			});
		}
	});

    $(document).ready(function () {
        $('body').on('mouseover', '[data-lang-entity] > :first-child', function() {
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
					return '<div class="description">' + langEntityData.description + '</div>' + (showMore && langEntity != 'tutorial' ? '<p class="more">Click now for more information...</p>' : '');
				},
				container: 'body',
				placement: function(tip, element) {
                    var $thisement, above, actualHeight, actualWidth, below, boundBottom, boundLeft, boundRight, boundTop, elementAbove, elementBelow, elementLeft, elementRight, isWithinBounds, left, pos, right;
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
            includeSelectAllOption: true,
            selectAllText: "(Un)check all",
            selectAllValue: 'all',
            dropRight: false,
            buttonText: function (checkedOptions) {
                var options = $(arguments[1][0]).find('option').map(function () {
                    return $(this).val();
                });
                $.each(options, function(i){
                    if(options[i] === 'all') {
                        options.splice(i,1);
                        return false;
                    }
                });
                $.each(checkedOptions, function(i){
                    if(checkedOptions[i] === 'all') {
                        checkedOptions.splice(i,1);
                        return false;
                    }
                });

				var $btn = $el.find('button.multiselect');
				if(checkedOptions.length == options.length) {
					$btn.removeClass('btn-warning');
					$btn.addClass('btn-primary');
				} else {
					$btn.addClass('btn-warning');
					$btn.addClass('btn-primary');
				}

                if (options.length == checkedOptions.length) return '<i class="fa fa-filter"></i><small>all visible</small>';
                else if (options.length - checkedOptions.length == 1) return '<i class="fa fa-filter"></i><small>1 excluded</small>';
                else if (checkedOptions.length == 0) return '<i class="fa fa-filter"></i><small>all excluded</small>';
                else return '<i class="fa fa-filter"></i><small>' + (options.length - checkedOptions.length) + ' excluded</small>';
            },
            onChange: function ($element, checked) {
                //$allLabel = $('.multiselect-container').find('li:first-child label');
                //$allLabel.contents().last()[0].textContent = $allLabel.text()[0] == 'U' ? 'Check all' : 'Uncheck all';
            }
        });
        $el.find('.multiselect-container [value=all]').prop('checked', true);
        
        $el.find('button.multiselect').attr('title', '');

        // copies the options value to <li>'s data-lang-entity-container attribute classes of the parent li element (for easier styling)
        $el.find('.multiselect-container input[value]').each(function () {
            $this = $(this);
            if($this.val() == 'all') return;
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
            output: $("#results"),
            data: window.searchData,
            stopWords: ["and", "be", "by", "do", "for", "he", "how", "if", "is", "it", "my", "not", "of", "or", "the", "to", "up", "what", "when"], // filtered out of query
            replaceWords: [ // words replaced in the query
                ["test_A", "test_B"]
            ],
            stemWords: [ // silently adds the stem if the corresponding word was found in the query
                ["javascript", "script"]
            ],
            langEntityGroups: [ // TODO create from window[langEntities] (try console.log(window[langEntities])
                ["grouped_typedef", "typedef"],
                ["global_typedef", "typedef"],
                ["member_typedef", "typedef"],
                ["grouped_tag", "tag"],
                ["global_variable", "variable"],
                ["local_variable", "variable"],
                ["member_variable", "variable"]
            ],
            callback: function($form, $results) {
            	if($form.find('input[type=search]').val().length == 0) {
            		$("html").removeClass('shows-results');
            	    $("#results").fadeOut();
            	    $el.find('.pre-action').slideDown();
                } else {
                	$("html").addClass('shows-results');
                    $("#results").fadeIn();
                    $el.find('.pre-action').slideUp();
                }
            }
        });
    }

    $(document).ready(function () {
        createFilterableSearch($('#search'));
        $('#search').fadeIn();
        
        try {
            // search immediately if query was passed within the url
    		var q = decodeURI((RegExp('q=' + '(.+?)(&|$)').exec(parent.location.search)||[,null])[1]);
    		if(q && q != 'null') {
    			$('#search [type=search]').val(q).change().focus();
    		}
        } catch(e) {
            // some browsers like Chrome don't allow this cross-frame access if using file://
        }
    });
    
    // hide form and results until they are pimped
	$('head').append('<style>#search, #results { display: none; }</style>');

})(jQuery);