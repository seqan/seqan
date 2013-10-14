/*
JSON Search Engine 0.1
Copyright 2013 Björn Kahlert, Manuel Holtgrewe, Freie Universität Berlin
JSON Search is released under the MIT License
and based on the Tipue Search, http://www.tipue.com
*/
(function ($) {

	/**
	 * Constructs a new JSON search based on a given <form>.
	 * Options are all optional and default to the values given below.
	 */
    $.fn.jsonsearch = function (options) {
        var $form = $(this);
 
        var settings = $.extend({

            numElementsPerPage: 7,
            target: '_self',
            showURL: true,
            raw: false,
            minimumLength: 3,
            descriptiveWords: 25,
            highlightTerms: true,
            highlightEveryTerm: false,
            data: [],
            stopWords: [],
            replaceWords: [],
            stemWords: [],
            langEntityGroups: [],
            searchOnKeyPress: false,
            queryInput: $form.find('input[type=text],input[type=search]'),
            queryLangEntityInput: $form.find('select'),
            button: $form.find('input[type=submit],input[type=button]'),
            output: $form.find('.results'),
            callback: function($form, $results) {}
// TODO: Statisch den Popover-Code rausschreiben
        }, options);

        if (settings.output.length == 0) settings.output = $('<div class="results"></div>').appendTo($form);
        
        // if a stop word was found it won't appear in the result
        function getNonStopWords(words, stopWords) {
        	var nonStopWords = [];
			$.each(words, function(i, word) {
				var keep = true;
                $.each(stopWords, function(j, stopWord) {
                	if (word == stopWord) {
                        keep = false;
                        return false; // break;
                    }
                });
                
                if (keep) {
                    nonStopWords.push(word);
                }
            });
            return nonStopWords;
        }
        
        // if a word was found it won't appear in the result but its replacement
        function getReplacedWords(words, replaceWords) {
        	var replacedWords = [];
            $.each(words, function(i, word) {
            	var replaced = false;
            	$.each(replaceWords, function(j, replaceWord) {
                    if (word == replaceWord[0]) {
                    	replacedWords.push(replaceWord[1]);
                    	replaced = true;
                    	return false; // break
                    }
                });
				
				if(!replaced) {
					replacedWords.push(word);
				}
			});
			return replacedWords;
        }
        
        // adds a stemmed word to the result if the original word was found
        function getStemmedWords(words, stemWords) {
        	var stemmedWords = [];
            $.each(words, function(i, word) {
            	stemmedWords.push(word);
            	$.each(stemWords, function(j, stemWord) {
                    if (word == stemWord[0]) {
                    	stemmedWords.push(stemWord[1]);
                    }
                });
			});
			return stemmedWords;
        }

        return this.each(function () {

            var data = $.extend({}, settings.data);
            var ankerTarget = settings.target ? ' target="' + settings.target + '"' : '';

            function getURLP(name) {
                return decodeURIComponent((new RegExp('[?|&]' + name + '=' + '([^&;]+?)(&|#|;|$)').exec(location.search) || [, ""])[1].replace(/\+/g, '%20')) || null;
            }
            if (getURLP('q')) {
                settings.queryInput.val(getURLP('q'));
                search(0, true);
            }

			// make the search button invoke the search
            settings.button.click(function () {
                search(0, true);
            }).submit(function () {
				return false;
            });
            $form.submit(function () {
				return false;
            });
            settings.queryInput.keyup(function (event) {
                if (event.keyCode == '13' || settings.searchOnKeyPress) {
                    search(0, true);
                }
            });
            settings.queryLangEntityInput.change(function (event) {
                search(0, true);
            });

			/**
			 * The search algorithm itself.
			 *
			 * @param start 0-based index of where to start to search within the data. A value of 5 for example would mean that the first 6 items are skipped.
			 * @param replace if true the replaceWords are applied.
			 *
			 * @return void; the output is directly written to the output element(s)
			 */
            function search(start, replace) {
                settings.output.hide();
                var out = '';
                var results = '';
                
                var langEntities = settings.queryLangEntityInput.val();
                if(!langEntities) langEntities = [];

                var words = $.trim(settings.queryInput.val().toLowerCase()).split(' ');
                var nonStopWords = getNonStopWords(words, settings.stopWords);
                var replacedWords;
                var stemmedWords;

                if (nonStopWords.join(" ").length < settings.minimumLength) {
                	if (words.length != nonStopWords.length) {
                        out += '<div class="warning_head">Nothing found</div><div class="warning">Common words are largely ignored</div>';
                    } else {
                        out += '<div class="warning_head">Search too short</div>';
                        if (settings.minimumLength == 1) {
                            out += '<div class="warning">Should be one character or more</div>';
                        } else {
                            out += '<div class="warning">Should be ' + settings.minimumLength + ' characters or more</div>';
                        }
                    }
                } else {
                    replacedWords = replace ? getReplacedWords(nonStopWords, settings.replaceWords) : nonStopWords;
                    stemmedWords = getStemmedWords(replacedWords, settings.stemWords);
                    
                    /**
                     * The search algorithm's core itself
                     *
                     * Preconditions
                     * - cleanedWords contains an array of search terms
                     *
                     * Invariant
                     * - this always resolves to the currently iterated page
                     * - one page is described by one line in search.data.js
                     *   e.g. this.title refers to the title of the currently iterated page
                     *
                     * Postconditions
                     * - found contains an unsorted array of result elements
                     * - each result element must contain the following fields:
                     *   - score (the lower the better)
                     *   - title
                     *   - location
                     *   - langEntity (e.g. class or variable)
                     */
                    var cleanedWords = stemmedWords;
                    var found = [];
                    $.each(data, function(i) {
                    	this.langEntity = getReplacedWords([this.langEntity], settings.langEntityGroups)[0];
                    	if($.inArray(this.langEntity, langEntities) < 0) return;

                        var score = 1000000000;
                         
                        for (var j = 0; j < cleanedWords.length; j++) {
                            var pattern = new RegExp(cleanedWords[j], 'i');
                            if (this.title.search(pattern) != -1) {
                                score -= (200000 - i);
                            }
                            if (this.text.search(pattern) != -1) {
                                score -= (150000 - i);
                            }

                            if (settings.highlightTerms) {
                            	var hightlightPattern = new RegExp('(' + cleanedWords[j] + ')', 'i');
                                if (settings.highlightEveryTerm) hightlightPattern = new RegExp('(' + cleanedWords[j] + ')', 'gi');
                                text = this.text.replace(hightlightPattern, "<b>$1</b>");
                            }
                            
                            if (this.tags.search(pattern) != -1) {
                                score -= (100000 - i);
                            }

                        }
                        
                        if (score < 1000000000) {
                            found.push({ score: score, title: this.title, text: this.text, location: this.loc, langEntity: this.langEntity });
                        }
                        
                    });

                    if (found.length == 0 && !settings.raw) {
                    	out += '<div class="warning_head">Nothing found</div>';
                    } else {
                    	if(!settings.raw) {
                        	if ($.grep(replacedWords,function(x) {return $.inArray(x, words) < 0}).length > 0) {
     	    	               	// evaluates to true if the cleanedWords are not all contained in words = words have been replaced
        	                	// stop words are not relevant since they never make this condition fail
                	            out += '<div class="warning_head">Showing results for ' + words.join(' ') + '</div>';
                    	        out += '<div class="warning">Search for <a href="javascript:void(0)" class="replaced">' + cleanedWords.join(' ') + '</a></div>';
                        	}
                        
                        	if (found.length == 1) {
                          		out += '<div class="results_count">1 result</div>';
                        	} else {
                          		c_c = found.length.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
                          		out += '<div class="results_count">' + c_c + ' results</div>';
                        	}
                        }

                        found.sort(function(r1, r2) { return r1.score - r2.score; });
                        var l_o = 0;
                        out += '<ol class="results">';
                        for (var i = 0; i < found.length; i++) {
                            if (settings.numElementsPerPage < 0 || (l_o >= start && l_o < settings.numElementsPerPage + start)) {
                            	var langEntity = found[i].langEntity;
                            	var langEntityEntry = window.langEntities[langEntity];
                            	if(!langEntityEntry) langEntityEntry = { name: 'UNKNOWN', ideogram: 'UNKNOWN', color: '#FF0000', description: 'Unknown language entity' };
                            	
                            	out += '<li class="result" data-lang-entity-container="' + langEntity + '">' +
                                       '<h2>' +
                                         '<span data-lang-entity="' + langEntity + '" data-pimped="true">' +
                                            '<a href="language_entities.html#' + langEntity + '">' + langEntityEntry.ideogram + '</a>' +
                                            '<a href="' + found[i].location + '"' + ankerTarget + '>' + found[i].title + '</a>' +
                                         '</span>' +
                                       '</h2>' +
                                       '<div>';

                                var t = found[i].text;
                                var t_d = '';
                                var t_w = t.split(' ');
                                if (t_w.length < settings.descriptiveWords) {
                                    t_d = t;
                                } else {
                                    for (var f = 0; f < settings.descriptiveWords; f++) {
                                        t_d += t_w[f] + ' ';
                                    }
                                }
                                t_d = $.trim(t_d);
                                if (t_d.charAt(t_d.length - 1) != '.') {
                                    t_d += ' ...';
                                }
                                out += '<div class="text">' + t_d + '</div>';

                                if (settings.showURL) {
                                    t_url = found[i].location;
                                    if (t_url.length > 45) {
                                        t_url = found[i].location.substr(0, 45) + ' ...';
                                    }
                                    out += '<div class="location"><a href="' + found[i].location + '"' + ankerTarget + '>' + t_url + '</a></div>';
                                }
                                out += '</div>';
                                out += '</li>';
                            }
                            l_o++;
                        }
                        out += '</ol>';

                        if (settings.numElementsPerPage > 0 && found.length > settings.numElementsPerPage) {
                            var pages = Math.ceil(found.length / settings.numElementsPerPage);
                            var page = (start / settings.numElementsPerPage);
                            out += '<div class="paging"><ul>';

                            if (start > 0) {
                                out += '<li><a href="javascript:void(0)" class="' + (start - settings.numElementsPerPage) + '_' + replace + '">Prev</a></li>';
                            }

                            if (page <= 2) {
                                var p_b = pages;
                                if (pages > 3) {
                                    p_b = 3;
                                }
                                for (var f = 0; f < p_b; f++) {
                                    if (f == page) {
                                        out += '<li class="current">' + (f + 1) + '</li>';
                                    } else {
                                        out += '<li><a href="javascript:void(0)" class="' + (f * settings.numElementsPerPage) + '_' + replace + '">' + (f + 1) + '</a></li>';
                                    }
                                }
                            } else {
                                var p_b = pages + 2;
                                if (p_b > pages) {
                                    p_b = pages;
                                }
                                for (var f = page; f < p_b; f++) {
                                    if (f == page) {
                                        out += '<li class="current">' + (f + 1) + '</li>';
                                    } else {
                                        out += '<li><a href="javascript:void(0)" class="' + (f * settings.numElementsPerPage) + '_' + replace + '">' + (f + 1) + '</a></li>';
                                    }
                                }
                            }

                            if (page + 1 != pages) {
                                out += '<li><a href="javascript:void(0)" class="' + (start + settings.numElementsPerPage) + '_' + replace + '">Next</a></li>';
                            }

                            out += '</ul></div>';
                        }
                    }
                }

                settings.output.html(out);
                settings.output.slideDown(200);

                $form.find('replaced').click(function () {
                    search(0, false);
                });

                $form.find('foot_box').click(function () {
                    var id_v = $(this).attr('id');
                    var id_a = id_v.split('_');

                    search(parseInt(id_a[0]), id_a[1]);
                });
                
                if(typeof settings.callback == 'function') {
                	settings.callback($form, settings.output);
                }
            }

        });
    };

})(jQuery);