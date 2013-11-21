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
            showURL: false,
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
            queryInput: $form.find('input[type=text],input[type=search]'),
            queryLangEntityInput: $form.find('select'),
            langEntityDefaultOrder: [
            	'concept', 'class', 'enum',
				'typedef', 'grouped_typedef', 'global_typedef', 'member_typedef',
				'metafunction', 'global_metafunction', 'interface_metafunction',
				'tag', 'grouped_tag', 'generic',
				'function', 'global_function', 'interface_function', 'member_function',
				'variable', 'global_variable', 'member_variable',
				'adaption', 'macro',
				'group', 'page',
				'unknown'],
			maxResultsPerGroup: 5,
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
        
        // function findAndScore(query, entry).
        //
        // Return -1 if nothing was found.
        function findAndScore(query, entry) {
            var BASE_SCORE = 1000000000;
            
            var RANK_FACTOR = 10000;
            var DAMPENING_TEXT = 0.99;  // dampening per text sequence
            var DAMPENING_QUERY = 0.99;  // dampening per query sequence
        
            var resultScore = -1;
        
            for (var j = 0; j < query.length; j++) {
                var pattern = new RegExp(query[j], 'i');
                for (var k = 0; k < entry.length; k++) {
                    var result = entry[k].match(pattern);
                    if (result === null || result.length == 0)
                        continue;  // No match.
                    for (var x = 0; x < result.length && x < 1; x++) {  // x < 1 reduces repeat effect
                        // Compute begin/end position by searching for the result.
                        var beginPos = entry[k].search(result[x]);
                        var endPos = beginPos + result[x].length;
        
                        // Rank points: (4) full match, (3) prefix match, (2) suffix match, (1) infix
                        // match.
                        var rankPoints = 1;
                        rankPoints += 2 * (beginPos == 0);
                        rankPoints += (endPos == entry[k].length);
        
                        // Ratio of matched text.
                        var ratio = 1.0 * (endPos - beginPos) / entry[k].length;
        
                        // Compute final score.
                        var scoreDelta = rankPoints * RANK_FACTOR * ratio;
                        for (var y = 0; y < j; ++y)
                            scoreDelta *= DAMPENING_QUERY;
                        for (var y = 0; y < k; ++y)
                            scoreDelta *= DAMPENING_TEXT;
                        var newScore = BASE_SCORE - scoreDelta;
                        if (resultScore == -1 || newScore < resultScore)
                            resultScore = newScore;
                    }
                    break;  // Found a match, done.
                }
        
                if (resultScore != -1)
                    break;  // Found a match, done.
            }
            
            return resultScore;
        }
        
        // Create match with highlighting.
        function highlightedMatch(score, obj, query) {
            var result = {score: score, title:obj.title, text:obj.text, location:obj.loc,
                          subentries:[], akas:null, langEntity:obj.langEntity,
                          hiTitle: false, hiAka: false, hiSubentry: false};
        
            function highlightString(str, pattern) {
                if (str.match(pattern))
                    return str.replace(pattern, "<b>$1</b>");
                else
                    return str;
            }
        
            var akas = obj.akas.split(',');
            var subentries = obj.subentries.split(',');
        
            for (var j = 0; j < query.length; j++) {
                var pattern = new RegExp('(' + query[j] + ')', 'i');
        
                if (!result.hiTitle && result.title.search(pattern) != -1)
                {
                    result.title = highlightString(result.title, pattern);
                    result.hiTitle = true;
                }
        
                if (!result.hiTitle && !result.hiAka)
                    for (var i = 0; i < akas.length; i++)
                        if (akas[i].search(pattern) != -1)
                        {
                            result.hiAka = true;
                            result.aka = highlightString(akas[i], pattern);
                            break;
                        }
        
                if (!result.hiTitle && !result.hiAka && result.subentries.length < 4)
                    for (var i = 0; i < subentries.length; ++i)
                    {
                        var xs = subentries[i].split(' ', 2);
                        var kind = xs[0];
                        var title = xs[1];
                        if (subentries[i].search(pattern) != -1)
                        {
                            result.hiSubentry = true;
                            result.subentries.push([kind, highlightString(title, pattern)]);
                        }
                    }
            }
        
            return result;
        }

        return this.each(function () {

            var data = $.extend({}, settings.data);
            var ankerTarget = settings.target ? ' target="' + settings.target + '"' : '';

            if ($.urlParam('q')) {
                settings.queryInput.val($.urlParam('q'));
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
            settings.queryInput.change(function() {
                search(0, true);
            });
            settings.queryInput.keyup(function(event) {
                if (event.keyCode == '13') {
                    try {
                        $firstResult = settings.output.find('.result a:nth-child(2)').first();
                        var target = $firstResult.attr('target') || 'main';
                        var targetFrame = window.parent.frames[target] || parent;
                        targetFrame.location.href = $firstResult.attr('href');
                    } catch(e) {
                        // security exceptions if using file://
                    }
                } else {
                    search(0, true);
                }
            });
            settings.queryLangEntityInput.change(function (event) {
                search(0, true);
            });
            $('#results').on('click', 'li.more:not(.result) a', function() {
                var i=0;
                var $more = $(this).parents('[data-lang-entity-container]').find('.more');
                $more.each(function() {
            		$this = $(this);
            		if($this.hasClass('result')) $this.slideDown();
            		else $this.slideUp();
            	});
            	return false;
            });

			/**
			 * The search algorithm itself.
			 *
			 * @param start 0-based index of where to start to search within the data. A value of 5 for example would mean that the first 6 items are skipped.
			 * @param replace if true the replaceWords are applied.
			 *
			 * @return void; the output is directly written to the output element(s)
			 */
            var lastQuery = false;
            var lastLangEntities = false;
            function search(start, replace) {
                var out = '';
                var results = '';
                
                var langEntities = settings.queryLangEntityInput.val();
                if(!langEntities) langEntities = [];
                
                if(settings.queryInput.val() != '') {
                	settings.queryInput.addClass('not-empty');
                } else {
                	settings.queryInput.removeClass('not-empty');
                }

                var words = $.trim(settings.queryInput.val().toLowerCase()).split(' ');
                var nonStopWords = getNonStopWords(words, settings.stopWords);
                var replacedWords;
                var stemmedWords;

                if (nonStopWords.join(" ").length < settings.minimumLength) {
                    lastQuery = false;
                    lastLangEntities = false;
                    if(!settings.raw) {
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
                    var query = cleanedWords;
                    if(query.join(' ') == lastQuery && langEntities.join(' ') == lastLangEntities) { return; }
                    lastQuery = query.join(' ');
                    lastLangEntities = langEntities.join(' ');
                    settings.output.hide();
                    
                    var found = [];
                    $.each(data, function(i) {
                    	this.langEntity = getReplacedWords([this.langEntity], settings.langEntityGroups)[0];
                    	if($.inArray(this.langEntity, langEntities) < 0) return;
                    	
                        var score = findAndScore(query, [this.title, this.text, this.akas, this.subentries]);
                        var result = highlightedMatch(score, this, query);

                        if (score != -1 && score < 1000000000) {
                            found.push(result);
                        }
                    });

                    if (found.length == 0) {
                    	if(settings.raw) out += '<ol class="results empty"><li>Nothing found.</li></ol>';
                    	else out += '<div class="warning_head">Nothing found</div>';
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
                        
                        found.sort(function(r1, r2) {
                        	var o1 = $.inArray(r1.langEntity, settings.langEntityDefaultOrder);
                        	var o2 = $.inArray(r2.langEntity, settings.langEntityDefaultOrder);
                        	
                        	if(o1 < 0) o1 = 9999;
                        	if(o2 < 0) o2 = 9999;
                        	
                        	if(o1 == o2) return r1.score - r2.score;
                        	return o1-o2;
                        });
                        
                        var l_o = 0;
                        var entriesInGroup;
                        var lastLangEntity = false;
                        out += '<ol class="results">';
                        for (var i = 0; i < found.length; i++) {
                            if (settings.numElementsPerPage < 0 || (l_o >= start && l_o < settings.numElementsPerPage + start)) {
                            	var langEntity = found[i].langEntity;
                            	var langEntityEntry = window.langEntities[langEntity];
                            	if(!langEntityEntry) langEntityEntry = { name: 'UNKNOWN', ideogram: 'UNKNOWN', color: '#FF0000', description: 'Unknown language entity' };
                            	
                            	// groups entries by their lang entity
                            	if(lastLangEntity != langEntity) {
                            		if(lastLangEntity) out += '</ol></li>';
                            		entriesInGroup = 1;
                            		out += '<li data-lang-entity-container="' + langEntity + '" data-pimped="true"><span data-lang-entity="' + langEntity + '"><a href="page_LanguageEntities.html#' + langEntity + '">' + langEntityEntry.ideogram + '</a><span>' + langEntityEntry.name + '</span></span><ol class="nav">';
                            		lastLangEntity = langEntity;
                            	} else {
                            		entriesInGroup++;
                            	}
                            	
                            	var isLastEntryInGroup = found.length > i+1 && found[i+1].langEntity != langEntity;
                            	var entryIsOneTooMuch = entriesInGroup == settings.maxResultsPerGroup+1;

                                // adds more link if there is more than one item to be displayed of the same group
                                if(entryIsOneTooMuch && !isLastEntryInGroup){
									out += '<li class="more"><a href="#">...</a></li>';
                                }
                                
                                out += '<li class="result' + (entriesInGroup > settings.maxResultsPerGroup && !(entryIsOneTooMuch && isLastEntryInGroup) ? ' more' : '') + '">\
                                        <h2>\
                                            <span data-lang-entity="' + langEntity + '" data-pimped="true">\
                                                <a href="page_LanguageEntities.html#' + langEntity + '">' + langEntityEntry.ideogram + '</a>\
                                                <a href="' + found[i].location + '"' + ankerTarget + '>' + found[i].title + '<div>';
                                        
                                if (found[i].aka) {
                                    out += '<div class="aka">' + found[i].aka + '</div>';
                                }
                                
                                if (found[i].subentries.length > 0) {
                                    out += '<ul class="subentries">';
                                    for (var j = 0; j < found[i].subentries.length; j++) {
                                        out += '<li>' + found[i].subentries[j][1] + '</li>';
                                    }
                                    out += '</ul>';
                                }
                                        
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
									out += '<div class="location">' + t_url + '</div>';
								}
								
								out += '        </div></a>\
                                            </span>\
                                        </h2>\
                                        <div>';
								out += '</li>';
                            }
                            l_o++;
                        }
                        
                        if(lastLangEntity) out += '</ol></li>';
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
