/**
 * Language Entity Labels
 */
(function ($) {

    var languageEntityLabels = {
        typedef: ['typedef', '<p>A typedef allows to assign a simple name to a complex type, e.g. <code>DnaString</code> is nothing but <code>String&lt;Dna&gt;</code>.</p>'],
        grouped_typedef: ['typedef', 'todo'],
        global_typedef: ['typedef', 'todo'],
        member_typedef: ['class { typedef }', 'todo'],

        concept: ['concept', '<p>Since SeqAn relies on template programming the majority of the classes\' member functions are technically global. A concept can be seen as an interface which allows ...</p><p>E.g. the so called interface function <code>size(c)</code> returns the size of every object that implements the <code>ContainerConcept</code>.'],
        'class': ['Class', '<p>A class in OOP typically wraps member functions and member variables. In SeqAn classes are mainly used to ...</p>'],
        'enum': ['enum', '<p>An enum is an enumerated type that consists of finite number of values of the same type.</p>'],

        meta_function: ['Fn<>', '<p>A metafunction is a function that is run at compile-time (in contrast to an ordinary function that is run at run-time).</p><p>SeqAn heavily relies on them for performance reasons.</p><p>E.g. ...'],
        interface_meta_function: ['#Fn<>', 'todo'],

        'function': ['fn()', '<p>Because SeqAn relies on template programming there are functions that used to be member functions but technically are global functions.</p><p>Such functions are called interface functions.</p>E.g. ...'],
        global_function: ['fn()', 'todo'],
        interface_function: ['#fn()', 'todo'],
        member_function: ['.fn()', 'todo'],

        grouped_tag: ['Tag', 'todo'],
        tag: ['Tag', 'todo'],

        variable: ['int x', 'todo'],
        global_variable: ['int x', 'todo'],
        local_variable: ['int x', 'todo'],
        member_variable: ['.x', 'todo'],

        adaption: ['foreign::', 'todo'],
        macro: ['#define', 'todo'],

        group: ['Group', 'todo'],
        page: ['Page', 'todo']
    }

    languageEntityLabels['metafunction'] = languageEntityLabels['meta_function'];
    
	$.fn.extend({
		pimpLanguageEntityLabels: function() {
			return this.each(function() {
				$(this).find('*[data-lang-entity]').each(function () {
					var $this = $(this);
					if($this.prop('pimped')) return true;
					$this.prop('pimped', true);
            
					var languageEntity = $this.attr('data-lang-entity');
					var oneLanguageEntity = 'a' + ($.inArray(languageEntity.substr(0, 1).toLowerCase(), ['a', 'e', 'i', 'o', 'u']) >= 0 ? 'n ' : ' ') + '<strong>' + languageEntity.replace('_', ' ') + '</strong>';

					var entry = languageEntityLabels[languageEntity];
					if (!entry) entry = ["undefined", "undefined"];

					$('<a href="language_entities.html#' + languageEntity + '">' + entry[0] + '</a>')
					    .prependTo(this)
					    .attr('title', 'What is ' + oneLanguageEntity + '?')
					    .popover({
					        html: true,
					        placement: 'right',
					        trigger: 'hover',
					        content: entry[1] + '<p>Click now for more information</p>',
					        container: 'body'
					    });
				});
			});
		}
	});

    $(document).ready(function () {
    	$('body').pimpLanguageEntityLabels();
    });

})(jQuery);




/**
 * Search Engine
 */




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
            minimumLength: 0,
            descriptiveWords: 25,
            highlightTerms: true,
            highlightEveryTerm: true,
            output: $("#results"),
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
            	$results.pimpLanguageEntityLabels();
            }
        });
    }

    $(document).ready(function () {
        createFilterableSearch($('#search'));
    });

})(jQuery);