window.langEntities = {
	typedef: {
		ideogram: 'typedef',
		description: '<p>A typedef allows to assign a simple name to a complex type, e.g. <code>DnaString</code> is nothing but <code>String&lt;Dna&gt;</code>.</p>',
		belongsTo: null
		},
	grouped_typedef: {
		ideogram: 'typedef',
		description: 'todo',
		belongsTo: 'typedef'
		},
	global_typedef: {
		ideogram: 'typedef',
		description: 'todo',
		belongsTo: 'typedef'
		},
	member_typedef: {
		ideogram: 'class { typedef }',
		description: 'todo',
		belongsTo: 'typedef'
		},
		
	concept: {
		ideogram: 'Concept',
		description: '<p>Since SeqAn relies on template programming the majority of the classes\' member functions are technically global. A concept can be seen as an interface which allows ...</p><p>E.g. the so called interface function <code>size(c)</code> returns the size of every object that implements the <code>ContainerConcept</code>.</p>',
		belongsTo: null
		},
	'class': {
		ideogram: 'Class',
		description: '<p>A class in OOP typically wraps member functions and member variables. In SeqAn classes are mainly used to ...</p>',
		belongsTo: null
		},
	'enum': {
		ideogram: 'typedef',
		description: '<p>An enum is an enumerated type that consists of finite number of values of the same type.</p>',
		belongsTo: null
		},
		
	metafunction: {
		ideogram: 'Fn<>',
		description: 'todo',
		belongsTo: null
		},
	interface_metafunction: {
		ideogram: '#Fn<>',
		description: 'todo',
		belongsTo: 'metafunction'
		},
		
	'function': {
		ideogram: 'fn()',
		description: 'todo',
		belongsTo: null
		},
	global_function: {
		ideogram: 'fn()',
		description: 'todo',
		belongsTo: 'function'
		},
	interface_function: {
		ideogram: '#fn()',
		description: 'todo',
		belongsTo: 'function'
		},
	member_function: {
		ideogram: '.fn()',
		description: 'todo',
		belongsTo: 'function'
		},
		
	tag: {
		ideogram: 'Tag',
		description: 'todo',
		belongsTo: null
		},
	grouped_tag: {
		ideogram: 'Tag',
		description: 'todo',
		belongsTo: 'tag'
		},

	variable: {
		ideogram: 'int x',
		description: 'todo',
		belongsTo: null
		},
	global_variable: {
		ideogram: 'int x',
		description: 'todo',
		belongsTo: 'variable'
		},
	local_variable: {
		ideogram: 'int x',
		description: 'todo',
		belongsTo: 'variable'
		},
	member_variable: {
		ideogram: '.x',
		description: 'todo',
		belongsTo: 'variable'
		},       

	adaption: {
		ideogram: 'foreign::',
		description: 'todo',
		belongsTo: null
		},
	macro: {
		ideogram: '#define',
		description: 'todo',
		belongsTo: null
		},
		
	group: {
		ideogram: 'Group',
		description: 'todo',
		belongsTo: null
		},
	page: {
		ideogram: 'Page',
		description: 'todo',
		belongsTo: null
   	}
};