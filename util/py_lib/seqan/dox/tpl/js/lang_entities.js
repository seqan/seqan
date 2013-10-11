window.langEntities = {
	typedef: {
		ideogram: 'typedef',
		description: 'C++ typedefs creates type aliases, e.g. alias complex template instantiations to a simple name.',
		belongsTo: null
		},
	grouped_typedef: {
		ideogram: 'Typedef in a semantic group.',
		description: 'todo',
		belongsTo: 'typedef'
		},
	global_typedef: {
		ideogram: 'Global typedef.',
		description: 'todo',
		belongsTo: 'typedef'
		},
	member_typedef: {
		ideogram: 'class { typedef }',
		description: 'Typedef withing a class/struct.',
		belongsTo: 'typedef'
		},
		
	concept: {
		ideogram: 'Concept',
		description: 'Informal <b>interface</b> for types.',
		belongsTo: null
		},
	'class': {
		ideogram: 'Class',
		description: 'C++ Data structure with both data and functions.',
		belongsTo: null
		},
	'enum': {
		ideogram: 'enum',
		description: 'Custom C++ type for a predetermined number of values, each has a name.',
		belongsTo: null
		},
		
	metafunction: {
		ideogram: 'Fn<>',
		description: 'Compile-time evaluated function that returns a type as a function of types or compile-time constants.  In C++, implemented as class templates.',
		belongsTo: null
		},
	interface_metafunction: {
		ideogram: '#Fn<>',
		description: 'A metafunction that is part of a type\'s global interface.',
		belongsTo: 'metafunction'
		},
		
	'function': {
		ideogram: 'fn()',
		description: 'C++ function.',
		belongsTo: null
		},
	global_function: {
		ideogram: 'fn()',
		description: 'Global C++ function',
		belongsTo: 'function'
		},
	interface_function: {
		ideogram: '#fn()',
		description: 'Global C++ function that is part of the interface of a class.',
		belongsTo: 'function'
		},
	member_function: {
		ideogram: '.fn()',
		description: 'A class\' or struct\'s member function.',
		belongsTo: 'function'
		},
		
	tag: {
		ideogram: 'Tag',
		description: 'Class that is only used for its type (e.g. in <em>tag dispatching<em>).',
		belongsTo: null
		},
	grouped_tag: {
		ideogram: 'Tag',
		description: 'Tag that belongs to a semantic group.',
		belongsTo: 'tag'
		},

	variable: {
		ideogram: 'int x',
		description: 'Variable.',
		belongsTo: null
		},
	global_variable: {
		ideogram: 'int x',
		description: 'Global variable.',
		belongsTo: 'variable'
		},
    // TODO(holtgrew): remove!
	local_variable: {
		ideogram: 'int x',
		description: 'Local variable',
		belongsTo: 'variable'
		},
	member_variable: {
		ideogram: '.x',
		description: 'Member variable of a class or struct.',
		belongsTo: 'variable'
		},       

	adaption: {
		ideogram: 'foreign::',
		description: 'Allows to use non-SeqAn types with SeqAn concept interfaces.',
		belongsTo: null
		},
	macro: {
		ideogram: '#define',
		description: 'C preprocessor macro',
		belongsTo: null
		},
		
	group: {
		ideogram: 'Group',
		description: 'Set of functions and/or tags that belong together.',
		belongsTo: null
		},
	page: {
		ideogram: 'Page',
		description: 'Static explanatory page.',
		belongsTo: null
   	}
};