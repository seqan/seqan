function darken(hex, amount) {
	var color = new less.tree.Color(hex.substr(1), 1.0);
	amount = new less.tree.Value(amount);
	return less.tree.functions.darken(color, amount).toCSS();
}

window.colors = {
	red: '#de5b5b',
	orange: '#debd5b',
	lgreen: '#9dde5b',
	rgreen: '#5bde7c',
	lblue: '#5bdede',
	rblue: '#5b7cde',
	purple: '#9d5bde',
	pink: '#de5bbd'
}

window.langEntities = {
	typedef: {
		name: 'Typedef',
		ideogram: 'typedef',
		color: window.colors.rgreen,
		description: 'C++ typedefs creates type aliases, e.g. alias complex template instantiations to a simple name.',
		belongsTo: null
		},
	grouped_typedef: {
		name: 'Grouped Typedef',
		ideogram: 'typedef',
		color: window.colors.rgreen,
		description: 'Typedef in a semantic group.',
		belongsTo: 'typedef'
		},
	global_typedef: {
		name: 'Global Typedef',
		ideogram: 'typedef',
		color: window.colors.rgreen,
		description: 'todo',
		belongsTo: 'typedef'
		},
	member_typedef: {
		name: 'Member Typedef',
		ideogram: 'class { typedef }',
		color: window.colors.rgreen,
		description: 'Typedef withing a class/struct.',
		belongsTo: 'typedef'
		},
		
	concept: {
		name: 'Concept',
		ideogram: 'Concept',
		color: darken(window.colors.red, 20),
		description: 'Informal <b>interface</b> for types.',
		belongsTo: null
		},
	'class': {
		name: 'Class',
		ideogram: 'Class',
		color: window.colors.red,
		description: 'C++ Data structure with both data and functions.',
		belongsTo: null
		},
	'enum': {
		name: 'Enum',
		ideogram: 'enum',
		color: darken(window.colors.red, -20),
		description: 'Custom C++ type for a predetermined number of values, each has a name.',
		belongsTo: null
		},
		
	metafunction: {
		name: 'Metafunction',
		ideogram: 'Fn<>',
		color: window.colors.rblue,
		description: 'Compile-time evaluated function that returns a type as a function of types or compile-time constants.  In C++, implemented as class templates.',
		belongsTo: null
		},
	interface_metafunction: {
		name: 'Interface Metafunction',
		ideogram: '#Fn<>',
		color: window.colors.rblue,
		description: 'A metafunction that is part of a type\'s global interface.',
		belongsTo: 'metafunction'
		},
		
	'function': {
		name: 'Function',
		ideogram: 'fn()',
		color: window.colors.lblue,
		description: 'C++ function.',
		belongsTo: null
		},
	global_function: {
		name: 'Global Function',
		ideogram: 'fn()',
		color: window.colors.lblue,
		description: 'Global C++ function',
		belongsTo: 'function'
		},
	interface_function: {
		name: 'Interface Function',
		ideogram: '#fn()',
		color: window.colors.lblue,
		description: 'Global C++ function that is part of the interface of a class.',
		belongsTo: 'function'
		},
	member_function: {
		name: 'Member Function',
		ideogram: '.fn()',
		color: window.colors.lblue,
		description: 'A class\' or struct\'s member function.',
		belongsTo: 'function'
		},
		
	tag: {
		name: 'Tag',
		ideogram: 'Tag',
		color: window.colors.purple,
		description: 'Class that is only used for its type (e.g. in <em>tag dispatching</em>).',
		belongsTo: null
		},
	grouped_tag: {
		name: 'Grouped Tag',
		ideogram: 'Tag',
		color: window.colors.purple,
		description: 'Tag that belongs to a semantic group.',
		belongsTo: 'tag'
		},

	variable: {
		name: 'Variable',
		ideogram: 'int x',
		color: window.colors.lgreen,
		description: 'Variable.',
		belongsTo: null
		},
	global_variable: {
		name: 'Global Variable',
		ideogram: 'int x',
		color: window.colors.lgreen,
		description: 'Global variable.',
		belongsTo: 'variable'
		},
    // TODO(holtgrew): remove!
	local_variable: {
		name: 'Local Variable',
		ideogram: 'int x',
		color: window.colors.lgreen,
		description: 'Local variable',
		belongsTo: 'variable'
		},
	member_variable: {
		name: 'Member Variable',
		ideogram: '.x',
		color: window.colors.lgreen,
		description: 'Member variable of a class or struct.',
		belongsTo: 'variable'
		},       

	adaption: {
		name: 'Adaption',
		ideogram: 'foreign::',
		color: darken(window.colors.orange, 10),
		description: 'Allows to use non-SeqAn types with SeqAn concept interfaces.',
		belongsTo: null
		},
	macro: {
		name: 'Macro',
		ideogram: '#define',
		color: darken(window.colors.orange, -10),
		description: 'C preprocessor macro',
		belongsTo: null
		},
		
	group: {
		name: 'Group',
		ideogram: 'Group',
		color: darken(window.colors.pink, -10),
		description: 'Set of functions and/or tags that belong together.',
		belongsTo: null
		},
	page: {
		name: 'Page',
		ideogram: 'Page',
		color: darken(window.colors.pink, 10),
		description: 'Static explanatory page.',
		belongsTo: null
   	},
   	
   	unknown: {
		name: 'Unknown Language Entity',
		ideogram: 'UNKNOWN',
		color: '#f00',
		description: 'This is an unknown language entity.',
		belongsTo: null
   	}
};