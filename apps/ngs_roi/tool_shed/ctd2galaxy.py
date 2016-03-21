#!/usr/bin/env python2
"""Conversion of the CTD format into Galaxy XML.

The CTD parser should be reusable but is not in its own module since it is
only used here at the moment.
"""

# TODO(holtgrew): Option lists do not work at the moment.

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

try:
    import argparse
except ImportError:
    import argparse26 as argparse
import operator
import sys
import xml.sax
import xml.sax.saxutils


# The suffix to identify file extension options (for '--arg-${NO}${SUFFIX}'
# and '--${PARAM_NAME}${SUFFIX}').
FILE_EXT_SUFFIX = '-file-ext'


class CTDFormatException(Exception):
    """Raised when there is a format error in CTD."""
    

class CLIElement(object):
    """Represents a <clielement> tag.

    :ivar option_identifier: with parameters (e.g. --param), empty if argument.
    :type option_identifier: str
    :ivar is_list: whether the element is a list.
    :type is_list: bool
    :ivar param_node: link to ParametersNode, set after parsing, None if unset
    :ivar is_list: w or not this element is a list.
    :type is_list: bool
    """

    def __init__(self, option_identifier='', mapping_path='', is_list=False):
        """Initialize object."""
        self.option_identifier = option_identifier
        self.param_node = None  # Link to ParametersNode, set after parsing.
        self.mapping_path = mapping_path
        self.is_list = is_list

    def __str__(self):
        """String representaiton of CLIElement."""
        t = (self.option_identifier, self.mapping_path, self.is_list)
        return 'CLIElement(%s, %s, %s)' % tuple(map(repr, list(t)))


class ParametersNode(object):
    """Represents a <NODE> tag inside the <PARAMETERS> tags.

    :ivar name: name attribute of the node
    :ivar description: text for description attribute of the node
    :ivar value: value attribute of the node
    :ivar type_: type attribute of the node
    :ivar tags: tags attribute of the node
    :ivar supported_formats: supported_format attribute of the node
    :ivar restrictions: restrictions attribute of the node
    :ivar path: the path to the node
    :ivar path: list of strings
    :ivar parent: link to the parent of the node
    :ivar children: children of the node
    :type children: dict with name to node mapping
    :ivar cli_element: CLIElement that this parameter is mapped to.
    :ivar required: Whether or not this parameter is required.
    :type required: bool
    """
    
    def __init__(self, kind='', name='', description='', value='', type_='', tags='',
                 restrictions='', supported_formats='', required=False):
        """Initialize the object."""
        self.kind = kind
        self.name = name
        self.description = description
        self.value = value
        self.type_ = type_
        self.tags = tags
        self.supported_formats = supported_formats
        self.restrictions = restrictions
        self.path = None  # root if is None
        self.parent = None  # not set, usually a list
        self.children = {}
        self.cli_element = None
        self.required = required

    def computePath(self, is_root=True, path=[]):
        """Compute path entry from parent links.

        :param is_root: whether or not this is the root node
        :type is_root: bool
        :param path: path to this node, excluding root
        :type path: list of strings
        """
        self.path = list(path)
        if not is_root:
            self.path.append(self.name)
        if not self.children:
            return  # nothing to do: early exit.
        for name, child in self.children.items():
            child.computePath(False, self.path)

    def applyFunc(self, f):
        """Apply f to self and all children."""
        f(self)
        for c in self.children.values():
            c.applyFunc(f)
            
    def find(self, path):
        """Return ParametersNode object at the path below the node."""
        if not path:
            return self
        if not self.children.get(path[0]):
            return None
        return self.children[path[0]].find(path[1:])

    def __str__(self):
        """Return string representation."""
        t = (self.name, self.description, self.value, self.type_, self.tags,
             self.supported_formats, self.children, self.path, self.required)
        return 'ParametersNode(%s, %s, %s, %s, %s, %s, %s, path=%s, %s)' % tuple(map(repr, t))

    def __repr__(self):
        """Return programmatic representation, same as __str__()."""
        return str(self)


class Tool(object):
    """Represents the top-level <tool> tag from a CTD file.

    :ivar name: name attribute value
    :type name: str
    :ivar executable_name: executableName attribute value
    :type executable_name: str
    :ivar version: version attribute value
    :type version: str
    :ivar description: description attribute value
    :type description: str
    :ivar manual: manual attribute value
    :type manual: str
    :ivar doc_url: docurl attribute value
    :type doc_url: str
    :ivar category: category attribute value
    :type category: str
    :ivar cli_elements: list of CLIElement objects
    :ivar parameters: root parameters node
    :type parameters: ParametersNode
    """

    def __init__(self, name='', executable_name='', version='',
                 description='', manual='', doc_url='',
                 category=''):
        self.name = name
        self.executable_name = executable_name
        self.version = version
        self.description = description
        self.manual = manual
        self.doc_url = doc_url
        self.category = category
        self.cli_elements = []
        self.parameters = None

    def parsingDone(self):
        """Called after parsing is done.

        The method will compute the paths of the parameter nodes and link the
        CLIElement objects in self.cli_elements to the ParameterNode objects.
        """
        self.parameters.computePath()
        for ce in self.cli_elements:
            if not ce.option_identifier:
                continue  # Skip arguments
            path = ce.mapping_path.split('.')
            node = self.parameters.find(path)
            if not node:
                raise CTDFormatException('Unknown parameter %s' % '.'.join(path))
            ce.param_node = node
            node.cli_element = ce

    def __str__(self):
        t = (self.name, self.executable_name, self.version, self.description,
             self.manual, self.doc_url, self.category)
        return 'Tool(%s, %s, %s, %s, %s, %s, %s)' % tuple(map(repr, list(t)))

        

class CTDHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        self.result = None
        # A stack of tag names that are currently open.
        self.stack = []
        # The current parameter to append nodes below.
        self.parameter_node = None

    def startElement(self, name, attrs):
        """Handle start of element."""
        # Maintain a stack of open tags.
        self.stack.append(name)
        # Handle the individual cases.  The innermost tag is self.stack[-1].
        if self.stack == ['tool']:
            # Create the top level Tool object.
            self.tool = Tool()
            self.result = self.tool
            if not attrs.get('name'):
                raise CTDFormatException('No attribute "name" in <tool> tag.')
            self.tool.name = attrs.get('name')
        elif self.stack == ['tool', 'cli', 'clielement']:
            # Create a new CLIElement object for a <clieelement> tag.
            if not attrs.get('isList'):
                raise CTDFormatException('No attribute isList in <clielement>.')
            if attrs.get('optionIdentifier') is None:
                raise CTDFormatException('no attribute optionIdentifier in <clielement>.')
            is_list = (attrs.get('isList') == 'true')
            option_identifier = attrs.get('optionIdentifier')
            self.tool.cli_elements.append(CLIElement(option_identifier=option_identifier, is_list=is_list))
        elif self.stack == ['tool', 'cli', 'clielement', 'mapping']:
            # Handle a <mapping> sub entry of a <clieelement> tag.
            if not attrs.get('referenceName'):
                raise CTDFormatException('no attribute referenceName in <mapping>')
            self.tool.cli_elements[-1].mapping_path = attrs['referenceName']
        elif self.stack == ['tool', 'PARAMETERS']:
            # Handle the <PARAMETERS> entry by creating a new top parameters node.
            self.tool.parameters = ParametersNode(kind='node', name='<root>')
            self.parameter_node = self.tool.parameters
        elif self.stack[:2] == ['tool', 'PARAMETERS'] and self.stack[-1] == 'NODE':
            # Create a new node ParametersNode for the <PARAMETERS> entry.
            if not attrs.get('name'):
                raise CTDFormatException('no attribute name in <NODE>')
            name = attrs.get('name')
            node = ParametersNode(kind='node', name=name)
            node.parent = self.parameter_node
            self.parameter_node.children[name] = node
            self.parameter_node = node
        elif self.stack[:2] == ['tool', 'PARAMETERS'] and self.stack[-1] in ['ITEM', 'ITEMLIST']:
            # Create a new item ParametersNode for the <ITEM>/<ITEMLIST> entry.
            if not attrs.get('name'):
                raise CTDFormatException('no attribute name in <ITEM>/<ITEMLIST>')
            name = attrs.get('name')
            value = attrs.get('value')
            type_ = attrs.get('type')
            tags = attrs.get('tags')
            description = attrs.get('description')
            restrictions = attrs.get('restrictions')
            required = attrs.get('required') == 'true'
            supported_formats = attrs.get('supported_formats', '')
            kind = {'ITEM': 'item', 'ITEMLIST': 'itemlist'}[self.stack[-1]]
            child = ParametersNode(
                kind=kind, name=name, description=description, value=value,
                type_=type_, tags=tags, supported_formats=supported_formats,
                restrictions=restrictions)
            self.parameter_node.children[name] = child

    def endElement(self, name):
        """Handle closing tag."""
        # Maintain stack.
        self.stack.pop()
        # Go up one node in the parameters tree if </NODE>
        if name == 'NODE':
            self.parameter_node = self.parameter_node.parent

    def characters(self, content):
        """Handle characters in XML file."""
        if self.stack == ['tool', 'executableName']:
            self.tool.executable_name += content
        elif self.stack == ['tool', 'version']:
            self.tool.version += content
        elif self.stack == ['tool', 'description']:
            self.tool.description += content
        elif self.stack == ['tool', 'manual']:
            self.tool.manual += content
        elif self.stack == ['tool', 'docurl']:
            self.tool.doc_url += content
        elif self.stack == ['tool', 'category']:
            self.tool.category += content


class CTDParser(object):
    """Parser for CTD files."""

    def __init__(self):
        self.handler = CTDHandler()

    def parse(self, path):
        # Parse XML into Tool object.
        parser = xml.sax.make_parser()
        parser.setContentHandler(self.handler)
        parser.parse(path)
        # Compute paths for tool's parameters.
        self.handler.result.parsingDone()
        return self.handler.result


class XMLWriter(object):
    """Base class for XML writers.


    :ivar result: list of strings that are joined for the final XML
    :ivar indent_level: int with the indentation level
    """

    def __init__(self):
        self.result = []
        self.indent_level = 0

    def indent(self):
        """Return indentation whitespace."""
        return '    ' * self.indent_level

    def appendTag(self, tag, text='', args={}):
        """Append a tag to self.result with text content only or no content at all."""
        e = xml.sax.saxutils.quoteattr
        args_str = ' '.join('%s=%s' % (key, e(str(value))) for key, value in args.items() if value is not None)
        if args_str:
            args_str = ' '+ args_str
        vals = {'indent': self.indent(),
                'tag': tag,
                'text': text.strip(),
                'args': args_str}
        if text:
            self.result.append('%(indent)s<%(tag)s%(args)s>%(text)s</%(tag)s>\n' % vals)
        else:
            self.result.append('%(indent)s<%(tag)s%(args)s />\n' % vals)

    def openTag(self, tag, args={}):
        """Append an opening tag to self.result."""
        e = xml.sax.saxutils.quoteattr
        args_str = ' '.join('%s=%s' % (key, e(str(value))) for key, value in args.items())
        if args_str:
            args_str = ' ' + args_str
        vals = {'indent': self.indent(),
                'tag': tag,
                'args': args_str}
        self.result.append('%(indent)s<%(tag)s%(args)s>\n' % vals)

    def closeTag(self, tag):
        """Append a closing tag to self.result."""
        vals = {'indent': self.indent(), 'tag': tag}
        self.result.append('%(indent)s</%(tag)s>\n' % vals)

    def handleParameters(self, node):
        """Recursion for appending tags for ParametersNode."""
        for pn in node.children.values():
            if pn.kind in ['item', 'itemlist']:
                args = {'name': pn.name,
                        'value': pn.value,
                        'type': pn.type_,
                        'description': pn.description,
                        'restrictions': pn.restrictions,
                        'tags': pn.tags}
                self.appendTag(pn.kind.upper(), args=args)
            else:  # node.kind == 'node'
                args = {'name': pn.name,
                        'description': pn.description}
                self.openTag('NODE', args=args)
                self.indent_level += 1
                self.handleParameters(pn)
                self.indent_level -= 1
                self.closeTag('NODE')


class CTDWriter(XMLWriter):
    """Write a Tool to CTD format."""

    def run(self, tool, f):
        """Write the given Tool to file f."""
        self.result.append('<?xml version="1.0" encoding="UTF-8"?>\n')
        self.openTag('tool')
        self.indent_level += 1
        self.appendTag('name', tool.name)
        self.appendTag('executableName', tool.executable_name)
        self.appendTag('version', tool.version)
        self.appendTag('description', tool.description)
        self.appendTag('manual', tool.manual)
        self.appendTag('docurl', tool.doc_url)
        self.appendTag('category', tool.category)
        # <cli> and <clielement> group
        self.openTag('cli')
        self.indent_level += 1
        for ce in tool.cli_elements:
            self.openTag('clielement', args={'optionIdentifier': ce.option_identifier,
                                             'isList': {True: 'true', False: 'false'}[ce.is_list]})
            self.indent_level += 1
            self.appendTag('mapping', args={'referenceName': ce.mapping_path})
            self.indent_level -= 1
            self.closeTag('clielement')
        self.indent_level -= 1
        self.closeTag('cli')
        # <PARAMETERS>, <NODE>, <ITEM>, <ITEMLIST> group
        self.openTag('PARAMETERS', args={'version': 1.4,
                                         'xsi:noNamespaceSchemaLocation': 'http://open-ms.sourceforge.net/schemas/Param_1_4.xsd',
                                         'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance'})
        self.indent_level += 1
        self.handleParameters(tool.parameters)
        self.indent_level -= 1
        self.closeTag('PARAMETERS')
        self.indent_level -= 1
        self.closeTag('tool')
        # Write result
        for x in self.result:
            f.write(x)


class GalaxyCommandSnippet(object):
    """Stores a snippet for the Galaxy <command> tag.

    Such a snippet consists of a list of text that will be concatenated using
    space as the separator.

    Optionally, a condition can be given which will be pasted verbatimly into
    an #if condition that will also be properly closed.  As a bonus, the
    snippet will be properly indented.
    """

    def __init__(self, words, condition=None, indent=4, level=1):
        self.words = words
        self.condition = condition
        self.indent = indent
        self.level = level

    def build(self):
        res = []
        if self.condition:
            res.append('#if %s' % self.condition)
        res.append(' '.join(self.words))
        if self.condition:
            res[-1] = ' ' * self.indent + res[-1]
            res.append('#end if')
        return '\n'.join([' ' * self.indent * self.level + l for l in res])


class GalaxyWriter(XMLWriter):
    """Write a Tool to the Galaxy format."""

    def run(self, tool, f):
        """Write the given Tool to file f."""
        self.result.append('<?xml version="1.0" encoding="UTF-8"?>\n')
        self.openTag('tool', {'id': tool.executable_name, 'name': tool.name})
        self.indent_level += 1
        self.addCommandTag(tool)
        self.appendTag('description', text=tool.description)
        self.openTag('inputs')
        self.indent_level += 1
        tool.parameters.applyFunc(lambda x: self.addInputParam(x))
        self.indent_level -= 1
        self.closeTag('inputs')
        self.openTag('outputs')
        self.indent_level += 1
        tool.parameters.applyFunc(lambda x: self.addOutputParam(x))
        self.indent_level -= 1
        self.closeTag('outputs')
        self.openTag('stdio')
        self.indent_level += 1
        self.appendTag('exit_code', args={'range': '1:', 'level': 'fatal'})
        self.appendTag('exit_code', args={'range': ':-1', 'level': 'fatal'})
        self.indent_level -= 1
        self.closeTag('stdio')
        self.indent_level -= 1
        self.closeTag('tool')
        # Write result
        for x in self.result:
            f.write(x)

    def addInputParam(self, param_node):
        """Add a ParametersNode object if it is to go to <inputs>."""
        if param_node.type_ == 'output-file':
            return  # Skip output files
        if param_node.kind not in ['item', 'itemlist']:
            return  # Skip if not item.
        if param_node.name.endswith('-file-ext'):
            return  # Skip if extension to override.
        args = {}
        if not param_node.required:
            args['optional'] = 'true'  # false would be default
        if param_node.type_ == 'input-file':
            args['type'] = 'data'
            args['format'] = ','.join([x.replace('*', '').replace('.', '')
                                       for x in param_node.supported_formats.split(',')])
            args['name'] = '_'.join(param_node.path).replace('-', '_').replace('.', '_')
            args['label'] = param_node.description
            args['type'] = 'data'
            self.appendTag('param', args=args)
        else:
            TYPE_MAP = {
                'string': 'text',
                'double': 'float',
                'int': 'integer'
            }
            args['type'] = TYPE_MAP[param_node.type_]
            args['name'] = '_'.join(param_node.path).replace('-', '_').replace('.', '_')
            args['label'] = param_node.description
            if param_node.type_ == 'string' and param_node.restrictions and \
               sorted(param_node.restrictions.split(',')) == ['false', 'true']:
                args['type'] = 'boolean'
                if param_node.value == 'true':
                    args['checked'] = 'true'
                args['truevalue'] = param_node.cli_element.option_identifier
                args['falsevalue'] = ''
                self.appendTag('param', args=args)
                return
            args['value'] = param_node.value
            if param_node.type_ == 'string' and param_node.restrictions:
                args['type'] = 'select'
                self.openTag('param', args=args)
                self.indent_level += 1
                for v in param_node.restrictions.split(','):
                    self.appendTag('option', v, {'value': v})
                self.indent_level -= 1
                self.closeTag('param')
            else:
                self.appendTag('param', args=args)
            
    def addOutputParam(self, param_node):
        """Add a ParametersNode object if it is to go to <inputs>."""
        if param_node.type_ != 'output-file':
            return  # Only add for output files.
        if param_node.name.endswith('-file-ext'):
            return  # Skip if extension to override.
        args = {}
        if '.' in param_node.supported_formats:
            args['format'] = param_node.supported_formats.split(',')[0].split('.')[-1]
        else:
            args['format'] = param_node.supported_formats.split(',')[0].split('*')[-1]
        args['name'] = '_'.join(param_node.path).replace('-', '_').replace('.', '_')
        args['label'] = param_node.description
        self.appendTag('data', args=args)
            
    def addCommandTag(self, tool):
        """Write <command> tag to self.result."""
        file_ext_elements = []
        # Process non-file-extension arguments.
        snippets = []
        for ce in tool.cli_elements:
            if ce.param_node.name.endswith('-file-ext'):
                file_ext_elements.append(ce)
                continue  # Skip -file-ext options.
            # The name of the variable that is used.
            var_name = '$' + ce.mapping_path.replace('-', '_').replace('.', '_')
            # Check whether it is optional.
            optional = not ce.param_node.required
            # Check whether it is a boolean.
            bool_param = False
            if ce.param_node.type_ == 'string' and ce.param_node.restrictions and \
               sorted(ce.param_node.restrictions.split(',')) == ['false', 'true']:
                bool_param = True
            # Get variable name.
            val = '"' + var_name + '"'
            # Build the snippet for the command.
            if bool_param:
                # The true value for boolean parameters is the argument itself.
                snippets.append(GalaxyCommandSnippet([var_name]))
            else:
                condition = {True: var_name, False: None}.get(optional)
                snippets.append(GalaxyCommandSnippet([ce.option_identifier, val],
                                                     condition=condition))
        # Process file extension arguments.
        ext_overrides = []
        for ce in file_ext_elements:
            if ce.option_identifier == '--write-ctd-file-ext':
                continue  # Skip special args.
            # The name of the variable that is used.
            var_name = ce.mapping_path[:-len('-file-ext')].replace('-', '_').replace('.', '_')
            snippets.append(GalaxyCommandSnippet([ce.option_identifier, '${%s.ext}' % var_name],
                                                 condition='$%s' % var_name))
        # Finalize building the command line.
        txt = GalaxyCommandSnippet([tool.executable_name]).build() + '\n' + '\n'.join([s.build() for s in snippets])
        self.appendTag('command', text=txt)


def main():
    """Main function."""
    # Setup argument parser.
    parser = argparse.ArgumentParser(description='Convert CTD to Galaxy XML')
    parser.add_argument('-i', '--in-file', metavar='FILE',
                        help='CTD file to read.', dest='in_file',
                        required=True)
    parser.add_argument('-o', '--out-file', metavar='FILE',
                        help='File to write. Output type depends on extension.',
                        dest='out_file', required=True)

    args = parser.parse_args()

    # Parse input.
    sys.stderr.write('Parsing %s...\n' % args.in_file)
    ctd_parser = CTDParser()
    tool = ctd_parser.parse(args.in_file)

    # Write output.
    sys.stderr.write('Writing to %s...\n' % args.out_file)
    if args.out_file.endswith('.ctd'):
        writer = CTDWriter()
    else:
        writer = GalaxyWriter()
    with open(args.out_file, 'wb') as f:
        writer.run(tool, f)

    return 0
        

if __name__ == '__main__':
    sys.exit(main())
