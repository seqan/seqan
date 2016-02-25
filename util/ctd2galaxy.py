#!/usr/bin/env python2

import argparse
import sys
import xml.sax


class CLIElement(object):
    """Represents a <clielement> tag.

    option_identifier -- str with parameters (e.g. --param), empty if argument.
    is_list -- bool whether the element is a list.
    """

    def __init__(self, option_identifier='', mapping_path='', is_list=False):
        self.option_identifier = option_identifier
        self.mapping = None  # Link to ParametersNode, set after parsing.
        self.mapping_path = mapping_path
        self.is_list = is_list

    def __str__(self):
        t = (self.option_identifier, self.mapping_path, self.is_list)
        return 'CLIElement(%s, %s, %s)' % tuple(map(repr, list(t)))


class ParametersNode(object):
    """Represents a <NODE> tag inside the <PARAMETERS> tags."""
    
    def __init__(self, kind='', name='', description='', value='', type_='', tags='',
                 restrictions='', supported_formats=''):
        self.name = name
        self.description = description
        self.value = value
        self.type_ = type_
        self.tags = tags
        self.supported_formats = supported_formats
        self.restrictions = None
        self.path = None  # root if is None
        self.parent = None  # not set, usually a list
        self.children = {}

    def computePath(self):
        """Compute path entry from parent links."""

    def __str__(self):
        t = (self.name, self.description, self.value, self.type_, self.tags,
             self.supported_formats, self.children)
        return 'ParametersNode(%s, %s, %s, %s, %s, %s, %s)' % tuple(map(repr, t))

    def __repr__(self):
        return str(self)


class Tool(object):
    """Represents the top-level <tool> tag from a CTD file."""

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

    def __str__(self):
        t = (self.name, self.executable_name, self.version, self.description,
             self.manual, self.doc_url, self.category)
        return 'Tool(%s, %s, %s, %s, %s, %s, %s)' % tuple(map(repr, list(t)))

        

class CTDFormatException(Exception):
    """Raised when there is a format error in CTD."""

    
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
        if self.stack == ['tool']:
            # Create the top level Tool object.
            self.tool = Tool()
            self.result = self.tool
        elif self.stack == ['tool', 'cli', 'clielement']:
            # Create a new CLIElement object for a <clieelement> tag.
            if not attrs.get('isList'):
                raise CTDFormatException('No attribute isList in <clielement>.')
            if attrs.get('optionIdentifier') is None:
                raise CTDFormatException('no attribute optionIdentifier in <clielement>.')
            is_list = (attrs.get('isList') == 'false')
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
            self.parameter_node = node.parent
        elif self.stack[:2] == ['tool', 'PARAMETERS'] and self.stack[-1] == 'ITEM':
            # Create a new item ParametersNode for the <ITEM> entry.
            if not attrs.get('name'):
                raise CTDFormatException('no attribute name in <ITEM>')
            name = attrs.get('name')
            value = attrs.get('value')
            type_ = attrs.get('type')
            tags = attrs.get('tags')
            description = attrs.get('description')
            restrictions = attrs.get('restrictions')
            supported_formats = attrs.get('supported_formats')
            child = ParametersNode(
                kind='item', name=name, description=description, value=value,
                type_=type_, tags=tags, supported_formats=supported_formats)
            self.parameter_node.children[name] = child

    def endElement(self, name):
        self.stack.pop()
        if name == 'NODE':
            self.parameter_node = self.parameter_node.parent

    def characters(self, content):
        if self.stack == ['tool', 'name']:
            self.tool.name += content
        elif self.stack == ['tool', 'executableName']:
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
        parser = xml.sax.make_parser()
        parser.setContentHandler(self.handler)
        parser.parse(path)
        return self.handler.result


def main():
    parser = argparse.ArgumentParser(description='Convert CTD to Galaxy XML')
    parser.add_argument('-i', '--in-file', metavar='FILE',
                        help='CTD file to read.', dest='in_file',
                        required=True)

    args = parser.parse_args()

    ctd_parser = CTDParser()
    tool = ctd_parser.parse(args.in_file)
    print tool
    for cli in tool.cli_elements:
        print '  %s' % cli
    print tool.parameters
        

if __name__ == '__main__':
    sys.exit(main())