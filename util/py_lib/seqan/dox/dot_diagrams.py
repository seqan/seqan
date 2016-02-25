#!/usr/bin/env python2
"""Representation of inheritance diagrams for using with dot/graphviz."""


class IDNode(object):
    """Node in an inheritance diagram.

    @ivar name: Name of the node.
    @ivar entity: Represented proc_doc.ProcEntry.
    """

    def __init__(self, name='', entity=None):
        self.name = name
        self.entity = entity


class InheritanceDiagram(object):
    """Represents an inheritance diagram.

    @ivar nodes: Dict of IDNode objects, indexed by name.
    @ivar edges: Set of (IDNode, IDNode) pairs.
    """

    def __init__(self):
        self.nodes = {}
        self.edges = set()

    def addNode(self, node):
        """Add the given class diagram node."""
        self.nodes[node.name] = node

    def addEdge(self, a, b):
        """Add an edge a -> b."""
        self.edges.add((a.name, b.name))


class InheritanceDiagramRenderer(object):
    def render(self, diagram):
        return None


class DiagramRenderResult(object):
    """Result for rendering a dot diagram.

    @ivar image_path: String with the path to the image.
    @ivar map_path: String with the path to the <map> file.
    """

    def __init__(self):
        self.image_path = None
        self.map_path = None
