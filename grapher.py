#!/usr/local/bin/python
"""
Create and work with POAGraphs

>>> p = POAGraph()
>>> v1 = p.add_vertex('A')
>>> v2 = p.add_vertex('C')
>>> e = p.add_edge(v1, v2)
>>> p.draw('output.png')

"""

import graph_tool as gt

from   graph_tool.all import sfdp_layout, planar_layout, graph_draw


class POAGraph(object):
    """
    Represents a POA graph
    """
    
    def __init__(self):

        # The Primary graph (Have a smaller name for readability)
        self.graph = gt.Graph(directed=True)
        self._g    = self.graph 

        # The vertex names
        self.v_names = self._g.new_vertex_property("string")

        # Sequences in a vertex
        self.v_seqs = self._g.new_vertex_property("vector<string>")

        # The edge names
        self.e_names = self._g.new_edge_property("string")

        # The edge weights
        # (Not used for now as currently separate edge for each connection
        #  Possibly improve to combine edges into one)
        self.e_weights = self._g.new_edge_property("double")


    def add_vertex(self, node_name):
        """
        Add a vertex to the graph
        """
        vertex = self._g.add_vertex()
        self.v_names[vertex] = node_name
        return vertex


    def add_edge(self, v1, v2):
        """
        Add an edge between vertices v1 and v2
        """
        edge = self._g.add_edge(v1, v2)
        self.e_names = str(v1) + "-" + str(v2)
        return edge


    def vertex(self, v_name):
        """
        returns a vertex corresponding to a given name
        """
        return self._g.vertex(v_name)


    def get_vertex_name(self, vertex):
        """
        returns a name corresponding to a given vertex
        """
        return self.v_names[vertex]


    def get_all_vertices(self):
        """
        returns all vertices in graph
        """
        return self._g.vertices()


    def get_all_edges(self):
        """
        returns all edges in graph
        """
        return self._g.edges()


    def draw(self, out_filename):
        """
        Draws the POA graph to given out_file (out_filename is a filename which
        must contain the .png extension)
        """
        # sfdp_layout seems to work for now. Try more later
        layout = sfdp_layout(self._g)

        graph_draw(self._g,
                   pos=layout,
                   vertex_text=self.v_names,
                   vertex_font_size=12,
                   output=out_filename)

