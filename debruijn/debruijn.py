#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
from operator import itemgetter
import random
import statistics
from collections import defaultdict
import networkx as nx
random.seed(9001)

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=5, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Generator that reads the given file.
      :Parameters:
          fastq_file: File to read
    """
    with open(fastq_file) as myfile:
        for line in myfile: # @EVA...
            yield next(myfile)[:-1] # Séquence
            next(myfile) # +
            next(myfile) # JJJ...


def cut_kmer(read, kmer_size):
    """Generator that extracts kmers from a sequence.
      :Parameters:
          read: Generator that gives sequences
          kmer_size: Size of a kmer
    """
    for i in range(0, len(read) - kmer_size +1): # remove the "\n" at the end
        yield read[i:i+kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    """Extract from a file a dict with a kmer as key and the number of occurences as value.
      :Parameters:
          fastq_file: File to read
          kmer_size: Size of a kmer
    """
    kmer_dict = defaultdict(lambda: 0)
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            kmer_dict[kmer] += 1
    return kmer_dict

def build_graph(kmer_dict):
    """Build a graph from the given dict;
      :Parameters:
          kmer_dict: Dict of kmer.
    """
    kmer_graph = nx.DiGraph()
    for key in kmer_dict:
        node1 = key[:-1]
        node2 = key[1:]
        kmer_graph.add_node(node1)
        kmer_graph.add_node(node2)
        kmer_graph.add_edge(node1,node2,weight=kmer_dict[key])
    return kmer_graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove the given paths from a given graph.
      :Parameters:
          graph: The graph from which we remove paths.
          path_list: The list of paths we want to remove.
          delete_entry_node: Boolean that tell if the entry nodes must be removed.
          delete_sink_node: Boolean that tell if the sink nodes must be removed.
    """
    for path in path_list:
        for node_number in range(0,len(path)-1) :
            graph.remove_edge(path[node_number], path[node_number+1])
            if node_number == 0:
                if delete_entry_node:
                    graph.remove_node(path[node_number])
            else:
                graph.remove_node(path[node_number])
        if delete_sink_node:
            graph.remove_node(path[len(path)-1])
    return graph

def std(data):
    """Return the standard deviation of list of values.
      :Parameters:
          data: List of values.
    """
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Return the graph without the unwanted paths.
      :Parameters:
          graph: The graph we want to clean.
          path_list: List of all the possible paths of the graph.
          path_length: List of the length of the paths.
          weight_avg_list: List of the average weight of the paths.
          delete_entry_node: Boolean that tell if the entry nodes must be removed.
          delete_sink_node: Boolean that tell if the sink nodes must be removed.
    """
    # determine if a path is more frequent than others
    max_weight = max(weight_avg_list)
    max_weight_indexes = [i for i, x in enumerate(weight_avg_list) if x == max_weight]
    if len(max_weight_indexes) == 1:
        path_to_keep = max_weight_indexes[0]
    else:
        # determine if a path is longer than others
        max_length = max(path_length)
        max_length_indexes = [i for i, x in enumerate(path_length) if x == max_length]
        if len(max_length_indexes) == 1:
            path_to_keep = max_length_indexes[0]
        else:
            # choose a path randomly
            path_to_keep = random.randint(0,len(path_list)-1)
    path_to_keep = path_list[path_to_keep]
    # make the list of paths to remove
    paths_to_remove = []
    for path in path_list:
        if path != path_to_keep:
            paths_to_remove.append(path)
    graph = remove_paths(graph,paths_to_remove,delete_entry_node,delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """Return the average weight of a path.
      :Parameters:
          graph: The graph which contains the path.
          path : The path we want to evaluate.
    """
    average_weight = []
    for i in range(0, len(path)-1):
        data = graph.get_edge_data(path[i],path[i+1])
        average_weight.append(data['weight'])
    return statistics.mean(average_weight)

def solve_bubble(graph, ancestor_node, descendant_node):
    """Solve a bubble from given ancestor and descendant nodes in a given graph.
      :Parameters:
          graph: The graph which contains the paths.
          ancestor_node: The origin node from the paths.
          descendant_node: The ending node from the paths.
    """
    path_list = list(nx.all_simple_paths(graph,ancestor_node, descendant_node))
    path_length = []
    path_weight = []
    for path in path_list:
        path_length.append(len(path))
        path_weight.append(path_average_weight(graph,path))
    graph = select_best_path(graph, path_list, path_length, path_weight)
    return graph

def simplify_bubbles(graph):
    """Remove all bubbles from a given graph.
      :Parameters:
          graph: The graph from which we want to remove bubbles.
    """
    nodes = graph.nodes()
    bubbles = []
    for node in nodes:
        predecessors = list(graph.predecessors(node))
        for predecessor in predecessors:
            if len(predecessors) > 1:
                p_0 = predecessors[0]
                p_1 = predecessors[1]
                common_ancestor =  nx.algorithms.lowest_common_ancestor(graph,p_0,p_1)
                for i in range(1,len(predecessors)-1):
                    pred = predecessors[i]
                    common_ancestor = nx.algorithms.lowest_common_ancestor(graph, common_ancestor, pred, default=common_ancestor)
                if common_ancestor != None:
                    #graph = solve_bubble(graph,common_ancestor,node)
                    bubbles.append((node,common_ancestor))
    for bubble in bubbles:
        graph = solve_bubble(graph,bubble[1],bubble[0])
    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove all entry tips from the given graph.
      :Parameters:
          graph: The graph from which we want to remove tips.
          starting_nodes: The entry nodes from the paths.
    """
    reverse_graph = nx.reverse(graph)
    solved_graph = solve_out_tips(reverse_graph,starting_nodes)
    graph = nx.reverse(solved_graph)
    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove all ending tips from the given graph.
      :Parameters:
          graph: The graph from which we want to remove tips.
          ending_nodes: The ending nodes from the paths.
    """
    if len(ending_nodes) <= 1:
        return graph
    common_ancestor = nx.algorithms.lowest_common_ancestor(graph, ending_nodes[0], ending_nodes[1])
    for i in range(2,len(ending_nodes)):
        pred = ending_nodes[i]
        if common_ancestor != None :
            common_ancestor = nx.algorithms.lowest_common_ancestor(graph, common_ancestor, pred, default=common_ancestor)
    if common_ancestor == None:
        return graph
    paths = []
    path_length = []
    path_weight = []
    for ending_node in ending_nodes:
        simple_paths = list(nx.all_simple_paths(graph,common_ancestor,ending_node))
        for simple_path in simple_paths:
            paths.append(simple_path)
            path_length.append(len(simple_path))
            path_weight.append(path_average_weight(graph,simple_path))
    graph = select_best_path(graph, paths, path_length, path_weight,delete_sink_node=True)
    return graph

def get_starting_nodes(graph):
    """Get all the starting nodes of a graph.
      :Parameters:
          graph: The graph from which we want the starting nodes.
    """
    starting_nodes = []
    for node in graph:
        if len(graph.in_edges(node)) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    """Get all the sink nodes of a graph.
      :Parameters:
          graph: The graph from which we want the sink nodes.
    """
    sink_nodes = []
    for node in graph:
        if len(graph.out_edges(node)) == 0:
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Get all the possible paths of the given graph, with the given starting and ending nodes.
      :Parameters:
          graph: The graph from which we want the paths.
          starting_nodes: Starting nodes of the graph.
          ending_nodes: Endind nodes of the graph.
    """
    contigs = []
    for source in starting_nodes:
        for target in ending_nodes:
            paths = nx.all_simple_paths(graph,source,target)
            for path in paths :
                final_path = path[0]
                for edge in path[1:]:
                    final_path += edge[-1]
                contigs.append((final_path,len(final_path)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    """Save the contigs in a file.
      :Parameters:
          contigs_list: List of contigs.
          output_file: File in which we save the contigs.
    """
    myfile = open(output_file,'w')
    contig_number = 0
    for contig, len_contig in contigs_list:
        my_string = '>contig_' + str(contig_number) + ' len=' + str(len_contig) + '\n'
        myfile.write(my_string + fill(contig) + '\n')
        contig_number += 1
    myfile.close()

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments and create graph
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    # Remove bubbles
    graph = simplify_bubbles(graph)
    # Remove entry and out tips
    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    sink_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph, sink_nodes)
    # Save contigs
    source = get_starting_nodes(graph)
    target = get_sink_nodes(graph)
    contigs = get_contigs(graph,source,target)
    save_contigs(contigs,args.output_file)

if __name__ == '__main__':
    main()
