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
from random import randint
import statistics
from collections import defaultdict
import networkx as nx
import matplotlib
random.seed(9001)
from matplotlib import pyplot as plt

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
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as myfile:
        for line in myfile: # @EVA...
            yield next(myfile)[:-1] # Séquence
            next(myfile) # +
            next(myfile) # JJJ...


def cut_kmer(read, kmer_size):
    for i in range(0, len(read) - kmer_size +1): # remove the "\n" at the end
        yield read[i:i+kmer_size]
    #for value in read[:-kmer_size]:
    #    yield read[i]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = defaultdict(lambda: 0)
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            kmer_dict[kmer] += 1
    return kmer_dict

def build_graph(kmer_dict):
    kmer_graph = nx.DiGraph()
    for key in kmer_dict:
        node1 = key[:-1]
        node2 = key[1:]
        kmer_graph.add_node(node1)
        kmer_graph.add_node(node2)
        kmer_graph.add_edge(node1,node2,weight=kmer_dict[key])
    return kmer_graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
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
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
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
    average_weight = []
    for i in range(0, len(path)-1):
        data = graph.get_edge_data(path[i],path[i+1])
        average_weight.append(data['weight'])
    return statistics.mean(average_weight)

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph,ancestor_node, descendant_node))
    path_length = []
    path_weight = []
    for path in path_list:
        path_length.append(len(path))
        path_weight.append(path_average_weight(graph,path))
    graph = select_best_path(graph, path_list, path_length, path_weight)
    return graph

def test_solve_bubble():
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 15),
                                     (4, 5, 15), (2, 10,10), (10, 5,10),
                                     (2, 8, 3), (8, 9, 3), (9, 5, 3),
                                     (5, 6, 10), (5, 7, 10)])
    #graph_1 = solve_bubble(graph_1, 2, 5)

    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 10),
                                     (4, 5, 10), (2, 10,10), (10, 5,10),
                                     (2, 8, 10), (8, 9, 10), (9, 5, 10),
                                     (5, 6, 10), (5, 7, 10)])
    #graph_2 = solve_bubble(graph_2, 2, 5)
    simplify_bubbles(graph_1)



def simplify_bubbles(graph):
    nodes = graph.nodes()
    bubbles = []
    for node in nodes:
        predecessors = list(graph.predecessors(node))
        print("Node "+str(node))
        for predecessor in predecessors:
            print(predecessor)
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
        print("Bubble : ")
        print(str(bubble[1]) + " " + str(bubble[0]))
        graph = solve_bubble(graph,bubble[1],bubble[0])

    return graph


def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph:
        if len(graph.in_edges(node)) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes = []
    for node in graph:
        if len(graph.out_edges(node)) == 0:
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for source in starting_nodes:
        for target in ending_nodes:
            paths = nx.all_simple_paths(graph,source,target)
            for path in paths :
                final_path = path[0]
                for edge in path[1:]:
                    final_path += edge[1]
                contigs.append((final_path,len(final_path)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    myfile = open(output_file,'w')
    contig_number = 0
    for contig, len_contig in contigs_list:
        my_string = '>contig_' + str(contig_number) + ' len=' + str(len_contig) + '\n'
        myfile.write(my_string + fill(contig) + '\n')
        contig_number += 1
    myfile.close()



# def draw_graph(graph, graphimg_file):
#     """Draw the graph
#     """
#     fig, ax = plt.subplots()
#     elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
#     #print(elarge)
#     esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
#     #print(elarge)
#     # Draw the graph with networkx
#     #pos=nx.spring_layout(graph)
#     pos = nx.random_layout(graph)
#     nx.draw_networkx_nodes(graph, pos, node_size=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
#                            edge_color='b', style='dashed')
#     #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
#     # save image
#     plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    mydict = build_kmer_dict('data/eva71_hundred_reads.fq', 21)
    my_graph = build_graph(mydict)
    #draw_graph(my_graph, 'graph.png')
    source = get_starting_nodes(my_graph)
    target = get_sink_nodes(my_graph)
    contigs = get_contigs(my_graph,source,target)
    #save_contigs(contigs,'tests/test.fna')
    test_solve_bubble()

if __name__ == '__main__':
    main()
