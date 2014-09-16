import sys
import networkx as nx
import matplotlib.pyplot as plt
import itertools as itl

if (len(sys.argv) != 2):
    print "Syntax: %s <cluster-file>" % (sys.argv[0])
    exit(0)

cls_file = sys.argv[1]

def read_cluster(file_name):
    cluster = []
    with open(file_name, "r") as f:
        for line in f:
            cluster.append(line.rstrip().split())
        f.close()
    return cluster

def cluster_to_map(cluster):
    cluster_id = 0
    cls_map = {}
    for cls in cluster:
        for element in cls:
            cls_map[element] = cluster_id
        cluster_id += 1
    return cls_map

def create_dot(graph, node_weights, edge_weights, dot_file):
    with open(dot_file, "w") as f:
        print >>f, "digraph G {"
#       for i in range(len(orig)):
#           print "    O%d [label=\"%d:%d\", fillcolor=\"yellow\", style=\"filled\"]" % (i, i, len(orig[i]))
        reduced = 0
        for x in graph.keys():
            for y in graph[x]:
                if len(graph[x]) == 1 and len(graph[y]) == 1:
                    reduced += 1
                    continue
                else :
                    print "%d, %d, %d, %d" % (x, y, edge_weights[(x,y)], edge_weights[(y,x)])
                    print >>f, "    %d.%d -> %d.%d [label=%d];" % (x, node_weights[x], y, node_weights[y], edge_weights[(x,y)])
        print >>f, "}"
        print "Reduced %d pairs of clusters" % (reduced/2)


def build_graph(cls_fname):
    cluster = read_cluster(cls_fname)
    cls_map = cluster_to_map(cluster)
    edge_weights = {}
    node_weights = {}
    graph = {}
    xcls = 0
    for cls in cluster:
        graph[xcls] = set([])
        node_weights[xcls] = len(cls)
        for x in cls:
            if x[-1] == '1':
                y = x[0:-1] + '2'
            else:
                y = x[0:-1] + '1'
            if y in cls_map:
                ycls = cls_map[y]
                if ycls in graph[xcls]:
                    edge_weights[(xcls, ycls)] += 1
                else :
                    graph[xcls].add(ycls)
                    edge_weights[(xcls, ycls)] = 1
        xcls += 1
    return (graph, node_weights, edge_weights)

(graph, node_weights, edge_weights) = build_graph(cls_file)
create_dot(graph, node_weights, edge_weights, "temp.dot")

