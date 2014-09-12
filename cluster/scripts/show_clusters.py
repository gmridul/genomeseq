import sys
import networkx as nx
import matplotlib.pyplot as plt
import itertools as itl

if (len(sys.argv) != 3):
    print "Syntax: %s <original-cluster-file> <new-cluster-file>" % (sys.argv[0])
    exit(0)

orig_cls_file = sys.argv[1]
new_cls_file = sys.argv[2]

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

def create_dot(orig, new, edges):
    print "digraph G {"
    for i in range(len(orig)):
        print "    O%d [label=\"%d:%d\", fillcolor=\"blue\", style=\"filled\"]" % (i, i, len(orig[i]))
    for i in range(len(new)):
        print "    N%d [label=\"%d:%d\" fillcolor=\"red\", style=\"filled\"]" % (i, i, len(new[i]))
    for (x,y) in edges:
        print "    O%d -> N%d [label=%d];" % (x, y, edges[(x,y)])
    print "}"


def build_graph(orig_fname, new_fname):
    orig_cluster = read_cluster(orig_fname)
    new_cluster = read_cluster(new_fname)
    new_cls_map = cluster_to_map(new_cluster)
    edges = {}
    cluster_id = 0
    for cls in orig_cluster:
        for y in cls:
            xcls = cluster_id
            ycls = new_cls_map[y]
            if (xcls, ycls) in edges:
                edges[(xcls, ycls)] += 1
            else :
                edges[(xcls, ycls)] = 1
        cluster_id += 1
    create_dot(orig_cluster, new_cluster, edges)

build_graph(orig_cls_file, new_cls_file)

