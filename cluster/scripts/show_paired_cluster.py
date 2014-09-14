import sys
import networkx as nx
import matplotlib.pyplot as plt
import itertools as itl

if (len(sys.argv) != 2):
    print "Syntax: %s <original-cluster-file> <new-cluster-file>" % (sys.argv[0])
    exit(0)

orig_cls_file = sys.argv[1]

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

def build_graph(orig_fname):
    orig_cluster = read_cluster(orig_fname)
    orig_cls_map = cluster_to_map(orig_cluster)
    edges = {}
    nodes = []
    cluster_id = 0
    for cls in orig_cluster:
        nodes.append(1)#len(cls))
        cluster_id += 1
        for x in cls:
            if x[-1]=='1':
                xcls = orig_cls_map[x]
                ycls = orig_cls_map[x[0:-1]+'2']
#            print x, y
#            print xcls, ycls
#            print
                if (xcls != ycls):
                    mn = min(xcls, ycls)
                    mx = max(xcls, ycls)
                    if (mn, mx) in edges:
                        edges[(mn, mx)] += 1
                    else :
                        edges[(mn, mx)] = 1
                elif(xcls==ycls):
                    print "pairs should not be in the same cluster"
                    exit(0)

    print edges
    print nodes
    return (nodes, edges)


def draw_graph(nodes, edges, labels=None, graph_layout='shell',
               node_size=800, node_color='blue', node_alpha=0.3,
               node_text_size=10,
               edge_color='black', edge_alpha=0.3, edge_tickness=5,
               edge_text_pos=0.5,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.Graph()


    # add edges
    for i in range(len(nodes)):
        G.add_node(i)

    # add edges
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos, node_size=[400+50*x for x in nodes],
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(edges))

    edge_labels = dict(zip(edges, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels,
                                 label_pos=edge_text_pos)

    #nx.draw_graphviz(G)

    # show graph
    plt.show()

    nx.write_dot(G,"file.dot")

graph = [(0, 1), (1, 5), (1, 7), (4, 5), (4, 8), (1, 6), (3, 7), (5, 9),
         (2, 4), (0, 4), (2, 5), (3, 6), (8, 9)]

# you may name your edge labels
labels = map(chr, range(65, 65+len(graph)))

#graph = build_graph(orig_cls_file, new_cls_file)
(nodes, edges) = build_graph(orig_cls_file)
# if edge labels is not specified, numeric labels (0, 1, 2...) will be used
draw_graph(nodes, edges.keys(), [str(x) for x in edges.values()], graph_layout='random')

