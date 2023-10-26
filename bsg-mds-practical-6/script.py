from treeNode import *
from tools import *
from lca import *
from measure import *

def buildRestrictedTree(targetRanks, nodeFilePath):
    root_node, _, treeMapping = build_taxonomy_tree(nodes_file_path)
    filter_taxonomy_tree(root_node, target_ranks, root_node)
    return treeMapping, root_node

def buildSkeletonTree(listOfNodes, root_node):
    ancestorsByNode = findAncestors(listOfNodes)
    skeletonTreeRoot = buildSkeletonTree(root_node, ancestorsByNode, listOfNodes)
    return skeletonTreeRoot

def findLca(listOfNodes):
    ancestorsByNode = findAncestors(listOfNodes)
    lowestCommonAncestor = findCommonAncestors(listOfNodes, ancestorsByNode)[-1]
    return lowestCommonAncestor
        
if __name__ == "__main__":
    nodes_file_path = "nodes.dmp"
    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    treeMapping, rootNode  = buildRestrictedTree(target_ranks, nodes_file_path)
    sequences = extractSequence("mapping.txt", treeMapping)
    # list to stock for the root of the skeleton tree for each sequence
    # where the root is also the youngest common ancestor.
    listOfCommonAncestors = []
    for sequenceName, listOfNodes in sequences.items():
        lca = treeMapping[findLca(listOfNodes)]
        optimalNode = calculateMaximumFMeasure(listOfNodes, lca)
        print(f'{sequenceName=}, taxid of node : {optimalNode[1].taxid}, rank of node : {optimalNode[1].rank}, best fmeasure : {optimalNode[0]}')
        
