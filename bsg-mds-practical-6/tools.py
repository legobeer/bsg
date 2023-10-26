from collections import defaultdict

def extractSequence(sequenceFile, taxonomy_tree):
    """Extract the different sequences with their name"""
    sequences = defaultdict(list)
    with open(sequenceFile, 'r') as f:
        for line in f:
            tmp = line.split()
            sequenceName = tmp[0]
            sequenceList = [int(sequence) for sequence in tmp[1:]]
            for e in sequenceList:
                sequences[sequenceName].append(taxonomy_tree[e])
    return sequences

def count_nodes(root_node):
        """Count the total number of node in the tree."""
        count = 1  
        for child in root_node.children:
            count += count_nodes(child)  
        return count