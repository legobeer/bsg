from collections import defaultdict

class TreeNode:
    def __init__(self, taxid, rank, parent=None):
        self.taxid = taxid
        self.rank = rank
        self.parent = parent
        self.children = []
    
    def __str__(self) -> str:
        string = f'taxid : {self.taxid}'
        return string

def build_taxonomy_tree(nodes_file_path):
    """Build the tree"""
    taxonomy_tree = {}
    countNode = set()
    coutNodeInTargetRank = set()
    with open(nodes_file_path, 'r') as nodes_file:
        for line in nodes_file:
            parts = line.strip().split('|')
            taxid = int(parts[0].strip())
            rank = parts[2].strip()
            countNode.add(taxid)
            current_node = TreeNode(taxid, rank)
            taxonomy_tree[taxid] = current_node
    with open(nodes_file_path, 'r') as nodes_file:
        for line in nodes_file:
            parts = line.strip().split('|')
            taxid = int(parts[0].strip())
            rank = parts[2].strip()
            parent_taxid = int(parts[1].strip())
            current_node = taxonomy_tree[taxid]
            if parent_taxid in taxonomy_tree and parent_taxid != taxid:
                parent_node = taxonomy_tree[parent_taxid]
                current_node.parent = parent_node
                parent_node.children.append(current_node)
            taxonomy_tree[taxid] = current_node
    root = taxonomy_tree.get(1, None)
    return (root, len(list(countNode)), taxonomy_tree)

def filter_taxonomy_tree(start_node, target_ranks, last_valid_parent):
    """
    Build the restricted tree depending if the rank of the node
    is in target_ranks.
    """
    if start_node.rank in target_ranks or start_node.parent is None:
        for child in start_node.children[:]:
            filter_taxonomy_tree(child, target_ranks, start_node)
        return 
    children = start_node.children
    last_valid_parent.children.remove(start_node)
    last_valid_parent.children.extend(children)
    for child in children:
        child.parent = last_valid_parent
        filter_taxonomy_tree(child, target_ranks, last_valid_parent)

def findAncestors(listOfNodes):
    """
    Find all the ancestors for each nodes
    present in the list listOfNodes.
    """
    ancestersByNode = dict()
    for node in listOfNodes:
        ancestors = []
        currentAncestor = node
        while currentAncestor:
            ancestors.append(currentAncestor)
            currentAncestor = currentAncestor.parent
        ancestersByNode[node.taxid] = ancestors
    return ancestersByNode


def findLineage(node):
    """Return the lineage of the node."""
    ancestors = findAncestors([node])
    lineage = []
    ancestor = next(iter(ancestors. values()))
    for i in range(len(ancestor)):
        lineage.append(ancestor[len(ancestor) - 1 - i].rank)
    return lineage[1:-1]  


def findCommonAncestors(listOfNode, ancestorsByNode):
    """
    Return all the common ancestors of the nodes
    which are in listOfNode.
    """
    commonAncestors = []
    isBreak = False
    for values in ancestorsByNode.values():
        tmp = values[-1].taxid
        for _values in ancestorsByNode.values():
            if not _values:
                isBreak = True
                break
            if tmp != _values[-1].taxid:
                isBreak = True
                break
            _values.pop()
        if isBreak:
            break
        commonAncestors.append(tmp)
    return commonAncestors

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


def pruneTree(root, nodeFromTheSequence):
    if root is None:
        return None

    root.children = [pruneTree(child, nodeFromTheSequence) for child in root.children]

    if len(root.children) == 1 and root.taxid not in nodeFromTheSequence:
        return root.children[0]

    return root

def buildSkeletonTree(lastCommonAncestor, allAncestors, listOfNodes):
    root = TreeNode(lastCommonAncestor.taxid, lastCommonAncestor.rank)
    ancestors = set()
    nodeFromTheSequence = set(listOfNodes)
    for a in allAncestors.values():
        ancestors = ancestors.union(set(a))
    # Iterate through the ancestors and copy nodes if they are in the ancestors set
    stack = [(root, lastCommonAncestor)]
    numberOfNode = len(list(ancestors))
    count = 0
    while stack:
        new_node, old_node = stack.pop()
        common_elements = set(old_node.children) & nodeFromTheSequence
        common_count = len(common_elements)
        for child in old_node.children:
            if child in ancestors:
                new_child = TreeNode(child.taxid, child.rank, new_node)
                count += 1
                if new_child:
                    stack.append((new_child, child))
                    new_node.children.append(new_child)
    # Post-process the tree to prune and promote nodes
    root = pruneTree(root, nodeFromTheSequence)
    return root




def count_nodes(root_node):
        """Count the total number of node in the tree."""
        count = 1  
        for child in root_node.children:
            count += count_nodes(child)  
        return count
    
        
if __name__ == "__main__":
    # To test
    # nodes_file_path = "test.dmp"
    nodes_file_path = "nodes.dmp"
    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    # build the tree
    root_node, numberOfNode, treeMapping  = build_taxonomy_tree(nodes_file_path)
    # filter the tree with the target rank
    filter_taxonomy_tree(root_node, target_ranks, root_node)
    # extract the sequences
    sequences = extractSequence("mapping.txt", treeMapping)
    # list to stock for the root of the skeleton tree for each sequence
    # where the root is also the youngest common ancestor.
    listOfCommonAncestors = []
    for sequenceName, listOfNodes in sequences.items():
        ancestorsByNode = findAncestors(listOfNodes)
        ancestorsByNodeCopy = {key: value for key, value in ancestorsByNode.items()}
        lastCommonAncestor = findCommonAncestors(listOfNodes, ancestorsByNodeCopy)[-1]
        # listOfCommonAncestors.append(treeMapping[lastCommonAncestor])
        ancestorsByNode = findAncestors(listOfNodes)
        skeletonTreeRoot = buildSkeletonTree(root_node, ancestorsByNode, listOfNodes)
        numberOfNode = count_nodes(skeletonTreeRoot)
        rank = treeMapping[lastCommonAncestor].rank
        print(f'{sequenceName=}, {lastCommonAncestor=}, {numberOfNode=}, {rank=}')
    # example to print a lineage
    print(findLineage(treeMapping[1190344]))
    for lastCommonAncestor in listOfCommonAncestors:
        # for each lastCommonAncestor
        # we delete his parent to get
        # the skeleton tree
        lastCommonAncestor.parent = None
