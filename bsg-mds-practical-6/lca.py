from treeNode import *

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



