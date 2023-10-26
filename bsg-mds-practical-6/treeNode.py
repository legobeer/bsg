class TreeNode:
    def __init__(self, taxid, rank, parent=None):
        self.taxid = taxid
        self.rank = rank
        self.parent = parent
        self.children = []
    
    def __str__(self) -> str:
        string = f'taxid : {self.taxid}'
        return string
    
    def getLeaves(self):
        leaves = set()  # Create an empty set to store the leaves

        # Use a stack for iterative DFS
        stack = [self]

        while stack:
            node = stack.pop()

            if not node.children:
                # If the node has no children, it's a leaf, so add it to the set
                leaves.add(node)

            # Add the children to the stack for further traversal
            stack.extend(node.children)

        return leaves


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

