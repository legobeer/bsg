from lca import *
from treeNode import *
from tools import *

def calculateFMeasure(read, rootNode):
    if (not rootNode.children):
        return 0
    leaves = rootNode.getLeaves()
    intersection = leaves & set(read)
    tp = len(intersection)
    if not tp:
        return 0
    fn = len(read) - tp
    fp = len(leaves) - tp
    p = abs(tp) / ((abs(tp) + abs(fp)))
    r = abs(tp) / ((abs(tp) + abs(fn)))
    return 2 * p * r / (p + r)

def calculateMaximumFMeasure(read, lca):
    fMeasure = calculateFMeasure(read, lca)
    result = (fMeasure, lca)
    stack = [lca]
    while stack:
        node = stack.pop()
        fMeasure = calculateFMeasure(read, node)
        if fMeasure >= result[0]:
            result = (fMeasure, node)
        for child in node.children:
            stack.append(child)
    return result
            
    

