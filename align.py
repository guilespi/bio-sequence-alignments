import re
import numpy
import ast

#read all sequences from fasta file
def read_sequences(file_name):
    f = open(file_name, 'r')
    lines = f.readlines()
    f.close()
    sequences = []
    current_sequence=None
    new_fasta_re = re.compile("^>")
    for l in lines:
        if new_fasta_re.match(l):
            print "Sequence found: %s" % l
            if current_sequence != None:
                sequences.append(current_sequence)
            current_sequence = ""
        else:
            current_sequence += l.rstrip()
    if current_sequence != None:
        sequences.append(current_sequence)
    else:
        raise NameError('Invalid fasta file empty sequence')
    return sequences

#read match matrix from file
def read_matrix(file_name):
    score_matrix = None
    with open(file_name, 'r') as f:
        s = f.read()
        score_matrix = ast.literal_eval(s)
    return score_matrix

#definitions for building paths
MATCH_PATH=2
INSERT_PATH=4
DELETE_PATH=8

#Needleman wunsch for two sequences a and b
def needleman_wunsch(sequence_a, sequence_b, score_matrix):
    gap_penalty = -3
    sequence_a = "0" + sequence_a
    sequence_b = "0" + sequence_b
    alignment_matrix = numpy.empty((len(sequence_a), len(sequence_b), ))
    path_matrix = numpy.empty((len(sequence_a), len(sequence_b), ))
    alignment_matrix[:] = numpy.NAN
    path_matrix[:] = numpy.NAN

    for i in range(len(sequence_a)):
        alignment_matrix[i, 0] = gap_penalty * i
    for j in range(len(sequence_b)):
        alignment_matrix[0, j] = gap_penalty * j

    for i in range(1, len(sequence_a)):
        for j in range(1, len(sequence_b)):
            a=sequence_a[i]
            b=sequence_b[j]
            match = score_matrix[(a, b)] if (a, b) in score_matrix else score_matrix[(b, a)]
            match = alignment_matrix[i-1, j-1] + match
            delete = alignment_matrix[i-1, j] + gap_penalty
            insert = alignment_matrix[i, j-1] + gap_penalty
            best_score = max(match, delete, insert)
            path = 0
            if best_score == match:
                path += MATCH_PATH
            if best_score == delete:
                path += DELETE_PATH
            if best_score == insert:
                path += INSERT_PATH
            alignment_matrix[i, j] = best_score
            path_matrix[i, j] = path
    return {'alignments' : alignment_matrix,
            'paths' : path_matrix}

#Create a new node for the traverse tree
def create_path_node():
    return  { 'match'  : None,
              'insert' : None,
              'delete' : None,
              'val-a'  : None,
              'val-b'  : None }

#Traverse the best score matrix and create
#a tree with the best alignments found
def recurse_path(path_matrix,
                 sequence_a,
                 sequence_b,
                 x, y,
                 parent):
    paths = path_matrix[x, y]
    paths = int(paths) if not numpy.isnan(paths) else 0
    print("Recursing %s, %s => %s" % (x,y, paths))
    if paths & MATCH_PATH:
        node = create_path_node()
        node['val-a'] = sequence_a[x-1]
        node['val-b'] = sequence_b[y-1]
        parent['match'] = node
        if paths > 0:
            recurse_path(path_matrix, sequence_a, sequence_b, x-1, y-1, node)
    if paths & INSERT_PATH or (paths == 0 and y > 0):
        node = create_path_node()
        node['val-a'] = '-'
        node['val-b'] = sequence_b[y-1]
        parent['insert'] = node
        if paths > 0:
            recurse_path(path_matrix, sequence_a, sequence_b, x, y-1, node)
    if paths & DELETE_PATH or (paths == 0 and x > 0):
        node = create_path_node()
        node['val-a'] = sequence_a[x-1]
        node['val-b'] = '-'
        parent['delete'] = node
        if paths > 0:
            recurse_path(path_matrix, sequence_a, sequence_b, x-1, y, node)
    return parent

#prints the solution found as a tree with al possible paths
def print_path_tree(t, prefix, operation, name):
    total_nodes = (1 if t['match'] else 0) + (1 if t['insert'] else 0) + (1 if t['delete'] else 0)
    def select_prefix(index):
        if index >= total_nodes:
            return "    "
        else:
            return "   |"
    if operation != None:
        print("%s %s %s" % (prefix, operation, name))
    if t['match'] != None:
        name = t['match']['val-a'] + '/' + t['match']['val-b']
        print_path_tree(t['match'], prefix + select_prefix(1), "=", name)
    if t['insert'] != None:
        name = t['insert']['val-a'] + '/' + t['insert']['val-b']
        print_path_tree(t['insert'], prefix + select_prefix(2), '+', name)
    if t['delete'] != None:
        name = t['delete']['val-a'] + '/' + t['delete']['val-b']
        print_path_tree(t['delete'], prefix + select_prefix(3), '-', name)

def main():
    sequences = read_sequences('sequences.fasta')
    score_matrix = read_matrix('blosum.txt')
    sequence_a = sequences[0]
    sequence_b = sequences[1]
    result = needleman_wunsch(sequence_a, sequence_b, score_matrix)
    print "Alignment Matrix:"
    print result['alignments']
    print "\nPaths Matrix:"
    print result['paths']
    #build result tree
    x = len(sequence_a)
    y = len(sequence_b)
    root = create_path_node()
    tree = recurse_path(result['paths'],
                        sequence_a,
                        sequence_b,
                        x, y,
                        root)
    print("\nBest score found = [%s]" % int(result['alignments'][x][y]))
    #print the tree
    print("\nSequences possible alignments tree:")
    print_path_tree(tree, "", "  *", "")

