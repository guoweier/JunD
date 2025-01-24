# Weier Guo Python
# Create: 01/15/2025
# Introduction: Use PRICE assembly for producing junction contigs
#### UPDATE: to add mismatches for contigs ####

## PACKAGES ##
import sys, math, os, time
import argparse
from collections import defaultdict, deque, Counter

## PREPROCESSING ##
# Remove duplicate reads to reduce redundancy.
def deduplicate_reads(seed):
    """
    Remove duplicated sequence reads in .fasta file.
    Args:
        seed (file): .fasta file carrying seeds.
    Returns:
        list: List of all unique sequence. 
    """
    seeds = []
    with open(seed, "r") as s:
        for line in s:
            if line[0] == ">":
                name = line[:-1]
            else:
                seq = line.split("\n")[0]
                if seq not in seeds:
                    seeds.append(seq)
                else:
                    continue 
    return seeds

# Remove subset short reads.
def find_subset(seq1, seq2):
    """
    Find subset maximum overlaps between seq1 and seq2. 
    Args:
        seq1 (str): Read sequence 1. 
        seq2 (str): Read sequence 2. 
    Returns:
        int: maximum overlaps between seq1 and seq2. 
    """
    max_overlap = 0
    # Iterate through possible overlap lengths
    for i in range(min(len(seq1), len(seq2))+1):
        if seq1[-i:] == seq2[:i]:
            max_overlap = i 
        elif seq2[-i:] == seq1[:i]:
            max_overlap = -i
    return max_overlap

def desubset_seq(seq1, seq2):
    """
    Remove subset sequence.
    Args:
        seq1 (str): Read sequence 1. 
        seq2 (str): Read sequence 2. 
    Returns: 
        list: If desubset successfully, only the long seq in list; if no desubset, both seqs in list.
    """
    subset_length = find_subset(seq1, seq2)
    if abs(subset_length) == len(seq1) and len(seq1) < len(seq2):
        return [seq2]
    elif abs(subset_length) == len(seq2) and len(seq2) < len(seq1):
        return [seq1]
    else:
        return [seq1, seq2]

def desub_reads(seeds):
    """
    Remove short reads that are subset of long reads.
    Args: 
        seeds (list): List of read sequence.
    Return: 
        list: List of read sequence with no subset short reads.
    """
    desub_seqs = []
    while seeds:
        seq = seeds.pop(0)
        desub = False 
        for i in range(len(seeds)):
            other_seq = seeds[i]
            desub_seq = desubset_seq(seq, other_seq)
            if len(desub_seq) == 1:
                # Desub successfully
                seeds[i] = "NA"
                desub = True 
                break 
        if not desub:
            if seq != "NA":
                desub_seqs.append(seq)
    uni_desub_seqs = list(set(desub_seqs))
    return uni_desub_seqs

## BUILD DE BRUIJN GRAPH ##
def calculate_overlap(node1, node2):
    """
    Calculates the maximum overlap between two strings (nodes).
    Args:
        node1 (str): The first string.
        node2 (str): The second string. 
    Returns:
        int: The length of the maximum overlap.
    """
    max_overlap = 0
    len1, len2 = len(node1), len(node2)
    # Check suffix of node1 against prefix of node2
    for i in range(1, min(len1, len2)+1):
        if node1[-i:] == node2[:i]:
            max_overlap = i 
    return max_overlap

def add_graph_overlap(graph):
    """
    Add overlap value for each edge in simplified graph.
    Args:
        graph (dict): Directed graph represented as an adjacency list. Keys are nodes, values are lists of neighbors.
    Returns: 
        graph (dict): Directed graph represented as an adjacency list. Keys are nodes, values are lists of (neighbor, overlap).
    """
    graph_overlap = {}
    for node in graph:
        for value in graph[node]:
            overlap = calculate_overlap(node, value)
            if node not in graph_overlap:
                graph_overlap[node] = []
            graph_overlap[node] += [(value, overlap)]
    return graph_overlap

def build_de_bruijn_graph(seeds, k):
    """
    Builds a De Bruijn graph from a set of DNA sequences using kmers as nodes.
    Args:
        seeds (list): List carrying reads from .fasta file.  
    Returns:
        dict: De Bruijn graph represented as an adjacency list. Keys are kmers (nodes), and values are lists of (neighbor, overlap).
    """
    de_bruijn_graph = {}
    for read in seeds:
        # Generate kmers from the read
        for i in range(len(read)-k+1):
            kmer_prefix = read[i:i+k-1] 
            kmer_suffix = read[i+1:i+k]
            # Add edge to graph
            if kmer_prefix not in de_bruijn_graph:
                de_bruijn_graph[kmer_prefix] = []
            if kmer_suffix in de_bruijn_graph[kmer_prefix]:
                continue 
            else:
                de_bruijn_graph[kmer_prefix] += [kmer_suffix]
    graph_overlap = add_graph_overlap(de_bruijn_graph)
    return dict(graph_overlap)

## SIMPLIFY GRAPH ##
# step2 (remove tips) functions
def remove_tips(graph, max_tip_length):
    """
    Remove tips (dead-end paths) in a De Bruijn graph.
    Args:
        graph (dict): The De Bruijn graph represented as an adjacency list. Keys are nodes, and values are lists of tuples (neighbor, overlap).
        max_tip_length (int): The maximum length for a tip. Tips longer than this length are not removed.
    Returns: 
        dict: The De Bruijn graph with tips removed.
    """
    def calculate_path_length(path):
        """Calculate the total length of a path."""
        total_length = 0
        for node, overlap in path:
            total_length += len(node)*2-overlap
        return total_length

    def is_startnode(node):
        """Define the node is the appropriate node (dead-end node or its parent is outgoing node) for tracing tips."""
        return sum(1 for edges in list(graph.values()) if (node in dict(edges) and len(edges) > 1)) > 0

    
    def find_tip(start_node):
        """
        Trace a tip starting from a given node.
        Args:
            start_node (str): The node to start tracing the tip from.
        Returns: 
            list: The nodes in the tip path.
        """
        tip_path = [(start_node, len(start_node))]
        current_node = start_node 
        while current_node in graph:
            if len(graph[current_node]) == 1:
                next_node, overlap = graph[current_node][0]
                if sum(1 for edges in list(graph.values()) if next_node in dict(edges)) == 1:
                    tip_path.append((next_node, overlap))
                    current_node = next_node
                elif sum(1 for edges in list(graph.values()) if next_node in dict(edges)) > 1:
                    tip_path = []
                    break
            elif len(graph[current_node]) > 1:
                tip_path = []
                break 
            elif len(graph[current_node]) == 0:
                break 
        return tip_path

    # Identify and remove tips
    nodes_to_remove = set()
    for node in list(graph.keys()):
        # Check for dead-end paths (no outgoing edges)
        if is_startnode(node):
            tip_path = find_tip(node)
            tip_length = calculate_path_length(tip_path)
            # Mark nodes in the tip for removal if it's shorter than the max_tip_length
            if tip_length <= max_tip_length:
                nodes_to_remove.update(tip_path)
    # Remove marked nodes from the graph
    for node, overlap in nodes_to_remove:
        graph.pop(node, None)
        for neighbors in list(graph.values()):
            neighbors[:] = [n for n in neighbors if n[0] != node]
    return graph


# step3 (resolve bubble) functions
# find reconvergeing nodes (breadth-first search)
def find_reconverge_nodes(outgoing_edges, graph):
    """
    Find nodes where paths from multiple outgoing edges of a given node reconverge. 
    Args:
        outgoing_edges (list): List of nodes connected by outgoing edges.
        graph (dict): Be Bruijn graph represented as an adjacency list.
    Returns:
        set: A set of nodes where paths from the outgoing eges reconverge. 
    """
    visited = defaultdict(set)
    reconverging_nodes = set()
    # BFS starting from each outgoing edge
    for edge in outgoing_edges:
        queue = deque([edge])
        while queue:
            current = queue.popleft()
            if current in visited and len(visited[current]) >= 1:
                reconverging_nodes.add(current)
            visited[current].add(edge)
            for neighbor in graph.get(current, []):
                if neighbor not in visited[current]:
                    queue.append(neighbor)
    return reconverging_nodes

# resolve bubble by path confidence
def resolve_bubble(node, reconverge_node, graph, path_confident_func=None):
    """
    Resolve a bubble in a De Bruijn graph by retaining the most confident path.
    Args:
        node (str): The starting node where the bubble diverges.
        reconverge_node (str): The node where the paths reconverge. 
        graph (dict): De Bruijn graph represented as an adjacency list. 
        path_confidence_func (function): Optional. A function to calculate path confidence. Defaults to path length if not provided. 
    Returns:
        None: Modifies the graph in-place to remove the less confident path. 
    """
    # nested function to find all paths from start to end
    def find_paths(start, end):
        paths = []
        queue = deque([[start]])
        while queue:
            path = queue.popleft()
            current = path[-1]
            if current == end:
                paths.append(path)
                continue 
            for neighbor in graph.get(current, []):
                if neighbor not in path:
                    queue.append(path+[neighbor])
        return paths 
    paths = find_paths(node, reconverge_node)
    # Default confidence function: path length (longer paths are less confident)
    if path_confident_func is None:
        path_confident_func = lambda path: -len(path)
    path_confidences = sorted([(path, path_confident_func(path)) for path in paths], key=lambda x: x[1], reverse=True)
    best_path = path_confidences[0][0]
    for path, _ in path_confidences[1:]:
        for i in range(len(path)-1):
            start, end = path[i], path[i+1]
            if end in graph.get(start, []):
                graph[start].remove(end)
                if not graph[start]:
                    del graph[start]

# step4 (remove small cycle) functions
# find cycles
def find_cycles(graph):
    """
    Detected all cycles in a directed graph.
    Args:
        graph (dict): Directed graph represented as an adjacency list. Keys are nodes, and values are lists of neighbors. 
    Returns:
        list: A list of cycles, where each cycle is represeneted as a list of nodes. 
    """
    def dfs(node, path, visited, stack, cycles):
        visited.add(node)
        stack.add(node)
        path.append(node)
        for neighbor in graph.get(node, []):
            if neighbor not in visited:
                # Recur if the neighbor hasn't been visited
                dfs(neighbor, path, visited, stack, cycles)
            elif neighbor in stack:
                # Cycle detected, extract the cycle
                cycle_start_index = path.index(neighbor)
                cycles.append(path[cycle_start_index:])
        # Remove the current node from the stack and path 
        stack.remove(node)
        path.pop()
    visited = set()
    stack = set()
    cycles = []
    for node in graph:
        if node not in visited:
            dfs(node, [], visited, stack, cycles)
    return cycles 

# remove cycles
def remove_cycle(cycle, graph):
    for i in range(len(cycle)):
        # Get the current node and the next node in the cycle
        current_node = cycle[i]
        next_node = cycle[(i+1)%len(cycle)]
        # Remove the edge from the graph
        if next_node in graph.get(current_node, []):
            graph[current_node].remove(next_node)
        # If current node has no more outgoing edges, remove it. 
        if not graph[current_node]:
            del graph[current_node]

# main function for simplify graph
def simplify_graph(graph_ini, threshold_length, max_tip_length):
    # 2. remove tips
    graph = remove_tips(graph_ini, max_tip_length)
    #print(graph)
    # 3. resolve bubbles
    for node in list(graph.keys()):
        outgoing_edges = graph[node]
        if len(outgoing_edges) > 1:
            reconverging_nodes = find_reconverge_nodes(outgoing_edges, graph)
            for reconverge_node in reconverging_nodes:
                resolve_bubble(node, reconverge_node, graph)
    # 4. remove small cycles
    for cycle in find_cycles(graph):
        if len(cycle) < threshold_length:
            remove_cycle(cycle, graph) 
    return graph

## ASSEMBLE CONTIGS ##
def traverse_graph(graph, min_overlap):
    visited_edges = set()
    paths = []
    def is_linear(node):
        """Check if the given node is part of a linear path."""
        neighbor, overlap = graph[node][0]
        return len(graph[node]) == 1 and (neighbor in list(graph.keys())) and len(graph[neighbor]) == 1

    def is_startnode(node):
        """Define the start node (dead-end)."""
        return sum(1 for n, edges in graph.items() if node in dict(edges)) == 0

    def merge_kmers(path):
        """Merge kmers in the path."""
        kmer1 = path.pop(0)
        sequence = kmer1
        while len(path) > 0:
            kmer2 = path.pop(0)
            overlap = find_overlap(sequence, kmer2, min_overlap)
            if overlap + 1 == len(kmer2):
                sequence += kmer2[overlap:]
        return sequence

    def extend_path(start_node):
        """
        Extends a path in a graph while respecting the minimum overlap constraint.
        Args:
            node (str): The starting node for the extension.
            path (list): The current path being extended.
            graph (dict): Directed graph represented as an adjacency list. Keys are nodes, values are lists of (neighbor, overlap).
            min_overlap (int): Minimum overlap required to extend the path.
            visited_edges (set): Set of visited edges (tuples of (node, neighbor)).
        Returns:
            list: The extended path.
        """
        path = [start_node]
        current_node = start_node
        while current_node in graph and is_linear(current_node):
            next_node, overlap = graph[current_node][0]
            if overlap >= min_overlap:
                path.append(next_node)
                visited_edges.add(current_node)
                current_node = next_node
        visited_edges.add(current_node)
        sequence = merge_kmers(path)
        return sequence 

    for node in graph:
        if is_startnode(node):
            sequence = extend_path(node)
            paths.append(sequence)
    return paths


## POST-ASSEMBLY ##
# correct contig error by aligning reads
def align_read_to_contig(read, contig):
    """
    Align a single read to a contig and return the aligned read and contig positions.
    Args:
        read (str): The read sequence.
        contig (str): The contig sequence.
    Returns: 
        list: List of aligned positions (index pairs) where the read matches the contig.
    """
    aligned_positions = []
    for i in range(len(contig)-len(read)+1):
        match_count = sum(1 for j in range(len(read)) if contig[i+j]==read[j])
        if match_count == len(read):
            aligned_positions.append(i)
    return aligned_positions

def correct_contig_with_reads(contig, reads, min_overlap):
    """
    Corrects a contig by aligning it with the reads and adjusting bases according to the majority vote. 
    Args:
        contig (str): The contig sequence. 
        reads (list): List of read sequences.
        min_overlap (int): Minimum overlap of a read to a contig for correction.
    Returns:
        str: The error-corrected contig.
    """
    corrected_contig = list(contig)
    # For each position in the contig, correct it using the aligned reads.
    for i in range(len(contig)):
        aligned_bases = []
        for read in reads:
            if i+len(read) <= len(contig):
                # Get aligned positions for this read
                alignments = align_read_to_contig(read, contig)
                for pos in alignments:
                    if pos <= i < pos+len(read):
                        aligned_bases.append(read[i-pos])
        if aligned_bases:
            # Use the majority base to correct the position
            base_count = Counter(aligned_bases)
            most_common_base, _ = base_count.most_common(1)[0]
            corrected_contig[i] = most_common_base
    return ''.join(corrected_contig)

def error_correction(contigs, reads, min_overlap):
    """
    Perform error correction on a list of contigs using the reads.
    Args:
        contigs (list): List of contig sequences.
        reads (list): List of read sequences. 
        min_overlap (int): Minimum overlap of a read to a contig for correction.
    Returns:
        list: List of error-corrected contigs. 
    """
    corrected_contigs = []
    for contig in contigs:
        corrected_contig = correct_contig_with_reads(contig, reads, min_overlap)
        corrected_contigs.append(corrected_contig)
    return corrected_contigs

# merge overlapping contigs
def find_overlap(contig1, contig2, min_overlap):
    """
    Finds the longest overlap between two contigs where the suffix of contig1 matches the prefix of contig2. 
    Args:
        contig1 (str): The first contig.
        contig2 (str): The second contig.
        min_overlap (int): The minimum overlap required to merge the contigs. 
    Returns:
        int: The length of the overlap, or 0 if no valid overlap exists.
    """
    max_overlap = 0
    # Iterate through possible overlap lengths
    for i in range(min_overlap, min(len(contig1), len(contig2))+1):
        if contig1[-i:] == contig2[:i]:
            max_overlap = i 
        elif contig2[-i:] == contig1[:i]:
            max_overlap = -i
    return max_overlap

def merge_two_contigs(contig1, contig2, min_overlap):
    """
    Merges two contigs if they overlap, otherwise return the two contigs as is. 
    Args:
        contig1 (str): The first contig.
        contig2 (str): The second contig.
        min_overlap (int): The minimum overlap required to merge the contigs.
    Returns:
        str: The merged contig.
    """
    overlap_length = find_overlap(contig1, contig2, min_overlap)
    if overlap_length > 0:
        # Merge the contigs by appending the non-overlapping portion of contig2 to contig1. 
        return contig1+contig2[overlap_length:]
    elif overlap_length < 0:
        # Merge the contigs by appending the non-overlapping portion of contig1 to contig2.
        overlap_length_m = -overlap_length
        return contig2+contig1[overlap_length_m:]
    else:
        # No overlap, return the contigs as a list of separate contigs
        return contig1, contig2 

def get_reverse_complement_seq(seq):
    """
    Get reverse complement sequence.
    Args:
        seq (str): read sequence.
    Returns:
        str: reverse complemented sequence.
    """
    def complement_rule(n):
        """Complement ATGC."""
        if n == "A":
            out = "T"
        elif n == "T":
            out = "A"
        elif n == "C":
            out = "G"
        elif n == "G":
            out = "C"
        return out 

    co_seq = []
    for n in seq:
        out = complement_rule(n)
        co_seq.append(out)
    reco_seq = list(reversed(co_seq))
    text = ''.join(reco_seq)
    return text

def merge_overlapping_contigs(contigs, min_overlap):
    """
    Merges overlapping contigs in a list of contigs.
    Args:
        contigs (list): List of error-corrected contigs (strings).
        min_overlap (int): Minimum overlap required to merge the contigs.
    Returns:
        list: List of merged contigs.
    """
    # merge kmers inside each extended path
    merged_contigs = []
    while contigs:
        contig = contigs.pop(0)
        merged = False 
        for i in range(len(contigs)):
            other_contig = contigs[i]
            merged_contig = merge_two_contigs(contig, other_contig, min_overlap)
            if isinstance(merged_contig, str):
                # Merged successfully
                contigs[i] = merged_contig
                merged = True 
                break 
        if not merged:
            # If no merge happened, keep the contig in the result 
            merged_contigs.append(contig)

    # add reverse compliment sequence
    reco_merged_contigs = []
    for seq in merged_contigs:
        reco_seq = get_reverse_complement_seq(seq)
        reco_merged_contigs.append(reco_seq)
    merged_contigs = merged_contigs+reco_merged_contigs

    # merge merged_contigs
    final_contigs = []
    while merged_contigs:
        merged_contig = merged_contigs.pop(0)
        merged = False 
        for j in range(len(merged_contigs)):
            other_merged_contig = merged_contigs[j]
            final_contig = merge_two_contigs(merged_contig, other_merged_contig, min_overlap)
            if isinstance(final_contig, str):
                # Merge merged_contig successfully
                merged_contigs[j] = final_contig
                merged = True 
                break 
        if not merged:
            # If no merge for merged_contigs, keep the contig in the result
            final_contigs.append(merged_contig)
    print(final_contigs)
    uni_final_contigs = desub_reads(final_contigs)
    return uni_final_contigs

## MAIN ##
def de_novo_assembly(seeds, k, threshold_length, min_overlap, max_tip_length):
    # 1. Generate kmers
    kmer_graph = build_de_bruijn_graph(seeds, k)
    # 2. simplify graph
    simplified_graph = simplify_graph(kmer_graph, threshold_length, max_tip_length)
    # 3. Assemble contigs
    contigs = traverse_graph(simplified_graph, min_overlap)
    # 4. Post-assembly
    error_corrected_contigs = error_correction(contigs, seeds, min_overlap)
    merged_contigs = merge_overlapping_contigs(error_corrected_contigs, min_overlap)

    return merged_contigs

def write_contigs(seedpath, k, threshold_length, min_overlap, max_tip_length):
    li = os.listdir(seedpath)
    seedfiles = list(filter(lambda x: x.endswith(".fasta") or x.endswith(".fa"), li))
    os.system(f"mkdir contigs")
    for seed in seedfiles:
        name = seed[6:].split(".")[0]
        print(name)

        with open(f"contigs/contigs_{name}.fasta", "w") as o:
            seeds = deduplicate_reads(seedpath+seed)
            merged_contigs = de_novo_assembly(seeds, k, threshold_length, min_overlap, max_tip_length)
            for i in range(len(merged_contigs)):
                o.write(f">Contig{i+1}\n")
                text = merged_contigs[i]
                o.write(f"{text}\n")


## CALL ARGUMENTS ##
def parse_arguments():
    parser = argparse.ArgumentParser(description="Contig assembly")
    parser.add_argument("--seedpath", type=str, default="seeds/", help="Path for seed files (default=seeds/).")
    parser.add_argument("--k", type=int, default=30, help="Kmer length (default=30).")
    parser.add_argument("--threshold_length", type=int, default=3, help="The threshold for removing cycles (default=3).")
    parser.add_argument("--min_overlap", type=int, default=25, help="Minimum kmer overlaps (default=25).")
    parser.add_argument("--max_tiplength", type=int, default=30, help="Maximum tip length (default=30).")
    return parser.parse_args()

if __name__ in "__main__":
    args = parse_arguments()
    write_contigs(args.seedpath, args.k, args.threshold_length, args.min_overlap, args.max_tiplength)
    


