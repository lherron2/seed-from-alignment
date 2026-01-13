import argparse
import sys
import heapq

def get_levenshtein_distance(s1, s2):
    """
    Calculates the Levenshtein distance between two strings.
    """
    if len(s1) < len(s2):
        return get_levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def sanitize_structure(structure):
    """
    Removes whitespace and alignment gap characters (-, ~, _).
    Preserves '.' as it represents unpaired bases.
    """
    s = structure.strip()
    # Remove alignment gaps
    s = s.replace('-', '').replace('~', '').replace('_', '')
    return s

def main():
    parser = argparse.ArgumentParser(
        description="Find the top K nearest dot-bracket structure matches in a database file."
    )
    parser.add_argument(
        '--key', 
        required=True, 
        help="Path to the file containing the query dot-bracket structure."
    )
    parser.add_argument(
        '--db', 
        required=True, 
        help="Path to the database file containing rows of structures."
    )
    parser.add_argument(
        '--top', '-k',
        type=int,
        default=10,
        help="Number of matches to return (default: 10)."
    )

    args = parser.parse_args()

    # 1. Load and Sanitize the Query Key
    query_structure = ""
    try:
        with open(args.key, 'r') as kf:
            for line in kf:
                clean_line = sanitize_structure(line)
                if clean_line:
                    query_structure = clean_line
                    break
        
        if not query_structure:
            print("Error: Key file appears to be empty or invalid.")
            sys.exit(1)
            
    except FileNotFoundError:
        print(f"Error: Key file '{args.key}' not found.")
        sys.exit(1)

    print(f"""Query Structure (sanitized): 
{query_structure}""")
    print(f"Length: {len(query_structure)}")
    print("=" * 60)

    # 2. Search the Database using a Heap
    # We use a min-heap to store the top K results.
    # To keep the 'smallest' distances, we store (-distance) in the heap.
    # Python's heap is a min-heap. If we push negative numbers, the "smallest"
    # number in the heap is actually the largest distance (worst match) in our top K.
    # When the heap grows > K, we pop the smallest number (the worst match), 
    # leaving us with the best matches.
    
    top_k_heap = []

    try:
        with open(args.db, 'r') as dbf:
            # Skip header
            next(dbf, None)

            for line_num, line in enumerate(dbf, start=2):
                target = sanitize_structure(line)
                
                if not target:
                    continue

                # Calculate distance
                dist = get_levenshtein_distance(query_structure, target)

                # Store in heap: (-distance, -line_num, content)
                # We use -line_num to handle tie-breaking consistently
                # (though strict order matters less here).
                entry = (-dist, line_num, line.strip())
                
                heapq.heappush(top_k_heap, entry)

                # Maintain only K items
                if len(top_k_heap) > args.top:
                    heapq.heappop(top_k_heap)

    except FileNotFoundError:
        print(f"Error: Database file '{args.db}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

    # 3. Process and Sort Results
    # Convert heap back to a sorted list of positive distances
    results = []
    while top_k_heap:
        neg_dist, line_num, content = heapq.heappop(top_k_heap)
        results.append((-neg_dist, line_num, content))
    
    # Sort by distance (ascending), then by line number
    results.sort(key=lambda x: (x[0], x[1]))

    # 4. Output Results
    if results:
        print(f"Top {len(results)} Matches:\n")
        for i, (dist, line_num, content) in enumerate(results):
            print(f"Rank {i+1} | Line: {line_num} | Edit Distance: {dist}")
            # Printing structure on a new line as requested
            print(content)
            print("-" * 40)
    else:
        print("No valid structures found in the database.")

if __name__ == "__main__":
    main()
