# =========== Task 1: Suffix array and binary search===========
def brute_force(text, p):
    positions = []
    if not p:
        return positions

    for i in range(len(text) - len(p) + 1):
        if text[i:i + len(p)] == p:
            positions.append(i)
    return positions


def suffix_array(text):
    suffixes_with_indices = []
    for i in range(len(text)):
        suffix = text[i:]
        suffixes_with_indices.append((i, suffix))

    def sort_key(item):
        # Receives a tuple (index, suffix), returns the suffix
        return item[1]

        # Sort the list of tuples, NOT range()

    suffixes_with_indices.sort(key=sort_key)

    sa = []
    for item in suffixes_with_indices:
        index = item[0]
        sa.append(index)

    return sa


def find_block(text, sa, p):
    #uses two binary searches
    plen = len(p)

    # binary searcg block a (first suffix >= p)
    lo, hi = 0, len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + plen] < p:
            lo = mid + 1
        else:
            hi = mid
    a = lo

    # binary search block b (first > p)
    lo, hi = 0, len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        if text[sa[mid]:sa[mid] + plen] <= p:
            lo = mid + 1
        else:
            hi = mid
    b = lo

    return a, b

# ============ Task 2: The Burrows-Wheeler transform and its inverse=======================
def bwt(text, sa):
    return "".join(text[i - 1] for i in sa)


def first_occurrence(letters):
    counts = {}
    for char in letters:
        counts[char] = counts.get(char, 0) + 1

    first_occ = {}
    total = 0
    for char in sorted(counts.keys()):
        first_occ[char] = total
        total += counts[char]

    return first_occ


def inverse_bwt(bwt_string):
    # Compute the occurrence number k for each letter in L (the BWT string)
    occ = []
    counts = {}
    for char in bwt_string:
        counts[char] = counts.get(char, 0) + 1
        occ.append(counts[char])

    # get the FirstOccurrence array
    first_occ = first_occurrence(bwt_string)

    # traverse backwards using the LF mapping
    row = 0
    text_rev = []
    for _ in range(len(bwt_string) - 1):
        char = bwt_string[row]
        k = occ[row]
        text_rev.append(char)
        row = first_occ[char] + k - 1

    # Reverse the collected characters and append the $ sign
    return "".join(reversed(text_rev)) + "$"

# ================== Task 3: The FMIndex class==================


#============ Self-checks =========
if __name__ == "__main__":
    # ================ TASK 1 ===========
    text = "AAGGCGGCTTCA$"
    p = "GGC"

    print("Text:", text)
    print("Pattern:", p)

    # Test brute_force
    print("\n=============TASK 1 ============")
    print("Output of brute_force:")
    print(brute_force(text, p))

    # Test suffix_array
    sa = suffix_array(text)
    print("Output of suffix_array (SA):")
    print(sa)

    # Test find_block
    a, b = find_block(text, sa, p)
    print("Output of find_block (start and end of the block):")
    print(a, b)

    # Print specific positions from the found block
    print("Pattern positions obtained using SA:")
    print(sa[a:b])

    print("\n=============TASK 2 ============")
    text = "CAT$"
    sa = [3, 1, 0, 2]

    print("Original text:", text)

    # Test bwt
    b = bwt(text, sa)
    print("BWT:", b)

    # Test first_occurrence
    f_occ = first_occurrence(b)
    print("First occurrence:", f_occ)

    # Test inverse_bwt
    recovered_text = inverse_bwt(b)
    print("Recovered text:", recovered_text)

    print("\n=============TASK 3 ============")