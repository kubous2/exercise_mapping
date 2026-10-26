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
class FMIndex:
    def __init__(self, text, s):
        sa = suffix_array(text)
        self.bwt = "".join(text[i - 1] for i in sa)
        self.first_occurrence = first_occurrence(self.bwt)
        self.s = s

        # self.counts[X][i] is the number of X in self.bwt[:i]
        alphabet = list(self.first_occurrence.keys())
        self.counts = {char: [0] * (len(self.bwt) + 1) for char in alphabet}

        for i, char in enumerate(self.bwt):
            for a in alphabet:
                self.counts[a][i + 1] = self.counts[a][i]
            self.counts[char][i + 1] += 1

        # Sampled suffix array
        self.sampled = {row: sa[row] for row in range(len(sa)) if sa[row] % s == 0}

    def lf(self, row):
        char = self.bwt[row]
        return self.first_occurrence[char] + self.counts[char][row]

    def count(self, p):
        top, bottom = 0, len(self.bwt)
        for char in reversed(p):
            if char not in self.first_occurrence:
                return top, top  # Pattern does not occur

            top = self.first_occurrence[char] + self.counts[char][top]
            bottom = self.first_occurrence[char] + self.counts[char][bottom]

            if top >= bottom:
                return top, top
        return top, bottom

    def locate(self, p):
        top, bottom = self.count(p)
        positions = []
        for r in range(top, bottom):
            curr_row = r
            steps = 0
            while curr_row not in self.sampled:
                curr_row = self.lf(curr_row)
                steps += 1
            positions.append(self.sampled[curr_row] + steps)
        return positions


# ================== Task 4: Mapping reads ==================
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))


def hamming(a, b):
    assert len(a) == len(b), "Strings must be of the same length"
    return sum(1 for x, y in zip(a, b) if x != y)


def seeds(read, n):
    pieces = []
    L = len(read)
    for i in range(n):
        start = i * L // n
        end = (i + 1) * L // n if i < n - 1 else L
        pieces.append((start, read[start:end]))
    return pieces


def map_read(read, fm, ref, d=2, n_seeds=None):
    if n_seeds is None:
        n_seeds = d + 1

    hits = []
    total_candidates = 0
    read_len = len(read)

    for strand, query in [("+", read), ("-", reverse_complement(read))]:
        candidates = set()
        for offset, seed in seeds(query, n_seeds):
            for p in fm.locate(seed):
                start = p - offset
                # Zabrání přesahu mimo referenční genom
                if 0 <= start <= len(ref) - read_len:
                    candidates.add(start)

        total_candidates += len(candidates)
        # Verifikace kandidátů pomocí Hammingovy vzdálenosti
        for start in candidates:
            mismatches = hamming(query, ref[start:start + read_len])
            if mismatches <= d:
                hits.append((strand, start, mismatches))

    return sorted(hits), total_candidates


def brute_map(read, ref, d=2):
    hits = []
    read_len = len(read)
    for strand, query in [("+", read), ("-", reverse_complement(read))]:
        for start in range(len(ref) - read_len + 1):
            mismatches = hamming(query, ref[start:start + read_len])
            if mismatches <= d:
                hits.append((strand, start, mismatches))
    return sorted(hits)

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
    text = "CAT$"
    s = 2
    idx = FMIndex(text, s)

    print("Text:", text)
    print("Sampling step:", s)

    print("BWT:", idx.bwt)

    print("Counts for 'T':", idx.counts["T"])

    print("Sampled suffix array:", idx.sampled)

    print("LF mapping for row 1:", idx.lf(1))

    print("Count interval for 'AT':", idx.count("AT"))

    print("Locate positions for 'AT':",idx.locate("AT"))

    print("\n=============TASK 4 ============")
    # Test reverse_complement
    print("reverse_complement('AACG'):", reverse_complement("AACG"))

    # Test hamming
    print("hamming('ACGT', 'AGGA'):", hamming("ACGT", "AGGA"))

    # Test seeds
    print("seeds('ACGTACGTAC', 3):", seeds("ACGTACGTAC", 3))

    # Malý test pro map_read a brute_map
    ref_test = "ACGTACGTACGTACGT"
    read_test = "ACGTACGT"
    fm_test = FMIndex(ref_test + "$", 2)

    hits, cands = map_read(read_test, fm_test, ref_test, d=0)
    print("\nmap_read (d=0) hits:", hits, "| candidates checked:", cands)

    print("brute_map (d=0) hits:", brute_map(read_test, ref_test, d=0))