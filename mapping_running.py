from Bio import SeqIO
import random
from mapping_exercise import FMIndex, map_read, brute_map, reverse_complement, hamming


def run_mapping_and_checks():
    try:
        # Load files
        ref = str(SeqIO.read("reference.fasta", "fasta").seq)
        reads = list(SeqIO.parse("reads.fasta", "fasta"))
    except FileNotFoundError:
        print("The reference.fasta or reads.fasta files were not found in the same folder.")
        return

    # FMIndex
    fm = FMIndex(ref + "$", 32)

    # Map all reads
    mapped = 0
    unmapped = []
    multi_hits = 0

    for i, record in enumerate(reads):
        read_seq = str(record.seq)
        hits, _ = map_read(read_seq, fm, ref, d=2)

        if not hits:
            unmapped.append(read_seq)
        else:
            mapped += 1
            if len(hits) > 1:
                multi_hits += 1

        # Print hits for the first 5 reads
        if i < 5:
            print(f"Read {record.id} hits:")
            for hit in hits:
                strand, position, mismatches = hit
                print(f"  Strand: {strand}, Position: {position}, Mismatches: {mismatches}")

    print(f"\nTotal mapped reads: {mapped}")
    print(f"Total unmapped reads: {len(unmapped)}")
    print(f"Reads with multiple hits: {multi_hits}")

    # Self-checks

    #Verify positions and mismatches for the first 20 reads
    for record in reads[:20]:
        read_seq = str(record.seq)
        hits, _ = map_read(read_seq, fm, ref, d=2)
        for strand, pos, mm in hits:
            assert 0 <= pos <= len(ref) - 100, f"Position {pos} is out of bounds"
            query = read_seq if strand == "+" else reverse_complement(read_seq)
            actual_mm = hamming(query, ref[pos:pos + 100])
            assert actual_mm == mm, f"Mismatch count error: {actual_mm} != {mm}"
            assert actual_mm <= 2, f"Too many mismatches: {actual_mm}"

    # 30 random reads match brute_map results
    random_reads = random.sample(reads, min(30, len(reads)))
    for record in random_reads:
        read_seq = str(record.seq)
        fm_hits, _ = map_read(read_seq, fm, ref, d=2)
        bf_hits = brute_map(read_seq, ref, d=2)
        assert fm_hits == bf_hits, f"brute_map mismatch for read {record.id}"

    # No unmapped read has a hit in brute_map
    for read_seq in unmapped:
        assert not brute_map(read_seq, ref, d=2), "Unmapped read found in brute_map!"

    # Reference ends
    end_read_1 = ref[:100]
    hits_1, _ = map_read(end_read_1, fm, ref, d=2)
    assert ("+", 0, 0) in hits_1, "Failed to map first 100bp of reference"

    end_read_2 = reverse_complement(ref[-100:])
    hits_2, _ = map_read(end_read_2, fm, ref, d=2)
    assert ("-", len(ref) - 100, 0) in hits_2, "Failed to map last 100bp of reference"

if __name__ == "__main__":
    run_mapping_and_checks()