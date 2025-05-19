import matplotlib.pyplot as plt


def readFastQ(filename):
    """
    This function will take in a fastq file and parse through through the lines saving each sequence read and
    return the list of sequences and quality scores in that order
    """
    sequences = []
    quality_scores = []

    with open(filename) as f:
        while True:
            f.readline()
            read = f.readline().rstrip()
            f.readline()
            quality_of_reads = f.readline().rstrip()
            if len(read) == 0:
                break
            sequences.append(read)
            quality_scores.append(quality_of_reads)
        return sequences, quality_scores


def naive_with_mismatches(p, t, max_mismatches=2):
    """
    Implements a naive pattern matching algorithm that allows for mismatches.

    This function searches for all occurrences of pattern 'p' within text 't'
    using a brute-force approach that allows up to 'max_mismatches' differences
    between the pattern and the text.

    Args:
        p (str): The pattern string to search for
        t (str): The text string to search within
        max_mismatches (int): Maximum number of mismatches allowed (default: 2)

    Returns:
        list: A list of starting positions (0-based indices) where the pattern occurs in the text
              with no more than the specified number of mismatches
    """
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            occurrences.append(i)
    return occurrences


def naive(p, t):
    """
    Implements a naive exact pattern matching algorithm.

    This function searches for all occurrences of pattern 'p' within text 't'
    using a brute-force approach that checks each possible position in the text.

    Args:
        p (str): The pattern string to search for
        t (str): The text string to search within

    Returns:
        list: A list of starting positions (0-based indices) where the pattern occurs in the text
    """
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i + j] != p[j]:  # Fixed the comparison operator from -- to !=
                match = False
                break
        if match:  # Fixed indentation - this should be outside the inner loop
            occurrences.append(i)
    return occurrences  # Added return statement


def histogram(quals):
    """
    Create a histogram of quality scores.
    Converts ASCII quality characters to numeric values before counting.
    Additionally saves the plot to a file called 'quality_histogram.png'.
    """
    quality = []
    if not quals or len(quals) == 0:
        return [], []
    hist = [0] * 50  # Assuming quality scores range from 0-49
    for qual in quals:
        for phred in qual:
            # Convert ASCII character to quality score
            q = phred33toQuality(phred)
            if 0 <= q < len(hist):  # Ensure the index is valid
                hist[q] += 1
            quality.append(q)

    plt.figure(figsize=(10, 6))
    plt.bar(range(len(hist)), hist)
    plt.title("Quality Score Distribution")
    plt.xlabel("Quality Score")
    plt.ylabel("Frequency")
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.savefig("quality_histogram.png", dpi=300, bbox_inches="tight")
    plt.close()

    return hist, quality


def phred33toQuality(qual):
    """This function will turn the ASCII in an FastQ file to a numeric value to interperet"""

    return ord(qual) - 33


def Qualitytophred33(qual):
    """This function will turn a numerical quality value in a fastq file to it's ascii representation"""
    return chr(qual + 33)


def ReverseComplement(sequence):
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        sequence (str): A DNA sequence string containing A, C, G, T bases

    Returns:
        str: The reverse complement of the input sequence
    """
    Bases = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    complement = ""
    sequence = sequence.upper()
    for base in sequence:
        try:
            complement += Bases[base]
        except KeyError:
            print(f"invalid nucleotide {base}")
    # Return the reversed complement
    return complement


def findGCByPos(reads):
    """
    Calculate the GC content at each position across all reads.
    Additionally creates and saves a plot of GC content by position to 'gc_content_by_position.png'.

    Args:
        reads: List of sequence reads

    Returns:
        Tuple of (gc_content_by_position, total_bases_by_position)
    """
    gc = [0] * len(reads[0])
    totals = [0] * len(reads[0])

    for read in reads:
        for i in range(len(read)):
            if read[i] == "C" or read[i] == "G":
                gc[i] += 1
            totals[i] += 1

    # Calculate GC content percentage
    gc_percent = [0] * len(gc)
    for i in range(len(gc)):
        if totals[i] > 0:
            gc_percent[i] = (gc[i] / float(totals[i])) * 100

    # Create and save the plot
    plt.figure(figsize=(12, 6))
    plt.plot(range(len(gc_percent)), gc_percent, "b-", linewidth=2)
    plt.title("GC Content by Position")
    plt.xlabel("Position in Read")
    plt.ylabel("GC Content (%)")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.axhline(y=50, color="r", linestyle="--", alpha=0.5, label="50% GC content")
    plt.legend()
    plt.ylim(0, 100)
    plt.savefig("gc_content_by_position.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Return the original values for compatibility
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])

    return gc, totals


if __name__ == "__main__":
    pass
