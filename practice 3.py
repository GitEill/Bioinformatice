from itertools import permutations

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

human = readGenome('Bioinformatic\chr1.GRCh38.excerpt.fasta')




def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def approxMatch(x,y):
    D= []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0  
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            dishor = D[i][j-1]+1
            disver = D[i-1][j]+1
            if x[i-1] == y[j-1]:
                disdiag = D[i-1][j-1]
            else:
                disdiag = D[i-1][j-1]+1
            D[i][j] = min(dishor, disver, disdiag)
    return min(D[-1])



dis1 = approxMatch('GCTGATCGATCGTACG', human)
dis2 = approxMatch('GATTTACCAGATTGAG', human)
dis1
dis2

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def naive_overlap_map(reads, k):
    olaps = {}
    for a, b in permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a, b)] = olen
    return olaps

def overlap_map(reads, k):
    kmer = {}
    for read in reads:
        for i in range(len(read)-k+1):
            km= read[i:i+k]
            if km not in kmer.keys():
                kmer[km] = {read}
            else:
                kmer[km].add(read)
    olaps = {}
    for i in reads:
        km = i[len(i)-k:len(i)]
        readlst = list(kmer[km])
        for j in readlst:
            if i!= j:
                olen = overlap(i, j, min_length=k)
                if olen > 0:
                    olaps[(i, j)] = olen
    return list(olaps.keys())
        
PHIseq = readFastq('Bioinformatic\ERR266411_1.for_asm.fastq')[0]
q3 = overlap_map(PHIseq,30)
len(q3)
