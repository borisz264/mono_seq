import os
import operator
import itertools
import gzip
import numpy as np
from scipy import stats

def old_chris_formula(R, k, read_len):
    """
     Implements the formula Chris burge derived
    """
    return (4 ** k -1) * ( R * (read_len - k + 1) - (read_len - k) ) / (4 ** k + read_len - k - (R *(read_len - k +1)))


def chris_formula(R, k, read_len):
    """
     Implements the formula Chris burge derived
    """
    return (4. ** k - 1) * (read_len - k +1 ) * (R - 1) / (4. ** k -R * (read_len -k +1))


def get_adjacent_kmers(kmer):
    """
    returns all the k+1 mers that contain the kmer
    """
    return ['A' + kmer, 'C' + kmer, 'G' + kmer, 'T' + kmer,
      kmer + 'A', kmer + 'C', kmer + 'G', kmer + 'T']


def iter_RNAfold_output(energy_file):
    """
    iterates through RNAfold input and returns an energy iterator
    """
    for l1, l2 in iterLinePairs(energy_file):
        yield float(l2.split(' (')[1].replace(')', ''))


def iterLinePairs(inFile):
    for l1, l2 in iterNlines(inFile, 2):
        yield l1, l2


def make_dir(dirname):
    """
    Makes the directory; doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            print 'The directory was made by another thread extremely recently.'


def file_exists(fname):
    """
    makes sure a given file exists
    """
    if not os.path.exists(fname):
        return False
    fstats = os.stat(fname)
    if not fstats[6]:
        return False
    if not os.access(fname, os.R_OK):
        raise ValueError('Input File %s cannot be read' % fname)
    return True


def getBinIndex_soft_upper(v, bins):
    for i in range(len(bins) - 1):
        if v > bins[i] and v <= bins[i + 1]:
            return i
    return -1


def get_barcode(line):
    """
    - Extracts the barcode from the first line of a fastq quartet
        - Assumes the first line is of the form:
            @D5FF8JN1:4:1101:1220:2099#ACTTGA/1
    """
    return line.split('#')[-1].split('/')[0]


def get_index_from_kmer(kmer):
    """
    returns the base10 version of the base 4 DNA representation
    """
    index = 0
    base2face = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(kmer):
        if not base in 'ACGT':
            return -1
        power = len(kmer) - 1 - i
        index += base2face[base] * (4 ** power)
    return index


def get_kmer_from_index(kmax, index):
    """
    takes a number (essentially base 4)
    and returns the kmer it corresponds to in alphabetical order
    eg.
    AAAA = 0*1
    CA = 4*4 + 0*1
    GC = 3*4 + 1 * 1
    """
    bases = 'ACGT'
    out = ''
    for k in range(kmax - 1, -1, -1):
        face, index = divmod(index, 4 ** k)
        out += bases[face]
    return out


def yield_kmers(k):
    """
    An iterater to all kmers of length k in alphabetical order
    """
    bases = 'ACGT'
    for kmer in itertools.product(bases, repeat=k):
        yield ''.join(kmer)


def aopen(file, mode='r'):
    if file[-3:] == '.gz':
        return gzip.open(file, mode + 'b')
    else:
        return open(file, mode)


def hamming_N(str1, str2):
    if not len(str1) == len(str2):
        raise(ValueError, 'lengths don\'t match')
    str1 = str1.upper()
    str2 = str2.upper()
    str1 = str1.replace('N', '#')
    return sum(itertools.imap(operator.ne, str1, str2))


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


# from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))


def iterNlines(inFile, N):
    assert N >= 1
    with aopen(inFile) as f:
        lines = [f.readline() for i in range(N)]
        while True:
            yield lines
            lines = [f.readline() for i in range(N)]
            if lines[0] == '':
                break


def save_fig(fig1, path, extensions=['png', 'pdf']):
    for ext in extensions:
        fig1.savefig(path + '.' + ext, transparent=True, dpi = 900)


def simpleaxis(sp):
    sp.spines['top'].set_visible(False)
    sp.spines['right'].set_visible(False)
    sp.get_xaxis().tick_bottom()
    sp.get_yaxis().tick_left()

def close_float_value(a, b, max_percent=1.0):
    if a == 0 and b == 0:
        return True
    if not (a > 0 and b > 0):
        return False
    ratio = float(max(a, b)) / float(min(a, b))
    percent_increase = (ratio - 1.0) * 100.0
    return percent_increase < max_percent


def significantly_enriched(xs, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
    xs = stats.zscore(xs)
    return [x > zthresh for x in xs]


def iter4Lines(inFile):
    return iterNlines(inFile, 4)

def getAllMismatchedSeqs(kmer, mismatchPositions):
    nucs = ['A', 'C', 'G', 'T']
    #generate tuples of allowed nucs at each mismatch position using a recursive algorithm
    allowedNucs = {}
    mismatchPositions = np.array(mismatchPositions)
    assert len(set(mismatchPositions)) == len(mismatchPositions)
    if len(mismatchPositions) == 0:
        yield kmer
    else:
        mismatchNucs = [] + nucs
        #print kmer
        #print mismatchPositions
        #print mismatchPositions[0]
        #print kmer[mismatchPositions[0]]
        mismatchNucs.remove(kmer[mismatchPositions[0]])
        downstreamMismatchSeqs = getAllMismatchedSeqs(kmer[mismatchPositions[0]+1:], mismatchPositions[1:]-(mismatchPositions[0]+1))
        for mismatchNuc in mismatchNucs:
            for downstreamMismatchSeq in downstreamMismatchSeqs:
                returnSeq = kmer[:mismatchPositions[0]] + mismatchNuc +downstreamMismatchSeq
                assert len(returnSeq) == len(kmer)
                yield returnSeq

def getPaddedMismatchedAdjacentKmers(kmerSequence, padding, numMismatches):
    '''
    Yield all sequences of length (len(kmerSequence)+padding )that contain the given kmer, with exactly the given number of mismatches.
    The order yielded is as follows:
        First mismatches are allowed at position 0 to (numMismatches-1)
            For each register:
                Iterate through all possible nucs at mismatch position in alphabetical order
                    Iterate through each nucleotide in padding positions in alphabetical order.
                Shift to next register
            Move most 3' mismatch position down by one, but not past the end of the kmerSequence if end of KmerSequence
            is reached, shift secondmost 3' mismatch 1 nt 3', and reset most 3' mismatch to 1nt 3' of that one
    '''

    # for troubleshooting, want to check that no repeats are generated, so will assert that size of this list and set
    # must be the same
    kmer_set = set()
    kmer_list =[]
    upper_to_combined = {}
    nucs = 'ACGT'
    #initialize mismatchPositions
    #print numMismatches
    if numMismatches == 0:
        for mismatchedKmer in [kmerSequence]:
            for shift in range(padding+1):
                #generate all possible mismatches to the kmer
                for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                    for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                        paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                        if paddedSeq not in kmer_set:
                            kmer_list.append(paddedSeq)
                        kmer_set.add(paddedSeq)
    else:
        mismatchPositionsList = itertools.combinations(range(len(kmerSequence)), numMismatches)
        for mismatchPositions in mismatchPositionsList:
            #print mismatchPositions
            for mismatchedKmer in getAllMismatchedSeqs(kmerSequence, mismatchPositions):
                for shift in range(padding+1):
                    #generate all possible mismatches to the kmer
                    for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                        for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                            paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                            paddedUpper = paddedSeq.upper()
                            if paddedUpper not in kmer_set:
                                kmer_list.append(paddedUpper)
                            kmer_set.add(paddedUpper)

    #print kmer_list
    #print kmer_set
    #print len(kmer_list), len(kmer_set)
    #assert len(kmer_list) == len(kmer_set)

    return kmer_list






