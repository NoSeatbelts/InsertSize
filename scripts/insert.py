import array
import operator
import sys

import pysam
import numpy as np
from itertools import chain
try:
    from cytoolz import frequencies as freq
except ImportError:
    try:
        from stl.util import freq
    except ImportError:
        from collections import Counter as freq


def get_pileup_pv(pileup_read):
    if pileup_read.alignment.is_reverse:
        return pileup_read.alignment.opt("PV")[
            pileup_read.alignment.qlen - pileup_read.query_position - 1]
    else:
        return pileup_read.alignment.opt("PV")[pileup_read.query_position]


def get_bc_stack(pileup):
    from collections import defaultdict
    ret = defaultdict(list)
    for i in pileup.pileups:
        bc = i.alignment.seq[i.query_position]
        ret[bc].append(i)
    return ret


def get_expectation_count(bc_stacks, basecall='T', minFM=2):
    """
    bc_stack is a dict with base calls as keys
    and lists of of pileup reads with the same base call.
    """
    pos_exp = sum(1 - 10**(-.1 * get_pileup_pv(pileup_read)) for
                  pileup_read in bc_stacks[basecall])
    # Expected number of false positives given observed other bases
    neg_exp = sum(10**(-.1 * get_pileup_pv(pileup_read)) / 3 for
                  pileup_read in chain.from_iterable(stack for key, stack in
                                                     bc_stacks.items()
                                                     if key != basecall))
    '''
    # Expected number given p values of other base calls
    Could use this, but unsure if it's valid or safe.
    '''
    print "pos_exp", pos_exp, "neg_exp", neg_exp
    if(sum(1 for pileup_read in bc_stacks[basecall] if
           pileup_read.alignment.opt("FM") >= minFM) == 0):
        print "No sufficient FM support. Abort!"
        return 0  # Filter this out
    skepticals = sum(1 for pileup_read in bc_stacks[basecall] if
                     pileup_read.alignment.opt("FM") < minFM)
    if neg_exp >= skepticals:
        return int(pos_exp - neg_exp + 0.5)
    else:
        return int(pos_exp - skepticals + 0.5)


class AAFObj(object):
    def __init__(self, varcount, allcount):
        self.frac = varcount * 1. / allcount
        self.varcount = varcount
        self.allcount = allcount

    def __str__(self):
        return "frac:%f,var:%i,all%i" % (self.frac, self.varcount,
                                         self.allcount)


class AAFStats(object):
    def __init__(self, pileup_column, minFM=2, ref='C', alt='T'):
        all_counts = freq(i.alignment.seq[i.query_position] for
                          i in pileup_column.pileups)
        fm_counts = freq(i.alignment.seq[i.query_position] for
                         i in pileup_column.pileups
                         if i.alignment.opt("FM") >= minFM)
        if alt not in all_counts:
            all_counts[alt] = 0
        if alt not in fm_counts:
            fm_counts[alt] = 0
        self.all = AAFObj(all_counts[alt], sum(all_counts.values()))
        self.fm = AAFObj(fm_counts[alt], sum(fm_counts.values()))
        self.all_discounted = AAFObj(fm_counts[alt], sum(all_counts.values()))
        bc_stacks = get_bc_stack(pileup_column)
        self.all_expected = AAFObj(get_expectation_count(bc_stacks, alt,
                                                         minFM),
                                   sum(all_counts.values()))

    def __str__(self):
        return "%s\t%s\t%s\t%s\t" % (self.all, self.fm, self.all_discounted,
                                     self.all_expected)


def print_aaf_rates(outfile, bamlist, contig, pos, ref='C', alt='T', minFM=2):
    fh = open(outfile, "w")
    fh.write("Bam name\tAll[aaf/aac/all_allele count]\t"
             "FM>=%i[aaf/aac/all_allele_count]\t"
             "DiscountSingles[aaf/aac/all_allele count]\t"
             "ExpectationCount[aaf/aac/all_allele count]\n" % minFM)
    for bam in bamlist:
        a = pysam.AlignmentFile(bam)
        b = a.pileup(contig, pos)
        c = b.next()
        while c.pos < pos:
            c = b.next()
        fh.write(bam + "\t")
        fh.write(str(AAFStats(c, minFM=minFM, ref=ref, alt=alt)))
        fh.write("\n")
    fh.close()


def aaf_ret(bam, contig, pos, ref='C', alt='T', minFM=1):
    """
    Returns the frequency of ref allele 'C' in pileup in a bam.
    """
    a = pysam.AlignmentFile(bam)
    b = a.pileup(contig, pos)
    c = b.next()
    while c.pos < pos:
        c = b.next()
    counts = freq(i.alignment.seq[i.query_position] for i in c.pileups
                  if i.alignment.opt("FM") >= minFM)
    try:
        return counts[alt] / counts[ref]
    except KeyError:
        if ref not in counts:
            return -137.
        return 0.


def aaf_pos(bam, contig, pos, ref='C', alt='T', minFM=1):
    """
    Returns the frequency of ref allele 'C' in pileup in a bam.
    """
    a = pysam.AlignmentFile(bam)
    b = a.pileup(contig, pos)
    c = b.next()
    while c.pos < pos:
        c = b.next()
    counts = freq(i.alignment.seq[i.query_position] for i in c.pileups
                  if i.alignment.opt("FM") >= minFM)
    try:
        return counts[alt] / counts[ref]
    except KeyError:
        if ref not in counts:
            return -137.
        return 0.


def get_tlens(pileup, allele, tlen_cutoff, fm_cutoff):
    if allele == "ALL":
        return np.abs({i.alignment.query_name: abs(
                       i.alignment.template_length) for
                       i in pileup.pileups if
                       i.alignment.is_proper_pair and
                       i.alignment.template_length and
                       not (i.is_del or i.is_refskip) and
                       i.alignment.opt("FM") >= fm_cutoff and
                       abs(i.alignment.template_length) < tlen_cutoff
                       }.values())
    return np.abs({i.alignment.query_name: abs(
                   i.alignment.template_length) for
                   i in pileup.pileups if
                   i.alignment.is_proper_pair and
                   i.alignment.template_length and
                   not (i.is_del or i.is_refskip) and
                   i.alignment.seq[i.query_position] == allele and
                   i.alignment.opt("FM") >= fm_cutoff and
                   abs(i.alignment.template_length) < tlen_cutoff
                   }.values())


class AlleleInfo(object):
    def __init__(self, pileup, allele, tlen_cutoff, fm_cutoff):
        tlens = get_tlens(pileup, allele, tlen_cutoff, 1)
        self.all_mean = np.mean(tlens)
        self.all_std = np.std(tlens)
        self.all_n = tlens.size
        tlens = get_tlens(pileup, allele, tlen_cutoff, fm_cutoff)
        self.fm_mean = np.mean(tlens)
        self.fm_std = np.std(tlens)
        self.fm_n = tlens.size


class AlleleInfoSummary(object):
    def __init__(self, pileup, ref, alt, tlen_cutoff, fm_cutoff):
        self.ref_info = AlleleInfo(pileup, ref, tlen_cutoff, fm_cutoff)
        self.alt_info = AlleleInfo(pileup, alt, tlen_cutoff, fm_cutoff)
        self.all_info = AlleleInfo(pileup, "ALL", tlen_cutoff, fm_cutoff)


class AlleleInsertSummary(object):
    """
    What do we want at the end?
    Mean/std insert size [tumor allele]
    Mean/std insert size [ref allele]
    Mean/std insert size [all alleles]

    And again with minFM 2
    Mean/std insert size [tumor allele]
    Mean/std insert size [ref allele]
    Mean/std insert size [all alleles]
    """
    def __init__(self, path, contig, pos, ref, alt,
                 tlen_cutoff=250, fm_cutoff=2):
        self.handle = pysam.AlignmentFile(path, "r")
        self.contig = contig
        self.pos = pos
        self.ref = ref
        tmp = self.handle.pileup(contig, pos)
        self.pileup = tmp.next()
        while self.pileup.pos < pos:
            self.pileup.pos = tmp.next()
        self.stds = np.array([], dtype=np.double)
        self.means = np.array([], dtype=np.double)
        for allele in (ref, alt):
            # Singletons
            pass

    def get_summary(self):
        ret = array.array('d', [0, 0, 0, 0])
        return ret

'''
Also add a histogram for each.

Counts for inserts at each size for each sample [bp at a time]
[Do the same thing with only the reference]
[Do the same thing with only the tumor]

Pick 3-4 normals and 3-4 tumors.

'''


def get_freqs(path, contig, pos, ref, alt, tlen_cutoff, fm_cutoff):
    it = pysam.AlignmentFile(path).pileup(contig, pos)
    column = it.next()
    while column.pos < pos:
        column = it.next()
    ref_counts = freq(get_tlens(column, ref, tlen_cutoff, fm_cutoff))
    alt_counts = freq(get_tlens(column, alt, tlen_cutoff, fm_cutoff))
    all_counts = freq(get_tlens(column, "ALL", tlen_cutoff, fm_cutoff))
    return {"ref": ref_counts, "alt": alt_counts, "all": all_counts}


def write_freqs(path, freqs, minFM):
    for key, counts in freqs.iteritems():
        handle = open("%s.%s.minFM%i.txt" % (".".join(path.split(".")[:-1]),
                                             key, minFM),
                      "w")
        handle.write("##%s\n" % path)
        handle.write("##%s\n" % key)
        handle.write("#Insert Size\tCount\n")
        [handle.write("%s\t%s\n" % (k, v)) for k, v in
         sorted(counts.iteritems(), key=operator.itemgetter(0))]


def write_all_freqs(path, contig, pos, ref, alt,
                    tlen_cutoff=250, fm_cutoff=2):
    write_freqs(path, get_freqs(path, contig, pos, ref, alt,
                                tlen_cutoff, 1), 1)
    write_freqs(path, get_freqs(path, contig, pos, ref, alt,
                                tlen_cutoff, fm_cutoff), fm_cutoff)


def get_global_insert_freqs(path, minFM=1, max_insert=250):
    ret = {}
    tlen = 0
    for read in pysam.AlignmentFile(path):
        if read.flag & 3980:
            # Read 2, secondary, supplementary, qcfail,
            # optical, unmapped, mate unmapped, duplicate, qc fail
            continue
        if read.flag & 2 and read.opt("FM") >= minFM and read.template_length:
            tlen = abs(read.template_length)
            if tlen < max_insert:
                try:
                    ret[tlen] += 1
                except KeyError:
                    ret[tlen] = 1
    return ret


def write_global_freqs(path, freqs, max_insert):
    with open("%s.%i.global.txt" % (".".join(path.split(".")[:-1]),
                                    max_insert),
              "w") as handle:
        handle.write("##%s\n" % path)
        handle.write("##max_insert:%i\n" % max_insert)
        handle.write("#Insert Size\tCount\n")
        handle.write("\n".join("%s\t%s" % (k, v) for k, v in
                               sorted(freqs.iteritems(),
                                      key=operator.itemgetter(0))))


def write_freqs_for_set(fns, max_insert):
    [write_global_freqs(fn, get_global_insert_freqs(fn, max_insert=max_insert),
                        max_insert) for fn in fns]


def write_global_freqs_for_set(fns, start=150, stop=550, step=50):
    for insert_size in xrange(start, stop, step):
        write_freqs_for_set(fns, insert_size)


def bam_read_generator(bampath, minFM=1):
    for read in pysam.AlignmentFile(bampath):
        if not read.flag & 2816 and read.opt("FM") >= minFM:
            yield read


def get_mean_insert(path, minFM=1, max_insert=250):
    count, total = 0, 0
    for read in pysam.AlignmentFile(path):
        if read.flag & 3980:
            # Read 2, secondary, supplementary, qcfail,
            # optical, unmapped, mate unmapped, duplicate, qc fail
            continue
        if read.flag & 2 and read.opt("FM") >= minFM:
            if read.template_length and read.template_length < max_insert:
                count += 1
                total += abs(read.template_length)
    ret = count * 1. / total
    sys.stderr.write("Mean insert size for %s: %f\n" % (path, ret))
    return ret


def mean_insert_dict(files, minFM=2, max_insert=250):
    return {fn: get_mean_insert(fn, minFM, max_insert) for fn in files}


def parse_samtools_insert_stats(path):
    import collections
    ret = collections.defaultdict(lambda: 0)
    for line in open(path):
        toks = line.split('\t')
        ret[int(toks[0])] = int(toks[1])
    return ret


def build_samtools_insert_dict(fns):
    return {fn.split('/')[-1].split('.')[0]:
            parse_samtools_insert_stats(fn)
            for fn in fns}


def write_samtools_insert_table(fns, outfile):
    ofh = open(outfile, "w")
    insert_data = build_samtools_insert_dict(fns)
    # Make the union of the insert size counts
    key_union = set()
    for d in insert_data.values():
        key_union |= set(d.keys())
    sorted_keys = sorted(key_union)
    del key_union
    ofh.write("#Insert Size")
    sorted_prefixes = sorted(insert_data.keys())
    for prefix in sorted_prefixes:
        ofh.write("\t%s" % prefix)
    ofh.write("\n")
    for key in sorted_keys:
        ofh.write("%i" % key)
        for prefix in sorted_prefixes:
            ofh.write("\t%i" % insert_data[prefix][key])
        ofh.write("\n")
