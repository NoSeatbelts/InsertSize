import array
import operator
import sys
from subprocess import check_output

import pysam
import numpy as np
try:
    from cytoolz import frequencies as freq
except ImportError:
    try:
        from stl.util import freq
    except ImportError:
        from collections import Counter as freq


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


def get_gene_depths(bampath, bedpath, minMQ=0):
    gene_footprints = {}
    for line in open(bedpath, "r"):
        if len(line) < 20:
            continue
        contig, start, stop, name = line.strip().split('\t')
        try:
            gene_footprints[name.split("_")[0]] += int(stop) - int(start)
        except KeyError:
            gene_footprints[name.split("_")[0]] = int(stop) - int(start)
    gene_depths = {}
    bedcov = check_output("samtools bedcov -Q%i %s %s" % \
                          (minMQ, bedpath, bampath),
                          shell=True)
    for line in bedcov.split('\n'):
        if len(line) < 1:
            continue
        print line
        toks = line.strip().split('\t')
        contig, start, stop, name, counts = toks[:5]
        try:
            gene_depths[name.split("_")[0]] += int(counts)
        except KeyError:
            gene_depths[name.split("_")[0]] = int(counts)
    ret = {}
    for key in gene_footprints.iterkeys():
        ret[key] = gene_depths[key] * 1. / gene_footprints[key]
    return {key: value * 1. / gene_footprints[key] for
            key, value in gene_depths.iteritems()}


def write_gene_depths(bampath, bedpath):
    with open(bampath + ".depth_by_gene.txt", "w") as f:
        f.write("#Filename\tminMQ0\tminMQ1\n")
        minMQ0 = get_gene_depths(bampath, bedpath, 0)
        minMQ1 = get_gene_depths(bampath, bedpath, 1)
        f.write("\n".join("%s\t%s\t%s\t" % (key, value, minMQ1[key]) for
                          key, value in minMQ0.iteritems()))

def batch_gene_depths(bedpath, l=[]):
    if len(l) == 0:
        l = check_output("ls *.rsq.bam", shell=True).split('\n')
    for fn in l:
        write_gene_depths(fn, bedpath)
