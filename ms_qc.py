from collections import defaultdict
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats
import subprocess
import os
import cPickle
import ms_utils
import numpy as np
import itertools
import pysam
import math
import gzip

class TPS_qc:
    def __init__(self, tpse, experiment_settings, threads):
        """
        Constructor for Library class
        """
        self.threads = threads
        self.tpse = tpse
        self.experiment_settings = experiment_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir
        ms_utils.make_dir(self.tpse.rdir_path('QC'))

    def write_mapping_summary(self, output_file):

        f = open(output_file, 'w')

        f.write('sample name\ttotal reads\tread1 adaptors\tread2 adaptors\ttoo short\tpass trimming filter\tmapping input\tunaligned pairs\tuniquely aligned pairs\tmultiply aligned pairs\ttotal alignment %\n')
        for lib_settings in self.experiment_settings.iter_lib_settings():
            mapping_input, paired_reads, unaligned_pairs, uniquely_aligned_pairs, multiply_aligned_pairs,\
            overall_alignment_percent = self.parse_paired_end_mapping_stats(lib_settings.get_pool_mapping_stats())
            total_reads, read1_adaptors, read2_adaptors, too_short, passing_filter = self.get_adaptor_trimming_stats(lib_settings.get_log())
            f.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n' % (lib_settings.sample_name, total_reads, read1_adaptors,
                                                      read2_adaptors, too_short, passing_filter, mapping_input,
                                                      unaligned_pairs, uniquely_aligned_pairs,
                                                      multiply_aligned_pairs, overall_alignment_percent))
            f.write('%s percents\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (lib_settings.sample_name,
                                                                               100*total_reads/float(total_reads),
                                                                               100*read1_adaptors/float(total_reads),
                                                                               100*read2_adaptors/float(total_reads),
                                                                               100*too_short/float(total_reads),
                                                                               100*passing_filter/float(total_reads),
                                                                               100*mapping_input/float(total_reads),
                                                                               100*unaligned_pairs/float(total_reads),
                                                                               100*uniquely_aligned_pairs/float(total_reads),
                                                                               100*multiply_aligned_pairs/float(total_reads),
                                                                               100*(uniquely_aligned_pairs+multiply_aligned_pairs)/float(total_reads)))
        f.close()

    def parse_paired_end_mapping_stats(self, alignment_summary_file):
        '''
        example alignment summary:

        random stuff here

        100000 reads; of these:
          100000 (100.00%) were paired; of these:
            37476 (37.48%) aligned concordantly 0 times
            60871 (60.87%) aligned concordantly exactly 1 time
            1653 (1.65%) aligned concordantly >1 times
            ----
            37476 pairs aligned 0 times concordantly or discordantly; of these:
              74952 mates make up the pairs; of these:
                66796 (89.12%) aligned 0 times
                1106 (1.48%) aligned exactly 1 time
                7050 (9.41%) aligned >1 times
        66.60% overall alignment rate

        more stuff here
        '''
        f = open(alignment_summary_file)
        for line in f:
            if line.strip().endswith('reads; of these:'):
                total_reads = int(line.strip().split()[0])
                line=f.next()
                paired_reads = int(line.strip().split()[0])
                line =f.next()
                unaligned_pairs = int(line.strip().split()[0])
                line =f.next()
                uniquely_aligned_pairs = int(line.strip().split()[0])
                line =f.next()
                multiply_aligned_pairs = int(line.strip().split()[0])
                line =f.next()
                overall_alignment_percent = float(line.strip().split()[0][:-1])
        f.close()
        return total_reads, paired_reads, unaligned_pairs, uniquely_aligned_pairs, multiply_aligned_pairs, overall_alignment_percent

    def get_adaptor_trimming_stats(self, log_file):
        '''
        example:
        This is cutadapt 1.10 with Python 2.7.10
        Command line parameters: -a GCTGCACGGTGACGTCTCNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --overlap 5 -u 7 -U 21 -q 10 --trim-n --minimum-length 30 --pair-filter=both -o /Users/boris/Gilbert_Lab/Book_4/4.146/80S_seq_pilot_pe/adaptor_removed/80S_1_1.fastq.gz -p /Users/boris/Gilbert_Lab/Book_4/4.146/80S_seq_pilot_pe/adaptor_removed/80S_1_2.fastq.gz /Users/boris/Gilbert_Lab/Book_4/4.146/truncated_pe_data/160519Gil_D16-4658_1_sequence.fastq.gz /Users/boris/Gilbert_Lab/Book_4/4.146/truncated_pe_data/160519Gil_D16-4658_2_sequence.fastq.gz
        Trimming 1 adapter with at most 10.0% errors in paired-end mode ...
        Finished in 10.02 s (40 us/read; 1.50 M reads/minute).

        === Summary ===

        Total read pairs processed:            250,000
          Read 1 with adapter:                  60,994 (24.4%)
          Read 2 with adapter:                       0 (0.0%)
        Pairs that were too short:              60,994 (24.4%)
        Pairs written (passing filters):       189,006 (75.6%)

        Total basepairs processed:    20,000,000 bp
          Read 1:    10,000,000 bp
          Read 2:    10,000,000 bp
        Quality-trimmed:                   3,840 bp (0.0%)
          Read 1:             7 bp
          Read 2:         3,833 bp
        Total written (filtered):      9,825,137 bp (49.1%)
          Read 1:     6,237,184 bp
          Read 2:     3,587,953 bp

        === First read: Adapter 1 ===

        Sequence: GCTGCACGGTGACGTCTCNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC; Type: regular 3'; Length: 54; Trimmed: 60994 times.

        No. of allowed errors:
        0-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-39 bp: 3; 40-49 bp: 4; 50-54 bp: 5

        Bases preceding removed adapters:
          A: 0.1%
          C: 0.2%
          G: 0.1%
          T: 0.2%
          none/other: 99.3%

        Overview of removed sequences
        length	count	expect	max.err	error counts
        5	21	244.1	0	21
        6	6	61.0	0	6
        7	18	15.3	0	18
        8	10	3.8	0	10
        9	11	1.0	0	11
        10	5	0.2	1	5
        11	8	0.1	1	7 1
        12	12	0.0	1	12
        13	8	0.0	1	8
        14	9	0.0	1	9
        15	9	0.0	1	8 1
        16	6	0.0	1	6
        17	9	0.0	1	8 1
        18	11	0.0	1	9 1 1
        19	13	0.0	1	13
        20	7	0.0	2	7
        21	15	0.0	2	14 0 1
        22	19	0.0	2	16 3
        23	10	0.0	2	8 1 1
        24	18	0.0	2	15 2 1
        25	23	0.0	2	18 2 3
        26	23	0.0	2	22 1
        27	32	0.0	2	29 3
        28	14	0.0	2	14
        29	21	0.0	2	20 1
        30	14	0.0	3	9 3 2
        31	31	0.0	3	12 1 16 2
        32	23	0.0	3	21 1 1
        33	60588	0.0	3	102 52134 7310 1042
        '''
        f = open(log_file)
        for line in f:
            if line.strip().startswith('Total read pairs processed:'):
                total_reads = int(line.strip().split()[-1].replace(',', ''))
                line = f.next()
                read1_adaptors = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                read2_adaptors = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                too_short = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                passing_filter = int(line.strip().split()[-2].replace(',', ''))
        f.close()
        return total_reads, read1_adaptors, read2_adaptors, too_short, passing_filter


    def get_collapsed_read_fractions(self, lib_settings):
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC','collapsed_fracs',
          '%(sample_name)s.collapsed_read_fractions.pkl' % {'sample_name': lib_settings.sample_name})
        if not ms_utils.file_exists(out_name) and not self.experiment_settings.get_property('force_recollapse'):
            collapsed_reads_file = lib_settings.get_collapsed_reads()
            read_counts = []
            f = gzip.open(collapsed_reads_file)
            for line in f:
                if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
                    if line.startswith('>'):#> marks the start of a new sequence
                        num_reads = int(line[1:].strip().split('-')[1])
                        read_counts.append(num_reads)
                    else:
                        continue
            f.close()
            read_fractions = np.array(read_counts)/float(sum(read_counts))
            bzUtils.makePickle(read_fractions, out_name)
        else:
            read_fractions = bzUtils.unPickle(out_name)

        return (lib_settings.sample_name, read_fractions)

    def get_library_enrichment_correlation(self, lib1, lib2):
        lib1_enrichments = []
        lib2_enrichments = []
        for sequence in lib1.pool_sequence_mappings:
            lib1_enrichments.append(lib1.pool_sequence_mappings[sequence].enrichment)
            lib2_enrichments.append(lib2.pool_sequence_mappings[sequence].enrichment)
        spearmanR, spearmanP = stats.spearmanr(lib1_enrichments, lib2_enrichments)
        pearsonR, pearsonP = stats.pearsonr(lib1_enrichments, lib2_enrichments)
        return pearsonR, spearmanR, pearsonP, spearmanP

    def get_library_count_correlation(self, lib1, lib2):
        lib1_counts = []
        lib2_counts = []
        for sequence in lib1.pool_sequence_mappings:
            lib1_counts.append(lib1.pool_sequence_mappings[sequence].total_passing_reads)
            lib2_counts.append(lib2.pool_sequence_mappings[sequence].total_passing_reads)
        spearmanR, spearmanP = stats.spearmanr(lib1_counts, lib2_counts)
        pearsonR, pearsonP = stats.pearsonr(lib1_counts, lib2_counts)
        return pearsonR, spearmanR, pearsonP, spearmanP

    def get_library_count_distribution(self, lib):
        return [lib.pool_sequence_mappings[sequence].total_passing_reads for sequence in lib.pool_sequence_mappings]

    def print_library_count_concordances(self):
        out_name =  os.path.join(self.experiment_settings.get_rdir(), 'QC',
          'count_concordances.txt')
        f = open(out_name, 'w')
        header = 'sample1\tsample2\tpearson r\t pearson p\t spearman r\t spearman p\n'
        f.write(header)
        for libi, libj in itertools.combinations(self.tpse.libs, 2):
            pearsonR, spearmanR, pearsonP, spearmanP = self.get_library_count_correlation(libi, libj)
            line = '%s\t%s\t%f\t%f\t%f\t%f\n' % (libi.get_sample_name(), libj.get_sample_name(),
                                                         pearsonR, pearsonP, spearmanR, spearmanP)
            f.write(line)
        f.close()


    def plot_average_read_positions(self):
        for lib in self.tpse.libs:
            self.plot_average_read_positions_one_lib(lib)

    def plot_average_read_positions_one_lib(self, lib, min_x = 0, max_x = 150):
        positions = np.array(range(min_x, max_x+1))
        averages = [np.average([pool_sequence_mapping.fraction_at_position(position) for pool_sequence_mapping in lib.pool_sequence_mappings.values() if pool_sequence_mapping.total_passing_reads>0]) for position in positions]

        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        plot.bar(positions , averages,color=bzUtils.rainbow[0], lw=0)
        plot.set_xticks(positions[::10]+0.5)
        plot.set_xticklabels(positions[::10])
        plot.set_xlabel("position of read 5' end from RNA end")
        plot.set_ylabel("average read fraction")
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.read_positions.pdf' % {'sample_name': lib.get_sample_name ()})
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()


    def plot_count_distributions(self):
        num_libs = len(self.tpse.libs)
        fig = plt.figure(figsize=(16,16))
        plot_index = 1
        cutoff = 100
        hbins = np.arange(0, 400, 10)
        hbins = np.append(hbins, 10000000)
        for lib in self.tpse.libs:
            plot = fig.add_subplot(math.sqrt(bzUtils.next_square_number(num_libs)), math.sqrt(bzUtils.next_square_number(num_libs)), plot_index)
            sample_name = lib.lib_settings.sample_name
            dist = self.get_library_count_distribution(lib)
            plot.hist(dist, bins = hbins, color=bzUtils.skyBlue, histtype='stepfilled', edgecolor = None, lw = 0)
            plot.set_xlabel("# reads", fontsize = 10)
            plot.set_ylabel("# genes (%d have >= %d reads)" % (bzUtils.number_passing_cutoff(dist, cutoff), cutoff), fontsize = 10)
            plot.set_xlim(0, 400)
            #plot.set_ylim(0,1)
            plot.axvline(cutoff, ls = 'dashed')
            plot.set_title(sample_name, fontsize = 8)
            plot_index += 1
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15, wspace=0.4, hspace=0.6)
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          'count_distributions.pdf')
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

        """
        def plot_insert_size_distributions(self):
            #plot distribution of insert sizes from cutadapt output
            TODO - need to parse log file to get this info
            num_libs = len(self.tpse.libs)
            fig = plt.figure(figsize=(16,16))
            plot_index = 1
            cutoff = 100
            hbins = np.arange(0, 51, 1)
            for lib in self.tpse.libs:
                plot = fig.add_subplot(math.sqrt(bzUtils.next_square_number(num_libs)), math.sqrt(bzUtils.next_square_number(num_libs)), plot_index)
                sample_name = lib.lib_settings.sample_name
                dist = self.get_insert_sizes(lib)
                plot.hist(dist, bins = hbins, color=bzUtils.skyBlue, histtype='stepfilled', edgecolor = None, lw = 0)
                plot.set_xlabel("insert size", fontsize = 10)
                plot.set_ylabel("fraction of reads" % (bzUtils.number_passing_cutoff(dist, cutoff), cutoff), fontsize = 10)
                plot.set_xlim(0, 400)
                #plot.set_ylim(0,1)
                plot.axvline(cutoff, ls = 'dashed')
                plot.set_title(sample_name, fontsize = 8)
                plot_index += 1
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15, wspace=0.4, hspace=0.6)
            out_name =  os.path.join(
              self.experiment_settings.get_rdir(),
              'QC',
              'insetrt_size_distributions.pdf')
            plt.savefig(out_name, transparent='True', format='pdf')
            plt.clf()
        """