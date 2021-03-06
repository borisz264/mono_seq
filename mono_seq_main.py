__author__ = 'boris zinshteyn'
"""
Intended for processing of 80s monosome-seq data from defined RNA pools
Based on Alex Robertson's original RBNS pipeline, available on github
"""
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import subprocess

import ms_settings
import ms_utils
import ms_lib
import ms_qc
import ms_plotting
from collections import defaultdict
import numpy as np
import scipy.stats as stats

class mse:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.remove_adaptor()
        self.trim_reference_pool_fasta()
        self.build_bowtie_index()
        self.map_reads()
        self.initialize_libs()

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        ms_utils.make_dir(self.rdir_path('sequence_counts'))
        self.libs = []

        ms_utils.parmap(lambda lib_settings: ms_lib.initialize_pool_sequence_mappings(self.settings, lib_settings),
                        self.settings.iter_lib_settings(), nprocs=self.threads)
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')
        self.monosome_libs = [self.find_lib_by_sample_name(sample_name) for
                              sample_name in self.settings.get_property('monosome_libraries')]
        self.mrnp_libs = [self.find_lib_by_sample_name(sample_name) for
                              sample_name in self.settings.get_property('mrnp_libraries')]
        self.total_libs = [self.find_lib_by_sample_name(sample_name) for
                          sample_name in self.settings.get_property('total_libraries')]
        self.input_libs = [self.find_lib_by_sample_name(sample_name) for
                           sample_name in self.settings.get_property('total_libraries')]


    def find_lib_by_sample_name(self, sample_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == sample_name:
                return lib
        assert False #if this triggers, your settings file is broken.


    def initialize_lib(self, lib_settings):
        lib = ms_lib.ms_Lib(self.settings, lib_settings)
        self.libs.append(lib)

    def needs_calculation(self, lib_settings, count_type, k):
        if self.settings.get_force_recount(count_type):
            return True
        return not lib_settings.counts_exist(count_type, k)


    def make_tables(self):
        ms_utils.make_dir(self.rdir_path('tables'))
        self.make_counts_table()
        self.make_counts_table(fractional=True)
        self.make_monosome_recruitment_table()
        self.write_sequence_subset(0, read_cutoff=128)
        self.write_sequence_subset(0.8, read_cutoff=128)
        self.write_sequence_subset(0.7, read_cutoff=128)
        for anno_filename in self.settings.get_property('matched_set_annotations'):
            self.make_matched_recruitment_change_table(anno_filename,
                                                       read_cutoff=self.settings.get_property('comparison_read_cutoff'))

    def make_plots(self):
        ms_utils.make_dir(self.rdir_path('plots'))

        #ms_plotting.all_library_rpm_scatter(self)
        #ms_plotting.monosome_over_mrnp_reproducibility(self)
        #ms_plotting.monosome_over_total_reproducibility(self)
        #ms_plotting.monosome_over_mrnp_plus_monosome_reproducibility(self)
        for anno_filename in self.settings.get_property('matched_set_annotations'):
            ms_plotting.plot_recruitment_violins(self, anno_filename,
                                                 read_cutoff=self.settings.get_property('comparison_read_cutoff'))
            '''
            ms_plotting.recruitment_change_rank_value_plot_static(self, anno_filename,
                                                 read_cutoff=self.settings.get_property('comparison_read_cutoff'))
            ms_plotting.reverse_recruitment_change_rank_value_plot_static(self, anno_filename,
                                                 read_cutoff=self.settings.get_property('comparison_read_cutoff'))

            if self.settings.get_property('make_interactive_plots'):
                ms_plotting.recruitment_change_rank_value_plot_interactive(self, anno_filename,
                                                                           read_cutoff=self.settings.get_property('comparison_read_cutoff'))

                ms_plotting.recruitment_fold_change_rank_value_plot_interactive(self, anno_filename,
                                                                               read_cutoff=self.settings.get_property(
                                                                                   'comparison_read_cutoff'))
            '''
    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            ms_utils.make_dir(self.rdir_path('adaptor_removed'))
            ms_utils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting), self.settings.iter_lib_settings(), nprocs = self.threads)

    def remove_adaptor_one_lib(self, lib_settings):
        lib_settings.write_to_log('adaptor trimming')
        """
        -a specifies the 3' adaptor to trim from the forawrd read (read1)
        -G specifies the 5' adaptor to trim from the reverse read (read2)
        -o is the read1 output file
        -p is the read2 output file
        """
        if not self.settings.get_property('read2_5p_adaptor_sequence').strip()=='':
            command_to_run = 'cutadapt -a %s -G %s --overlap 5 -u %d -U %d -q %d --trim-n --minimum-length %d --pair-filter=both -o %s -p %s %s %s 1>>%s 2>>%s' % (
                self.settings.get_property('read1_3p_adaptor_sequence'), self.settings.get_property('read2_5p_adaptor_sequence'),
                self.settings.get_property('read1_5p_bases_to_trim'), self.settings.get_property('read2_5p_bases_to_trim'),
                self.settings.get_property('quality_cutoff'), self.settings.get_property('min_post_adaptor_length'),
                lib_settings.get_adaptor_trimmed_reads()[0], lib_settings.get_adaptor_trimmed_reads()[1],
                lib_settings.get_paired_fastq_gz_files()[0], lib_settings.get_paired_fastq_gz_files()[1],
                lib_settings.get_log(), lib_settings.get_log())
        else:
            command_to_run = 'cutadapt -a %s --overlap 5 -u %d -U %d -q %d --trim-n --minimum-length %d --pair-filter=both -o %s -p %s %s %s 1>>%s 2>>%s' % (
                self.settings.get_property('read1_3p_adaptor_sequence'),
                self.settings.get_property('read1_5p_bases_to_trim'), self.settings.get_property('read2_5p_bases_to_trim'),
                self.settings.get_property('quality_cutoff'), self.settings.get_property('min_post_adaptor_length'),
                lib_settings.get_adaptor_trimmed_reads()[0], lib_settings.get_adaptor_trimmed_reads()[1],
                lib_settings.get_paired_fastq_gz_files()[0], lib_settings.get_paired_fastq_gz_files()[1],
                lib_settings.get_log(), lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def build_bowtie_index(self):
        """
        builds a bowtie 2 index from the input fasta file
        recommend including barcode+PCR sequences just in case of some no-insert amplicons
        """
        self.settings.write_to_log('building bowtie index')
        if self.settings.get_property('force_index_rebuild') or not self.settings.bowtie_index_exists():
            ms_utils.make_dir(self.rdir_path('bowtie_indices'))
            index = self.settings.get_bowtie_index()
            subprocess.Popen('bowtie2-build -f --offrate 0 %s %s 1>>%s 2>>%s' % (self.settings.get_trimmed_pool_fasta(),
                                                                      self.settings.get_bowtie_index(), self.settings.get_log()+'.bwt',
                                                                      self.settings.get_log()+'.bwt'), shell=True).wait()
        self.settings.write_to_log('building bowtie index complete')

    def trim_reference_pool_fasta(self):
        '''
        Trims the reference sequences to the length of the trimmed reads + a buffer
        '''
        trim_5p = self.settings.get_property('pool_5p_bases_to_trim') #nucleotides to cut from 5' end
        trim_3p = self.settings.get_property('pool_3p_bases_to_trim') #nucleotides to cut from 3' end
        f = open(self.settings.get_property('pool_fasta'))
        g = open(self.settings.get_trimmed_pool_fasta(), 'w')
        for line in f:
            if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
                if line.startswith('>'):#> marks the start of a new sequence
                    g.write(line)
                else:
                    g.write(self.settings.get_property('pool_prepend')+line.strip()[trim_5p:len(line.strip())-trim_3p]+self.settings.get_property('pool_append')+'\n')
        f.close()
        g.close()

    def map_reads(self):
        """
        map all reads using bowtie
        :return:
        """
        self.settings.write_to_log('mapping reads')
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.mapped_reads_exist():
                    break
            else:
                return
        ms_utils.make_dir(self.rdir_path('mapped_reads'))
        ms_utils.make_dir(self.rdir_path('mapping_stats'))
        ms_utils.make_dir(self.rdir_path('unmapped_reads'))

        ms_utils.parmap(lambda lib_setting: self.map_one_library(lib_setting), self.settings.iter_lib_settings(),
                       nprocs = self.threads)
        self.settings.write_to_log( 'finished mapping reads')

    def map_one_library(self, lib_settings):
        lib_settings.write_to_log('mapping_reads')
        subprocess.Popen('bowtie2 -q --very-sensitive-local --norc --no-mixed --no-overlap --no-discordant -t -x %s -p %d -1 %s -2 %s --un-conc-gz %s -S %s 1>> %s 2>>%s' % (self.settings.get_bowtie_index(), self.threads,
                                                                                                   lib_settings.get_adaptor_trimmed_reads()[0], lib_settings.get_adaptor_trimmed_reads()[1], lib_settings.get_unmappable_reads_prefix(), lib_settings.get_mapped_reads_sam(),
                                                                                                                      lib_settings.get_log(), lib_settings.get_pool_mapping_stats()), shell=True).wait()
        #subprocess.Popen('samtools view -b -h -o %s %s 1>> %s 2>> %s' % (lib_settings.get_mapped_reads(), lib_settings.get_mapped_reads_sam(), lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        #also, sort bam file, and make an index

        #samtools view -uS myfile.sam | samtools sort - myfile.sorted
        subprocess.Popen('samtools view -uS %s | samtools sort - %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_mapped_reads_sam(), lib_settings.get_mapped_reads_sam(),
                                                                          lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()


        #subprocess.Popen('samtools sort %s %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_mapped_reads_sam(), lib_settings.get_mapped_reads_sam(),
        #                                                                  lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        subprocess.Popen('mv %s.temp_sorted.bam %s' % (lib_settings.get_mapped_reads_sam(),
                                                                          lib_settings.get_mapped_reads()), shell = True).wait()
        subprocess.Popen('samtools index %s' % (lib_settings.get_mapped_reads()), shell = True).wait()
        subprocess.Popen('rm %s' % (lib_settings.get_mapped_reads_sam()), shell = True).wait()
        lib_settings.write_to_log('mapping_reads done')
    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        out_path = self.rdir_path(*args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return ms_utils.aopen(out_path, 'w')

    def perform_qc(self):
        qc_engine = ms_qc.ms_qc(self, self.settings, self.threads)
        qc_engine.write_mapping_summary(self.settings.get_overall_mapping_summary())
        qc_engine.print_library_count_concordances()
        qc_engine.plot_average_read_positions()
        qc_engine.plot_fragment_length_distributions()
        qc_engine.plot_count_distributions()
        qc_engine.read_cutoff_choice_plot()

    def make_counts_table(self, fractional=False):
        """
        write out number of fragments mapping to each TL in each dataset
        :param fractional: if True, replace raw counts with library fraction in reads per million
        :return:
        """
        if fractional:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'rpm.txt'), 'w')
        else:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'raw_counts.txt'), 'w')

        header = 'sequence name\t' + '\t'.join([lib.lib_settings.sample_name for lib in self.libs]) + '\n'
        summary_file.write(header)
        if fractional:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' % ((10**6)*lib.pool_sequence_mappings[sequence_name].fragment_count/float(lib.total_mapped_fragments)) for lib in self.libs]))
                summary_file.write(out_line)
        else:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' %
                                                    lib.pool_sequence_mappings[sequence_name].fragment_count
                                                    for lib in self.libs]))
                summary_file.write(out_line)
        summary_file.close()


    def make_monosome_recruitment_table(self, read_cutoff=128):
        """
        write out 80S recruitment metric for each TL in each replicate
        :param read_cutoff: require this many read between mRNP and monosome to include this TL.
        :return:
        """
        output_file = open(os.path.join(
            self.rdir_path('tables'),
            'monosome_recruitment.txt'), 'w')
        trimmed_sequences = ms_utils.convertFastaToDict(self.settings.get_trimmed_pool_fasta())
        header = 'sequence name\tsequence\t' + '\t'.join(['%s/(%s+%s)' % (self.monosome_libs[i].lib_settings.sample_name,
                                                                self.monosome_libs[i].lib_settings.sample_name,
                                                                self.mrnp_libs[i].lib_settings.sample_name)
                                                for i in range(len(self.monosome_libs))]) + '\n'
        output_file.write(header)
        for sequence_name in self.monosome_libs[0].pool_sequence_mappings:
            out_line = '%s\t%s\t%s\n' % (sequence_name, trimmed_sequences[sequence_name],
                                     '\t'.join(['%f' %
                                                (self.monosome_libs[i].get_rpm(sequence_name)/
                                                (self.monosome_libs[i].get_rpm(sequence_name)+
                                                 self.mrnp_libs[i].get_rpm(sequence_name)))
                                                if (self.monosome_libs[i].get_counts(sequence_name) +
                                                 self.mrnp_libs[i].get_counts(sequence_name)) >= read_cutoff else ''
                                                for i in range(len(self.monosome_libs)) ]))
            output_file.write(out_line)
        output_file.close()


    def write_sequence_subset(self, recruitment_cutoff, read_cutoff=128, as_RNA=True):
        """
        write out fasta of all sequences that pass a certain recruitment cutoff in all libraries
        :return:
        """
        output_file = open(os.path.join(
            self.rdir_path('tables'),
            'recruitment_above_%f.fasta' % recruitment_cutoff), 'w')
        trimmed_sequences = ms_utils.convertFastaToDict(self.settings.get_trimmed_pool_fasta())

        for sequence_name in self.monosome_libs[0].pool_sequence_mappings:
            rec_scores = [(self.monosome_libs[i].get_rpm(sequence_name) / (self.monosome_libs[i].get_rpm(sequence_name) +
                                                                           self.mrnp_libs[i].get_rpm(sequence_name)))
                          for i in range(len(self.monosome_libs)) if (self.monosome_libs[i].get_counts(sequence_name) +
                              self.mrnp_libs[i].get_counts(sequence_name)) >= read_cutoff]
            if (len(rec_scores) == len(self.monosome_libs)):
                average_score = np.average(rec_scores)
                if average_score >=recruitment_cutoff:
                    output_file.write('>%s_rec_%f\n' % (sequence_name, average_score))
                    seq = trimmed_sequences[sequence_name]
                    if as_RNA:
                        seq = ms_utils.rna(seq)
                    output_file.write('%s\n' % (seq))
        output_file.close()

    def make_matched_recruitment_change_table(self, annotation_file, read_cutoff=128):
        """
        write out number of fragments mapping to each TL in each dataset
        :param read_cutoff: require this many read between mRNP and monosome to include this TL.
        :return:
        """
        set_name1, set_name2, matched_set = self.parse_matched_set_annotation(annotation_file)

        output_file = open(os.path.join(
            self.rdir_path('tables'),
            '%s_%s_matched_monosome_recruitment_change.txt' % (set_name1, set_name2)), 'w')

        header = '%s\t%s\t' % (set_name1, set_name2) + '\t'.join(['%s %s-%s recruitment score' % (self.monosome_libs[i].lib_settings.sample_name,
                                                                                               set_name1, set_name2)
                                                                for i in range(len(self.monosome_libs))]) + '\t'+\
                 '\t'.join(['%s %s/%s recruitment score' % (self.monosome_libs[i].lib_settings.sample_name,
                                                                                               set_name1, set_name2)
                                                                for i in range(len(self.monosome_libs))])+'\tttest p\n'
        output_file.write(header)
        for matched_pool_seqs in matched_set:
            set1_scores = []
            set2_scores = []
            for i in range(len(self.monosome_libs)):
                set_1_counts = self.monosome_libs[i].get_counts(matched_pool_seqs[0])\
                               + self.mrnp_libs[i].get_counts(matched_pool_seqs[0])
                set_2_counts = self.monosome_libs[i].get_counts(matched_pool_seqs[1]) \
                               + self.mrnp_libs[i].get_counts(matched_pool_seqs[1])
                # include only comparisons where the average number of reads is high enough
                if set_1_counts >= read_cutoff and set_2_counts >= read_cutoff:
                    set1_score = self.monosome_libs[i].get_rpm(matched_pool_seqs[0]) / \
                                 (self.monosome_libs[i].get_rpm(matched_pool_seqs[0]) +
                                  self.mrnp_libs[i].get_rpm(matched_pool_seqs[0]))
                    set2_score = self.monosome_libs[i].get_rpm(matched_pool_seqs[1]) / \
                                 (self.monosome_libs[i].get_rpm(matched_pool_seqs[1]) +
                                  self.mrnp_libs[i].get_rpm(matched_pool_seqs[1]))
                else:
                    set1_score = float('nan')
                    set2_score = float('nan')
                set1_scores.append(set1_score)
                set2_scores.append(set2_score)
            recruitment_changes = np.array(set1_scores)-np.array(set2_scores)
            recruitment_fold_changes = np.array(set1_scores)/np.array(set2_scores)
            scores_1_filtered, scores_2_filtered = ms_utils.filter_x_y_pairs(set1_scores, set2_scores)
            if len(scores_1_filtered)>0 and len(scores_2_filtered)>0:
                t, p = stats.ttest_ind(scores_1_filtered, scores_2_filtered)
            else:
                p = float('nan')
            out_line = '%s\t%s\t%s\t%s\t%f\n' % (matched_pool_seqs[0], matched_pool_seqs[1],
                                             '\t'.join(['%f' % score_change for score_change in recruitment_changes]),
                                             '\t'.join(['%f' % score_change for score_change in recruitment_fold_changes]),
                                             p)
            output_file.write(out_line)
        output_file.close()

    def parse_matched_set_annotation(self, filename):
        matched_set = set()# a set of tuples matching a sequence name to one matched to it. sequences cana ppear multiple times, but the pairs ought to be unique

        f = open(filename)
        lines = f.readlines()
        header = lines[0]
        set_name1, set_name2 = header.strip().split('\t')
        for line in lines[1:]:
            seq_name1, seq_name2 = line.strip().split('\t')
            matched_set.add((seq_name1, seq_name2))
        f.close()
        return set_name1, set_name2, matched_set

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--perform-qc",
                        help="performs quality control analysis.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--comparisons",
                        help="Does comparisons to other experiments",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables, folding and comparisons",
                        action='store_true')
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 8)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = ms_settings.ms_settings(args.settings_file)
    ms_experiment = mse(settings, args.threads)
    print 'mse ready'
    if args.perform_qc or args.all_tasks:
        print 'QC'
        settings.write_to_log('performing QC')
        ms_experiment.perform_qc()
        settings.write_to_log('done performing QC')
    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        ms_experiment.make_tables()
        settings.write_to_log('done making tables')

    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        ms_experiment.make_plots()
        settings.write_to_log('done making plots')
    '''
    if args.comparisons or args.all_tasks:
        settings.write_to_log('doing comparisons')
        ms_experiment.compare_all_other_experiments()
    '''

main()