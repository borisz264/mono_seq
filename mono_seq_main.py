import operator
import aColors

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

    def make_plots(self):
        ms_utils.make_dir(self.rdir_path('plots'))
        self.plot_AUG_reads()
        self.plot_AUG_reads(unique_only = True,)
        self.plot_last_AUG_reads()
        self.plot_last_AUG_reads(unique_only = True,)
        self.plot_AUG_reads(which_AUG = 2, unique_only = True)
        self.plot_AUG_reads(which_AUG = 2)


    def make_table_header(self, of):
        """
        takes a file handle and writes a good header for it such that
        each lane is a column.
        """
        of.write('#')
        for lib in self.libs:
            of.write('\t' + lib.get_barcode())
        of.write('\n[%s]' % self.settings.get_property('protein_name'))
        for lib in self.libs:
            of.write('\t%s' % lib.get_conc())
        of.write('\nwashes')
        for lib in self.libs:
            of.write('\t%i' % lib.get_washes())
        of.write('\nT (C)')
        for lib in self.libs:
            of.write('\t%s' % lib.get_temperature())
        of.write('\n')

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
        subprocess.Popen('bowtie2 -q --very-sensitive-local --norc --no-mixed --dovetail --no-discordant -t -x %s -p %d -1 %s -2 %s --un-conc-gz %s -S %s 1>> %s 2>>%s' % (self.settings.get_bowtie_index(), self.threads,
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
        #qc_engine.write_mapping_summary(self.settings.get_overall_mapping_summary())
        #qc_engine.print_library_count_concordances()
        qc_engine.plot_average_read_positions()
        qc_engine.plot_fragment_length_distributions()
        qc_engine.plot_count_distributions()

    def make_counts_table(self, fractional=False):
        """
        write out number of fragments mapping to each TL in each dataset
        :param fractional: if True, replace raw counts with library fraction
        :return:
        """
        #TODO - finish this method!
        if fractional:
            summary_file = os.path.join(
                self.get_rdir(),
                'tables',
                'fractional_counts.txt')
        else:
            summary_file = os.path.join(
                self.get_rdir(),
                'tables',
                'raw_counts.txt')

        header = ['sequence name\t'] + '\t'.join([lib.settings.sample_name for lib in self.libs]) + ['\n']
        summary_file.write(header)
        for sequence_name in self.libs[0].pool_sequence_mappings:

            out_line = '%s\t%s\n' % (sequence_name, '\t'.join(['%f' % lib.pool_sequence_mappings[sequence_name].fragment_count]))

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
    '''
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

    if args.comparisons or args.all_tasks:
        settings.write_to_log('doing comparisons')
        ms_experiment.compare_all_other_experiments()
    '''

main()