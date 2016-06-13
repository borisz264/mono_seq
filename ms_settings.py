import os
import ConfigParser
import simplejson
import itertools
import shutil
import datetime

import ms_utils

class ms_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)
    def get_force_recount(self, count_type):
        return self.settings['force_%s_recount' % count_type]
    def get_settings_file(self):
        return self.settings_file
    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)
    def get_rdir(self):
        ms_utils.make_dir(self.rdir)
        return self.rdir
    def get_wdir(self):
        ms_utils.make_dir(self.wdir)
        return self.wdir
    def get_input_barcode(self):
        return self.settings['library_seq_barcode']

    def iter_lib_settings(self):
        for i in range(len(self.sample_names)):
            yield ms_lib_settings(self,
                                  self.sample_names[i],
                                  self.read1_fastq_gz_file_handles[i], self.read2_fastq_gz_file_handles[i])

    def process_settings(self, settings_file):
        """
        - reads the settings file and converts str to float, list, etc.
        - stores result in self.settings as a dict()
        """
        int_keys = [ 'read1_5p_bases_to_trim', 'read2_5p_bases_to_trim', 'minimum_reads_for_inclusion',
                     'pool_5p_bases_to_trim', 'pool_3p_bases_to_trim', 'min_post_adaptor_length', 'quality_cutoff']
        #float_keys = []
        str_keys = ['read1_suffix', 'read2_suffix', 'read1_3p_adaptor_sequence', 'read2_5p_adaptor_sequence', 'rrna_index', 'genome_index', 'pool_append',
                    'pool_prepend']
        boolean_keys = ['force_remapping', 'force_recount', 'force_index_rebuild', 'force_retrim', 'trim_adaptor']
        list_str_keys = ['fastq_gz_prefixes', 'sample_names', 'mrnp_libraries', 'monosome_libraries', 'total_libraries', 'input_pools', 'matched_set_annotations']
        #list_float_keys = ['concentrations', 'input_rna']
        extant_files = ['pool_fasta',]
        config = ConfigParser.ConfigParser()
        config.read(settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True
        for k in int_keys:
            settings[k] = int(settings[k])
        for k in str_keys:
            settings[k] = settings[k]
        #for k in float_keys:
        #    settings[k] = float(settings[k])
        for k in boolean_keys:
            if not settings[k].lower() in ['true', 'false']:
                raise ValueError(
                  'Boolean value %s must be "true" or "false"' % k)
            settings[k] = settings[k].lower() == 'true'
        #for k in list_float_keys:
        #    settings[k] = map(float, simplejson.loads(settings[k]))
        #for k in list_int_keys:
        #    settings[k] = map(int, simplejson.loads(settings[k]))
        for k in list_str_keys:
            settings[k] = simplejson.loads(settings[k])
        self.fqdir = settings['fastq_dir']
        self.sample_names = settings['sample_names']
        #for paired end reads, there are now 2 fastq files per sample, at least that's how MIT provides the data.
        self.fastq_gz_read1_files = [fastq_gz_prefix+settings['read1_suffix'] for fastq_gz_prefix in
                                      settings['fastq_gz_prefixes']]
        self.fastq_gz_read2_files = [fastq_gz_prefix + settings['read2_suffix'] for fastq_gz_prefix in
                                     settings['fastq_gz_prefixes']]
        self.fastq_gz_files = self.fastq_gz_read1_files + self.fastq_gz_read2_files
        self.read1_fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file)
                                            for fastq_gz_file in self.fastq_gz_read1_files]
        self.read2_fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file)
                                            for fastq_gz_file in self.fastq_gz_read2_files]

        self.fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file) for fastq_gz_file in self.fastq_gz_files]
        for file_handle in self.fastq_gz_file_handles:
            assert ms_utils.file_exists(file_handle)
        for k in extant_files:
            assert ms_utils.file_exists(settings[k])
        self.settings = settings
        self.wdir = settings['working_dir']
        self.rdir = settings['results_dir']
        ms_utils.make_dir(self.rdir)
        shutil.copy(settings_file, self.rdir)

    def check_barcode_lens(self):
        """
        verifies that all the barcodes are the same length
        """
        barcode_lens = set(map(len, self.settings['barcodes']))
        if 1 != len(barcode_lens):
            raise ValueError('all barcodes must be the same length')
        self.barcode_len = barcode_lens.pop()
        self.settings['barcode_len'] = self.barcode_len

    def check_barcodes_are_separated(self):
        """
        makes sure the barcodes are all totally distinguishable
        """
        for b1, b2 in itertools.combinations(self.settings['barcodes'], 2):
            hamming_dist = ms_utils.hamming_distance(b1, b2)
            if hamming_dist < 2:
                raise ValueError('The barcodes supplied are not well '
                  'separated: %s-%s' % (b1, b2))

    def get_bowtie_index(self):
        index = os.path.join(
          self.get_rdir(),
          'bowtie_indices',
          'pool_index')
        return index

    def get_rRNA_bowtie_index(self):
        index = index = self.get_property('rrna_index')
        return index

    def get_genome_bowtie_index(self):
        index = self.get_property('genome_index')
        return index
    def bowtie_index_exists(self):
        return ms_utils.file_exists(self.get_bowtie_index() + '.1.bt2')

    def get_log(self):
        log = os.path.join(
          self.get_rdir(),
          'log.txt')
        return log
    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()
    def get_trimmed_pool_fasta(self):
        log = os.path.join(
          self.get_rdir(),
          'trimmed_pool_seqs.fasta')
        return log
    def get_overall_mapping_summary(self):
        summary_file = os.path.join(
          self.get_rdir(),
          'QC',
          'mapping_summary.txt')
        return summary_file

class ms_lib_settings:
    def __init__(self, experiment_settings, sample_name, read1_fastq_gz_filehandle, read2_fastq_gz_filehandle):
        self.experiment_settings = experiment_settings
        self.sample_name = sample_name
        self.read1_fastq_gz_filehandle = read1_fastq_gz_filehandle
        self.read2_fastq_gz_filehandle = read2_fastq_gz_filehandle

    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_log(self):
        ms_utils.make_dir(os.path.join(self.experiment_settings.get_rdir(), 'logs'))
        log = os.path.join(
          self.experiment_settings.get_rdir(),
          'logs',
          '%(sample_name)s.log' %
           {'sample_name': self.sample_name})
        return log
    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    def get_paired_fastq_gz_files(self):
        return self.read1_fastq_gz_filehandle, self.read2_fastq_gz_filehandle

    def get_collapsed_reads(self):
        collapsed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'collapsed_reads',
          '%(sample_name)s.fasta.gz' %
           {'sample_name': self.sample_name})
        return collapsed_reads

    def get_adaptor_trimmed_reads(self):
        read1_trimmed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'adaptor_removed',
          '%(sample_name)s_1.fastq.gz' %
           {'sample_name': self.sample_name})
        read2_trimmed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'adaptor_removed',
          '%(sample_name)s_2.fastq.gz' %
           {'sample_name': self.sample_name})
        return read1_trimmed_reads, read2_trimmed_reads

    def get_pool_mapping_stats(self):
        pool_mapping_stats = os.path.join(self.experiment_settings.get_rdir(), 'mapping_stats', '%(sample_name)s.pool.txt' % {'sample_name': self.sample_name})
        return pool_mapping_stats

    def get_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)s.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_mapped_reads_sam(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)s.sam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_unmappable_reads_prefix(self):
        unmapped_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'unmapped_reads',
          '%(sample_name)s.unmappable.fastq.gz' %
           {'sample_name': self.sample_name})
        return unmapped_reads

    def get_unmappable_reads(self):
        unmapped_reads_1 = os.path.join(
          self.experiment_settings.get_rdir(),
          'unmapped_reads',
          '%(sample_name)s.unmappable.fastq.1.gz' %
           {'sample_name': self.sample_name})
        unmapped_reads_2 = os.path.join(
          self.experiment_settings.get_rdir(),
          'unmapped_reads',
          '%(sample_name)s.unmappable.fastq.1.gz' %
           {'sample_name': self.sample_name})
        return unmapped_reads_1, unmapped_reads_2


    def get_trimmed_reads(self):
        trimmed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'trimmed_reads',
          '%(sample_name)s.trimmed.fasta.gz' %
           {'sample_name': self.sample_name})
        return trimmed_reads

    def get_sequence_counts(self):
        sequence_counts = os.path.join(
          self.experiment_settings.get_rdir(),
          'sequence_counts',
          '%(sample_name)s.counts.pkl' %
           {'sample_name': self.sample_name})
        return sequence_counts

    def get_overall_contamination_summary(self):
        summary_file = os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.contamination_summary.txt' %
           {'sample_name': self.sample_name})
        return summary_file

    def collapsed_reads_exist(self):
        collapsed_reads = self.get_collapsed_reads()
        return ms_utils.file_exists(collapsed_reads)

    def adaptorless_reads_exist(self):
        adaptorless_reads = self.get_adaptor_trimmed_reads()
        return ms_utils.file_exists(adaptorless_reads[0]) and ms_utils.file_exists(adaptorless_reads[1])

    def trimmed_reads_exist(self):
        trimmed_reads = self.get_trimmed_reads()
        return ms_utils.file_exists(trimmed_reads[0]) and ms_utils.file_exists(trimmed_reads[1])

    def mapped_reads_exist(self):
        mapped_reads = self.get_mapped_reads()
        return ms_utils.file_exists(mapped_reads)

    def sequence_counts_exist(self):
        sequence_counts = self.get_sequence_counts()
        return ms_utils.file_exists(sequence_counts)