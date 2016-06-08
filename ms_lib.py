from collections import defaultdict
import matplotlib.pyplot as plt
import re
import scipy.stats
import subprocess
import os
import cPickle
import ms_utils
import numpy as np
from collections import Counter

import pysam
import ms_utils

class ms_Lib:
    def __init__(self, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir
        self.pool_sequence_mappings = {}
        self.initialize_pool_sequence_mappings()
        self.enrichment_sorted_mappings = None


    def initialize_pool_sequence_mappings(self):
        if self.get_property('force_recount') or not self.lib_settings.sequence_counts_exist():
            print "counting"
            gene_names = []
            trimmed_sequences = ms_utils.convertFastaToDict(self.experiment_settings.get_trimmed_pool_fasta())
            samfile = pysam.Samfile(self.lib_settings.get_mapped_reads(), "rb" )
            for sequence_name in trimmed_sequences:
                gene_name = sequence_name.split('_')[0] #Sequence names are assumed to be of type <gene_name>_TL_description.
                                                        # For example: YDL112W_-41_WT or YDL112W_-41_mut_24-32
                gene_names.append(gene_name)
                self.pool_sequence_mappings[sequence_name] = pool_sequence_mapping(sequence_name,
                                                                                   trimmed_sequences[sequence_name],
                                                                                   samfile)

            samfile.close()
            #self.compute_lib_fractions()
            #gene_counts = Counter(gene_names)
            ms_utils.makePickle(self.pool_sequence_mappings, self.lib_settings.get_sequence_counts())
        else:
            self.pool_sequence_mappings = ms_utils.unPickle(self.lib_settings.get_sequence_counts())


    def get_single_TL_mappings(self, names_only = False):
        single_TL_mappings = set()
        single_TL_names = set()
        for mapping_name in self.pool_sequence_mappings:
            if self.pool_sequence_mappings[mapping_name].is_only_tl:
                single_TL_mappings.add(self.pool_sequence_mappings[mapping_name])
                single_TL_names.add(mapping_name)

        if names_only:
            return single_TL_names
        else:
            return single_TL_mappings

    def compute_lib_fractions(self):
        total_passing_library_reads = float(sum([pool_sequence_mapping.total_passing_reads for pool_sequence_mapping in self.pool_sequence_mappings.values()]))
        for pool_sequence_mapping in self.pool_sequence_mappings.values():
            pool_sequence_mapping.lib_fraction = pool_sequence_mapping.total_passing_reads/total_passing_library_reads

    def calculate_enrichments(self, input_lib):
        for pool_sequence_mapping in self.pool_sequence_mappings.values():
            input_pool_sequence_mapping = input_lib.get_pool_sequence_mapping(pool_sequence_mapping.sequence_name)
            assert input_pool_sequence_mapping != None
            try:
                pool_sequence_mapping.enrichment = pool_sequence_mapping.lib_fraction/input_pool_sequence_mapping.lib_fraction
            except:
                pool_sequence_mapping.enrichment = 0

    def get_pool_sequence_mapping(self, sequence_name):
        if sequence_name in self.pool_sequence_mappings:
            return self.pool_sequence_mappings[sequence_name]
        else:
            return None
    def get_conc(self):
        return self.lib_settings.get_conc()

    def get_washes(self):
        return self.lib_settings.washes

    def get_poly_ic(self):
        return self.lib_settings.poly_ic_conc

    def get_temperature(self):
        return self.lib_settings.temperature

    def get_rna_conc(self):
        return self.lib_settings.input_rna

    def get_sample_name(self):
        return self.lib_settings.sample_name

    def plot_pcr_bias(self):
        collapsed_reads_file = self.lib_settings.get_collapsed_reads()
        read_counts = np.array()
        f = open(collapsed_reads_file)
        for line in f:
            if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
                if line.startswith('>'):#> marks the start of a new sequence
                    num_reads = int(line[1:].strip().split('-')[1])
                    read_counts.append(num_reads)
                else:
                    continue
        f.close()
        read_fractions = read_counts/float(sum(read_counts))
        read_fractions

    def get_counts(self, sequence_name):
        return self.pool_sequence_mappings[sequence_name].total_passing_reads

    def get_mappings_with_minimum_reads(self, minimum_reads, names_only = False):
        passing_mappings = set()
        for mapping in self.pool_sequence_mappings.values():
            if mapping.get_number_rt_stops() >= minimum_reads:
                passing_mappings.add(mapping)

        if names_only:
            return set([passing_mapping.sequence_name for passing_mapping in passing_mappings])
        else:
            return passing_mappings


class pool_sequence_mapping:
    """
    Represents a single sequence from the input pool
    Stores
        The original RNA sequence used in the pool (No adaptor)
        The Trimmed sequence used for mapping
        The positions of all reads mapping to this sequence
        Total number of reads mapping to this sequence
        Fraction of library reads mapping here
        Enrichment relative to input library

    """
    def __init__(self, sequence_name, full_sequence, sam_file):
        self.sequence_name = sequence_name
        self.full_sequence = full_sequence
        self.total_passing_reads = 0
        self.fragment_5p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_3p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_lengths_at_position = defaultdict(list)  # will map position to a list of fragment lengths with 5' ends at that position
        self.fragment_lengths = []
        self.paired_end_mapping_tags = defaultdict(int)
        self.assign_read_ends_from_sam(sam_file)
        self.fragment_count = len(self.fragment_lengths)

    def assign_read_ends_from_sam(self, sam_file):
        all_mapping_reads = sam_file.fetch(reference = self.sequence_name)
        for read in all_mapping_reads:
            # only need to look at the forward mapping read in each pair, since the necessary mate info is there
            if read.is_read1:
                #read1 should be on the forawrd strans since --norc should be specified
                #this alignment should be the primary one. IF this throws erros, I will need to write more logic.
                assert (not read.is_reverse) and (not read.is_secondary)
                pe_mapping_tag = read.get_tag('YT') #this should always return 'CP' for a concordantly-mapped pair
                self.paired_end_mapping_tags[pe_mapping_tag]+=1
                fragment_start = read.reference_start #0-based start of fragment
                fragment_length = read.template_length
                fragment_end = fragment_start + fragment_length
                self.fragment_5p_ends_at_position[fragment_start] += 1
                self.fragment_3p_ends_at_position[fragment_end] += 1
                self.fragment_lengths_at_position[fragment_start].append(fragment_length)
                self.fragment_lengths.append(fragment_length)

    def contains_subsequence(self, subsequence):
        if subsequence in self.full_sequence:
            return True
        else:
            return False

    def positions_of_subsequence(self, subsequence):
        #this regex will NOT return overlapping sequences
        return [m.start() for m in re.finditer(subsequence, self.full_sequence)]

    def fraction_at_position(self, position):
        if position < 0 or position > len(self.full_sequence)-1:
            return None
        else:
            #return self.reads_at_position[position]/float(self.total_passing_reads)
            if self.get_number_rt_stops() == 0:
                return 0
            else:
                return self.fragment_5p_ends_at_position[position] / self.get_number_rt_stops()