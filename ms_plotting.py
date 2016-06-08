import ms_utils
import mono_seq_main
import numpy
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
from collections import defaultdict

def plot_AUG_reads(self, which_AUG=1, unique_only=False, min_x=-30, max_x=30, read_cutoff=100):
    # 1 is for the first AUG, 2 the 2nd and so on. Only TLs with enough AUGs are counted
    assert which_AUG > 0
    positions = numpy.array(range(min_x, max_x + 1))
    mappings_passing_cutoff_in_all_libs = self.libs[0].get_mappings_with_minimum_reads(read_cutoff, names_only=True)
    for lib in self.libs[1:]:
        mappings_passing_cutoff_in_all_libs = \
            mappings_passing_cutoff_in_all_libs.intersection(lib.get_mappings_with_minimum_reads(read_cutoff,
                                                                                                 names_only=True))

    if unique_only:
        out_name = os.path.join(
            self.settings.get_rdir(),
            'plots',
            'unique_AUG%d_density.pdf' % (which_AUG))
        mapping_names = mappings_passing_cutoff_in_all_libs.intersection(lib.get_single_TL_mappings(names_only=True))
    else:
        out_name = os.path.join(
            self.settings.get_rdir(),
            'plots',
            'AUG%d_density.pdf' % (which_AUG))
        mapping_names = mappings_passing_cutoff_in_all_libs

    fig = plt.figure(figsize=(8, 8))
    plot = fig.add_subplot(111)
    positions = range(min_x, max_x)
    color_index = 0
    genes_plotted = set()
    for lib in self.libs:
        offset_sum = defaultdict(float)
        offset_counts = defaultdict(int)
        num_genes_counted = 0
        for mapping_name in mapping_names:
            mapping = lib.pool_sequence_mappings[mapping_name]
            AUG_positions = mapping.positions_of_subsequence('ATG')
            if len(AUG_positions) >= which_AUG:
                genes_plotted.add(mapping_name)
                num_genes_counted += 1
                alignment_position = AUG_positions[which_AUG - 1]
                for position in positions:
                    AUG_relative_position = alignment_position - position
                    read_fraction_at_position = mapping.fraction_at_position(AUG_relative_position)
                    if read_fraction_at_position != None:
                        offset_sum[position] += read_fraction_at_position
                        offset_counts[position] += 1
        offset_averages = {}
        for position in positions:
            # print position, offset_sum[position], float(offset_counts[position])
            offset_averages[position] = offset_sum[position] / float(offset_counts[position])
        offset_average_array = [offset_averages[position] for position in positions]
        plot.plot(positions, offset_average_array, color=ms_utils.rainbow[color_index], lw=2,
                  label='%s (%d)' % (lib.get_sample_name(), num_genes_counted))
        color_index += 1
    plot.axvline(16, ls='--')
    plot.axvline(19, ls='--')

    lg = plt.legend(loc=2, prop={'size': 10}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xticks(positions[::3])
    plot.set_xticklabels(positions[::3])
    plot.set_xlabel("position of read 5' end from AUG %d" % (which_AUG))
    plot.set_ylabel("average read fraction")
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()
    print genes_plotted
    for mapping_name in genes_plotted:
        self.plot_single_sequence_read_distributions(mapping_name)


def plot_last_AUG_reads(self, unique_only=False, min_x=-30, max_x=30, read_cutoff=100):
    # 1 is for the first AUG, 2 the 2nd and so on. Only TLs with enough AUGs are counted
    positions = numpy.array(range(min_x, max_x + 1))
    mappings_passing_cutoff_in_all_libs = self.libs[0].get_mappings_with_minimum_reads(read_cutoff, names_only=True)
    for lib in self.libs[1:]:
        mappings_passing_cutoff_in_all_libs = \
            mappings_passing_cutoff_in_all_libs.intersection(lib.get_mappings_with_minimum_reads(read_cutoff,
                                                                                                 names_only=True))

    if unique_only:
        out_name = os.path.join(
            self.settings.get_rdir(),
            'plots',
            'unique_last_AUG_density.pdf')
        mapping_names = mappings_passing_cutoff_in_all_libs.intersection(lib.get_single_TL_mappings(names_only=True))
    else:
        out_name = os.path.join(
            self.settings.get_rdir(),
            'plots',
            'last_AUG_density.pdf')
        mapping_names = mappings_passing_cutoff_in_all_libs

    fig = plt.figure(figsize=(8, 8))
    plot = fig.add_subplot(111)
    positions = range(min_x, max_x)
    color_index = 0
    for lib in self.libs:
        offset_sum = defaultdict(float)
        offset_counts = defaultdict(int)
        num_genes_counted = 0
        for mapping_name in mapping_names:
            mapping = lib.pool_sequence_mappings[mapping_name]
            AUG_positions = mapping.positions_of_subsequence('ATG')
            if len(AUG_positions) >= 1:
                num_genes_counted += 1
                alignment_position = AUG_positions[-1]
                for position in positions:
                    AUG_relative_position = alignment_position - position
                    read_fraction_at_position = mapping.fraction_at_position(AUG_relative_position)
                    if read_fraction_at_position != None:
                        offset_sum[position] += read_fraction_at_position
                        offset_counts[position] += 1
        offset_averages = {}
        for position in positions:
            # print position, offset_sum[position], float(offset_counts[position])
            offset_averages[position] = offset_sum[position] / float(offset_counts[position])
        offset_average_array = [offset_averages[position] for position in positions]
        plot.plot(positions, offset_average_array, color=ms_utils.rainbow[color_index], lw=2,
                  label='%s (%d)' % (lib.get_sample_name(), num_genes_counted))
        color_index += 1
    plot.axvline(16, ls='--')
    plot.axvline(19, ls='--')
    lg = plt.legend(loc=2, prop={'size': 10}, labelspacing=0.2)
    lg.draw_frame(False)
    plot.set_xticks(positions[::3])
    plot.set_xticklabels(positions[::3])
    plot.set_xlabel("position of read 5' end from last AUG")
    plot.set_ylabel("average read fraction")
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()


def plot_single_sequence_read_distributions(self, sequence_name):
    fig = plt.figure(figsize=(8, 8))
    plot = fig.add_subplot(111)
    colorIndex = 0
    for lib in self.libs:
        mapping = lib.pool_sequence_mappings[sequence_name]
        positions = numpy.array(range(0, len(mapping.full_sequence)))
        fractions = [mapping.fraction_at_position(position) for position in positions]
        plot.plot(positions, fractions, color=ms_utils.rainbow[colorIndex], lw=1, label=lib.lib_settings.sample_name)
        colorIndex += 1
    for AUG_pos in mapping.positions_of_subsequence('ATG'):
        plot.axvline(AUG_pos + 16, ls='--')
        plot.axvline(AUG_pos + 19, ls='--')

    plot.set_xticks(positions[::10])
    plot.set_xticklabels(positions[::10])
    plot.set_xlim(-1, len(mapping.full_sequence))
    plot.set_xlabel("position of read 5' end from RNA end (--expected AUG toeprints)")
    plot.set_ylabel("read fraction")
    lg = plt.legend(loc=2, prop={'size': 10}, labelspacing=0.2)
    lg.draw_frame(False)
    out_name = os.path.join(
        self.settings.get_rdir(),
        'plots',
        '%(sequence_name)s.read_positions.pdf' % {'sequence_name': sequence_name})
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()


def write_counts_table()