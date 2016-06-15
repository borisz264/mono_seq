import ms_utils
import numpy
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import scipy.stats as stats
import operator

def all_library_rpm_scatter(mse):
    output_file = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        'all_scatter_plots.pdf')
    num_libs = len(mse.libs)
    num_plots_wide = num_libs-1
    num_plots_high = num_libs-1
    fig = plt.figure(figsize=(24,24))

    for i in range(len(mse.libs)):
        for j in range(i+1, len(mse.libs)):
            plot_index = (j-1)*(num_plots_wide)+(i+1)
            plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
            if j == num_plots_high:
                plot.set_xlabel("%s RPM" % (mse.libs[i].lib_settings.sample_name))
            if i == 0:
                plot.set_ylabel("%s RPM" % (mse.libs[j].lib_settings.sample_name))
            plot.set_xscale('symlog', linthreshx=0.1)
            plot.set_yscale('symlog', linthreshy=0.1)
            x = mse.libs[i].name_sorted_rpms()
            y = mse.libs[j].name_sorted_rpms()
            plot.scatter(x, y, color=ms_utils.black, s=3)
            plot.plot(numpy.arange(0,1000000,1), numpy.arange(0,1000000,1), color=ms_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)

            plot.annotate('rho=%.3f' % (rho), xy=(0, 0.8), xytext=(0, 0.8), textcoords='axes fraction')
            plot.set_xlim(0, 1000000)
            plot.set_ylim(0, 1000000)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.2, hspace=0.2)
    plt.savefig(output_file, transparent='True', format='pdf')

def monosome_over_mrnp_reproducibility(mse):
    output_file = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        'mono_over_mRNP_plots.pdf')
    num_libs = len(mse.monosome_libs)
    num_plots_wide = num_libs-1
    num_plots_high = num_libs-1
    fig = plt.figure(figsize=(8,8))

    for i in range(len(mse.monosome_libs)):
        for j in range(i+1, len(mse.monosome_libs)):
            plot_index = (j-1)*(num_plots_wide)+(i+1)
            plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
            if j == num_plots_high:
                plot.set_xlabel("%s / %s RPM" % (mse.monosome_libs[i].lib_settings.sample_name, mse.mrnp_libs[i].lib_settings.sample_name))
            if i == 0:
                plot.set_ylabel("%s / %s RPM" % (mse.monosome_libs[j].lib_settings.sample_name, mse.mrnp_libs[j].lib_settings.sample_name))
            plot.set_xscale('symlog', linthreshx=0.01)
            plot.set_yscale('symlog', linthreshy=0.01)
            x = mse.monosome_libs[i].name_sorted_rpms()/mse.mrnp_libs[i].name_sorted_rpms()
            y = mse.monosome_libs[j].name_sorted_rpms()/mse.mrnp_libs[j].name_sorted_rpms()
            plot.scatter(x, y, color=ms_utils.black, s=3)
            plot.plot(numpy.arange(0,1000,1), numpy.arange(0,1000,1), color=ms_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)
            fx, fy = ms_utils.filter_x_y_pairs(x, y)
            r, p = stats.pearsonr(fx, fy)
            plot.annotate('rho,r=%.3f,%.3f' % (rho, r), xy=(0, 0.9), xytext=(0, 0.9), textcoords='axes fraction')
            plot.set_xlim(0, 1000)
            plot.set_ylim(0, 1000)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
    plt.savefig(output_file, transparent='True', format='pdf')

def monosome_over_total_reproducibility(mse):
    output_file = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        'mono_over_total_plots.pdf')
    num_libs = len(mse.monosome_libs)
    num_plots_wide = num_libs-1
    num_plots_high = num_libs-1
    fig = plt.figure(figsize=(8,8))

    for i in range(len(mse.monosome_libs)):
        for j in range(i+1, len(mse.monosome_libs)):
            plot_index = (j-1)*(num_plots_wide)+(i+1)
            plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
            if j == num_plots_high:
                plot.set_xlabel("%s / %s RPM" % (mse.monosome_libs[i].lib_settings.sample_name, mse.total_libs[i].lib_settings.sample_name))
            if i == 0:
                plot.set_ylabel("%s / %s RPM" % (mse.monosome_libs[j].lib_settings.sample_name, mse.total_libs[j].lib_settings.sample_name))
            plot.set_xscale('symlog', linthreshx=0.01)
            plot.set_yscale('symlog', linthreshy=0.01)
            x = mse.monosome_libs[i].name_sorted_rpms()/mse.total_libs[i].name_sorted_rpms()
            y = mse.monosome_libs[j].name_sorted_rpms()/mse.total_libs[j].name_sorted_rpms()
            plot.scatter(x, y, color=ms_utils.black, s=3)
            plot.plot(numpy.arange(0,1000,1), numpy.arange(0,1000,1), color=ms_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)
            fx, fy = ms_utils.filter_x_y_pairs(x, y)
            r, p = stats.pearsonr(fx, fy)
            plot.annotate('rho,r=%.3f,%.3f' % (rho, r), xy=(0, 0.9), xytext=(0, 0.9), textcoords='axes fraction')
            plot.set_xlim(0, 1000)
            plot.set_ylim(0, 1000)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
    plt.savefig(output_file, transparent='True', format='pdf')

def monosome_over_mrnp_plus_monosome_reproducibility(mse):
    output_file = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        'mono_over_mRNP_plus_mono_plots.pdf')
    num_libs = len(mse.monosome_libs)
    num_plots_wide = num_libs-1
    num_plots_high = num_libs-1
    fig = plt.figure(figsize=(8,8))

    for i in range(len(mse.monosome_libs)):
        for j in range(i+1, len(mse.monosome_libs)):
            plot_index = (j-1)*(num_plots_wide)+(i+1)
            plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
            if j == num_plots_high:
                plot.set_xlabel("%s / (%s+%s) RPM" % (mse.monosome_libs[i].lib_settings.sample_name,
                                                      mse.monosome_libs[i].lib_settings.sample_name,
                                                      mse.mrnp_libs[i].lib_settings.sample_name))
            if i == 0:
                plot.set_ylabel("%s / (%s+%s) RPM" % (mse.monosome_libs[j].lib_settings.sample_name,
                                                      mse.monosome_libs[j].lib_settings.sample_name,
                                                      mse.mrnp_libs[j].lib_settings.sample_name))
            #plot.set_xscale('symlog', linthreshx=0.01)
            #plot.set_yscale('symlog', linthreshy=0.01)
            x = mse.monosome_libs[i].name_sorted_rpms()/(mse.mrnp_libs[i].name_sorted_rpms()+mse.monosome_libs[i].name_sorted_rpms())
            y = mse.monosome_libs[j].name_sorted_rpms()/(mse.mrnp_libs[j].name_sorted_rpms()+mse.monosome_libs[j].name_sorted_rpms())
            plot.scatter(x, y, color=ms_utils.black, s=3)
            plot.plot(numpy.arange(0,1000,1), numpy.arange(0,1000,1), color=ms_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)
            fx, fy = ms_utils.filter_x_y_pairs(x, y)
            r, p = stats.pearsonr(fx, fy)
            plot.annotate('rho,r=%.3f,%.3f' % (rho, r), xy=(0, 0.9), xytext=(0, 0.9), textcoords='axes fraction')
            plot.set_xlim(0, 1)
            plot.set_ylim(0, 1)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
    plt.savefig(output_file, transparent='True', format='pdf')

def recruitment_change_rank_value_plot_interactive(mse, annotation_file, read_cutoff = 128, corrected_p_cutoff = 0.05):
    from bokeh.plotting import figure, output_file, show, save, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict
    set_name1, set_name2, matched_set = mse.parse_matched_set_annotation(annotation_file)

    # output to static HTML file
    output_file_name = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        '%s_%s_recruitment_change_rank_value.html' % (set_name1, set_name2))
    output_file(output_file_name)

    all_change_scores = {}
    all_change_score_means = {}
    all_p_values = {}
    all_annotations = {}
    for matched_pool_seqs in matched_set:
        set1_scores = []
        set2_scores = []
        for i in range(len(mse.monosome_libs)):
            set_1_counts = mse.monosome_libs[i].get_counts(matched_pool_seqs[0]) \
                           + mse.mrnp_libs[i].get_counts(matched_pool_seqs[0])
            set_2_counts = mse.monosome_libs[i].get_counts(matched_pool_seqs[1]) \
                           + mse.mrnp_libs[i].get_counts(matched_pool_seqs[1])
            # include only comparisons where the average number of reads is high enough
            if set_1_counts >= read_cutoff and set_2_counts >= read_cutoff:
                set1_score = mse.monosome_libs[i].get_rpm(matched_pool_seqs[0]) / \
                             (mse.monosome_libs[i].get_rpm(matched_pool_seqs[0]) +
                              mse.mrnp_libs[i].get_rpm(matched_pool_seqs[0]))
                set2_score = mse.monosome_libs[i].get_rpm(matched_pool_seqs[1]) / \
                             (mse.monosome_libs[i].get_rpm(matched_pool_seqs[1]) +
                              mse.mrnp_libs[i].get_rpm(matched_pool_seqs[1]))
            else:
                set1_score = float('nan')
                set2_score = float('nan')
            set1_scores.append(set1_score)
            set2_scores.append(set2_score)
        scores_1_filtered, scores_2_filtered = ms_utils.filter_x_y_pairs(set1_scores, set2_scores)
        recruitment_changes = numpy.array(scores_1_filtered) - numpy.array(scores_2_filtered)
        if len(scores_1_filtered) > 0 and len(scores_2_filtered) > 0:
            comparison = (matched_pool_seqs[0], matched_pool_seqs[1])
            t, p = stats.ttest_ind(scores_1_filtered, scores_2_filtered)
            all_change_scores[comparison] = recruitment_changes
            all_p_values[comparison] = p
            average = numpy.average(recruitment_changes)
            all_change_score_means[comparison] = average
            all_annotations[comparison] = '%s-%s=%.3f, p=%f' % (matched_pool_seqs[0], matched_pool_seqs[1], average, p)
    bh_corrected_p_values = ms_utils.bonferroniCorrection(all_p_values)
    sorted_means = sorted(all_change_score_means.iteritems(), key=operator.itemgetter(1))
    sig_ranks = [] #will store rank values for passing p values
    sig_means = []
    sig_com1 = []
    sig_com2 = []
    sig_p = []
    sig_n = []

    insig_ranks = [] #will store rank values for failing p values
    insig_means = []
    insig_anno = []
    insig_com1 = []
    insig_com2 = []
    insig_p = []
    insig_n = []
    for rank in range(len(sorted_means)):
        comparison, mean = sorted_means[rank]
        if bh_corrected_p_values[comparison]<corrected_p_cutoff:
            sig_ranks.append(rank)
            sig_means.append(mean)
            sig_com1.append(comparison[0])
            sig_com2.append(comparison[1])
            sig_p.append(bh_corrected_p_values[comparison])
            sig_n.append(len(all_change_scores[comparison]))
        else:
            insig_ranks.append(rank)
            insig_means.append(mean)
            insig_anno.append(all_annotations[comparison])
            insig_com1.append(comparison[0])
            insig_com2.append(comparison[1])
            insig_p.append(bh_corrected_p_values[comparison])
            insig_n.append(len(all_change_scores[comparison]))
    all_ranks = range(len(sorted_means))
    all_max = [max(all_change_scores[sorted_means[rank][0]]) for rank in all_ranks]
    all_min = [min(all_change_scores[sorted_means[rank][0]]) for rank in all_ranks]
    source = ColumnDataSource(data=dict(x=insig_ranks, y=insig_means, com1=insig_com1, com2=insig_com2, p=insig_p, n=insig_n, value=insig_means))
    sig_source = ColumnDataSource(data=dict(x=sig_ranks, y=sig_means, com1=sig_com1, com2=sig_com2, p=sig_p, n=sig_n, value=sig_means))
    max_source = ColumnDataSource(data=dict(x=all_ranks, y=all_max))
    min_source = ColumnDataSource(data=dict(x=all_ranks, y=all_min))
    hover = HoverTool(names=['insig', 'sig'])
    TOOLS = "pan,wheel_zoom,reset,save"
    PlotFig = figure(x_axis_label="rank", y_axis_label="%s-%s recruitment change" % (set_name1, set_name2),
                     tools=[TOOLS,hover], toolbar_location="right")
    PlotFig.circle("x", "y", size=5, source=source, color=ms_utils.bokeh_black, name = 'insig')
    PlotFig.circle("x", "y", size=5, source=sig_source, color=ms_utils.bokeh_vermillion, name = 'sig')
    # adjust what information you get when you hover over it

    hover.tooltips = OrderedDict([("%s" % set_name1, "@com1"),
                                  ("%s" % set_name2, "@com2"),
                                  ("mean", "@value"),
                                  ("Bonf. p", "@p"),
                                  ("n", "@n")])

    PlotFig.line("x", "y", line_width=1, source=min_source, color=ms_utils.bokeh_skyBlue)
    PlotFig.line("x", "y", line_width=1, source=max_source, color=ms_utils.bokeh_skyBlue)
    PlotFig.x_range = Range1d(start=-1, end=len(sorted_means))
    PlotFig.y_range = Range1d(start=-1, end=1)
    save(PlotFig)

def recruitment_fold_change_rank_value_plot_interactive(mse, annotation_file, read_cutoff = 128, corrected_p_cutoff = 0.05):
    from bokeh.plotting import figure, output_file, show, save, ColumnDataSource, gridplot
    from bokeh.models import Range1d
    from bokeh.models import HoverTool
    from collections import OrderedDict
    set_name1, set_name2, matched_set = mse.parse_matched_set_annotation(annotation_file)

    # output to static HTML file
    output_file_name = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        '%s_%s_recruitment_fold_change_rank_value.html' % (set_name1, set_name2))
    output_file(output_file_name)

    all_change_scores = {}
    all_change_score_means = {}
    all_p_values = {}
    all_annotations = {}
    for matched_pool_seqs in matched_set:
        set1_scores = []
        set2_scores = []
        for i in range(len(mse.monosome_libs)):
            set_1_counts = mse.monosome_libs[i].get_counts(matched_pool_seqs[0]) \
                           + mse.mrnp_libs[i].get_counts(matched_pool_seqs[0])
            set_2_counts = mse.monosome_libs[i].get_counts(matched_pool_seqs[1]) \
                           + mse.mrnp_libs[i].get_counts(matched_pool_seqs[1])
            # include only comparisons where the average number of reads is high enough
            if set_1_counts >= read_cutoff and set_2_counts >= read_cutoff:
                set1_score = mse.monosome_libs[i].get_rpm(matched_pool_seqs[0]) / \
                             (mse.monosome_libs[i].get_rpm(matched_pool_seqs[0]) +
                              mse.mrnp_libs[i].get_rpm(matched_pool_seqs[0]))
                set2_score = mse.monosome_libs[i].get_rpm(matched_pool_seqs[1]) / \
                             (mse.monosome_libs[i].get_rpm(matched_pool_seqs[1]) +
                              mse.mrnp_libs[i].get_rpm(matched_pool_seqs[1]))
            else:
                set1_score = float('nan')
                set2_score = float('nan')
            set1_scores.append(set1_score)
            set2_scores.append(set2_score)
        scores_1_filtered, scores_2_filtered = ms_utils.filter_x_y_pairs(set1_scores, set2_scores)
        recruitment_changes = numpy.array(scores_1_filtered) / numpy.array(scores_2_filtered)
        if len(scores_1_filtered) > 0 and len(scores_2_filtered) > 0:
            comparison = (matched_pool_seqs[0], matched_pool_seqs[1])
            t, p = stats.ttest_ind(scores_1_filtered, scores_2_filtered)
            all_change_scores[comparison] = recruitment_changes
            all_p_values[comparison] = p
            average = numpy.average(recruitment_changes)
            all_change_score_means[comparison] = average
            all_annotations[comparison] = '%s-%s=%.3f, p=%f' % (matched_pool_seqs[0], matched_pool_seqs[1], average, p)
    bh_corrected_p_values = ms_utils.bonferroniCorrection(all_p_values)
    sorted_means = sorted(all_change_score_means.iteritems(), key=operator.itemgetter(1))
    sig_ranks = [] #will store rank values for passing p values
    sig_means = []
    sig_com1 = []
    sig_com2 = []
    sig_p = []
    sig_n = []

    insig_ranks = [] #will store rank values for failing p values
    insig_means = []
    insig_anno = []
    insig_com1 = []
    insig_com2 = []
    insig_p = []
    insig_n = []
    for rank in range(len(sorted_means)):
        comparison, mean = sorted_means[rank]
        if bh_corrected_p_values[comparison]<corrected_p_cutoff:
            sig_ranks.append(rank)
            sig_means.append(mean)
            sig_com1.append(comparison[0])
            sig_com2.append(comparison[1])
            sig_p.append(bh_corrected_p_values[comparison])
            sig_n.append(len(all_change_scores[comparison]))
        else:
            insig_ranks.append(rank)
            insig_means.append(mean)
            insig_anno.append(all_annotations[comparison])
            insig_com1.append(comparison[0])
            insig_com2.append(comparison[1])
            insig_p.append(bh_corrected_p_values[comparison])
            insig_n.append(len(all_change_scores[comparison]))
    all_ranks = range(len(sorted_means))
    all_max = [max(all_change_scores[sorted_means[rank][0]]) for rank in all_ranks]
    all_min = [min(all_change_scores[sorted_means[rank][0]]) for rank in all_ranks]
    source = ColumnDataSource(data=dict(x=insig_ranks, y=insig_means, com1=insig_com1, com2=insig_com2, p=insig_p, n=insig_n, value=insig_means))
    sig_source = ColumnDataSource(data=dict(x=sig_ranks, y=sig_means, com1=sig_com1, com2=sig_com2, p=sig_p, n=sig_n, value=sig_means))
    max_source = ColumnDataSource(data=dict(x=all_ranks, y=all_max))
    min_source = ColumnDataSource(data=dict(x=all_ranks, y=all_min))
    hover = HoverTool(names=['insig', 'sig'])
    TOOLS = "pan,wheel_zoom,reset,save"
    PlotFig = figure(x_axis_label="rank", y_axis_label="%s/%s fold recruitment change" % (set_name1, set_name2),
                     tools=[TOOLS,hover], toolbar_location="right", y_axis_type="log")
    PlotFig.circle("x", "y", size=5, source=source, color=ms_utils.bokeh_black, name = 'insig')
    PlotFig.circle("x", "y", size=5, source=sig_source, color=ms_utils.bokeh_vermillion, name = 'sig')
    # adjust what information you get when you hover over it

    hover.tooltips = OrderedDict([("%s" % set_name1, "@com1"),
                                  ("%s" % set_name2, "@com2"),
                                  ("mean", "@value"),
                                  ("Bonf. p", "@p"),
                                  ("n", "@n")])

    PlotFig.line("x", "y", line_width=1, source=min_source, color=ms_utils.bokeh_skyBlue)
    PlotFig.line("x", "y", line_width=1, source=max_source, color=ms_utils.bokeh_skyBlue)
    PlotFig.x_range = Range1d(start=-1, end=len(sorted_means))
    PlotFig.y_range = Range1d(start=.01, end=100)
    save(PlotFig)

def plot_recruitment_violins(mse, annotation_file, read_cutoff = 128):

    #Makes violin plots of recruitment scores
    set_name1, set_name2, matched_set = mse.parse_matched_set_annotation(annotation_file)

    # output to static HTML file
    output_file_name = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        '%s_%s_violin.pdf' % (set_name1, set_name2))

    legends = []
    data = []
    set1_seqs = [pair[0] for pair in matched_set]
    set2_seqs = [pair[1] for pair in matched_set]
    all_seqs = mse.monosome_libs[0].sorted_names()
    set_members = {set_name1:set1_seqs, set_name2:set2_seqs, 'all':all_seqs}

    for lib_index in range(len(mse.monosome_libs)):
        for set_type in [set_name1, set_name2]:
            legends.append('%s %s' % (mse.monosome_libs[lib_index].lib_settings.sample_name, set_type))
            scores = []
            for seq_name in set_members[set_type]:
                counts = mse.monosome_libs[lib_index].get_counts(seq_name) \
                         + mse.mrnp_libs[lib_index].get_counts(seq_name)
                if counts >= read_cutoff:
                    recruitment_score = mse.monosome_libs[lib_index].get_rpm(seq_name) / \
                                 (mse.monosome_libs[lib_index].get_rpm(seq_name) +
                                  mse.mrnp_libs[lib_index].get_rpm(seq_name))
                    scores.append(recruitment_score)
            data.append(scores)
    p_file_name = os.path.join(
        mse.settings.get_rdir(),
        'plots',
        '%s_%s_violin.ks_p.txt' % (set_name1, set_name2))
    g=open(p_file_name, 'w')
    g.write('datset1\tdataset2\tKS d\tKS p\tt-test t\tt-test p\n')
    for i in range(len(legends)):
        for j in range(len(legends)):
            d, p = stats.ks_2samp(data[i], data[j])
            tind, pind = stats.ttest_ind(data[i], data[j])
            g.write('%s\t%s\t%.3f\t%.3e\t%.3f\t%.3e\t\n' % (legends[i], legends[j], d, p, tind, pind))
    g.close()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # Hide the grid behind plot objects
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_axisbelow(True)
    ax1.set_ylabel('monosome recruitment score')
    #ax1.set_xlabel(ylabel)
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)

    pos = range(1,len(data)+1)  # starts at 1 to play nice with boxplot
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        d = [float(dm) for dm in d]
        k = stats.gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = numpy.arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        #print 'v=',v
        v = v/v.max()*w #scaling the violin to the available space
        if 'all' in legends[p-1]:
            color = (0, 0, 0)
        elif set_name1 in legends[p-1]:
            color = (0/255., 159/255., 115/255)
        elif set_name2 in legends[p-1]:
            color = (213/255., 94/255., 0)
        else:
            print legends[p-1]
        plt.fill_betweenx(x,p,v+p,facecolor=color,alpha=0.3)
        plt.fill_betweenx(x,p,-v+p,facecolor=color,alpha=0.3)
    if True:
        bplot = plt.boxplot(data,notch=1)
        plt.setp(bplot['boxes'], color='black')
        plt.setp(bplot['whiskers'], color='black')
        plt.setp(bplot['fliers'], color='red', marker='.')

    per50s = []
    i = 1
    for datum in data:
        #per50s.append(stats.scoreatpercentile(datum, 50))
        t = stats.scoreatpercentile(datum, 50)

        per50s.append(t)
        ax1.annotate(str(round(t,3)), xy=(i+0.1, t), xycoords='data', arrowprops=None, fontsize='small', color='black')
        i+= 1
    #ax1.set_xticks([0.0, 0.5, 1.0, 1.5])
    ax1.set_ylim(0, 1)
    xtickNames = plt.setp(ax1, xticklabels=legends)
    plt.setp(xtickNames, rotation=90, fontsize=6)
    plt.savefig(output_file_name, transparent='True', format='pdf')
