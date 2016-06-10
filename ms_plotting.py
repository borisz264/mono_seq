import ms_utils
import numpy
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import scipy.stats as stats

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