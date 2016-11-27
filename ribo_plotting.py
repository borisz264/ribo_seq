import ribo_utils
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import scipy.stats as stats
import uniform_colormaps
import math

'''
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
            plot.scatter(x, y, color=ribo_utils.black, s=3)
            plot.plot(numpy.arange(0,1000000,1), numpy.arange(0,1000000,1), color=ribo_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)

            plot.annotate('rho=%.3f' % (rho), xy=(0, 0.8), xytext=(0, 0.8), textcoords='axes fraction')
            plot.set_xlim(0, 1000000)
            plot.set_ylim(0, 1000000)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.2, hspace=0.2)
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
            plot.scatter(x, y, color=ribo_utils.black, s=3)
            plot.plot(numpy.arange(0,1000,1), numpy.arange(0,1000,1), color=ribo_utils.vermillion, lw = 1, linestyle='dashed')
            rho, pval = stats.spearmanr(x, y)
            fx, fy = ribo_utils.filter_x_y_pairs(x, y)
            r, p = stats.pearsonr(fx, fy)
            plot.annotate('rho,r=%.3f,%.3f' % (rho, r), xy=(0, 0.9), xytext=(0, 0.9), textcoords='axes fraction')
            plot.set_xlim(0, 1000)
            plot.set_ylim(0, 1000)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
    plt.savefig(output_file, transparent='True', format='pdf')
'''

def plot_fragment_length_distributions(experiment, min_x=15):
    fig = plt.figure(figsize=(16, 16))
    num_libs = len(experiment.libs)
    num_plots_wide = math.sqrt(ribo_utils.next_square_number(num_libs))
    num_plots_high = math.sqrt(ribo_utils.next_square_number(num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index+1)
        sample_name = lib.lib_settings.sample_name
        fragment_length_counts = lib.get_all_fragment_length_counts()
        bins = range(0, max(fragment_length_counts.keys()))
        fragment_length_fractions = np.array([fragment_length_counts[length] for length in bins]) / float(
            lib.total_mapped_fragments)
        # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
        plot.plot(bins, fragment_length_fractions,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot_index += 1
        plot.set_xlabel("fragment length", fontsize=8)
        plot.set_ylabel("fraction of fragments", fontsize=8)
        plot.set_title(sample_name,  fontsize=8)
        plot.set_xlim(min_x, max(bins))
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'fragment_length_distributions.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_readthrough_box(experiment):
    fig = plt.figure(figsize=(8, 8))
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = 1
    colormap = uniform_colormaps.viridis
    plot_index = 0
    plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
    data = []
    legends = []
    boxprops = dict(linewidth=2, color=ribo_utils.black)
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        readthroughs = [tx.compute_readthrough_ratio(16, read_end='3p', read_lengths='all', cds_cutoff=128) for
                        tx in lib.transcripts.values() if not tx.compute_readthrough_ratio(16, read_end='3p',
                                                                                           read_lengths='all',
                                                                                           cds_cutoff=128) == None ]
        data.append(readthroughs)
        legends.append('%s (%d)' % (sample_name, len(readthroughs)))
        # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
    plot.boxplot(data, notch=True, boxprops=boxprops)
    plot_index += 1
    #plot.set_xlabel("fragment length", fontsize=8)
    plot.set_ylabel("log10 readthrough fraction", fontsize=8)
    plot.set_title(sample_name,  fontsize=8)
    plot.set_xticklabels(legends, rotation=40, ha='right')
    #plot.set_xlim(min_x, max(bins))
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'readthrough_box.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_frame_distributions(experiment, read_lengths = ['all', [28], [29], [30]], read_ends = ['5p', '3p']):
    num_libs = len(experiment.libs)
    num_plots_wide = len(read_lengths) * len(read_ends)
    num_plots_high = num_libs
    fig = plt.figure(figsize=(num_plots_wide, num_plots_high))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        for read_end in read_ends:
            for read_length in read_lengths:
                plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index+1)
                sample_name = lib.lib_settings.sample_name
                frame_counts = np.zeros(3)
                offsets = {'5p':-15, '3p': 15}
                offset = offsets[read_end]
                for transcript in lib.transcripts.values():
                    frame_counts = frame_counts + transcript.get_read_frame_counts(transcript.cds_start+offset, transcript.cds_end+offset, read_end=read_end, read_lengths=read_length)

                # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
                bar_corners = (np.arange(3)+0.25)*.5
                bar_width = 0.5*.5
                bar_centers = bar_corners + bar_width/2.0
                plot.bar(bar_corners, frame_counts/sum(frame_counts), width=bar_width, label=sample_name, lw=0,
                         color=ribo_utils.black)
                plot.set_ylim(0, 1)
                plot.set_xticks(bar_centers)
                plot.set_xticklabels([str(n) for n in range(3)])
                #plot.set_xlabel("fragment length", fontsize=8)
                plot.set_title('%s, %s' % (read_end, str(read_length)), fontsize=8)
                if plot_index % num_plots_wide == 0:
                    plot.set_ylabel(sample_name, fontsize=8)
                plt.setp(plot.get_xticklabels(), fontsize=7)
                plt.setp(plot.get_yticklabels(), fontsize=7)
                #plot.set_xlim(min_x, max(bins))
                plot_index += 1
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    #plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'frame_ distributions.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_start_codon_average(experiment, up = 100, down = 500, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if read_end == '5p':
                start_offset = -15
                stop_offset = -12
            elif read_end == '3p':
                start_offset = 14
                stop_offset = 18
            cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
            if cds_reads >= min_cds_reads:
                tx_count, tx_inclusion = transcript.get_read_counts_array(transcript.cds_start, -1*up, down, read_end=read_end,
                                                                read_lengths=read_lengths)
                normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS start", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1 * up, down)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'start_codon_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_stop_codon_average(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if read_end == '5p':
                start_offset = -15
                stop_offset = -12
            elif read_end == '3p':
                start_offset = 14
                stop_offset = 18
            cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
            if cds_reads >= min_cds_reads:
                tx_count, tx_inclusion = transcript.get_read_counts_array(transcript.cds_end, -1*up, down, read_end=read_end,
                                                                read_lengths=read_lengths)
                normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'stop_codon_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_first_exon_average(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if len(transcript.exon_starts)>1:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    tx_count, tx_inclusion = transcript.get_read_counts_array(transcript.exon_starts[1], -1*up, down, read_end=read_end,
                                                                    read_lengths=read_lengths)
                    normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                    inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to second exon start", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'first_ej_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_stop_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if read_end == '5p':
                start_offset = -15
                stop_offset = -12
            elif read_end == '3p':
                start_offset = 14
                stop_offset = 18
            cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
            if cds_reads >= min_cds_reads:
                length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.cds_end, -1*up, down,
                                                                                       read_end=read_end)
                length_sum += length_sum_array
                count_sum += counts_array
        nt_positions = np.arange(-1*up, down+1)
        plot.plot(nt_positions, length_sum/count_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        plot.set_ylim(25, 35)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'stop_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_start_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if read_end == '5p':
                start_offset = -15
                stop_offset = -12
            elif read_end == '3p':
                start_offset = 14
                stop_offset = 18
            cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
            if cds_reads >= min_cds_reads:
                length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.cds_start, -1*up, down,
                                                                                       read_end=read_end)
                length_sum += length_sum_array
                count_sum += counts_array
        nt_positions = np.arange(-1*up, down+1)
        plot.plot(nt_positions, length_sum/count_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        plot.set_ylim(25, 35)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'start_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_first_exon_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if len(transcript.exon_starts)>1:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.exon_starts[1], -1*up, down,
                                                                                           read_end=read_end)
                    length_sum += length_sum_array
                    count_sum += counts_array
        nt_positions = np.arange(-1*up, down+1)
        plot.plot(nt_positions, length_sum/count_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        plot.set_ylim(25, 35)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'first_ej_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()