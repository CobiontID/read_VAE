# %% [markdown]
# # Interactive visualisation of read composition
# 
# ## License
# 
# <details>
# <summary>MIT License</summary>
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# </details>
# 
# 
# Copyright (c) 2022-2023 claudia-c-weber
# 
# More info: https://cobiontid.github.io/


# %% [markdown]
# Load dependencies

# %%
import argparse

import subprocess
import yaml
import os
import sys

import pandas as pd
import multiprocessing as mp
import datashader as ds
import panel as pn

import param
from holoviews.operation.datashader import datashade, dynspread
import holoviews as hv


import colorcet as cc
from matplotlib import cm, colors

import numpy as np
np.random.seed(42)

from io import StringIO

from itertools import chain, combinations

from collections import defaultdict

hv.extension('bokeh', logo=False)

#%%
# Skip this cell to specify config manually
parser = argparse.ArgumentParser(
    description='Dashboard to visualise 2D read embeddings.')

parser.add_argument("--config", help="Path to config file", required=True)

args = parser.parse_args()
print(args)
config_file = args.config

# %%
pn.extension(template='bootstrap', sizing_mode="scale_both", defer_load=True, loading_indicator=True, loading_spinner='arc', loading_color='blue')
#dash = pn.Column()
#dash.servable("Read VAE exploration dashboard")

#placeholder = pn.Row(
#            pn.Column(pn.indicators.LoadingSpinner(value=True, width=50, height=50)),
#        )
#placeholder.servable()
#%%
#config_file = "ilApaMono1.yml"

# %%
#!pip install Biopython
#from Bio.Blast.NCBIWWW import qblast
#from Bio.Blast import NCBIXML


# %% [markdown]
# # Set up plotting functions

# %%
# Utils to set up dataframe

def drop_samples(ix, max_reads):
    """Return indices of random samples from a batch of reads"""
    if len(ix) > max_reads:
        ix = np.random.choice(ix, size=max_reads, replace=False)
    return ix

def downsample(a, max_reads=50, nbins=250):
    """Drop samples from dataset in regions of the 2D histogram that exceed a threshold (max_reads)"""
    #nbins = 500
    h = np.histogram2d(a[:,0], a[:,1], bins=nbins)
    bins_y = np.searchsorted(h[2], a[:,1])
    bins_x = np.searchsorted(h[1], a[:,0])
    d = defaultdict(list)
    for i in range(len(a)):
        d[(bins_x[i], bins_y[i])].append(i)
    return  np.array(list(chain.from_iterable(drop_samples(v,max_reads) for v in d.values())))

def remove_ids(df, rm_list):
    return df[~df['id'].isin(rm_list)]

def get_category_labels(read_ids, label_files):
    """Assign labels to reads based on files with lists. Reads that do not appear in a list are labelled 0.
    All others are assigned integer labels corresponding to the order in which the files are listed - unless
     they are present in more than one set, in which case they are assigned to an additional category."""
    mask = np.array(len(read_ids)*[0], dtype="int32")
    if len(label_files) == 0:
        print("Nothing to label.")
        return mask
    # iterate over lists of classified reads and assign integer labels
    all_sets = []
    for i, cat in enumerate(label_files):
        seq_set = {j.strip("\n") for j in open(cat)}
        all_sets.append(seq_set)
        np.put(mask, np.where([seqid in seq_set for seqid in read_ids]), [i+1])
    #check if lists overlap
    print(f"{i+1} class(es).")
    nt = lambda a, b: all_sets[a].intersection(all_sets[b])
    if in_multiple := set().union(
        *[nt(*j) for j in combinations(range(i + 1), 2)]
    ):
        np.put(mask, np.where([seqid in in_multiple for seqid in read_ids]), [i+2])
        print(f"Adding extra bin containing intersect of sets: {i + 2}")
    return mask

def load_df(data, samples_bin=50, max_reads=50000000, keep_labelled=True):
    """Determine whether to downsample the dataset and fetch indices, then load into dataframe"""
    if len(data['vae']) > max_reads:
        print("Downsampling data.")
        idxs = downsample(data["vae"], samples_bin)
        # Keep labelled reads?
        if keep_labelled == True:
            labelled_idxs = np.where(data['classes'] > 0)[0]
            idxs = np.union1d(idxs, labelled_idxs)
        print("Downsampled.")
    else:
        idxs = np.array(range(len(data['vae'])))

    df = pd.DataFrame(data=data['vae'][idxs,:], columns=['x', 'y'])
    df['id'] = pd.Series(data["reads"][idxs])
    df['hex'] = pd.Series(data['annot'][idxs])
    df['fastk'] = pd.Series(data['fastk'][idxs])
    df['bin'] = 99
    df['classes'] = pd.Categorical(data['classes'][idxs], ordered=True)
    #print(df)
    return df

# From https://personal.sron.nl/~pault/ (included for Bokeh 2.x compatibility)
clrs = ['#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4',
                '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0',
                '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265',
                '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141',
                '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518',
                '#DC050C', '#A5170E', '#72190E', '#42150A']
indexes = [[9], [9, 25], [9, 17, 25], [9, 14, 17, 25], [9, 13, 14, 17,
            25], [9, 13, 14, 16, 17, 25], [8, 9, 13, 14, 16, 17, 25], [8,
            9, 13, 14, 16, 17, 22, 25], [8, 9, 13, 14, 16, 17, 22, 25, 27],
            [8, 9, 13, 14, 16, 17, 20, 23, 25, 27], [8, 9, 11, 13, 14, 16,
            17, 20, 23, 25, 27], [2, 5, 8, 9, 11, 13, 14, 16, 17, 20, 23,
            25], [2, 5, 8, 9, 11, 13, 14, 15, 16, 17, 20, 23, 25], [2, 5,
            8, 9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25], [2, 5, 8, 9, 11,
            13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 8, 9, 11,
            13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 7, 8, 9, 11,
            13, 14, 15, 16, 17, 19, 21, 23, 25, 27], [2, 4, 6, 7, 8, 9, 11,
            13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27], [1, 3, 4, 6, 7, 8,
            9, 11, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26, 27], [1, 3, 4,
            6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 19, 21, 23, 25, 26,
            27], [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20,
            22, 24, 25, 26, 27], [1, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15,
            16, 17, 18, 20, 22, 24, 25, 26, 27, 28], [0, 1, 3, 4, 6, 7, 8,
            9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28]]


# Interactive plotting

class Scatter(param.Parameterized):
    """Build scatterplot for reads. Set up widgets to control display parameters, and 
    the number of discrete bins for the coding density annotation (num_bins) and the 
    coverage range displayed (upper, lower)."""

    min_alpha = param.Integer(50,  bounds=(10, 255), doc="Set the minimum alpha value for points.", label="Minimum alpha")
    num_bins = param.Integer(5, bounds=(1, 10), doc="Select number of quantile bins for coding density (hex).", label="Number of bins")
    upper = param.Integer(32767, doc="Maximum k-mer coverage to display.", label="Max k-mer coverage")
    lower = param.Integer(0, doc="Minimum k-mer coverage to display.", label="Min k-mer coverage")
    bg = param.Selector(["white", "grey", "black"], doc="Select the background colour for the plot.", label="Background colour")
    show_class = param.ListSelector(default=[], objects=[], label='Select classes')
    color_cat = param.Selector(["glasbey_hv", "colorblind_bokeh", "tol_rainbow"], label='Categorical colour scheme')
    action = param.Action(lambda x: x.param.trigger('action'), label='Update histogram for current selection')
    reverse_colours = param.Boolean(doc="Reverse colour map for binned annotations", label="Reverse colours")

    pn.config.throttled = True
    
    def __init__(self, df_complete, sample_id, **kwargs):
        super(Scatter, self).__init__(**kwargs)
        
            ####
        #TODO: Column scaling

        # Initialise dataframe
        self.df_complete = df_complete
        self.sample_id = sample_id
        self.df = self.make_bins(self.num_bins)

        # Get points, selection box, and summary t6able
        self.points = hv.Points(data=self.df, kdims=['x','y'],vdims=['id'])
        self.box = hv.streams.BoundsXY(source=self.points, bounds=(-0.5, -0.5, 0.5, 0.5))
        #self.bounds, self.dmap = self.selections()
        self.bounds = hv.DynamicMap(lambda bounds: hv.Bounds(bounds), streams=[self.box])
        self.dmap = hv.DynamicMap(lambda bounds:  hv.Table(self.df[(self.df['x'] > bounds[0]) & (self.df['x'] < bounds[2]) & \
            (self.df['y'] > bounds[1]) & (self.df['y'] < bounds[3])].head(n=5000).round(2)).opts(editable=True, width=600), \
                             streams=[self.box])
        
        
        # Set up colours for classes
        self.n_classes = df_complete['classes'].nunique()
        class_list = list(range(self.n_classes))
        self.param.show_class.objects = class_list
        self.show_class = class_list

        self.cat_maps = {'glasbey_hv': cc.b_glasbey_hv,
         'colorblind_bokeh': list(hv.Cycle.default_cycles["Colorblind"]),
          }

        # Drop maps that are too short
        if self.n_classes > 23:
            self.param.color_cat.objects = ["glasbey_hv"]
        else:
            self.cat_maps["tol_rainbow"] =  [clrs[i] for i in indexes[self.n_classes]]
            if self.n_classes > 9:
                self.param.color_cat.objects = ["glasbey_hv", "tol_rainbow"]

    @pn.depends('action')
    def hist_coverage(self):
        if () in self.dmap.data:
            max_val = self.dmap.data[()]["fastk"].max()
            min_val = self.dmap.data[()]["fastk"].min()
            h_counts, h_bins = np.histogram(self.dmap.data[()]["fastk"], bins=min(50,  max_val-min_val), range=(max(1, min_val),min(10000, max_val)))
        else:
            h_counts, h_bins = np.histogram(self.df["fastk"], bins=50, range=(1,10000))
        return hv.Histogram((np.log1p(h_counts), h_bins)).opts(width=600, height=150, shared_axes=False, ylabel="log(Frequency)", xlabel="fastk")
    
    @pn.depends('action')
    def hist_hexamer(self):
        if () in self.dmap.data:
            #max_val = self.dmap.data[()]["hex"].max()
            #min_val = self.dmap.data[()]["hex"].min()
            h_counts, h_bins = np.histogram(self.dmap.data[()]["hex"], bins=50)
            #h_counts_base, h_bins_base = np.histogram(self.df["hex"], bins=50)

            return hv.Histogram((np.log1p(h_counts), h_bins)).opts(width=600, height=150, shared_axes=False, ylabel="log(Frequency)", xlabel="hexamer")
        else:
            h_counts, h_bins = np.histogram(self.df["hex"], bins=50)
            return hv.Histogram((np.log1p(h_counts), h_bins)).opts(width=600, height=150, shared_axes=False, ylabel="log(Frequency)", xlabel="hexamer")

    #@pn.depends('num_bins', watch=True)
    def make_bins(self, num_bins):
        """Bin coding density into quantiles, where num_bins is the number of bins. Then, call function to filter rows by coverage"""
        # If no data provided, fill with 0
        if self.df_complete['hex'].min() == self.df_complete['hex'].max():
            self.df_complete['bin'] = 0
        else:
            bins = pd.Categorical(pd.qcut(self.df_complete['hex'], num_bins, labels=False, duplicates='drop'))
            self.df_complete['bin'] = bins
        self.df = self.filter_df()
        return self.df

    def filter_df(self):
        """Get rows within specified coverage range, lower <= coverage <= upper. Filter selected classes"""
        if self.lower <= self.upper:
            self.df = self.df_complete.loc[(self.df_complete['fastk'] <= self.upper) & (self.df_complete['fastk'] >= self.lower)]
        else:
            self.df = self.df_complete
        if (self.df_complete['classes'].nunique() > 1) & (len(self.show_class) > 0):
            self.df = self.df[self.df['classes'].isin(self.show_class)]
        return self.df

    def colormap(self):
        """Set up colourmap using viridis, generate legend labels"""
        v = cm.get_cmap('viridis')
        colors_hex = [colors.rgb2hex(i)
                      for i in v(np.linspace(0, 1, self.num_bins))]
        if self.reverse_colours is True:
            colors_hex = colors_hex[::-1]
        legend_labels = dict(
            zip(
                list(range(self.num_bins)),
                [f"bin {i}" for i in range(self.num_bins)],
            )
        )
        return colors_hex, legend_labels
    
    def colormap_classes(self):
        """Set up categorical colourmap based on selection, generate legend labels"""
        colors_hex = self.cat_maps[self.color_cat][:self.n_classes]
        legend_labels = dict(
            zip(
                list(range(self.n_classes)),
                [f"class {i}" for i in range(self.n_classes)],
            )
        )
        return colors_hex, legend_labels
          
    @pn.depends('num_bins', 'upper', 'lower', 'show_class')
    def update_points(self):
        self.df = self.make_bins(self.num_bins)
        self.points = hv.Points(data=self.df, kdims=['x','y'],vdims=['id', 'bin', 'classes'])
        return 0
    
    @pn.depends('min_alpha', 'num_bins', 'upper', 'lower', 'bg', 'show_class', 'reverse_colours')
    def draw_scatter_table(self):
        self.update_points()
        colors_hex, legend_labels = self.colormap()
        return self.draw_shaded(colors_hex, legend_labels, "bin")

    @pn.depends('min_alpha', 'num_bins', 'upper', 'lower', 'bg', 'show_class', 'color_cat')
    def draw_scatter_table_classes(self):
        self.update_points()
        colors_hex, legend_labels = self.colormap_classes()
        return self.draw_shaded(colors_hex, legend_labels, "classes")

    def draw_shaded(self, colors_hex, legend_labels, col):
        shaded = datashade(
            self.points,
            aggregator=ds.count_cat(col),
            color_key=colors_hex,
            min_alpha=self.min_alpha,
        ).opts(
            bgcolor=self.bg,
            width=700,
            height=700,
            show_grid=True,
            tools=["box_select"],
            default_tools=[],
            legend_labels=legend_labels,
            legend_position='top_right',
            legend_offset=(0, 0),
            title=f'Reads for {self.sample_id}',
        )
        return hv.Overlay([dynspread(shaded, threshold=0.7)]).collate() * self.bounds.clone()


# %%

def get_seq_file(identifier, fasta, seq_preview):

    if indexed == True:
        seq_cmd = f'{samtools_path} faidx {fasta} "{identifier}" > temp.fa'
        #print(seq_cmd)
    else:
        if fasta.split(".")[-1] == "gz":
            grep = "zgrep"
        else:
            grep = "grep"
        seq_cmd = f"timeout 300s {grep} -A1 -m1 \"{identifier}\" {fasta} > temp.fa"
    seq_preview.object = "...fetching sequence..."
    seq = subprocess.run(seq_cmd, capture_output=True, shell=True)
    if seq.returncode != 0:
        return 1
    seq_preview.object = f'{open("temp.fa", "r").read()}'
    return 0


def box_selected_data_dl(box, df):
    """Select rows in table corresponding to selected points"""
    return df[(df['x'] > box.bounds[0]) & (df['x'] < box.bounds[2]) &
              (df['y'] > box.bounds[1]) & (df['y'] < box.bounds[3])]


def blast_function():
    """ Use local blast server (default).
    This function sends the retrieved fasta record to send the query to a local server, and returns the first five lines of the output.
    To use a different setup, simply modify "blast_cmd" in the line below.
    """
    blast_cmd = "timeout 300s curl -T temp.fa http://172.27.25.136:35227 | head -n5"
    blast = subprocess.run(blast_cmd, capture_output=True, shell=True)
    if blast.returncode == 0:
        return '{0}'.format(blast.stdout.decode('utf-8'))
    else:
        return f"Non-zero return code {blast} {blast.stderr.decode('utf-8')}"

def blast_remote():
    """This function will be called if the configuration specifies that the remote blast option should be used."""
    seq = "{}".format(open("temp.fa", "r").read())
    result_handle = qblast("blastn", "nt", seq, megablast=True)
    blast_record = NCBIXML.read(result_handle)
    return blast_record.alignments[0].title


def make_panel(scatter, fasta):

    def button_readid_click(event):
        text_readid.value = f'{scatter.dmap.data[()]["id"][0]}'

    def button_click_seq(event):
            if text_readid.value == "...":
                text_readid.value = f'{scatter.dmap.data[()]["id"][0]}'
            if get_seq_file(text_readid.value, fasta, seq_preview) == 0:
                seq = f'{open("temp.fa", "r").read()}'

    def button_click_blast(event):
        if text_readid.value == "...":
            text_readid.value = f'{scatter.dmap.data[()]["id"][0]}'
        blast_pane.object = "...running..."
        idle.value = True

        if get_seq_file(text_readid.value, fasta, seq_preview) == 0:
            seq = f'{open("temp.fa", "r").read()}'
            blast_pane.object = blast_function()
        else:
            blast_pane.object == 'Failed to get fasta'
        idle.value = False

    def find_read(event):
        if text_readid.value == "...":
            text_readid.value = f'{scatter.dmap.data[()]["id"][0]}'
            coord.value = ""
        locate = scatter.df.loc[scatter.df['id'] == text_readid.value]
        coord.value = "X: {:.4g}, Y: {:.4g}".format(
            locate['x'].values[0], locate['y'].values[0])

    def download_csv():
        reads = box_selected_data_dl(scatter.box, scatter.df)
        sio = StringIO()
        reads.to_csv(sio)
        sio.seek(0)
        sio.flush()
        return sio

    # Define buttons, panes, and actions
    button_readid = pn.widgets.Button(
        name='Get first sequence in selection', button_type='primary')
    text_readid = pn.widgets.TextInput(value='...')
    button_readid.on_click(button_readid_click)

    idle = pn.indicators.LoadingSpinner(value=False, width=50, height=50)

    button_blast = pn.widgets.Button(
        name='blastn selected sequence', button_type='primary')
    
    button_seq = pn.widgets.Button(name="Get sequence", button_type="primary")

    blast_pane = pn.pane.HTML("""Do megablast""", style={'background-color': '#fcfcfc', 'border': '1px solid black',
                                                         'padding': '5px', 'overflow': 'scroll', 'width': '310px', 'height': '100px'})

    button_blast.on_click(button_click_blast)
    button_seq.on_click(button_click_seq)

    coord = pn.widgets.StaticText(value="")
    button_find = pn.widgets.Button(
        name='Find read coordinates', button_type='primary')
    button_find.on_click(find_read)

    seq_preview = pn.pane.HTML(""" """, style={'background-color': '#fcfcfc', 'border': '1px solid black',
                                               'padding': '5px', 'overflow': 'scroll', 'width': '700px', 'height': '50px'})

    file_download = pn.widgets.FileDownload(callback=download_csv, filename='reads.txt',
                                 label='Download selected reads', button_type='success', auto=True, embed=False)

    def lay_out_elements():
        # Show class tab if multiple classes present
        if scatter.n_classes > 1:
            tabs = pn.Tabs(('Hexamer', scatter.draw_scatter_table),
                        ('Classes', scatter.draw_scatter_table_classes), dynamic=True)
        else:
            tabs = pn.Tabs(('Hexamer', scatter.draw_scatter_table))
        
        # Widgets for blast, displaying sequence 
        widgets_read_selection = pn.WidgetBox(pn.Row(text_readid), pn.Row(button_readid), pn.Row(pn.widgets.StaticText(
            name='Note', value='Click to update after drawing new selection')), pn.Row(button_seq), pn.Row(button_find), pn.Row(coord))

        widgets_blast = pn.WidgetBox(
            pn.Row(blast_pane), pn.Row(button_blast), idle,)

        # Configure class selection widget
        if len(scatter.param.show_class.objects) <= 11:
            multi_select = pn.widgets.CheckButtonGroup.from_param(
                scatter.param.show_class)
        else:
            multi_select = pn.widgets.MultiChoice.from_param(
                scatter.param.show_class)

        # Lay out widgets for parameter settings
        param_layout = pn.WidgetBox(scatter.param.min_alpha, scatter.param.num_bins, scatter.param.reverse_colours, scatter.param.upper, scatter.param.lower,
                                    'Background colour', pn.widgets.RadioButtonGroup.from_param(
                                        scatter.param.bg, name="Background colour"),
                                    )
        
        # Add elements to right column, depending on available annotations
        right_col = pn.Column(scatter.dmap)
        # Add histograms
        if plain_scatter != "True":
            for el in [scatter.hist_coverage, scatter.hist_hexamer, scatter.param.action]:
                right_col.append(el)
        else:
            scatter.param.num_bins.constant = True
            scatter.param.upper.constant = True
            scatter.param.lower.constant = True

        # Add widgets to filter by class
        if scatter.n_classes > 1:
            param_layout_classes = pn.WidgetBox(
                    "Filter classified sequences", scatter.param.color_cat, multi_select, width=600)
            right_col.append(param_layout_classes)


        filtered_view = pn.Row(
            pn.Column(param_layout,
                      widgets_read_selection,
                      widgets_blast,
                      pn.Row(file_download)),
            pn.Column(pn.panel(tabs), pn.Row(seq_preview)),
            right_col
        )
        return filtered_view

    return lay_out_elements()


# %% [markdown]
# # Load data
# 

# %%
#!wget https://github.com/CobiontID/CobiontID.github.io/raw/main/examples/ilCarKade1_204_downsampled.npz

# Load configuration

print("Loading configuration")
with open(config_file, 'r') as file:
    cfg = yaml.safe_load(file)


def ids_width(reads):
    """ Get max length for read ids to prevent truncation by np.loadtxt() """
    wc = subprocess.run(["wc", "-L", reads], capture_output=True)
    if wc.returncode == 0:
        return int(wc.stdout.decode('utf-8').split()[0])
    else:
        print(f"Couldn't get width for read identifiers in {reads}, check configuration")
        sys.exit(1)

tolid = cfg['tolid']
basedir = cfg["basedir"]

if "k" in cfg:
    k = cfg["k"]
else:
    k = 31

# Override default blast function
if "remote_blast" in cfg:
    if cfg["remote_blast"] is True:
        from Bio.Blast.NCBIWWW import qblast
        from Bio.Blast import NCBIXML
        blast_function = blast_remote
    
class_lists = cfg["class_lists"]
if class_lists is None:
    class_lists = []
else:
    if len(class_lists) != len(set(class_lists)):
        print("Warning: Check your class lists for duplicates!")

fasta = cfg["fasta"]

# Check whether to use an indexed fasta (for large files)
if ("samtools_path" in cfg) & ("indexed" in cfg):
    samtools_path = cfg["samtools_path"]
    indexed = cfg["indexed"]
    os.path.isfile(f'{fasta}.fai')
else:
    indexed = False

# Defaults are based on the directory structure produced by the Snakemake workflow for reads
default_path_dict = {
    "vae_path": f"{basedir}/vae/{tolid}/{tolid}.vae.out.2d.0",
    "fastk_path": f"{basedir}/fastk/reads/k_{k}/medians/{tolid}.median_{k}mer.txt",
    "hexamer_path": f"{basedir}/density/reads/{tolid}.reads.hexsum",
    "read_ids_path": f"{basedir}/kmer_counts/reads/k_4/ids/{tolid}.reads.ids.txt"
    }

# Override default path if specified
if "vae_path" in cfg:
    default_path_dict["vae_path"] = cfg["vae_path"]
if "fastk_path" in cfg:
    default_path_dict["fastk_path"] = cfg["fastk_path"]
if "hexamer_path" in cfg:
    default_path_dict["hexamer_path"] = cfg["hexamer_path"]
if "read_ids_path" in cfg:
    default_path_dict["read_ids_path"] = cfg["read_ids_path"]

width = ids_width(default_path_dict["read_ids_path"])

print(cfg)
# %%
plain_scatter = "False"
if "no_annot" in cfg:
    if cfg["no_annot"] is True:
        plain_scatter = "True"

def load_data_dict(default_path_dict, width, plain_scatter):
    """Retrieve data dictionary to pass to dataframe loader"""
    if plain_scatter == "False":
        data_dict = {
                    "fastk": np.loadtxt(default_path_dict["fastk_path"], dtype="int64"),
                    "annot": np.loadtxt(default_path_dict["hexamer_path"], dtype="float32"),
                    "reads": np.loadtxt(default_path_dict["read_ids_path"], dtype="U{}".format(width)),
                    }
    else:
        data_dict = {"reads": np.loadtxt(default_path_dict["read_ids_path"], dtype="U{}".format(width))}
        data_dict["annot"] = np.array((len(data_dict["reads"])-2)*[0] + [0.01, 1.], dtype="float32")
        data_dict["fastk"] = np.array(len(data_dict["reads"])*[1], dtype="int32")

    if ".npy" in default_path_dict["vae_path"]:
        data_dict['vae'] = np.load(default_path_dict["vae_path"])
    else:
        data_dict['vae'] = np.loadtxt(default_path_dict["vae_path"], dtype="float32")
    
    data_dict["classes"] = get_category_labels(data_dict["reads"], class_lists)

    return data_dict

# %%

# %%
assert os.path.isfile(fasta)

# %% [markdown]
# # Display data

# %%
#@title

sample_id = tolid

def load_dash_wrapper():
    print("Loading data...")
    data_dict = load_data_dict(default_path_dict, width, plain_scatter)
    # Show up to 15M points with 250 samples per bin
    xy = load_df(data_dict, 250, 15000000)
    print("Loading dashboard elements...")
    scatter = Scatter(xy, sample_id)
    view = make_panel(scatter, fasta)
    return(view)

#%%


# %% [markdown]
# ## How to read this plot
# 
# The x and y coordinates in the scatterplot below represent tetranucleotide counts for long sequencing reads, which have been reduced to two dimensions using a Variational Autoencoder. This allows sequences with different composition to be separated visually.
# 
# In the first tab ("Hexamer") each read is coloured by its estimated coding density (bins with smaller numbers, shown in purple, contain a smaller fraction of protein coding sequence than bins with larger colours, shown in yellow). If class labels are provided, the second tab shows reads coloured by class.

# %%
#pn.extension(template='bootstrap')
#hv.extension("bokeh")
#pn.serve(view)
pn.Column(load_dash_wrapper).servable("Read VAE dashboard")
