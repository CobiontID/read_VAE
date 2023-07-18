#!/usr/bin/env python
# cw21@sanger.ac.uk

import argparse

##################################
parser = argparse.ArgumentParser(
    description='Produce 2D representation of a multidimensional numpy array')

parser.add_argument(
    "--zfile", help="2D array with reduced counts", required=True)
parser.add_argument("--npy", help="Specify if input is a .npy file. Defaults to F (load from text)",
                    default="F", choices=["T", "F"])
parser.add_argument("--verbose", help="print extra model info",
                    default="F", choices=["T", "F"])
#TODO
parser.add_argument("--labels", help="add file with labels", default="")
parser.add_argument(
    "--fignames", help="add file with figure names", default="figure")
parser.add_argument("--superimpose", help="superimpose points",
                    default="F", choices=["T", "F"])
parser.add_argument("--categorical", help="Use contrasting categorical colours (defaults to gradient)",
                    default="F", choices=["T", "F"])
parser.add_argument("--outdir", help="Output directory", default="")
parser.add_argument(
    "--cmap", help="Custom colour map as RGB, e.g. \"#d60000 #8c3bff #018700\"", default="")
parser.add_argument("--remove", help="Remove all classified points (requires --labels to be provided)",
                    default="F", choices=["T", "F"])
parser.add_argument("--canvas_size", help="Width/height of output image (square canvas)",
                    type=int, default=1000, choices=[500, 1000, 1500])
parser.add_argument("--edges", help="Bin edges", default="", required=False)
parser.add_argument("--legend_y_label", help="Label for legend (applies if edges supplied)", default="", required=False)


args = parser.parse_args()
print(args)

zfile = args.zfile
categorical = args.categorical
verbose = args.verbose
labels = args.labels
figures = args.fignames
superimpose = args.superimpose
outdir = args.outdir
custom_cmap = args.cmap
remove_classified = args.remove
is_npy = args.npy
edges = args.edges
legend_y_label = args.legend_y_label

canvas_size = args.canvas_size

if custom_cmap != "":
    custom_cmap = custom_cmap.split()
    categorical = "T"
    print(custom_cmap)
    print(len(custom_cmap))

########################################

import sys
import pandas as pd
import datashader as ds

from holoviews.operation.datashader import datashade, dynspread
import holoviews as hv

import numpy as np

hv.extension('bokeh')
import colorcet as cc
from matplotlib import cm, colors

###########################################
## cmap defaults
default_cont_map = "viridis"
default_cat_map = list(hv.Cycle.default_cycles["Colorblind"])
longer_cat_map = cc.b_glasbey_hv

##########################################


def load_labels_filter(df, labels, remove_classified):
    """Load labels, and check length matches coordinates. Filter out all rows with labels that 
    do not equal zero if requested"""
    labels_int = [int(i.strip()) for i in open(labels)]
    if len(labels_int) == len(df):
        df['label'] = pd.Series(labels_int)
        if remove_classified == "T":
            df = df[df['label'] == 0]
        return df
    else:
        sys.exit(f"Length of label vector ({len(labels_int)}) doesn't match length of coordinates ({len(df)}). Exiting.")



def load_coords(is_npy, zfile, labels, remove_classified):
    """Load coordinates into dataframe from .npy or text"""
    if is_npy == "F":
        x_test_encoded = np.loadtxt(zfile, dtype="float32")
    else:
        x_test_encoded = np.load(zfile)

    df = pd.DataFrame(data=x_test_encoded, columns=['x', 'y'])

    if labels != "":
        df = load_labels_filter(df, labels, remove_classified)

    return df


def set_colormap(n_labels):
    """Set up colour maps, depending on whether a continuous or categorical map is needed.
    If a custom map is provided, check that the number of colours matches the number of labels"""
    if categorical == "F":
        v = cm.get_cmap(default_cont_map)
        return [colors.rgb2hex(i) for i in v(np.linspace(0, 1, n_labels))]
    elif (custom_cmap != "") & (len(custom_cmap) == n_labels):
        #TODO: Validity check?
        return custom_cmap
    elif n_labels <= len(default_cat_map):
        return default_cat_map[:n_labels]
    else:
        print(
            f"{n_labels} categories present, using longer colourmap. Consider merging or omitting some classes."
        )
        return longer_cat_map[:n_labels]

def get_legend_labels(n_labels):
    """Add text to legend labels (returns a dictionary)"""
    return dict(zip(list(range(n_labels)), [f"bin {i}" for i in range(n_labels)]))

def hook(plot, element):
    plot.state.yaxis.major_tick_line_color = None        # turn off x-axis major ticks
    plot.state.yaxis.minor_tick_line_color = None        # turn off x-axis minor ticks
    plot.state.yaxis.axis_line_color = None              # hide y-axis line
    plot.state.yaxis.axis_label_text_font = "sans serif"

def draw_legend_edges(edges, n_labels, colors_hex, legend_y_label, canvas_size):
    """Draw colourbar legend with bin edges as y ticks, using heatmap as workaround"""
    edges = [float(i.strip()) for i in open(edges)]
    assert len(edges) == n_labels+1
    hm_legend = hv.HeatMap(([0]*n_labels, range(n_labels), range(n_labels))).opts(
                cmap=colors_hex,
                ylabel=legend_y_label,
                xlabel="",
                xaxis=None,
                toolbar=None,
                active_tools=[],
                show_frame=False,
                border=2,
                frame_width=20,
                yaxis="right",
                frame_height=canvas_size,
                margin=10,
                yticks=[(i,j) for i,j in zip(
                    np.linspace(-0.49,n_labels+0.49, n_labels+2),
                    [f"{np.round(i, decimals=2)}" for i in edges]
                    )],
                shared_axes=False,
                hooks=[hook], alpha=0.9,
                padding=(0.0, 0.4),
            )
    return hm_legend

#######################################################
# Load data and render

df = load_coords(is_npy, zfile, labels, remove_classified)

# If no labels are provided, render using fire colormap
if labels == "":
    points = hv.Points(data=df, kdims=['x', 'y'])
    hv.save(
        dynspread(
            datashade(points, cmap=cc.fire).opts(
                bgcolor='black', width=canvas_size, height=canvas_size
            )
        ),
        f'{outdir}{figures}.2d_plot.png',
        fmt='png',
    )
else:
    n_labels = np.max(df['label'])+1
    colors_hex = set_colormap(n_labels)
    points = hv.Points(data=df, kdims=['x', 'y'], vdims=['label'])

    # If superimpose is set to False (default), render and shade all points together.
    if superimpose == "F":
        shaded = dynspread(
                datashade(
                    points,
                    aggregator=ds.count_cat('label'),
                    color_key=colors_hex,
                ).opts(
                    frame_width=canvas_size,
                    frame_height=canvas_size,
                    legend_labels=get_legend_labels(n_labels),
                )
            )
        if edges != "":
            hm_legend = draw_legend_edges(edges, n_labels, colors_hex, legend_y_label, canvas_size)
            shaded.opts(show_legend=False)
            shaded = shaded + hm_legend
        # Save plot
        hv.save(
            shaded,
            f'{outdir}{figures}.2d_plot_labelled.png',
            backend="bokeh",
            fmt='png',
        )
    else:
        points2 = hv.Points(data=df[df["label"] > 0], kdims=["x", "y"], vdims=["label"]).sort(
            "label").options(color="label", cmap=dict(zip(list(range(1, n_labels)), colors_hex)))

        hv.save(
            dynspread(
                datashade(points, cmap="grey").opts(
                    bgcolor='white', width=canvas_size, height=canvas_size
                )
            )
            * points2,
            f'{outdir}{figures}.2d_plot_superimpose.png',
            fmt='png',
        )
