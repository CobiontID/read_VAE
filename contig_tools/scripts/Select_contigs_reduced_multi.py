#!/usr/bin/env python
# Generates HTML file of decomposed tetranucleotide plots with binned annotations
# coding: utf-8
import sys
import argparse

##
parser = argparse.ArgumentParser(
    description='Plot decomposed contigs')
parser.add_argument("--infile", help="input .npy file of counts", required=True)
parser.add_argument("--outfile", help="Path to save html output file to.", default="contig_selection.html")
parser.add_argument("--seqidfile", help="input list of seqids", required=True)
parser.add_argument("--annotfiles", help="input hexsum, FastK results, coverage...", required=True)
parser.add_argument("--annotnames", help="list of names for dropdown menu", required=True)
parser.add_argument("--bins", help="input number of bins", default=5, type=int)

parser.add_argument("--speciesname", help="species name", default="sample_name")
parser.add_argument("--seqtype", help="Type of sequence (for plot label), defaults to contigs", default = "contigs")

parser.add_argument("--discretize", help="discretize by quantile or linear", default="quantile")
parser.add_argument("--pca", help="Do PCA instead of UMAP", default="F")
parser.add_argument("--keepbins", help="Use existing integer labels from file instead of binning (uses discrete colormap)", default = "F")

# Save coordinates
parser.add_argument("--save_coords", help="Path to file to store coordinates (.npy format)", default=None)
parser.add_argument("--save_tsv", help="Path to file to save tsv file", default=None)

# Plotting options
# Contig sizes
parser.add_argument("--ignore_sizes", help="Don't use contig size to scale points", default="F", choices=["F", "T"])
parser.add_argument("--rescale_max", help="Factor by which to divide contig size", default=15, type=int)
parser.add_argument("--rescale_min", help="Factor by which to divide contig size", default=150, type=int)
# Override axes
parser.add_argument("--override_x", help="Supply an alternative file with x coordinates", default="")
parser.add_argument("--override_y", help="Supply an alternative file with y coordinates", default="")
# Scale axes
parser.add_argument("--y_scale", help="Y axis scale", default="linear", choices=["linear", "log"])
parser.add_argument("--x_scale", help="X axis scale", default="linear", choices=["linear", "log"])


args = parser.parse_args()
print(args)

infile = args.infile
outfile = args.outfile
seqidfile = args.seqidfile
spname = args.speciesname
annotfiles = args.annotfiles.split()
annotnames = args.annotnames.split()
use_pca = args.pca
seqtype = args.seqtype
keep_bins = args.keepbins
save_coords = args.save_coords
save_tsv = args.save_tsv

ignore_sizes = args.ignore_sizes

RESCALE_MAX = args.rescale_max
RESCALE_MIN = args.rescale_min

override_x = args.override_x
override_y = args.override_y

y_axis_type = args.y_scale
x_axis_type = args.x_scale

# Override p_ctg label
if seqtype == "p_ctg":
    seqtype = "contigs"

try:
    assert len(annotfiles) == len(annotnames)
except:
    print("Error, annotation labels and annotation file list must have same length", file=sys.stderr)
    sys.exit(1)

n = args.bins

discretize = args.discretize

if None in [infile, outfile, seqidfile, annotfiles, annotnames]:
    print("--infile, --outfile, --annotfiles, --annotnames and --seqidfile must be provided", file=sys.stderr)
    sys.exit(1)


###
#imports
##

import numpy as np
np.random.seed(42)

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA



#from ipywidgets import interact
from bokeh.io import push_notebook, show, output_notebook, output_file, save
from bokeh.plotting import figure
from bokeh.models import Legend, LegendItem
from bokeh.models import ColumnDataSource, HoverTool, DataTable, TableColumn, CustomJS, Button, Select
from bokeh.events import ButtonClick
from bokeh.layouts import gridplot
from bokeh.models.glyphs import Circle

import matplotlib.colors
from matplotlib import cm


# ## Define helper functions

# In[2]:


def load_counts_norm(infile):
    counts = np.load(infile)
    sum_all = counts.sum(axis=1)[:,None]
    counts_n = counts/sum_all
    # Handle nan
    counts_n = np.nan_to_num(counts_n)
    return counts_n, sum_all



def do_umap(counts, n_components=2, y=None):
    scaler = StandardScaler()
    counts_t = scaler.fit_transform(counts)
    transform_umap = umap.UMAP(n_components=n_components, random_state=42).fit(counts_t, y=y)
    return transform_umap, scaler



def transform_umap(counts, transformer, scaler):
    return( transformer.transform(scaler.transform(counts)))


def do_pca(counts, n_components=2):
    scaler = StandardScaler()
    counts_t = scaler.fit_transform(counts)
    pca = PCA(n_components=2)
    return pca.fit_transform(counts_t)

# Plotting functions
def make_dict_multi(n_inputs):
    dict_string = "x=x,\ny=y,\nidentifier=contig_ids,\n"
    for i in range(n_inputs):
        dict_string += "color_{0}=[color_key[int(c)][1] for c in bins[{0}]],\nlabel_{0}=conts[{0}],\nmask_{0}=bins[{0}],\n".format(i)
    dict_string = "dict(\n{}\n)".format(dict_string[:-2])
    #print(dict_string)
    #dict_gen = ColumnDataSource(eval(dict_string))
    dict_gen = dict_string
    return dict_gen
    
def make_update_callback(cont_names):
    c_string = ""
    for i, c_name in enumerate(cont_names):
        c_string += "if (sel == '{0}') {{\
             for (var i=0;i<d1['color_{1}'].length; i++) {{\
               d2['color'].push(d1['color_{1}'][i]);\
               d2['label'].push(d1['label_{1}'][i]);\
               d2['mask'].push(d1['mask_{1}'][i]);\
             }}\
            }};\n\n".format(c_name, i)
    return c_string

def scale_points_axes(x, y, sizes, max_scale, min_scale):
    """If contig sizes are not disregarded, set sizes of points for largest and smallest sequences according to the coordinate ranges.
    The axis with the smaller range is considered.
    "max_scale": Roughly corresponds to how many times the largest contig should fit along the axis with the smaller range
    "min_scale": How many times the smallest contig should fit along the axis with the smaller range"""
    x_range = np.max(x) - np.min(x)
    y_range = np.max(y) - np.min(y)
    min_range = min(x_range, y_range)
    sizes_scaled = sizes/np.max(sizes)*(min_range/max_scale) + min_range/min_scale
    sizes_scaled = sizes_scaled/2
    return sizes_scaled

##
def draw_bokeh_multi(x, y, contig_ids, conts, bins, sizes, cont_names, spname, outfile, decomp):

    #color_mapper = LinearColorMapper(palette="Viridis256", low=0, high=cont1.max())
    hover = HoverTool(tooltips=[
        ("index", "$index"),
        ("(x,y)", "(@x, @y)"),
        ('label', '@label'),
        ('identifier', '@identifier'),
        ('size', '@sizes')
    ])

    # Caution: Hard-coded tetranucs
    p = figure(title="{} {} {} reduced with {}".format(spname, seqtype[:-1], "tetranucleotides", decomp),
       tools=[hover,'pan', 'box_zoom','wheel_zoom', 'reset','lasso_select','save'],
       x_axis_label='1', y_axis_label='2',  plot_width=600, plot_height=500,
       y_axis_type=y_axis_type, x_axis_type=x_axis_type
    )
    
    s = ColumnDataSource(eval(make_dict_multi(len(conts))))
    print(len(conts))
    
    sizes_scaled = scale_points_axes(x, y, sizes, RESCALE_MAX, RESCALE_MIN)
    
    s1 = ColumnDataSource(data=dict(x=x, y=y, identifier=contig_ids, color=[color_key[int(c)][1] for c in bins[0]],\
                                 label=conts[0], mask=bins[0], sizes_scaled=sizes_scaled, sizes=sizes)) 

    s2 = ColumnDataSource(data=dict(x=[], y=[], identifier=[], color=[], label=[], sizes=[]))

    if ignore_sizes == "F":
        p.circle( x='x', y='y', radius="sizes_scaled", alpha=0.5, color='color', source=s1)
    else:
        p.circle( x='x', y='y', radius=0.1, alpha=0.5, color='color', source=s1, size=0.1)
    
    line = {}
    for i in range(n):
        line_glyph = Circle(x=0, y=0, fill_color=color_key[int(i)][1])
        line[i] = p.add_glyph(line_glyph)
        line[i].visible = False
        
    legend = Legend(
        items=[("{}".format(i), [line[i]]) for i in range(n)],
        location="top_right", orientation='vertical',
    )
    p.add_layout(legend, 'right')
    
    
    columns = [
            TableColumn(field="identifier", title="ID"),
            TableColumn(field="mask", title="bin"),
            TableColumn(field="label", title="annot"),
            TableColumn(field="x", title="x"),
            TableColumn(field="y", title="y"),
            TableColumn(field="sizes", title="size")
        ]

    table = DataTable(source=s2, columns=columns, width=400, height=280, editable=True)

    
    # Callback for selection
    callback1 = CustomJS(args=dict(s1=s1, s2=s2, table=table), code="""
        var inds = cb_obj.indices;
        var d1 = s1.data;
        var d2 = s2.data;
        d2['x'] = []
        d2['y'] = []
        d2['identifier'] = []
        d2['color'] = []
        d2['label'] = []
        d2['mask'] = []
        d2['sizes'] = []
        for (var i = 0; i < inds.length; i++) {
            d2['x'].push(d1['x'][inds[i]])
            d2['y'].push(d1['y'][inds[i]])
            d2['identifier'].push(d1['identifier'][inds[i]])
            d2['color'].push(d1['color'][inds[i]])
            d2['label'].push(d1['label'][inds[i]])
            d2['mask'].push(d1['mask'][inds[i]])
            d2['sizes'].push(d1['sizes'][inds[i]])
        }
        s2.change.emit();
        table.change.emit();
    """)
    
    s1.selected.js_on_change('indices', callback1)
    
    # Dropdown 
    select = Select(title="Colour by:", value=cont_names[0], options=cont_names)
    
    # Callback for dropdown
    callback2_code = """
    var sel = cb_obj.value;
    var inds = s1.selected.indices;
    var d1 = s.data;
    var d2 = s1.data;
    var d3 = s2.data;
                
    d2['color'] = []
    d2['label'] = []
    d2['mask'] = []
                
    d3['label'] = []
    d3['mask'] = []
    d3['color'] = []

                
    {}
                
    for (var i = 0; i < inds.length; i++) {{
    d3['x'].push(d2['x'][inds[i]])
    d3['y'].push(d2['y'][inds[i]])
    d3['identifier'].push(d2['identifier'][inds[i]])
    d3['color'].push(d2['color'][inds[i]])
    d3['label'].push(d2['label'][inds[i]])
    d3['mask'].push(d2['mask'][inds[i]])
    }}
                
    s1.change.emit();
    s2.change.emit();
    table.change.emit();""".format(make_update_callback(cont_names))
    
    callback2 = CustomJS(args=dict(s=s, s1=s1, s2=s2, table=table), code=callback2_code)
    
    select.js_on_change("value", callback2)

    
    savebutton = Button(label="Save selection (tab delimited)", button_type="success")
    savebutton.js_on_event(ButtonClick, CustomJS(
        args=dict(source_data=s1),
        code="""
            var inds = source_data.selected.indices;
            var data = source_data.data;
            var out = "identifier\tlabel\tsize\tx\ty\\n";
            for (var i = 0; i < inds.length; i++) {
                out += data['identifier'][inds[i]] + "\t" + data['label'][inds[i]] + "\t" + data['sizes'][inds[i]] + "\t" + data['x'][inds[i]] + "\t" + data['y'][inds[i]] + "\\n";
            }
            var file = new Blob([out], {type: 'text/plain'});
            var elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(file);
            elem.download = 'selected-data.txt';
            document.body.appendChild(elem);
            elem.click();
            document.body.removeChild(elem);
            """
            )
                         )

    layout = gridplot([[select], [p, table], [savebutton]])
    output_file(outfile, mode='inline', title="{} {}".format(seqtype.capitalize(), spname))
    save(layout)


# ## Load counts

# In[6]:


counts_tetra_n, sum_counts_tetra =  load_counts_norm(infile)
contig_ids = [line.strip() for line in open(seqidfile)]

# ## Label bins and check annotations and counts have matching lengths

# In[9]:

continuous = []
digitized = []

print(annotfiles)

for j, contfile in enumerate(annotfiles):
    cont = np.array([float(i) for i in open(contfile).read().split()])
    try:
        assert len(cont) == len(contig_ids) == len(counts_tetra_n)
    except:
        print("Error, the length of the annotation vector ({}) must be equal to the number of contigs ({}) and identifiers ({}).\
            ".format(len(cont), len(counts_tetra_n), len(contig_ids)), file=sys.stderr)
        sys.exit(1)

    continuous.append(cont)
    if set(cont) == {0.0, 1.0}:
        print("not binning, already bool")
        n = 2
        digitized.append([int(i) if i == 0. else n-1 for i in cont])
    elif keep_bins == "T":
        print("Keeping input labels")
        digitized.append(cont.astype(int))
        # Get maximum number of labels
        n = len(set(cont))
        #n_labels = len(set(cont))
        #if n < n_labels:
        #    n = n_labels
        print(n)
    else:
        if discretize == "quantile":
            bins = [np.quantile(cont, i) for i in np.linspace(0, 1, n+1)]
            digitized.append(np.digitize(cont, bins[:-1]) - 1)
        else:
            bins = np.linspace(np.min(cont), np.max(cont), n+1)
            digitized.append(np.digitize(cont, bins[:-1]) - 1)

#import colorcet as cc
if keep_bins == "F":
    v = cm.get_cmap('viridis')
    color_key = list(enumerate([matplotlib.colors.rgb2hex(i) for i in v(np.linspace(0,1,n))]))
else:
    v = cm.get_cmap('tab20')
    color_key = list(enumerate(['#c0c0c0'] + [matplotlib.colors.rgb2hex(i) for i in v(np.linspace(0,1,n-1))]))

# Transform counts
if use_pca == "F":
    import umap #Don't load until input passes checks
    decomp = "UMAP"
    transformer_tetra, scaler_tetra = do_umap(counts_tetra_n[:,:])
    contig_tetra_transformed = transform_umap(counts_tetra_n[:,:], transformer_tetra, scaler_tetra)
else:
    contig_tetra_transformed = do_pca(counts_tetra_n)
    decomp = "PCA"


# In[10]:

#TODO: Add option to export "plain" html
if override_x == "":
    x_coord = contig_tetra_transformed[:,0]
else:
    x_coord = np.loadtxt(override_x)

if (override_y == ""):
    y_coord = contig_tetra_transformed[:,1]
else:
    y_coord = np.loadtxt(override_y)

output_notebook()
draw_bokeh_multi(x_coord, y_coord, contig_ids, continuous, digitized, sum_counts_tetra, annotnames, spname, outfile, decomp)

if save_coords is not None:
    print("Writing coordinates to {}".format(save_coords))
    np.save(save_coords, contig_tetra_transformed)

if save_tsv is not None:
    print("Saving tsv file to {}".format(save_tsv))
    np.savetxt(save_tsv, contig_tetra_transformed, fmt="%.6e", delimiter="\t", header="x_{0}\ty_{0}".format(decomp))

