#!/usr/bin/env python
# <cw21@sanger.ac.uk>
import argparse
import os

##################################
parser = argparse.ArgumentParser(
    description='Produce 2D representation of a multidimensional numpy array')

parser.add_argument(
    "--countfile", help="2D numpy array with tetranucleotide counts", required=True)
parser.add_argument("--verbose", help="print extra model info", default="F", choices=["T", "F"])
parser.add_argument("--fignames", help="add file with figure names", default="sample")
parser.add_argument("--kl", help="kl weight, corresponds to beta in a beta-VAE (beta == 1 corresponds to a standard VAE)", default=0.0025, type=float)
parser.add_argument("--extra_columns", help="Additional input feature data to load into matrix", default=None)
parser.add_argument("--batch_size", help="Batch size", default=256, type=int)
parser.add_argument("--epochs", help="Number of epochs", default=10, type=int)
parser.add_argument("--outdir", help="Output directory", default="")
parser.add_argument("--model_out", help="Path to save model weights to", default="")
parser.add_argument("--intermediate", help="Draw embeddings for intermediate steps?", default="F", choices=["T", "F"])
parser.add_argument("--plot_latents", help="Plot the latent space after training is complete (enabled by default)", default="T", choices=["T", "F"])
parser.add_argument("--remove", help="File specifying reads to remove before decomposition (0: keep; other integer value: drop; one integer per line)", default=None)
parser.add_argument("--latent_dim", help="Number of latent dimensions", default=2, type=int)
parser.add_argument("--reduce_beta", help="Reduce weight on KL loss (beta) if variance indicates posterior collapse (enabled by default)", default="T", choices=["T", "F"])
parser.add_argument("--learning_rate", help="Set learning rate", type=float, default=0.001)
parser.add_argument("--latent_format", help="Format to store latent reoresentations in (defaults to .txt, but .npy preferable for large datasets)", default="txt", choices=["txt", "npy"])
parser.add_argument("--cb", help="Apply CB correction (enabled by default)", default="T", choices=["T", "F"])

args = parser.parse_args()
print(args)

countfile = args.countfile
latent_dim = args.latent_dim

modelpath = args.model_out
verbose = args.verbose
figures = args.fignames
extra_columns = args.extra_columns
outdir = os.path.join(args.outdir, '')
kl = args.kl
epochs = args.epochs
batch_size =  args.batch_size
draw_intermediate = args.intermediate
remove = args.remove
reduce_beta = args.reduce_beta
learning_rate = args.learning_rate
plot_latents = args.plot_latents
latent_format = args.latent_format
cb = args.cb

#TODO: Check output dirs exist

print("Loading dependencies")
#################################
# imports
import gc

from sklearn.preprocessing import MinMaxScaler
from keras import backend as K
from keras import regularizers
from tensorflow.keras import layers
from tensorflow import keras
import tensorflow as tf
from tqdm.keras import TqdmCallback
import numpy as np

# Fix random seed
np.random.seed(42)
tf.random.set_seed(42)

# Load plotting extensions if needed
if (plot_latents == "T") | (draw_intermediate == "T"):
    import pandas as pd
    import datashader as ds
    from holoviews.operation.datashader import datashade, dynspread
    import holoviews as hv
    hv.extension('bokeh')

#from memory_profiler import profile

#################################
# Load data

#@profile
def load_counts_norm(infile, remove):
    """Load raw counts from .npy, divide by sums of rows to account for read length"""
    counts = np.load(infile).astype("float32")
    if remove is not None:
        counts = filter_rows(counts, remove)
    #sum_all = counts.sum(axis=1)[:, None]
    counts /= counts.sum(axis=1)[:, None] # Specify float32 to avoid issues with tf
    return counts

#@profile
def scale_counts(counts):
    """MinMax scale array"""
    scaler = MinMaxScaler(copy=False) # Scale in-place without copying
    counts = scaler.fit_transform(counts)
    return counts

def filter_rows(counts, remove):
    """Filter rows if requested"""
    keep_lines = np.array([int(i.strip()) for i in open(remove)])
    rows_to_keep = np.where(keep_lines == 0)[0]
    counts = counts[rows_to_keep]
    print("Filtered reads using {}, kept {} rows".format(remove, counts.shape[0]))
    return counts

#@profile
def load_all_scale(countfile, remove, extra_columns):
    """Open normalised counts, add extra columns if needed, scale data"""
    count_array = load_counts_norm(countfile, remove)
    if extra_columns is not None:
        new_cols = np.hstack([np.array(open(a).read().split(), dtype="float32").T[:, None] for a in extra_columns.split()])
        count_array = np.hstack((count_array, new_cols))
    count_array = scale_counts(count_array) # check memory increment
    return count_array

print("Loading data")
transformed_counts = load_all_scale(countfile, remove, extra_columns)

gc.collect()

##################################
# Config
original_dim = transformed_counts.shape[1]
n_rows = transformed_counts.shape[0]
intermediate_dim = 24
initializer = tf.keras.initializers.GlorotUniform(seed=42) #HeNormal()

print("Input has {} columns and {} rows".format(original_dim, n_rows))

##################################
# Use tf.dataset (passing the numpy array to model.fit() appears to cause increased memory usage in some instances)
shuffle = np.arange(n_rows)
np.random.shuffle(shuffle)
transformed_counts = transformed_counts[shuffle,:]
transformed_counts =  tf.data.Dataset.from_tensor_slices((transformed_counts)).batch(batch_size)

# shuffle in separate step to preserve order of sequences in output
# This was a hack to accommodate larger datasets
# Utils to fetch training and validation data
# x % 5 corresponds to val split 0.2
def is_test(x, _):
    return x % 5 == 0

def is_train(x, y):
    return not is_test(x, y)

recover = lambda x, y: y

# Split the dataset for training.
test_dataset = transformed_counts.enumerate() \
    .filter(is_test) \
    .map(recover)

# Split the dataset for testing/validation.
train_dataset = transformed_counts.enumerate() \
    .filter(is_train) \
    .map(recover)


##################################
# VAE functions


class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


# Encoder network
#@profile
def build_encoder(latent_dim, original_dim, intermediate_dim):
    """Build encoder network"""
    encoder_inputs = keras.Input(shape=(original_dim, ))
    x = layers.Dense(intermediate_dim, activation="relu", kernel_initializer=initializer)(encoder_inputs)
    x = keras.layers.Dropout(.2)(x)
    x = keras.layers.BatchNormalization()(x)
    x = layers.Dense(intermediate_dim/2, activation="relu", kernel_initializer=initializer)(x)
    x = layers.Dense(intermediate_dim/4, activation="relu",
                     activity_regularizer=regularizers.l1(10e-4), kernel_initializer=initializer)(x)
    z_mean = layers.Dense(latent_dim, name="z_mean")(x)
    z_log_var = layers.Dense(latent_dim, name="z_log_var")(x)
    z = Sampling()([z_mean, z_log_var])
    return keras.Model(
        encoder_inputs, [z_mean, z_log_var, z], name="encoder")


encoder = build_encoder(latent_dim, original_dim, intermediate_dim)


######
# Decoder

#@profile
def build_decoder(latent_dim, original_dim, intermediate_dim):
    """build decoder network"""
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(intermediate_dim/4 , activation="relu", kernel_initializer=initializer)(latent_inputs)
    x = layers.Dense(intermediate_dim/2 , activation="relu", kernel_initializer=initializer)(x)
    x = keras.layers.Dropout(.2)(x)
    x = keras.layers.BatchNormalization()(x)
    x = layers.Dense(intermediate_dim, activation="relu", kernel_initializer=initializer)(x)
    decoder_outputs = layers.Dense(original_dim, activation="sigmoid")(x)
    return keras.Model(latent_inputs, decoder_outputs, name="decoder")


decoder = build_decoder(latent_dim, original_dim, intermediate_dim)

if (verbose == "T"):
    print(encoder.summary())
    print(decoder.summary())

####
# Loss function customisation
# "lam" is lambda, i.e. BC(lambda)

# TODO: Fix docstring

def cont_bern_log_norm(lam, l_lim=0.49, u_lim=0.51):
    """Calculates correction term for CB BCE loss
    adapted from https://github.com/cunningham-lab/cb_and_cc
    """
    lam = tf.clip_by_value(lam, 1e-4, 1-1e-4)  # moved clipping code in here
    # computes the log normalizing constant of a continuous Bernoulli distribution in a numerically stable way.
    # returns the log normalizing constant for lam in (0, l_lim) U (u_lim, 1) and a Taylor approximation in
    # [l_lim, u_lim].
    # cut_y below might appear useless, but it is important to not evaluate log_norm near 0.5 as tf.where evaluates
    # both options, regardless of the value of the condition.
    cut_lam = tf.where(tf.logical_or(tf.less(lam, l_lim), tf.greater(
        lam, u_lim)), lam, l_lim * tf.ones_like(lam))
    log_norm = tf.math.log(tf.abs(
        2.0 * tf.math.atanh(1 - 2.0 * cut_lam))) - tf.math.log(tf.abs(1 - 2.0 * cut_lam))
    taylor = tf.math.log(2.0) + 4.0 / 3.0 * tf.pow(lam -
                                                   0.5, 2) + 104.0 / 45.0 * tf.pow(lam - 0.5, 4)
    return tf.where(tf.logical_or(tf.less(lam, l_lim), tf.greater(lam, u_lim)), log_norm, taylor)


def cb_loss(data, reconstructed):
    """Calculates CB-corrected BCE loss"""
    return tf.reduce_mean(
        K.binary_crossentropy(data, reconstructed)
        + cont_bern_log_norm(reconstructed), 1)

def ce_loss(data, reconstructed):
    """Calculates CE loss without correction"""
    return tf.reduce_mean(K.binary_crossentropy(data, reconstructed), 1)


def reconstruction_bce_kl(data, reconstruction, z_mean, z_log_var):
    """Calculate binary cross entropy with CB correction term and KL loss"""
    if cb == "T":
        reconstruction_loss = tf.reduce_mean(cb_loss(data, reconstruction))
    else:
        reconstruction_loss = tf.reduce_mean(ce_loss(data, reconstruction))
    reconstruction_loss *= original_dim
    kl_loss = 1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var)
    kl_loss = tf.reduce_mean(kl_loss)
    kl_loss *= -0.5
    return reconstruction_loss, kl_loss

#### Callbacks

# Reduce learning rate if validation loss does not improve
reduce_lr = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.2,
                                                 patience=3, min_lr=0.00001)

# Stop if validation loss fails to improve
early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=4, restore_best_weights=True)

## Draw intermediate steps
def scatter(epoch):
    """Draw scatter plot of latent space"""
    if (draw_intermediate != "F") & (latent_dim == 2):
        x_test_encoded = vae.encoder.predict(transformed_counts)
        df = pd.DataFrame(data=x_test_encoded[0],columns=['x', 'y'])
        points = hv.Points(data=df, kdims=['x','y'])
        plot = datashade(points).opts(width=700, height=700)
        hv.save(plot, f"{outdir}{figures}_epoch_{epoch}.png")
    return 0


#tf.keras.callbacks.TerminateOnNaN()
class ScatterCallback(tf.keras.callbacks.Callback):
    """Callback to plot intermediate steps after each epoch"""
    def on_epoch_end(self, epoch, logs=None):
        scatter(epoch)
    def on_train_begin(self, epoch, logs=None):
        scatter("start")

# Callback not in use? 
class MeanVar(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs=None):
        """At the end of each epoch, check for latent collapse.
        If variance is high, drop the weight on the KL loss by a specified multiplier (beta)."""
        x_test_encoded = vae.encoder.predict(transformed_counts)
        mu = tf.math.reduce_mean(x_test_encoded[0], axis=0)
        v = tf.math.reduce_mean(tf.math.exp(x_test_encoded[1]), axis=0)
        print("mu: {}, v: {}".format(mu, v))
        print(keras.backend.get_value(v))
        if np.any(v > 0.99):
            tf.print("Dropping KL weight *0.1 due to high variance in Z")
            weight = keras.backend.get_value(vae.kl_weight)
            tf.keras.backend.set_value(vae.kl_weight, weight*0.1)
        elif np.any(v > 0.75):
            tf.print("Dropping KL weight due to high variance in Z")
            weight = keras.backend.get_value(vae.kl_weight)
            tf.keras.backend.set_value(vae.kl_weight, weight*0.5)
        elif np.any(v > 0.5):
            tf.print("Dropping KL weight *0.1 due to high variance in Z")
            weight = keras.backend.get_value(vae.kl_weight)
            tf.keras.backend.set_value(vae.kl_weight, weight*0.75)

class VAE(keras.Model):
    """Fits VAE with CB-corrected BCE and KL loss"""
    #@profile
    def __init__(self, encoder, decoder, original_dim, kl_weight, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        original_dim = original_dim
        self.kl_weight = tf.Variable(
            kl_weight, trainable=False, name='kl_weight', dtype=tf.float32)
    
    #@profile
    def train_step(self, data):
        if isinstance(data, tuple):
            data = data[0]
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = encoder(data)
            reconstruction = decoder(z)
            reconstruction_loss, kl_loss = reconstruction_bce_kl(
                data, reconstruction, z_mean, z_log_var)
            total_loss = reconstruction_loss + kl_loss * self.kl_weight
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        gc.collect()
        return {
            "loss": total_loss,
            "reconstruction_loss": reconstruction_loss,
            "kl_loss": kl_loss,
            "kl_weight": self.kl_weight,
        }

    def test_step(self, data):
        if isinstance(data, tuple):
            data = data[0]
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = encoder(data, training=False)
            reconstruction = decoder(z, training=False)
            reconstruction_loss, kl_loss = reconstruction_bce_kl(
                data, reconstruction, z_mean, z_log_var)
            total_loss = reconstruction_loss + kl_loss * self.kl_weight
            gc.collect()
        return {
            "loss": total_loss,
            "reconstruction_loss": reconstruction_loss,
            "kl_loss": kl_loss,
            "kl_weight": self.kl_weight,
        }



#@profile
def compile_model():
    vae = VAE(encoder, decoder, original_dim, kl)
    vae.compile(optimizer=keras.optimizers.Adam(lr=learning_rate))
    return vae

vae = compile_model()

# TODO keep logs of training if requested
#@profile
callbacks = [reduce_lr, early_stop, TqdmCallback(verbose=0)]
if draw_intermediate == "T":
    callbacks.append(ScatterCallback())
if reduce_beta == "T":
    callbacks.append(MeanVar())

def train():
    print("Start training")
    history = vae.fit(train_dataset, validation_data=test_dataset, batch_size=batch_size, epochs=epochs, callbacks=callbacks, verbose=0, shuffle=True)
    return history

history = train()
print(history.history)

#FIXME
if (modelpath != ""):
    print("Save model weights")
    encoder.save(f'{modelpath.strip("/")}_encoder')
    decoder.save(f'{modelpath.strip("/")}_decoder')

print("Collect results")

x_test_encoded = encoder.predict(transformed_counts)

# Output embeddings
print("Save latent embeddings")
for i in range(3):
    if latent_format == "txt":
        np.savetxt(f"{outdir}{figures}.vae.out.{latent_dim}d.{i}", x_test_encoded[i][np.argsort(shuffle),:], delimiter="\t")
    elif latent_format == "npy":
        np.save(f"{outdir}{figures}.vae.out.{latent_dim}d.{i}.npy", x_test_encoded[i][np.argsort(shuffle),:])


# Plot embeddings if requested, and if latent dimensions == 2
if (latent_dim == 2) & (plot_latents == "T"):
    print("Plotting 2D representations.")
    import colorcet as cc
    for i in range(3):
        x_test_encoded[i]
        df = pd.DataFrame(data=x_test_encoded[i], columns=['x', 'y'])
        points = hv.Points(data=df, kdims=['x', 'y'])
        layout = dynspread(datashade(points, cmap=cc.fire).opts(bgcolor='black', width=1500, height=1500))
        hv.save(layout,f'{outdir}2d_plot_{figures}_{i}.png', fmt='png')
