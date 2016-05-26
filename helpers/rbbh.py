# ex02.py
#
# Functions and data useful in exercise 2 of the Plant and Pathogen
# Bioinformatics course at EMBL

from matplotlib.colors import LogNorm

from Bio import SeqIO

import matplotlib  # to get version
import matplotlib.pyplot as plt

import pandas as pd

import os

# GLOBALS
datadir = "pseudomonas_blastp/rbbh_output"

# Function to read BLAST data
def read_data(filename):
    """Reads BLASTP tabular output data from the file corresponding
    to passed accessions, returning a Pandas dataframe.
    """
    # We specify no index for these dataframes
    df = pd.DataFrame.from_csv(filename,
                               sep="\t", index_col=None)
    df.columns=['query_id', 'subject_id', 'query_length', 'subject_length', 
                'alignment_length', 'identical_sites', 'identity',
                'qcovs', 'Evalue', 'bitscore']
    df = calculate_coverage(df)
    return df
        
# Function to calculate query and subject coverage for a BLAST dataset
def calculate_coverage(d):
    """For the passed dataframe, calculates query and subject coverage,
    and returns the original dataframe with two new columns,
    query_coverage and subject_coverage.
    """
    d['query_coverage'] = 100 * d.alignment_length/d.query_length
    d['subject_coverage'] = 100 * d.alignment_length/d.subject_length
    return d

# Function to filter dataframe on percentage identity and coverage
def filter_cov_id(pid, cov, *df):
    """Filters the passed dataframe, returning only rows that meet
    passed percentage identity and coverage criteria.
    """
    return tuple([d[(d.identity > pid) & (d.query_coverage > cov) &
                    (d.subject_coverage > cov)] 
                    for d in df])
    
    
# Function to filter dataframe to best HSP match only, for any query
def filter_matches(*df):
    """Filters rows duplicated by query_id with a 
    Pandas DataFrame method: drop_duplicates. By default, this
        keeps the first encountered row, which is also the "best"
        HSP in BLAST tabular output.
    """
    return tuple([d.drop_duplicates(subset='query_id') for d in df])
    
# Function to return a dataframe of RBBH from two dataframes of BLAST matches
def find_rbbh(df1, df2, pid=0, cov=0):
    """Takes the accessions for two organisms, and 
    1. parses the appropriate BLASTP output into two dataframes, and 
       calculates coverage, discarding all but the "best" HSP, on the
       basis of bitscore, where there are multiple HSPs for a hit.
    2. cuts out rows that do not meet minimum percentage identity and
       coverage criteria.
    3. identifies RBBH from the remaining BLAST matches.
    """
    assert 0 <= pid <= 100, \
        "Percentage identity should be in [0,100], got %s" % str(pid)
    assert 0 <= cov <= 100, \
        "Minimum coverage should be in [0,100], got %s" % str(cov)
    # Filter to best HSP match only
    df1, df2 = filter_matches(df1, df2)
    # Filter datasets on minimum percentage identity and minimum coverage
    df1, df2 = filter_cov_id(pid, cov, df1, df2)
    # Identify RBBH from the two datasets on the basis of matching query and 
    # subject sequence ID
    matches = df1.merge(df2, left_on=("query_id", "subject_id"),
                        right_on=("subject_id", "query_id"))
    rbbh = matches[["query_id_x", "subject_id_x", "identity_x",
                    "identity_y", "query_coverage_x",
                    "query_coverage_y", "subject_coverage_x",
                    "subject_coverage_y", "bitscore_x", "bitscore_y",
                    "Evalue_x", "Evalue_y"]]
    #rbbh.to_csv(os.path.join(datadir, "rbbh_%s_vs_%s.tab" % (s1, s2)),
    #            sep='\t')
    return df1, df2, rbbh

# Function to plot 2D histogram, with colorbar scale, from
# two dataframe columns.
def plot_hist2d(col1, col2, xlab="", ylab="", header="", bins=100):
    """Plots a 2D histogram of the two passed dataframe columns."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # How we render 2D histograms/heatmaps varies, dependent on the 
    # Matplotlib version
    if float('.'.join(matplotlib.__version__.split('.')[:2])) > 1.2:
        c, x, y, f = ax.hist2d(list(col1), list(col2),
                               bins=bins, norm=LogNorm())
        fig.colorbar(f)
    else:
        h, x, y = histogram2d(list(col2), list(col1), bins=bins)
        im = ax.imshow(h, extent=[x[0], x[-1], y[0], y[-1]],
                       aspect="equal", origin="lower",
                       interpolation="nearest",
                       norm=LogNorm(vmin=1, vmax=len(col1)))
        fig.colorbar(im)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(header)
    
# Function for scatterplot from two dataframe columns.
def plot_scatter(col1, col2, xlab="", ylab="", header="",
                 origin_zero=False):
    """Plots a scatterplot of the two passed dataframe columns,
    colouring points by y-axis value.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    f = ax.scatter(col1, col2, c=col2,
                   cmap=plt.get_cmap("hot"))
    if origin_zero:
        ax.set_xlim([0, col1.max()])
        ax.set_ylim([0, col2.max()])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(header)


# Function to process GenBank files into a dictionary of CDS features,
# keyed by protein ID, where the values are a tuple of (source, start, end, 
# strand) information.
def read_genbank(*filenames):
    """Returns a dictionary of CDS annotations, where the dictionary keys
    are a tuple of CDS protein ID accession numbers (which may not be unique
    to a single genome, in RefSeq!) with the source filename stem, and the
    values are (source, start, end, strand, id) information about the CDS
    location on the chromosome.
        
    - *filenames, the organism's GenBank annotation files
    """
    ft_dict = {}
    for filename in filenames:
        filestem = os.path.splitext(os.path.split(filename)[-1])[0]
        with open(filename, 'rU') as fh:
            record = SeqIO.read(fh, 'genbank')
            # Reconstruct the name in the corresponding .fna file
            record_name = '|'.join(["gi", record.annotations['gi'],
                                    "ref", record.id])
            for ft in [f for f in record.features if f.type == "CDS" and 
                       "protein_id" in f.qualifiers]:
                ft_dict[(ft.qualifiers['protein_id'][0], filestem)] = \
                    (record_name, int(ft.location.start), 
                     int(ft.location.end), ft.location.strand)
    print("Loaded %d features" % len(ft_dict))
    return ft_dict

# Function to split a full sequence reference ID into only the last value
def split_seqid(seqid):
    """Split NCBI sequence ID to get last value."""
    if '|' not in seqid:
        return seqid
    return seqid.split('|')[-2]

# Function to write a single line of a Pandas RBBH dataframe to file
def write_line(line, features, fh, fwd, rev):
    """Write a single RBBH dataframe line to a file."""
    try:
        ft1 = features[(split_seqid(line['query_id_x']), fwd)]
    except:
        print("Could not find feature for %s (skipping)" %
              line['query_id_x'])
        return
    try:
        ft2 = features[(split_seqid(line['subject_id_x']), rev)]
    except:
        print("Could not find feature for %s (skipping)" %
              line['subject_id_x'])
        return
    fh.write(' '.join([str(line['bitscore_x']),
                       str(line['identity_x']),
                       str(ft1[2]) if ft1[3] < 0 else str(ft1[1]),
                       str(ft1[1]) if ft1[3] < 0 else str(ft1[2]),
                       str(ft1[0]),
                       str(ft2[2]) if ft2[3] < 0 else str(ft2[1]),
                       str(ft2[1]) if ft2[3] < 0 else str(ft2[2]),
                       str(ft2[0])
                   ]) + '\n')

# Function to write .crunch output for ACT visualisation, from the
# RBBH identified above
def write_crunch(rbbh, features, fwd, rev, outdir=".", filename="rbbh.crunch"):
    """Writes .crunch output in outdir, for those RBBH stored in the 
    passed dataframe

    fwd - the filestem for the 'forward' sequence
    rev - the filestem for the 'reverse' sequence
    """
    with open(os.path.join(outdir, filename), 'w') as fh:
        rbbh.apply(write_line, axis=1, args=(features, fh, fwd, rev))
    print("Wrote file to %s" % os.path.join(outdir, filename))
