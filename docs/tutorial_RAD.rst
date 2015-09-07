

Example *de novo* RADseq assembly using *pyRAD*
===============================================

pyrad v. 3.1.0 (9/6/2015)
~~~~~~~~~~~~~~~~~~~~~~~~~

About this notebook
~~~~~~~~~~~~~~~~~~~

All cell block in this document should be executed in your command line
(bash) shell, as indicated by the %%bash indicator.

Getting example data
~~~~~~~~~~~~~~~~~~~~

Execute the command below to download an example simulated RADseq data
set and unarchive it into your current directory.


.. code:: python

    %%bash
    wget -q dereneaton.com/downloads/test_rad.zip
    unzip simRADs.zip

The two necessary files below should now be located in your current directory.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  simRADs\_R1\_.fastq.gz : Illumina fastQ formatted reads (gzip
   compressed)
-  simRADs\_barcodes.txt : barcode map file

--------------

The first thing to do is to call *pyrad* from the command line.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be done in one of two ways:

If you installed pyrad following the instructions on the `download
page <https://github.com/dereneaton/pyrad>`__ then you should be able to
call pyrad using the command *pyrad* (lower case) from anywhere on your
machine:

.. code:: python

    %%bash
    pyrad --version

However, if you have problems installing for any reason you can call
pyrad without installation by entering *python* and the full path to
*pyRAD.py*, like this:

.. code:: python

    %%bash
    python /home/deren/pyrad-3.0.52/pyRAD --version

The first step of an analysis is to create a params file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    pyrad -n

The params file lists on each line one parameter followed by a **##**
mark, after which any comments can be left. In the comments section
there is a description of the parameter and in parentheses the step of
the analysis affected by the parameter. Lines 1-12 are required, the
remaining lines are optional. The params.txt file is further described
in the general tutorial.

Let's take a look at the default settings.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    cat params.txt

To change parameters you can edit params.txt in any text editor. Here to automate things I use the script below.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    sed -i '/## 7. /c\12                  ## 7. N processors... ' params.txt
    sed -i '/## 10. /c\.85                ## 10. lowered clust thresh... ' params.txt
    sed -i '/## 14. /c\c85m4p3            ## 14. outprefix... ' params.txt
    sed -i '/## 24./c\8                   ## 24. maxH raised ... ' params.txt
    sed -i '/## 30./c\*                   ## 30. all output formats... ' params.txt

Let's have a look at the changes:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    cat params.txt

--------------

**Let's take a look at what the raw data look like.**

Your input data will be in fastQ format, usually ending in .fq or
.fastq. Your data could be split among multiple files, or all within a
single file (de-multiplexing goes much faster if they happen to be split
into multiple files). The file/s may be compressed with gzip so that
they have a .gz ending, but they do not need to be. The location of
these files should be entered on line 2 of the params file. Below are
the first three reads in the example file.

.. code:: python

    %%bash
    less testrad_R1_.fastq.gz | head -n 12 | cut -c 1-90

--------------

Each read takes four lines. The first is the name of the read (its
location on the plate). The second line contains the sequence data. The
third line is a spacer. And the fourth line the quality scores for the
base calls. In this case arbitrarily high since the data were simulated.

These are 100 bp single-end reads prepared as RADseq. The first six
bases form the barcode and the next five bases (TGCAG) the restriction
site overhang. All following bases make up the sequence data.

--------------

Step 1: de-multiplexing
-----------------------

This step uses information in the barcodes file to sort data into a
separate file for each sample. Below is the barcodes file, with sample
names and their barcodes each on a separate line with a tab between
them.

.. code:: python

    %%bash
    cat testrad_barcodes.txt

Step 1 writes the de-multiplexed data to a new file for each sample in a
new directory created within the working directory called fastq/.

.. code:: python

    %%bash
    pyrad -p params.txt -s 1

You can see that this created a new file for each sample in the
directory 'fastq/'

.. code:: python

    %%bash
    ls fastq/

The statistics for step 1
^^^^^^^^^^^^^^^^^^^^^^^^^

A new directory called stats will also have been created. Each step of
the *pyRAD* analysis will create a new stats output file in this
directory. The stats output for step 1 is below:

.. code:: python

    %%bash
    cat stats/s1.sorting.txt

Step 2: quality filtering
~~~~~~~~~~~~~~~~~~~~~~~~~

This step filters reads based on quality scores, and can be used to
detect Illumina adapters in your reads, which is sometimes a problem
with homebrew type library preparations. Here the filter is set to the
default value of 0, meaning it filters only based on quality scores of
base calls. The filtered files are written to a new directory called
edits/.

.. code:: python

    %%bash
    pyrad -p params.txt -s 2

.. code:: python

    %%bash
    ls edits/

The filtered data are written in fasta format (quality scores removed)
into a new directory called edits/. Below I show a preview of the file
which you can view most easily using the ``less`` command (I use
``head`` here to make it fit in the text window better).

.. code:: python

    %%bash
    head -n 10 edits/1A_0.edit | cut -c 1-80

Step 3: clustering within-samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 3 de-replicates and then clusters reads within each sample by the
set clustering threshold and writes the clusters to new files in a
directory called clust.xx

.. code:: python

    %%bash
    pyrad -p params.txt -s 3

Once again, I recommend you use the unix command 'less' to look at the
clustS files. These contain each cluster separated by "//". For the
first few clusters below you can see that there is one or two alleles in
the cluster and one or a few reads that contained a (simulated)
sequencing error.

.. code:: python

    %%bash
    less clust.85/1A_0.clustS.gz | head -n 26 | cut -c 1-80

--------------

The stats output tells you how many clusters were found, and their mean
depth of coverage. It also tells you how many pass your minimum depth
setting. You can use this information to decide if you wish to increase
or decrease the mindepth before it is applied for making consensus base
calls in steps 4 & 5.

.. code:: python

    %%bash
    head -n 55 stats/s3.clusters.txt

Steps 4 & 5: Call consensus sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 4 jointly infers the error-rate and heterozygosity across samples.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    pyrad -p params.txt -s 4

.. code:: python

    %%bash
    less stats/Pi_E_estimate.txt

Step 5 calls consensus sequences using the parameters inferred above, and filters for paralogs.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    pyrad -p params.txt -s 5

The stats output for step 5
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    %%bash
    less stats/s5.consens.txt

Step 6: Cluster across samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 6 clusters consensus sequences across samples. It will print its
progress to the screen. This uses 6 threads by default. If you enter 0
for param 37 it will use all available processors.

.. code:: python

    %%bash
    pyrad -p params.txt -s 6 

Step 7: Assemble final data sets
--------------------------------

The final step is to output data only for the loci that you want to have
included in your data set. This filters once again for potential
paralogs or highly repetitive regions, and includes options to minimize
the amount of missing data in the output.

.. code:: python

    %%bash
    pyrad -p params.txt -s 7

Final stats output
~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    less stats/c85m4p3.stats

--------------

Output formats
--------------

We created 11 output files from our analysis. The standard two (.loci
and .excluded\_loci), as well as the 9 additional ones listed in the
params file. These are all shown below.

.. code:: python

    %%bash 
    ls outfiles/

Loci format
~~~~~~~~~~~

The ".loci" file contains each locus listed in a fasta-like format that
also shows which sites are variable below each locus. Autapomorphies are
listed as '-' and shared SNPs as '\*'. This is a custom format that is
human readable and also used as input to perform D-statistic tests in
pyRAD. This is the easiest way to visualize your results. I recommend
viewing the file with the command ``less``. Below I use a head and cut
to make it easy to view in this window.

.. code:: python

    %%bash 
    head -n 39 outfiles/c85m4p3.loci | cut -c 1-85

PHY format
~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.phy | cut -c 1-85

NEX format
~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.nex | cut -c 1-85

Alleles format
~~~~~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.alleles| cut -c 1-85

STRUCTURE (.str) format
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.str | cut -c 1-20

GENO (.geno) format (used in *Admixture*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 40 outfiles/c85m4p3.usnps.geno 

SNPs format
~~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.snps | cut -c 1-85

UNLINKED\_SNPs format
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash 
    head -n 50 outfiles/c85m4p3.unlinked_snps | cut -c 1-85

OTHER FORMATS
-------------

You may also produce some more complicated formatting options that
involve pooling individuals into groups or populations. This can be done
for the "treemix" and "migrate" outputs, which are formatted for input
into the programs *TreeMix* and *migrate-n*, respectively. Grouping
individuals into populations is done with the final lines of the params
file as shown below, and similar to the assignment of individuals into
clades for hierarchical clustering (see full tutorial).

Each line designates a group, and has three arguments that are separated
by space or tab. The first is the group name, the second is the minimum
number of individuals that must have data in that group for a locus to
be included in the output, and the third is a list of the members of
that group. Lists of taxa can include comma-separated names and wildcard
selectors, like below. Example:

.. code:: python

    %%bash 
    ## append group designations to the params file
    echo "pop1 4 1A0,1B0,1C0,1D0 " >> params.txt
    echo "pop2 4 2E0,2F0,2G0,2H0 " >> params.txt
    echo "pop3 4 3* " >> params.txt
    
    ## view params file
    cat params.txt

Creating population output files
--------------------------------

Now if we run *pyRAD* with the 'm' (migrate) or 't' (treemix) output
options, it will create their output files.

.. code:: python

    %%bash 
    pyrad -p params.txt -s 7

TREEMIX format
--------------

.. code:: python

    %%bash 
    less outfiles/c85m4p3.treemix.gz | head -n 30

MIGRATE-n FORMAT
----------------

.. code:: python

    %%bash 
    head -n 40 outfiles/c85m4p3.migrate | cut -c 1-85
