# Singularity Example Scientific Workflow

This repository provides a SCI-F (Standard Container Integration Format) Apps example to run a scientific workflow. This analysis was originally used to compare Singularity vs. Docker on cloud providers via [this container](https://github.com/vsoch/singularity-scientific-example) with [interesting results](https://vsoch.github.io/singularity-scientific-example/results/)! If you are interested in a Docker vs. Singularity implementation, see that project. For this small example, we want to give rationale for taking a SCI-F apps approach over a traditional Singularity image. We compare the following equivalent (but different!) implementations:

- [Singularity (without SCI-F)](Singularity.noscif)
- [Singularity (with SCI-F)](Singularity)

The original image is served on Singularity Hub:

```
singularity pull shub://vsoch/singularity-scientific-example
```

And the image with SCI-F app supports must be built with (not yet released) Singularity 2.4:

```
sudo singularity build scif.img Singularity
```

## Evaluation
For complete methods and results, see the written section of the corresponding paper (TBA). We can briefly state here that we reproduced results across the Docker and Singularity runs.  The aim here is to qualitatively evaluate SCI-F on its ability to expose container metadata, and information about the pipeline and executables inside. As each evaluation is scoped to the goal's of the container, for this example we focus on the purpose of deploying a set of steps that encompass a pipeline to perform variant calling. 

First, the goals of SCI-F:

- Containers are **consistent** to allow for comparison. I am able to easily discover relevant software and data for one or more applications defined by the container creator.
- Containers are transparent. If i discover a container and do not have any prior knowledge or metadata, the important executables and metadata are revealed to me.
- Container contents are easily available for **introspection** because they are programmatically parseable. I can run a function over a container, and know exactly the functions available to me, ask for help, and know where to interact with inputs and outputs.
 - Container internal infrastructure is **modular**. Given a set of SCI-F apps from different sources, I can import different install routines and have assurance that environment variables defined for each are sourced correctly for each, and that associated content does not overwrite previous content. Each software and data module must carry, minimally, a unique name and install location in the system.

For each of the above, let's define some basic tests.

### 1. Development Evaluation
For this use case, we are a container developer, and we are writing a singularity build recipe.

#### Can I easily define multiple entrypoints?

**Standard**
Singularity standard defaults to a single runscript per container. If I need to expose multiple functions, I either need to write a more complicated single entrypoint to direct the call, or I need to write detailed docs on different exec commands (and executables inside) for the user to run. For this real world use case, at the time when this runscript was written, SCI-F was not yet conceptualized. Given the sheer number of tools in the container, the runscript served to list a long list of executables, written to a text file, added to the path:

```
%runscript

    if [ $# -eq 0 ]; then
        echo "\nThe following software is installed in this image:"
        column -t /Software/.info | sort -u --ignore-case        
        echo "\Note that some Anaconda in the list are modules and note executables."
        echo "Example usage: analysis.img [command] [args] [options]"  
    else
        exec "$@"
    fi
```

Given that this container was intended to run a scientific workflow, this list doesn't help to make its usage transparent. It would be useful for an individual familiar with the core software, perhaps developing a specific workflow with it. Arguably, this complete listing should be provided, perhaps not as a main entrypoint, but a separate app to print metadata, or the software names and versions added as labels to the app(s) where they are relevant. It's important to note that this first container did not include any logic for performing the analysis, this was provided by the included [scripts](scripts). If the scripts are separated from the container, reproducing the analysis is no longer possible.

**SCI-F**
SCI-F has the advantage here because I can give names to different apps, and write a different executable runscript for each. I can still provide a main runscript, and it could either give the user a listing of possible options (below) or run the entire workflow the container provides. Here is a main runscript that instructs the user on how to use the apps, dynamically generating the list:

```
%runscript

    if [ $# -eq 0 ]; then
        echo "\nThe following software is installed in this image:"
        ls /scif/apps | sort -u --ignore-case        
        echo "Example usage: singularity --app <name> <container> [command] [args] [options]"  
    else
        exec "$@"
    fi

```

And here is an example app, specifically to downlad a component of the data:

```
%apprun download-reference
    mkdir -p $REF_DIR
    wget -P $REF_DIR ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.transcripts.fa.gz
    gzip -d $REF_DIR/gencode.v25.transcripts.fa.gz
    wget -P $REF_DIR ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gzip -d $REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Importantly, the scripts that were previously in [scripts](scripts) are now cleanly organized into sections in the build recipe. This modular organization and easy accessibility would have been very challening given the first container organization. The runscript would have needed to be long and complicated to infer what the user wanted to do, or the same functionality achieved by executing different scripts (inside the container), which is a non-trivial (or minimally more annoying to write) than simply writing lines into sections.


#### Can I easily install content known to modules?
Given that I have two functions for my container to perform, foo and bar, can I generate an install routine that will allow for both shared dependencies (i.e. global) and distinct dependencies?

**Standard**
Singularity standard has one mode to install global dependencies, everything goes into the `%post` script and any organization of required data, files, and libraries is scattered around the image. Other than coming up with a manual organization and documenting it, there is no way to cleanly define boundaries that will be easily discovered by the user. If you take a look at the [Standard Singularity](Singularity.noscif) recipe, you will see this reflected in one huge, single install procedure. As an example, for this container the tools `bwa` and `samtools` were generally used for an alignment step, and there is no way of knowing this. They are installed globally:

```
cd /Software
su -c 'git clone https://github.com/Linuxbrew/brew.git' singularity
su -c '/Software/brew/bin/brew install bsdmainutils parallel util-linux' singularity
su -c '/Software/brew/bin/brew tap homebrew/science' singularity
su -c '/Software/brew/bin/brew install art bwa samtools' singularity
su -c 'rm -r $(/Software/brew/bin/brew --cache)' singularity
``` 

In fact, the `art` tools are installed with the same manager (`brew`), but they belong to an entirely different step. If a research scientist (or user) were parsing this build recipe, or using an NLP algorithm that took distance into account, there would be no strong signal about how these software were used or related to one another.

**SCI-F**
With SCI-F, by simply defining an environment, labels, install, or runscript to be in the context of an app, the modularity is automatically generated. When I add a list of files to an app `foo`, I know they are added to the container's predictable location for `foo`. If I add a file to `bin` I know it goes into foo's bin, and is added to the path when `foo` is run. If I add a library to `lib`, I know it is added to `LD_LIBRARY_PATH` when foo is run. I don't need to worry about equivalently named files under different apps getting mixed up, or being called incorrectly because both are on the path.  For example, in writing these sections, a developer can make it clear that `bwa` and `samtools` are used together for alignment:

```
# =======================
# bwa index and align
# =======================

%appinstall bwa-index-align
    git clone https://github.com/lh3/bwa.git build
    cd build && git checkout v0.7.15 && make
    mv -t ../bin bwa bwakit

    apt-get install -y liblzma-dev
    cd .. && wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
    tar -xvjf samtools-1.5.tar.bz2
    cd samtools-1.5 && ./configure --prefix=${SINGULARITY_APPROOT}
    make && make install

%apprun bwa-index-align
    mkdir -p $DATADIR/Bam
    bwa index -a bwtsw $DATADIR/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    bwa mem -t $NUMCORES $DATADIR/Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa $DATADIR/Fastq/dna_1.fq.gz $DATADIR/Fastq/dna_2.fq.gz | samtools view -bhS - > $DATADIR/Bam/container.bam  

%applabels bwa-index-align
    bwa-version v0.7.15
    samtools-version v1.5
```

and that `art` is used to simulate reads:

```
# =======================
# simulate reads
# =======================

%apphelp simulate-reads
    Optionally set any of the following environment variables (defaults shown)
    READS (100000000)
    READ_LEN (150)
    GENOME_SIZE (3400000000)

%appenv simulate-reads
    READS=${READS:-100000000}
    READ_LEN=${READ_LEN:-150}
    GENOME_SIZE=${GENOME_SIZE:-3400000000}
    export GENOME_SIZE READ_LEN READS

%appinstall simulate-reads   
    wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
    tar -xzvf artbinmountrainier20160605linux64tgz.tgz 
    mv art_bin_MountRainier/* bin/
    chmod u+x bin/art_*

%apprun simulate-reads
    GENOME="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    FOLD_COVERAGE=$(python -c "print($READS*$READ_LEN/$GENOME_SIZE)")
    art_illumina --rndSeed 1 --in $GENOME --paired --len 75 --fcov $FOLD_COVERAGE --seqSys HS25 --mflen 500 --sdev 20 --noALN --out $FASTQ_DIR/dna_ && gzip $FASTQ_DIR/dna_1.fq && gzip $FASTQ_DIR/dna_2.fq

```

Whether I glanced at the recipe, or did some kind of text processing, I could easily see the relationships and purpose of the software in the container.

#### Can I associate environment and metadata with modules?
Given two different functions for my container to perform, foo and bar, can I define environment variables and labels (metadata) that I know will be sourced (environment) or exposed (inspect labels) in the context of the app?

**Standard**
Singularity standard also has one global shared `%environment`, and `%labels` section. If two functions in the container share the same environment variable and the value is different, this must be resolved manually. For this example, the first container didn't have any labels or environment, however in practice these global sections are usually used for high level, global variables like author and version (of the container). When I run the container, regardless of if different contexts or variables are needed for executables inside, I get the same environment.

**SCI-F**
With SCI-F, I simply write the different variables to their sections, and have confidence that they will be sourced (environment) or inspected (labels) with clear association to the app.

```

%appenv run-rtg
    MEM=${MEM:-4g}
    THREADS=${THREADS:2}
    export MEM THREADS

%applabel run-rtg
    rtg-version 3.6.2
```

#### Do I need to know standard locations in advance?
Given that a container has conformance to SCI-F, do I need to know how it works to use it?

**Standard**
With a standard location, we would be relying on Linux File System conventions (e.g., installation of content to `/usr/local` or intuitively infer that a folder called `/Software` (as with this scientific example) or `/code` is likely where the creator put important content.

**SCI-F**
Instead of requiring the user to know that an app's base might be at `/scif/apps/foo`, we instead expose environment variables (e.g., `SINGULARITY_APPBASE`) that can be referenced at runtime to refer to different apps. This is especially important if, for example, I need to reference the output of one app as input for another, or simply know it's install location. Regardless of which app is running, the container will also expose the top level folder for all apps installations, and data, `SINGULARITY_DATA` and `SINGULARITY_APPS` at `/scif/data` and `/scif/apps`, respectively.


### 2. Production Evaluation
For this use case, we operate under the scenario that we are familiar with Singularity and the commands to use SCi-F, but we know nothing about the containers. We are a user.

#### Do I know what the container does?
The most natural thing to do with a Singularity container, knowing that it is possible to execute, is to do exactly that. For this evaluation, we want to assess how well executing the container reveals the intentions of the creator. 

**Standard**
From the runscript we evaluated earlier, we are presented with a list of software and versions installed in the image, without detail to where or for what purpose. While this listing is comprehensive, it's most appropriate for a developer than a scientific workflow. In this listing, it isn't clear how the software is used or combine in the analysis. We are reliant on some external script or controller that drives the container. We don't have any hints about possible analysis steps the container can serve.

```
The following software is installed in this image:
alabaster           0.7.9     Anaconda
anaconda            4.3.0     Anaconda
anaconda-client     1.6.0     Anaconda
anaconda-navigator  1.4.3     Anaconda
argcomplete         1.0.0     Anaconda
art                 20160605  Homebrew
astroid             1.4.9     Anaconda
astropy             1.3       Anaconda
babel               2.3.4     Anaconda
backports           1.0       Anaconda
beautifulsoup4      4.5.3     Anaconda
bitarray            0.8.1     Anaconda
blaze               0.10.1    Anaconda
bokeh               0.12.4    Anaconda
boto                2.45.0    Anaconda
bottleneck          1.2.0     Anaconda
bsdmainutils        9.0.10    Homebrew
bwa                 0.7.15    Homebrew
cairo               1.14.8    Anaconda
cffi                1.9.1     Anaconda
chardet             2.3.0     Anaconda
chest               0.2.3     Anaconda
click               6.7       Anaconda
cloudpickle         0.2.2     Anaconda
clyent              1.2.2     Anaconda
colorama            0.3.7     Anaconda
conda               4.3.8     Anaconda
conda-build         1.21.3    Anaconda
conda-env           2.6.0     Anaconda
configobj           5.0.6     Anaconda
contextlib2         0.5.4     Anaconda
cryptography        1.7.1     Anaconda
curl                7.52.1    Anaconda
cycler              0.10.0    Anaconda
cython              0.25.2    Anaconda
cytoolz             0.8.2     Anaconda
dask                0.13.0    Anaconda
datashape           0.5.4     Anaconda
dbus                1.10.10   Anaconda
decorator           4.0.11    Anaconda
dill                0.2.5     Anaconda
docutils            0.13.1    Anaconda
dynd-python         0.7.2     Anaconda
entrypoints         0.2.2     Anaconda
et_xmlfile          1.0.1     Anaconda
expat               2.1.0     Anaconda
fastcache           1.0.2     Anaconda
flask               0.12      Anaconda
flask-cors          3.0.2     Anaconda
fontconfig          2.12.1    Anaconda
freetype            2.5.5     Anaconda
get_terminal_size   1.0.0     Anaconda
gevent              1.2.1     Anaconda
glib                2.50.2    Anaconda
greenlet            0.4.11    Anaconda
gsl                 2.3       Homebrew
gst-plugins-base    1.8.0     Anaconda
gstreamer           1.8.0     Anaconda
h5py                2.6.0     Anaconda
harfbuzz            0.9.39    Anaconda
hdf5                1.8.17    Anaconda
heapdict            1.0.0     Anaconda
htslib              1.3.1     Homebrew
icu                 54.1      Anaconda
idna                2.2       Anaconda
imagesize           0.7.1     Anaconda
ipykernel           4.5.2     Anaconda
ipython             5.1.0     Anaconda
ipython_genutils    0.1.0     Anaconda
ipywidgets          5.2.2     Anaconda
isort               4.2.5     Anaconda
itsdangerous        0.24      Anaconda
jbig                2.1       Anaconda
jdcal               1.3       Anaconda
jedi                0.9.0     Anaconda
jinja2              2.9.4     Anaconda
jpeg                9b        Anaconda
jsonschema          2.5.1     Anaconda
jupyter             1.0.0     Anaconda
jupyter_client      4.4.0     Anaconda
jupyter_console     5.0.0     Anaconda
jupyter_core        4.2.1     Anaconda
kallisto            0.43.0    Anaconda
lazy-object-proxy   1.2.2     Anaconda
libbsd              0.8.3     Homebrew
libdynd             0.7.2     Anaconda
libffi              3.2.1     Anaconda
libgcc              4.8.5     Anaconda
libgfortran         3.0.0     Anaconda
libiconv            1.14      Anaconda
libpng              1.6.27    Anaconda
libsodium           1.0.10    Anaconda
libtiff             4.0.6     Anaconda
libxcb              1.12      Anaconda
libxml2             2.9.4     Anaconda
libxslt             1.1.29    Anaconda
_license            1.1       Anaconda
llvmlite            0.15.0    Anaconda
locket              0.2.0     Anaconda
lxml                3.7.2     Anaconda
markupsafe          0.23      Anaconda
matplotlib          2.0.0     Anaconda
mistune             0.7.3     Anaconda
mkl                 2017.0.1  Anaconda
mkl-service         1.1.2     Anaconda
mpmath              0.19      Anaconda
multipledispatch    0.4.9     Anaconda
nb_anacondacloud    1.1.0     Anaconda
nb_conda            1.1.0     Anaconda
nb_conda_kernels    1.0.3     Anaconda
nbconvert           4.2.0     Anaconda
_nb_ext_conf        0.2.0     Anaconda
nbformat            4.2.0     Anaconda
nbpresent           3.0.2     Anaconda
ncurses             6.0_2     Homebrew
networkx            1.11      Anaconda
nltk                3.2.2     Anaconda
nose                1.3.7     Anaconda
notebook            4.3.1     Anaconda
numba               0.30.1    Anaconda
numexpr             2.6.1     Anaconda
numpy               1.11.3    Anaconda
numpydoc            0.6.0     Anaconda
odo                 0.5.0     Anaconda
openpyxl            2.4.1     Anaconda
openssl             1.0.2k    Anaconda
pandas              0.19.2    Anaconda
parallel            20170122  Homebrew
partd               0.3.7     Anaconda
patchelf            0.9_1     Homebrew
patchelf            0.9       Anaconda
pathlib2            2.2.0     Anaconda
path.py             10.0      Anaconda
patsy               0.4.1     Anaconda
pcre                8.39      Anaconda
pep8                1.7.0     Anaconda
pexpect             4.2.1     Anaconda
pickleshare         0.7.4     Anaconda
pillow              4.0.0     Anaconda
pip                 9.0.1     Anaconda
pixman              0.34.0    Anaconda
pkg-config          0.29.1_2  Homebrew
ply                 3.9       Anaconda
prompt_toolkit      1.0.9     Anaconda
psutil              5.0.1     Anaconda
ptyprocess          0.5.1     Anaconda
py                  1.4.32    Anaconda
pyasn1              0.1.9     Anaconda
pycosat             0.6.1     Anaconda
pycparser           2.17      Anaconda
pycrypto            2.6.1     Anaconda
pycurl              7.43.0    Anaconda
pyflakes            1.5.0     Anaconda
pygments            2.1.3     Anaconda
pylint              1.6.4     Anaconda
pyopenssl           16.2.0    Anaconda
pyparsing           2.1.4     Anaconda
pyqt                5.6.0     Anaconda
pytables            3.3.0     Anaconda
pytest              3.0.5     Anaconda
python              3.5.2     Anaconda
python-dateutil     2.6.0     Anaconda
pytz                2016.10   Anaconda
pyyaml              3.12      Anaconda
pyzmq               16.0.2    Anaconda
qt                  5.6.2     Anaconda
qtawesome           0.4.3     Anaconda
qtconsole           4.2.1     Anaconda
qtpy                1.2.1     Anaconda
readline            6.2       Anaconda
redis               3.2.0     Anaconda
redis-py            2.10.5    Anaconda
requests            2.12.4    Anaconda
rope                0.9.4     Anaconda
RTG                 3.6.2     User_Install
ruamel_yaml         0.11.14   Anaconda
samtools            1.3.1     Homebrew
scikit-image        0.12.3    Anaconda
scikit-learn        0.18.1    Anaconda
scipy               0.18.1    Anaconda
seaborn             0.7.1     Anaconda
setuptools          27.2.0    Anaconda
simplegeneric       0.8.1     Anaconda
singledispatch      3.4.0.3   Anaconda
sip                 4.18      Anaconda
six                 1.10.0    Anaconda
snowballstemmer     1.2.1     Anaconda
sockjs-tornado      1.0.3     Anaconda
sphinx              1.5.1     Anaconda
sphinx_rtd_theme    0.1.9     Anaconda
spyder              3.1.2     Anaconda
sqlalchemy          1.1.5     Anaconda
sqlite              3.13.0    Anaconda
statsmodels         0.6.1     Anaconda
sympy               1.0       Anaconda
terminado           0.6       Anaconda
tk                  8.5.18    Anaconda
toolz               0.8.2     Anaconda
tornado             4.4.2     Anaconda
traitlets           4.3.1     Anaconda
unicodecsv          0.14.1    Anaconda
util-linux          2.29      Homebrew
wcwidth             0.1.7     Anaconda
werkzeug            0.11.15   Anaconda
wheel               0.29.0    Anaconda
widgetsnbextension  1.2.6     Anaconda
wrapt               1.10.8    Anaconda
xlrd                1.0.0     Anaconda
xlsxwriter          0.9.6     Anaconda
xlwt                1.2.0     Anaconda
xz                  5.2.2     Anaconda
xz                  5.2.3     Homebrew
yaml                0.1.6     Anaconda
zeromq              4.1.5     Anaconda
zlib                1.2.8     Anaconda
\Note that some Anaconda in the list are modules and note executables.
Example usage: analysis.img [command] [args] [options]
```

**SCI-F**
When we run the SCi-F container, a different kind of information is presented. By listing the container contents at `/scif/apps`, the user knows what pipeline steps are included in the analysis.

```
singularity run scif.img

```

and in fact, this basic listing is generally useful for apps, so it's provided as it's own command, if the user doesn't write a runscript at all:

```
singularity apps scif.img
```

The runscript also hints that I can direct the `%help` command to better understand an app. While both SCiF and standard Singularity allows for specification of a global `%help` section, providing documentation on the level of the application is more focused and offers specificity to the user. Both also allow for global `%labels` that might point to documentation, code repositories, or other informative information.

#### Does moduarity come naturally?
For this metric, we want to know if using different functions of the container (modules) is intuitive. As we have seen above, the definition of a "module" could be anything from a series of script calls to perform a pipeline step (alignment using bwa and samtools), a single script call (such as a module just for bwa), or even the same function applied to a specific set of content (e.g., download "A" v.s. download "B"). 

**Standard**
Singularity proper represents a module on the level of the container. The entire container is intended for one entrypoint, and any deviation from that requires customization of that entrypoint, or special knowledge about other executables in the container to call with `exec`. The container is modular only in context of being a single step in a pipeline.

**SCI-F**
SCI-F, in that the software inside is defined and installed in a modular fashion, makes it easy to find three different modules:

 - download-fastq 
 - download-rtg
 - download-reference

and without looking further, infer that likely these three downloads can be run in parallel. While this doesn't represent any kind of statement or assurance of this, it allows for the container have a natural modularity. Consider steps that are named according to an order:

 - 1-download 
 - 2-preprocess
 - 3-analysis

The creator of the container, in choosing a careful naming, can deliver a huge amount of information about different modules.

#### Do I know what executables are important in the container?
Without much effort, I should have a high level understanding of the different functions that the container performs, as defined by the creator. For example, a container intended for development of variant calling will expose low level tools (e.g, bwa, seqtk) while a container that is intended will expose a pipeline (e.g., mapping).

**Standard**
For Singularity standard, if the container performs one function (and one only), then a single runscript / entrypoint is sufficient. Having multiple functions is completely reliant on the programming ability of the creator. If the creator is able to write a more complex runscript that explains different use cases, the container can satisfy this goal. If not, the container is a black box. 

**SCI-F**
For SCI-F, without any significant programming ability, the different apps are easily exposed to the user.


#### Can I easily get help for an executable?
I should be able to run a command to get help or a summary of what the container does (introspection).

**Standard**
For Singularity standard, you **can** ask a container for help, but it's a single global help.

**SCI-F**
For SCI-F apps, you can define help sections for all of the functions that your container serves, along with a global help.


### 3. Research Evaluation
An important attribute of having modular software apps is that it allows for separation of files and executables for research, and those that belong to the base system. From a machine learning standpoint, it provides labels for some subset of content in the container that might be used to better understand how different software relates to a pipeline. Minimally, it separates important content from the base, allowing, for example, a recursive tree generated at `/scif` to capture a large majority of additions to the container. Or simple parsing of the build recipe to see what software (possibly outside of this location) was intended for each app. Equally important, having container software installed at a global at `%post` also says important things about it - that it perhaps is important for more than one software module, or is more of a system library.


## Conclusion
SCI-F clearly has the advantage when it comes to a container creator embedding his or her work with implied metadata about software and container contents. SCI-F also makes it easier to package different run scripts with the container, and expose them easily to the user. However, this does not mean that the standard approach of using a container as a general toolbox and distributing it with a series of external callers is bad or wrong. The choice to use (or not use) SCI-F apps is largely dependent on the goals of the creator, and the intended users.
