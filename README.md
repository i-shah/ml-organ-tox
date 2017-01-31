# ml-organ-tox

This repository is in support of a research manuscript (currently in peer-review) that describes a machine learning approach to predicting chemical-induced target-organ toxicity. This repository provides the source code, data and supporting information for reproducing our results. The machine learning workflow is given below:

[Workflow](figs/algorithm.png)

The entire analysis is implemented using open source tools. Either [clone](https://help.github.com/articles/cloning-a-repository/)  or download this repository to get started ...  

# Installation

## Requirements

### Operating system 
This system has only been tested under linux: [Red Hat Enterprise Linux v6/v7](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and [Ubuntu Server 16.04 LTS](http://www.ubuntu.com). It should be possible to run it under macOS and possibly Windows (if you can install the following requirements). 

### MongoDB
Unless already available, install the latest version of [mongodb](http://www.mongodb.com). MongoDB is a document-oriented database that is quite handy for storing complex data structures [here's](https://docs.mongodb.com/getting-started/shell/introduction/) an introduction to MongoDB. Read the [tutorial on mongodb security](https://docs.mongodb.com/manual/tutorial) and confirm that authentication switched on (edit the /etc/mongod.conf add "authorization:enabled" under security, then restart mongodb).

## Create the mongo database 

Login to mongodb as the administrative user and create a database by the name "organtox_v1". [Create the user](https://docs.mongodb.com/manual/reference/method/db.createUser/) devel (password: devel) that will be used by the code to access the organtox_v1 database. 

## Download the mongodb files
The data required for this analysis is available as a mongodump file and it must be downloaded via [ftp](ftp://newftp.epa.gov/comptox/staff/ishah/mongodb_organtox_v1.tbz). Please note his is a large file and may take a long time to download. Untar this file using the following commmand (will need bunzip2 decompression):-

```
unix> tar jxvf mongodb_organtox_v1.tbz
```

This will create the directory: organtox_v1 with the following contents:
```
organtox_v1/tox_fp.bson
organtox_v1/bio_fp.bson
organtox_v1/tox_fp.metadata.json
organtox_v1/ml_run_v1.bson
organtox_v1/ml_run_v1.metadata.json
organtox_v1/bio_fp.metadata.json
organtox_v1/ml_summary_v1.metadata.json
organtox_v1/ml_summary_v1.bson
organtox_v1/ml_lr_v1.metadata.json
organtox_v1/chm_fp.metadata.json
organtox_v1/chm_fp.bson
organtox_v1/ml_lr_v1.bson
```

## Restore the organtox_v1 database from the downloaded files


The following unix command (see documentation on [mongorestore](https://docs.mongodb.com/manual/reference/program/mongorestore/)) loads the contents of the data/mongodump directory into the htsdb database (by user=devel with password = devel):-

```
unix> mongorestore -u devel -p devel -d organtox_v1 organtox_v1
```

Test out the mongodb installation by connecting to the db using the mongo command line client.

# Python

This code has been tested on [Python 2.7](http://python.org). The easiest way to install Python 2.7 is via [Anaconda](https://www.continuum.io/downloads). Install the python packages given in lib/requirements.txt as follows:

```
unix> pip install -r lib/requirements.txt
```

## Jupyter
[Jupyter notebook](http://jupyter.org/) is an interactive computing environment based and is freely available. It can be installed using Anaconda. After jupyter is installed read the [quickstart instructions](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) to create a configuration file. Set the following variables in your jupyter_notebook_config.py:

```python
c.NotebookApp.port = 7777
```

Enter the notebooks directort and start the notebook server from the command line using (I always run this in a [screen](https://www.gnu.org/software/screen/manual/screen.html)):

```
unix> jupyter notebook
```


# Testing the system

After you have completed the above steps open [this page](http://localhost:7777) in your browser to run different steps of the analysis. The machine learning analysis given in "organ-tox-ml.ipynb", the impact of different factors on F1 performance is given in "organ-tox-stats.ipynb", and the code for reproducing all figures / tables is given in "organ-tox-figs.ipynb".


# Help

Feel free to contact me if you have any questions.
