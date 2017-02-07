# ml-organ-tox

Can we use a [data-driven approach](https://en.wikipedia.org/wiki/Data_science) to predict the target-organ [toxicity](https://en.wikipedia.org/wiki/Toxicity) of chemicals in repeat dose [animal testing](https://en.wikipedia.org/wiki/Animal_testing) studies ? This  project frame chemical toxicity prediction as a supervised [machine learning (ML)](https://en.wikipedia.org/wiki/Machine_learning)  problem. Toxicity is assumed to be a categorical outcome, which is defined by histopathology data from legacy repeat-dose animal testing experiments. Chemicals that cause/do not cause histopathological effects in an organ are treated as positive/negative examples of each class. 

We use mulitple machine learning algorithms to build classifiers for each target organ class using different types of descriptors for chemicals including: (i) chemical - these are [standardised representations derived from the molecular graph](https://en.wikipedia.org/wiki/Molecular_descriptor), (ii) bioactivity descriptors - these are derived from [high-throughput screening](https://en.wikipedia.org/wiki/High-throughput_screening) experimental data, and (iii) hybrid descriptors - formed by a combination of chemical and bioactivity descriptors. Using cross-validation testing to evaluate performance, the impact of descriptor type, number of descriptors, number of +/- examples, and machine learning algorithm on predicting target organ outcomes was systematically evaluated. The algorithm for the workflow is given below:-

![](figs/algorithm.png)

# Installation
The entire analysis is implemented using open source tools. Either [clone](https://help.github.com/articles/cloning-a-repository/)  or download this repository to get started ...


## Operating system 
This system has only been tested under linux: [Red Hat Enterprise Linux v6/v7](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and [Ubuntu Server 16.04 LTS](http://www.ubuntu.com). It should be possible to run it under macOS and possibly Windows (if you can install the following requirements). 

## MongoDB
Unless already available, install the latest version of [mongodb](http://www.mongodb.com). MongoDB is a document-oriented database that is quite handy for storing complex data structures [here's](https://docs.mongodb.com/getting-started/shell/introduction/) an introduction to MongoDB. Read the [tutorial on mongodb security](https://docs.mongodb.com/manual/tutorial) and confirm that authentication switched on (edit the /etc/mongod.conf add "authorization:enabled" under security, then restart mongodb).

### Create the mongo database 

Login to mongodb as the administrative user and create a database by the name "organtox_v1". [Create an account](https://docs.mongodb.com/manual/reference/method/db.createUser/) that will be used by the code to access the organtox_v1 database. Currently, the username/password are set to devel/devel but these can be changed (make sure the code database connection code is changed in the jupyter notebooks - see below).

```
python> DB = openMongo(host='pb.epa.gov',user='devel',passwd='devel',db='organtox_v1')
```

### Download the mongodb files
The data required for this analysis is available as a mongodump file and it must be downloaded via [ftp](https://tinyurl.com/z4nfto6). Please note: this is a large file and it will take a long time to download. Untar this file using the following commmand (will need bunzip2 decompression):-

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

### Restore the organtox_v1 database from the downloaded files


The following unix command (see documentation on [mongorestore](https://docs.mongodb.com/manual/reference/program/mongorestore/)) loads the contents of the data/mongodump directory into the htsdb database (by user=devel with password = devel):-

```
unix> mongorestore -u devel -p devel -d organtox_v1 organtox_v1
```

Test out the mongodb installation by connecting to the db using the mongo command line client.
```
unix> mongo -u devel -p devel localhost/organtox_v1
```

## Python

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

After you have completed the above steps open [this page (if jupyter is running on port 7777)](http://localhost:7777) in your browser to run different steps of the analysis. The machine learning analysis given in notebooks/organ-tox-ml.ipynb, the analysis of different factors in machine learning on F1 performance is given in notebooks/organ-tox-stats.ipynb, and the code for reproducing all figures / tables is given in notebooks/organ-tox-figs.ipynb.


