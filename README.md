# PolyBoost

## Description
**PolyBoost** is a post-analysis tool for the batch processing output of [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2/) (I am not affiliated with the PolyPhen-2 group) that replaces the [naive Bayes classifier](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) with an extreme gradient boosting [XGBoost]([https://github.com/dmlc/xgboost) classifier.
## Citation

**PolyBoost: An enhanced genomic variant classifier** is currently under peer review. If, for some reason, you need to cite this in the meantime, please contact me.

## Usage

You will need to install PolyBoost (see below) and obtain the **batch mode output** from PolyPhen2 predictions (http://genetics.bwh.harvard.edu/pph2/bgi.shtml). An example of the batch mode output is found in polyphen2-example.txt in this repository.

After installation, run PolyBoost as follows:

    python3 -m polyboost [PolyPhen2 Output File] [Classifier]

Where [PolyPhen2 Output File] is the path to the **batch mode** output from PolyPhen-2 and Classifier is either **humdiv** or **humvar**. If you don't know which one to use, use **humvar**.   

## Installation

### QuickStart

Install PolyBoost with into Python 3.7 using:

    pip install polyboost
    
 If you are using Windows, XGBoost probably cannot be installed from PyPI. If you get an error message, follow the instructions for installation of XGBoost below and try this command again. 

### Requirements
PolyBoost requires Python 3.7, xgboost, numpy and scipy. Numpy and SciPy should be installed automatically as dependencies of PolyBoost and XGBoost. ***XGBoost will **not** be automatically installed because installation from PyPi does not work reliably on Windows at time of release. Install XGBoost before installing PolyBoost.***

Use of a Python [virtualenv](https://docs.python.org/3/library/venv.html) is recommended, but not required.


#### Python 3.7

You will need to install Python 3.7 through standard methods. 
#### XGBoost

Installation of XGBoost is ostensibly as easy as:

    pip3 install xgboost

I found, however, that I could not download this from the PyPi repository on Windows. In this case, you can download the XGBoost python wheel from [this location](https://www.lfd.uci.edu/~gohlke/pythonlibs/#xgboost) and install it like:

    pip3 install xgboost-{version}-{pythonversion}-{architecture}.whl

Example: I installed XGBoost 0.90 for Python 3.7 (32-bit), so I used:

    pip3 install xgboost-0.90-cp37-cp37m-win32.whl

I used the 32-bit module even though I am using 64-bit Windows because Python 3.7 was installed in 32 bit mode on my computer.

If you have difficulty, detailed instructions for installing XGBoost for your platform can be found here: https://xgboost.readthedocs.io/en/latest/build.html

#### PolyBoost

Install PolyBoost with:

    pip install polyboost
    
#### Numpy and Scipy

You should **not** need to install numpy and scipy manually, but you can do so with:

    pip install numpy scipy

## Options

### Number of Threads (--threads)
You can specify the number of threads to run predictions. You must choose between 1 and 16 threads. If you make no selection, the default is to use 4 threads.

Example using 8 threads:

    python3 -m polyboost polyphen2output.txt humvar --threads 8

### Threshold (--threshold)
You can manually choose a threshold between binary classification of "benign" and "damaging". The default choices are 0.504057 for HumVar and 0.4250919 for HumDiv. These defaults were determined during classifier development by maximizing the Youden index (sensitivity + specificity - 1) of the receiver operating characteristic (ROC) curve.

Example using a threshold value of 0.25:

    python3 -m polyboost polyphen2output.txt humvar --threshold 0.25 

### Output (--out)
By default, PolyBoost outputs to the console (standard output). You can optionally output to a file using --out.

Example redirecting to output.txt

    python3 -m polyboost polyphen2output.txt humvar --out output.txt

### 

## Questions?

Please e-mail me with questions. I will do my best to respond.