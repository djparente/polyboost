# PolyBoost

## Description
**PolyBoost** is a post-analysis tool for the batch processing output of [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2/) that replaces the [naive Bayes classifier](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) with an extreme gradient boosting [XGBoost]([https://github.com/dmlc/xgboost) classifier.

(Note: I am not affiliated with the PolyPhen-2 group).
## Citation

**PolyBoost: An enhanced genomic variant classifier using extreme gradient boosting** is currently under peer review. If, for some reason, you need to cite this in the meantime, please contact me.

## QuickStart

You will need to install PolyBoost and obtain the **batch mode output** from PolyPhen2 predictions (http://genetics.bwh.harvard.edu/pph2/bgi.shtml). An example of the batch mode output is found in polyphen2-example.txt in this repository.

Install PolyBoost with into Python 3.7 using:

    pip3 install polyboost
    
If you are using Windows, xgboost (a dependency) probably cannot be installed from PyPI. If you get an error message, follow the instructions for installation of XGBoost below and try this command again.

After installation, run PolyBoost as follows:

    python -m polyboost.polyboost [PolyPhen2 Output File] [Classifier]

Where [PolyPhen2 Output File] is the path to the **batch mode** output from PolyPhen-2 and Classifier is either **humdiv** or **humvar**. If you don't know which one to use, use **humvar**. Make sure the PolyPhen2 input file is in your working directory (i.e. the directory you are running that command from).

Example:

    python -m polyboost.polyboost polyphen2-example.txt humvar


On some systems with multiple python distributions, you may need to use python3 (or python3.7) instead of "python" to use the correct version of python.

## Installation

### Requirements
PolyBoost requires Python 3.7, xgboost, numpy and scipy. PolyBoost will attempt to install all of these dependencies automatically. ***However, XGBoost may **not** be automatically installed because installation from PyPi does not work reliably on Windows at time of release. If an error occurs, install XGBoost (see below) before installing PolyBoost.***

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

    pip3 install polyboost
    
#### Numpy and Scipy

You should **not** need to install numpy and scipy manually, but you can do so with:

    pip3 install numpy scipy

## Options

### Number of Threads (--threads)
You can specify the number of threads to run predictions. You must choose between 1 and 16 threads. If you make no selection, the default is to use 4 threads.

Example using 8 threads:

    python -m polyboost.polyboost polyphen2-example.txt humvar --threads 8

### Threshold (--threshold)
You can manually choose a threshold between binary classification of "benign" and "damaging". The default choices are 0.504057 for HumVar and 0.4250919 for HumDiv. These defaults were determined during classifier development by maximizing the Youden index (sensitivity + specificity - 1) of the receiver operating characteristic (ROC) curve.

Example using a threshold value of 0.25:

    python3 -m polyboost.polyboost polyphen2-example.txt humvar --threshold 0.25 

### Output (--out)
By default, PolyBoost outputs to the console (standard output). You can optionally output to a file using --out.

Example redirecting to output.txt

    python -m polyboost.polyboost polyphen2-example.txt humvar --out output.txt

## Output Example
    o_acc   o_pos   o_aa1   o_aa2   polyboost_probability   polyboost_prediction
    P26439  186     P       L       0.35185128              benign
    P26439  205     L       P       0.09412336              benign
    P26439  213     S       G       0.37042004              benign
    P26439  216     K       E       0.60328233              damaging
    P26439  222     P       H       0.06907171              benign
    P26439  222     P       Q       0.39627028              benign
    P26439  222     P       T       0.20633507              benign
    P26439  236     L       S       0.7706197               damaging
    P26439  245     A       P       0.17939752              benign
    P26439  253     Y       N       0.044733346             benign
    P26439  254     Y       D       0.2756629               benign
    P26439  259     T       M       0.027224064             benign

## Questions?

Please e-mail me with questions. I will do my best to respond.