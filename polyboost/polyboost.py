# Daniel J. Parente, MD PhD
# University of Kansas Medical Center
# 2019-10-17

import xgboost, argparse, sys, io
import numpy as np
import pkg_resources

# Get the fields from a line that is read, split by tabs and deleting newlines and spaces
def getFields(line):
    fields = line.replace("\r\n","").replace("\n","").replace(" ","").split("\t")
    return fields

# Parse a PolyPhen file to extract the features
def parsePolyPhen(inputfile):
    data = inputfile.readlines()

    # Header line
    headerFields = getFields(data[0][1:]) #Parse the header, ignoring the initial #
    headers = headerFields[0:4] + ['polyboost_probability','polyboost_prediction']

    headerToIndex = { x : i for i, x in enumerate(headerFields) }
    vectors = []
    identifiers = []
    for line in data[1:]:
        # Ignored commented lines
        if line.replace(" ","").startswith('#'):
            continue

        # Gets all the fields
        fields = getFields(line)
        identifiers.append(fields[0:4])

        # Compatibility with training files: Training files with "missing" features for the last few features
        # of a row omit the tabs; put these back in as missing data
        if len(fields) < len(headerFields):
            fields = fields + [''] * (len(headerFields) - len(fields))

        # Get the data points of interest
        v = [
            fields[headerToIndex['Score1']], #Float 0
            fields[headerToIndex['dScore']], #Float 1
            fields[headerToIndex['Nobs']], # Int 2
            fields[headerToIndex['dVol']], # Float 3
            fields[headerToIndex['PfamHit']], # Int 4 (Will need special handling)
            fields[headerToIndex['IdPmax']], # Float 5
            fields[headerToIndex['IdQmin']], # Float 6
            fields[headerToIndex['CpG']], # Int 7
            fields[headerToIndex['NormASA']], #Float 8
            fields[headerToIndex['B-fact']] ,#Float 9
            fields[headerToIndex['dProp']]  # Float 10
            ]

        # Convert the floats:
        for i in [0,1,3,5,6,8,9,10]:
            if v[i] == '?' or v[i] == "":
                v[i] = float('nan')
            else:
                v[i] = float(v[i])

        # Convert the ints
        for i in [2,7]:
            if v[i] == '?' or v[i] == "":   # PolyPhen-2 produces '?' but training files have whitespace-only fields (stripped out in proprocessing)
                v[i] = None
            else:
                v[i] = int(v[i])

        # Convert PfamHit into a boolean: When absent this results in "NO" being output by PolyPhen2, which we
        # will interpret as 0 (False). Any other string except ? we will interpret as 1 (True). I am not sure
        # if it will ever output "?", but we will treat this as missing data.
        # Convert the ints
        for i in [4]: # Yes, I know the loop over 1 element is inefficient coding, but it matches the formatting above; readability!
            if v[i] == '?' or v[i] == "":
                v[i] = None
            else:
                v[i] = 0 if v[i] == "NO" else 1
        vectors.append(v)

    # Convert into an xgboost DMatrix (via NumPy array)
    npa = np.array(vectors)
    matrix = xgboost.DMatrix(npa)

    # Return the matrix and the identifiers (Protein and Variant specificiation)
    return matrix, identifiers, headers

def main():
    # Initialize an argument parsers
    parser = argparse.ArgumentParser(description="PolyBoost applies an extreme gradient boosted classifier to the PolyPhen-2 feaature set using a classifier based on either the HumVar or HumDiv training datasets.")

    parser.add_argument('--threads', help="Number of threads for the classifier", type=int, choices=range(1,17))
    parser.add_argument('--threshold', help="Custom threshold value for binary prediction", type=float)
    parser.add_argument('inputfile', type=argparse.FileType('r'), help="PolyPhen-2 Batch Processing Output File")
    parser.add_argument('classifier', metavar="classifier", choices=["humvar", "humdiv"],
                        help="Select the classifier to use; either humvar or humdiv")
    parser.add_argument('--out', help="Output file", type=argparse.FileType('w'), default=sys.stdout)

    # Read the arguments
    args = parser.parse_args()

    # Extract an XGBoost DMatrix and identifiers (protein and variant descriptions) from the file
    matrix, identifiers, header = parsePolyPhen(args.inputfile)

    # Initialize a (multithreaded) booster
    threads = 4
    if args.threads is not None:
        threads = args.threads

    # Initialize a classifier
    bst = xgboost.Booster({'nthread': threads})

    # Load the appropriate model and threshold
    threshold = None
    if args.classifier == "humvar":
        modelpath = pkg_resources.resource_filename("polyboost.models", "humvar-final.model")
        bst.load_model(modelpath)
        threshold = 0.5085766
    elif args.classifier == "humdiv":
        modelpath = pkg_resources.resource_filename("polyboost.models", "humdiv-final.model")
        bst.load_model(modelpath)
        threshold = 0.4356405
    else:
        raise("Invalid classifier")

    # Set a custom threshold if any
    if args.threshold is not None:
        threshold = args.threshold

    # Perform predictions to obtain probabilities of damaging
    pred = 1-bst.predict(matrix)    # The probability of damaging is actually actually 1- the model's probability due to the way it was trained

    # Apply a threshold to produce a binary prediction
    binary_predictions = [ "damaging" if x > threshold else "benign" for x in pred ]

    # Write the header
    print("\t".join(header), file=args.out)

    # Write out the identifiers (protein and variant specificiation), numerical probability of damaging variant and binary prediction
    for x,y,z in zip(identifiers, pred, binary_predictions):
        print("\t".join(x + [str(y), z]), file=args.out)

if __name__ == '__main__':
    main()