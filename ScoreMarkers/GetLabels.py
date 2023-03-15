import functions.functions as fun
import argparse
from functions.functions import DefineLabel
parser = argparse.ArgumentParser(description="Function to score and label cells considering 2 set of markers.")

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument("-i", "--input_file",
                           help="anndata file with cells to be processed")
requiredNamed.add_argument("-c", "--csv",
                           help="CSV file with labels and markers, NOTE:\n"
                                "First column must contains the labels\n"
                                "Negative and positive markers must be in separated lines\n"
                                "The header need this format: 'Label' for the labels, "
                                "'PoN' for the column that specify if the marker is positive or negative\n"
                                "The other columns  containing markers can be unnamed")
parser.add_argument("-o", "--output",
                    help="Output file name, default is input_processed.h5ad",
                    default="")
parser.add_argument("-a", "--alpha",
                    help="Positive markers bias, correct the weight of positive markers. Default=1",
                    default=1)
parser.add_argument("-b", "--beta",
                    help="Negative markers bias, correct the weight of negative markers. Default=1",
                    default=1)
parser.add_argument("-l", "--label",
                    help="Name of the new variable in .obs, default = 'new_label'",
                    default="new_label")
parser.add_argument("-t", "--threshold_label",
                    help="Name of the label to use when cells scores doesn't pass the threshold, default = 'Other'",
                    default="Others")
parser.add_argument("-v", "--threshold_value",
                    help="Value of the percentile to set the threshold scores has pass to label the cells, default=10",
                    default=10)
parser.add_argument("-p", "--plot",
                    help="Save or not scores plots (True/False), default=True",
                    default=True)
args = parser.parse_args()

define = fun.DefineLabel(adata=args.input_file,
                         markers=args.csv)

define.get_label(alpha=args.alpha,
                 beta=args.beta,
                 newlabel=args.label,
                 newfile=args.output,
                 thresholdlab=args.threshold_label,
                 thresholdvalue=args.threshold_value,
                 plot=args.plot)
