# ScoreMarkers

Function to score the expression of positive and negative markers in single cell RNA seq data.

The project aim is to score and label the cells from an anndata file based on the expression of positive and negative markers. The scoring function work as follow:<br/>
- Data normalization and scaling (using scanpy built-in functions)
- Evaluate mean expression (not median because frequently is 0)
- For each label the score is given by the formula below:<br/>

$$
Score=\sum{(\alpha\times\frac{\sum{_{i=1}^{n_p}}{f(x_p)}}{n_p}-\beta\times\frac{\sum{ _{i=1} ^{n_n}}{f(x_n)}}{n_n})}
$$

With:
- $\alpha$ = weight for positive markers
- $\beta$ = weight for negative markers
- $n_p$ = number of positive markers
- $n_n$ = number of negative markers
- $f(x_p)$ = function to score positive markers
- $f(x_n)$ = function to score negative markers


Where $f(x_p)$ is:

$$
f(x_p)=
\begin{cases}
x=2 & \quad \text{when exp$(x_p)$ > mean(total expression)}\\ 
x=1 & \quad \text{when mean	&le; exp$(x_p)$ > 0}\\ 
x=0 & \quad \text{otherwise}
\end{cases}
$$

and $f(x_n)$ is:

$$
f(x_n)=
\begin{cases}
x=-1 & \quad \text{when exp$(x_n)$ > 0}\\ 
x=0 & \quad \text{otherwise}
\end{cases}
$$

With:
- exp $(x_p)$ = expression of positive marker
- exp $(x_n)$ = expression of negative marker

For each label is calculated the score in every cell, once all the scores are computed, all the absolute values of the scores are sorted and the confidence threshold is calculated as the n<sup>th</sup> percentile. By default the is considered the 10th percentile. Then the label (including the threshold label) with the highest score is assigned to the cell:

$$
Label = max\quad s\quad |\quad s \in Score
$$

# Installation

## Obtain code:

Clone the repository:
```
git clone https://github.com/VergaJU/ScoreMarkers.git
```
and enter into the directory

## Prepare environment

To avoid conflicts with the global environmnent is suggested to create and activate a new virtual environment:
```
python3 -m venv scoremarkers
source ./scoremarkers/bin/activate
```

And then install the required packages:

```
pip install requirements
```

## Usage

The script that run all the scoring functions and assign labels is `GetLabels.py`
```

usage: GetLabels.py [-h] [-i INPUT_FILE] [-c CSV] [-o OUTPUT] [-a ALPHA] [-b BETA] [-l LABEL] [-t THRESHOLD_LABEL] [-v THRESHOLD_VALUE]

Function to score and label cells considering 2 set of markers.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file name, default is input_processed.h5ad
  -a ALPHA, --alpha ALPHA
                        Positive markers bias, correct the weight of positive markers. Default=1
  -b BETA, --beta BETA  Negative markers bias, correct the weight of negative markers. Default=1
  -l LABEL, --label LABEL
                        Name of the new variable in .obs, default = 'new_label'
  -t THRESHOLD_LABEL, --threshold_label THRESHOLD_LABEL
                        Name of the label to use when cells scores doesn't pass the threshold, default = 'Other'
  -v THRESHOLD_VALUE, --threshold_value THRESHOLD_VALUE
                        Value of the percentile to set the threshold scores has pass to label the cells, default=10

required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        anndata file with cells to be processed
  -c CSV, --csv CSV     CSV file with labels and markers, NOTE: First column must contains the labels Negative and positive markers must be in separated lines
                        The header need this format: 'Label' for the labels, 'PoN' for the column that specify if the marker is positive or negative The other
                        columns containing markers can be unnamed
```

The required arguments are the anndata input, a single cell experiment dataset obtained with scanpy and a csv file containing the labels and markers with the following structure:
```
Label,PoN,,,,,
label1,pos,marker1,marker2,marker3,marker4,marker5
label1,neg,marker1,marker2,marker3,,
```

The optional arguments include:
- output: the output name, if not added it will be the input iname appended with `_processed.h5ad`
- alpha: Bias for the positive markers weight, increasing it increase also the weight assigned to them
- beta: Bias for the negative markers weight. Since the negative markers "must" not be present, the default parameters assigned an higher value to beta than to alpha
- label: the name of the .obs variable to assign the new labels
- threshold label: label to add to the cells whose score doesn't pass the threshold.
- threshold value: percentile of the scores distribution used to set the threshold, cells without scores above that value will be labelled as "threshold value" has to pass to assign the label to the cell. Before setting the threshold, the function distribute all the absolute values of the scores and then set the threshold as the selected percentile.

The outputs are:
- updated anndata file with labels
- csv file with cells and scores for each label

## Example

Get help page:
```
python GetLabels.py -h
```

Minimum example command:
```
python GetLabels.py -i /path/to/anndata/file.h5ad -c /path/to/marker/file.csv -o /path/for/output/file.h5ad
```