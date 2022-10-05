# ScoreMarkers

Function to score the expression of positive and negative markers in single cell RNA seq data.

The project aim is to score and label the cells from an anndata file based on the expression of positive and negative markers. The scoring function work as follow:<br/>
For each label the score is given by the formula below:<br/>

$$
\sum{(\alpha\times\frac{\sum_{i=1}^{n_p}{f(x_p)}}{n_p}+\beta\times\frac{\sum_{i=1}^{n_n}{f(x_n)}}{n_n})}
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
x=1 & \quad \text{when exp$(x_p)$ > 0}\\ 
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

For each label is calculated the score in every cell, once all the scores are computed, the label with the highest score is assigned to the cell.

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
python GetLabels.py -h
usage: GetLabels.py [-h] [-i INPUT_FILE] [-c CSV] [-o OUTPUT] [-a ALPHA] [-b BETA] [-l LABEL]

Function to score and label cells considering 2 set of markers (positive and negative)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file name, default is input_processed.h5ad
  -a ALPHA, --alpha ALPHA
                        Positive markers bias, correct the weight of positive markers. Default=0.8
  -b BETA, --beta BETA  Negative markers bias, correct the weight of negative markers. Default=1
  -l LABEL, --label LABEL
                        Name of the new variable in .obs, default = 'new_label'

required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        anndata file with cells to be processed
  -c CSV, --csv CSV     CSV file with labels and markers, NOTE: First column must contains the labels Negative and positive markers must be in separated lines The
                        header need this format: 'Label' for the labels, 'PoN' for the column that specify if the marker is positive or negativeThe other columns
                        containing markers can be unnamed

```