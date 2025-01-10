### README
## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Step-by-Step Guide in Jupyter-Notebook](#step-by-step-guide-in-jupyter-notebook)
- [Citation](#citation)

## Introduction 
This repository provides a Python implementation for calculating concordance scores for clinical indicators based on 
patient activity data. The **concordance score** evaluates the proportion of time clinical guidelines were followed for 
specific indicators (e.g., lab tests, vital sign measurements) during a defined evaluation period.

The **ratio model** is used to calculate concordance as follows:

$$ 
\text{Concordance Score} = \frac{\text{Number of concordant days in evaluation period}} {\text{Total number of days in evaluation period}}
$$


## Installation
1. Install Python:
   - Download and install Python from the [official Python website](https://www.python.org/).
   - Ensure `pip` is installed (it's included with Python versions 3.4 and above).

2. Clone this repository:
   ```bash
   git clone https://github.com/stephanievdl/ratio_code_test
   cd ratio_code_test
   ```


## Step-by-Step Guide in Jupyter-Notebook
To help you get started, we have included a Jupyter Notebook (`step_by_step_notebook.ipynb`) that demonstrates how to 
calculate concordance scores for selected indicators using 
1) A sample dataset. 
2) Your own dataset.

### Steps to Use the Notebook
1. **Launch JupyterLab**:
   Open a terminal or command prompt and run:
   ```bash
   jupyter-lab
   ```

2. **Open the Notebook**: 
   In JupyterLab, navigate to the location of `step_by_step_notebook.ipynb` and open the file.

3. **Follow the instructions**:
   The notebook provides step-by-step guidance to calculate concordance scores for selected indicators.


## Citation
If you use this code in your research or projects, please kindly cite the [following paper](https://):
Stephanie CC van der Lubbe, Lay Hoon Goh, Evangelos Kontopantelis, Wilson Wai San Tam, Jose M Valderas,
Measuring Guideline Concordance via Electronic Health Records: A New Model for Estimating Concordance Scores, 
*to be added*. DOI: 

