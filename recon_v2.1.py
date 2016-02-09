# -*- coding: utf-8 -*-

from sys import argv, exit
from string import split, join, strip, rsplit
from collections import defaultdict
from os.path import expanduser, isfile, isdir, dirname, abspath, exists, isabs
from os import listdir, system, getcwd
from commands import getoutput
from math import isnan, log, ceil, floor, exp, factorial
from scipy.optimize.lbfgsb import _minimize_lbfgsb as minimize
from datetime import datetime
from numpy import sqrt, array, interp, mean, arange
from argparse import ArgumentParser, REMAINDER, RawDescriptionHelpFormatter
from time import time
from ast import literal_eval
import warnings

# BEGIN PARSING COMMAND-LINE ARGUMENTS ====================

parser = ArgumentParser(description = """
        
Recon: Reconstruction of Estimated Communities from Observed Numbers

Joseph Kaplinsky, PhD
Ramy Arnaout, MD, DPhil



0. Contents
-----------

1. What is Recon?
3. Requirements
2. Terminology
4. Running Recon
5. Diversity measures
6. Error bars
7. Power calculations
8. Resampling a fit
9. Contact information
10. License



1. What is Recon?
-----------------

Recon is an algorithm for generating a description of an overall
population from a sample (see Section 2, Terminology).

More precisely, Recon uses the distribution of species counts in a
sample to estimate the distribution of species counts in the
population from which the sample was drawn.

Recon:

- assumes nothing about the shape of the frequency distribution of
  species in the population (e.g. exponential, power, etc.),

- avoids (over)fitting of sampling noise,

- scans many starting points in an attempt to find a global best fit,
  and

- outputs 95% confidence intervals for the number of species of a
  given size in the population.

The code provided in this directory:

- determines the species frequency distribution of the overall
  population using a modified maximum-likelihood approach that stops
  short of overfitting, given species frequencies observed in a sample
  (see Section 4, Running Recon);

- outputs diversity measurements for the overall population (see
  Section 5, Diversity measures);

- outputs error bars (see Section 6, Error bars), and

- performs power calculations of the sample size required to be able
  to detect a given difference between two overall populations (see
  Section 7, Power calculations).



2. Requirements
---------------

recon.py requires:

- Python 2.7 (https://www.python.org)
- the SciPy Python library (http://scipy.org)

Plotting additionally requires:

- d3.js (http://d3js.org)
- wkhtmltopdf (http://wkhtmltopdf.org)
- cpdf (http://community.coherentpdf.com)

It has been tested on Macintosh OS X (v10.8-11) and several other Unix
systems.



3. Terminology
--------------

A "species" is group made up of one or more "individuals." The number
of individuals of a given species is that species' "size."  A "sample"
is a set of individuals from one or more species that is drawn
randomly from an "overall" or "parent" population. Species represented
in the population that are not represented in the sample are called
"missing species." (For historical reasons, in the code and elsewhere
in this readme, "clone" is used interchangeably for "species.")

"Diversity" refers to any of a set of measures of the frequency
distribution in the population. These measures can be thought of as
effective numbers of species in the population. The "Hill numbers"
("qD," pronounced "D-q") are a family of diversity measures defined by
the parameter q, which determines the degree to which diversity
measures are weighted toward larger species. For example, 0D is
"species richness," a diversity measure that weights all species
equally (and so is just a count of the number of different
species). infD is the reciprocal of the Berger-Parker Index, the
effective number of species if all species were the size of the
largest species. Simple mathematical transformations of Shannon
entropy (q=1), Gini-Simpson Index (q=2), and other common diversity
measures correspond to different Hill numbers.



4. Running Recon (-R, --run_recon)
----------------------------------

Description:

Given a set of observations of species frequencies in a sample as
input, the -R option outputs a set of parameters that describe the
modified maximum-likelihood species frequency distribution in the
parent population, without overfitting.

Usage:

python recon.py -R [options...] -o OUTPUT_FILE INPUT_FILE

Example:

python recon.py -R -t 30 -c -o test_sample_1_fitfile.txt test_sample_1.txt

([0.88205137420509439, 0.11794862579490561], [0.306026605847127,
1.0723167141789034], 5224621, {1: 1833459, 2: 405423, 3: 86822, 4:
18467, 5: 3694, 6: 626, 7: 128, 8: 20, 9: 1}, None,
-1572523.5668443954, 2.2577288150787354)

Output:

The above command returns a tuple with the following elements, in
order:

a list of weights that, with means, describes the reconstructed parent
distribution. Each weight w_i is the fraction of all species in the
parent distribution that each contribute a mean m_i number of
individuals to the sample;

a list of means that, with weights, describes the reconstructed parent
distribution. Eacn mean m_i is the mean number of individuals a
species of this size contributes to the sample. Means are Poisson
parameters;

an integer of the number of missing species;

a dictionary of the the species-size distribution in the sample (the
observed distribution), where each key is a species size, and the
corresponding value is the number of species of that size. If the -c
option is given, the keys and values should correspond to the left and
right columns of input data (a useful check that your data was read in
successfully);

if given (the -A option), the true number of species in the parent
population;

a float of the log-likelihood of this fit; and

a float of the time in seconds for the fit.

It also writes to FILE_OUT a summary of the fit. We refer to files
of this type as "fitfiles;" fitfiles are used as input for other Recon
functions (below). The last block in the fitfile that is offset by
multiple equal signs ("=======") contains the final weights and means
(as a single list of weights-then-means; "fitted parameters") and
missing species ("estimated n0"). These are the same as in the output
described above.

Required praameters:

OUTPUT_FILE (preceded by -o): the filename to be used for output. Note
that if this file exists it will not be overwritten; instead, recon
will exit with an error message.

INPUT_FILE: a text file containing the number of individuals of each
species seen in the sample (i.e., the species sizes in the sample).

The default format is a tab-delimited file with lines of
the form:

species name <tab> species size

with a newline character delimiting lines. Species size is an
integer. Sample data in this format can be found in the file
test_sample_data_4.txt; in it, species 9_0 has a size of 9
individuals, individual, species 8_0 has 8 individuals, species
1_1833458 has 1 individual, etc.

Note that files in this format can be long and therefore take a few
seconds to read (test_sample_data_4.txt is over two million lines
long). The alternative species distribution format (-c below) is much
more compact and therefore faster to read. For example,
test_sample_data_1.txt contains the same information as
test_sample_data_4.txt, only in this more compact format, and runs in
~2 seconds vs. ~7 seconds for test_sample_data_4.txt.

Main command-line options:

-a, --aicc_multiple 

Sets the multiple of the observed number of datapoints that Recon
considers observations. E.g., if only singletons, doublets, and
triplets are observed but user believes not seeing quadruplets is
evidence of absence (as opposed to absence of evidence), -a 1.3 will
tell Recon to consider this as four observations for purposes of
calculating AICc.

-c, --clone_distribution_in_file
 
allows recon to read an alternative tab-delimited format with
lines of the form

species size <tab> number of species of this size

where both values are integers. test_sample_1.txt, test_sample_2.txt,
and test_sample_3.txt are all in this format. As mentioned above,
test_sample_4.txt contains the same data as test_sample_1.txt, only in
the default extended format.

-d, --bin_size BIN_SIZE (default: 1)

Average number of observations per individual. In many circumstances,
each individual in the sample will be observed and counted
once. However, there are cases where each individual in the sample
will be observed and counted multiple times. BIN_SIZE allows for this
possibility.

-l, --parameter_limit PARAMETER_LIMIT (default: 20)

The maximum number of parameters that the algorithm will use to fit
the data (default, 20). Recon will continue adding parameters until
the AICc indicates that additional parameters are not justified. In
practice the limit of 20 is essentially never reached.

-t --threshold THRESHOLD (default: 30)

The -t option allows you to modify THRESHOLD, the size above which
sampling error is considered small, which means that Recon will assume
the frequency of species of this size or greater in the population is
the same as the frequency in the sample. It defaults to 30, which
usually gives good results: this is because if, in a sample from a
well mixed population, species A is seen 30 times in a sample it is
very unlikely that there is another species B which is the same size
as A in the parent population but is missing or very poorly
reperesented in the sample. Smaller values will give faster run times
but less accurate results.



5. Diversity measures (-D, --make_table_of_D_numbers)
-----------------------------------------------------

Description:

Given a fitfile (either the output from Section 4 or any population
description in that format), outputs a table of diversity measures as
Hill numbers (see Section 3, Terminology). Note that measures for any
Hill number are obtainable, but the appropriate Hill number must have
been asked for when making the ERROR_BAR_PARAMETER_FILE (see section
6).

Usage:

python recon.py -D -b ERROR_BAR_PARAMETER_FILE [options...] -o
OUTPUT_FILE INPUT_FILE [INPUT_FILE_2 INPUT_FILE_3 ...]

Example:

python recon_v2.1.py -D -Q 0 1 inf -b error_bar_parameters.txt -o
test_D_number_table.txt test_sample_1_fitfile.txt
test_sample_2_fitfile.txt test_sample_3_fitfile.txt

Output:

A table of Hill numbers for the reconstructed distribution, one row
for each input file (see INPUT_FILE below). Columns prefixed "obs_"
show the Hill numbers from the observed data in the sample (sample
diversities). Columns prefixed "recon_" show the Hill numbers Recon
has estimated for the population (overall diversities). The difference
between recon_0D and obs_0D is n_0, the estimated number of missing
species. Columns with the suffixes "+" and "-" indicate upper and
lower error-bar limits, respectively.

Required parameters:

-D, make_table_of_D_numbers 

Tells Recon to output a table of diversity measures.

-b ERROR_BAR_PARAMETERS

A file that contains parameters for constructing error bars on
fits. The supplied file error_bar_params.txt can be
used. Alternatively, Recon can generate an error bar parameter file
from a set of gold standard fits (see section 6 below).

-o, --file_out OUTPUT_FILE

The desired name of the output file.

INPUT_FILE [INPUT_FILE_2 INPUT_FILE_3 ...]

The input file(s). A fit file output from the fit in section 4,
Running Recon. If additional input files are listed, each will
generate one row in the output table.

-Q HILL_NUMBER [HILL_NUMBER_2 HILL_NUMBER_3 ...]

Hill number parameter(s) for table. Note that the error-bar parameters
file must be run for whatever Hill numbers are desired for the power
table.



6. Error bars (-e, --make_error_bars)
-------------------------------------

Description:

Generates an error bar parameter file from a set of fits on data for
which the number of missing species is known (i.e., validation
datasets). Needed for D-number tables and power tables.

Usage:

python recon.py -e -o OUTPUT_FILE ERROR_BAR_FIT_DIRECTORY

Example:

python recon.py -e error_bar_parameters.txt Test_dir

Output:

An error bar parameter file.

Required parameters:

-e, --make_error_bars

Tells recon to make an error bar parameter file.

-o, --file_out OUTPUT_FILE

The name of the new error bar parameter file

ERROR_BAR_FIT_DIRECTORY

The name of a directory that contains the fits with known missing
species. The known missing species are encoded in the weights and
means of the population.



7. Power calculations (-p, --make_power_table)
----------------------------------------------

Description:

Generates a power table with the minimum sample size required to be
able to detect differences of a given magnitude in a given Hill number
for two populations. That is, if you have two populations, and want to
be able detect a difference in 1D of x%, the table tells you how big
your samples have to be.

Usage:

python recon.py -p [-q HILL_NUMBER -m MIN_NUMBER_OF_DOUBLETS]
-o FILE_OUT ERROR_BAR_PARAMETERS

Example:

python recon_v2.1.py -p -q 0 -o test_power_table.txt error_bar_parameters.txt

Output:

The rows and columns of the output table show the minimum fold
differences that can be detected for different population
sizes.

Required parameters:

-o, --file_out OUTPUT_FILE

The desired name of the output file.

ERROR_BAR_PARAMETERS

A file that contains parameters for constructing error
bars on fits. The supplied file error_bar_params.txt can
be used. Alternatively, Recon can generate an error bar
parameter file from a set of gold standard fits.

Optional parameters:

-C LIST_OF_NUMBER_OF_CLONES [LIST_OF_NUMBER_OF_CLONES ...]

The rough number of species in the overall population; same as the
columns of the power table.

-F LIST_OF_FOLD_DIFFERENCES [LIST_OF_FOLD_DIFFERENCES ...]

The rows of the table.

-q, --q HILL_NUMBER <default : 0>

The Hill number for which the power calculation is carried
out.

-m, --min_number_of_doublets MIN_NUMBER_OF_DOUBLETS

An additional statistical minimum required for good results. The
default of 100 should be good for most purposes.


8. Resampling a fit (-r, --resample)
------------------------------------

Description:

This allows resampling of a fit (model). The output gives the maximum
likelihood observed species size distribution of samples from the
model. This is the distribution that Recon attempts to make as close
as possible to the observed species size distribution. The closeness
of the fits can be compared to measure the goodness of fit.

Usage:

python recon.py -r -o OUTPUT_FILE INPUT_FILE

Example:

python recon_v2.1.py -r -o test_sample_1_resample.txt
test_sample_1_fitfile.txt

Output:

A list of species sizes up to the threshold that was used in the
original fit is outputted, together with a count of species for
each size. Output is written both to standard output and to
OUTPUT_FILE.

Required parameters:

OUTPUT_FILE

The desired name of the output file.

INPUT_FILE

A file that contains model fitted parameters as output from a
previously completely Recon fit.



9. Plotting resample against observed (-x, --make_resample_plot)
----------------------------------------------------------------

Example:

python recon_v2.1.py -x --x_max 30 -o test_sample_1_plotfile.txt -b error_bar_parameters.txt test_sample_1_fitfile.txt

Outputs a plotfile with the specified name and a .pdf file with the
same prefix as the input file. Requires d3.js (in resource_path, along
with style.css and plot_clone_size_distribution_ref.js), wkhtmltopdf
(somewhere in PATH; just a regular install usually takes care of
this), and cpdf (ditto) to run.



10. Contact information
-----------------------

Correspondence should be addressed to Ramy Arnaout at
rarnaout@gmail.com.


11. License
----------
PLEASE READ THIS AGREEMENT.  ANY USE OF THE SOFTWARE OR ANY OF THE SOURCE 
CODE, OR REPRODUCTION OR MODIFICATION OF THE SOFTWARE OR SOURCE CODE 
INDICATES YOUR ACCEPTANCE OF THE TERMS OF THIS AGREEMENT.  IF YOU DO NOT 
AGREE, DO NOT USE, COPY, OR MODIFY THE SOFTWARE.  YOU MAY PRINT THIS 
AGREEMENT FOR YOUR RECORDS.

THIS SOURCE CODE SOFTWARE LICENSE (the “Agreement”) is between Beth Israel 
Deaconess Medical Center, Inc. (“BIDMC”) and you (“Licensee”) as of the 
date that you accept these terms by clicking “Agree” below (“Effective 
Date”). You agree as follows:

1. Definitions
(a) “Derivative” means any translation, adaptation, alteration, 
transformation, or modification, including
inclusion as part of another software program or product, of the Software.
(b) “Embedded Terms” means any terms and conditions of this Agreement 
embedded in the Source Code.
(c) “Intellectual Property Rights” means all patents, patent rights, patent 
applications, copyrights, copyright registrations, trade secrets, 
trademarks and service marks (including, where applicable, all derivative 
works of the foregoing).
(d) “Object Code” means computer programs assembled, compiled, or converted 
to magnetic or electronic binary form, which are readable and useable by 
computer equipment.
(e) “Software” means the software program known as, “Recon: Reconstruction 
of Estimated Communities from Observed Numbers”, including its Object Code 
and Source Code.
(f) “Source Code” means computer programs written in higher-level 
programming languages and readable by humans.

2. License. Subject to the terms and conditions of this Agreement, BIDMC 
grants to Licensee a no cost, personal, non-exclusive, non-transferable, 
limited license (without the right to sublicense) to download the Software, 
and to copy, make Derivatives and use the Software for Licensee’s internal 
academic and research purposes during the term of this Agreement. Any 
rights not expressly granted in this Agreement are expressly reserved.
(a) Derivatives. Licensee agrees that from time to time, or upon request by 
BIDMC, Licensee will deliver all Derivatives to BIDMC and hereby grants to 
BIDMC a no cost, personal, perpetual, irrevocable, non-exclusive, non- 
transferable, limited license (with the right to sublicense) to download 
the Derivatives, to copy, distribute and make Derivatives of the 
Derivative, and to use the Derivative for BIDMC’s internal academic and 
research purposes.
(b) Commercial Restrictions on Use of the Software. The following are 
prohibited without obtaining a commercial license from the Office of 
Technology Ventures at BIDMC:
(i) Using the Software or any Derivative to produce any commercial product 
or to provide any commercial service.
(ii) Charging a fee for use of the Software or any Derivative for any 
purpose.
(iii) Distributing the Software or any Derivative to any other party, 
unless the distribution is made subject to this Agreement and the recipient 
“Agrees” to this Agreement through the website: 
https://github.com/ArnaoutLab/Recon
(c) Patents. BIDMC does not grant through this Agreement any licenses under 
any BIDMC patent or patent application.
(d) Intellectual Property Rights Notices; Embedded Terms. Licensee is 
prohibited from removing or altering any of the Intellectual Property 
Rights notice(s) and any Embedded Terms embedded in the Software. Licensee 
must reproduce the unaltered Intellectual Property Rights notice(s) and the 
Embedded Terms in any full or partial copies of the Source Code that 
Licensee makes.
(e) Export Compliance. Licensee acknowledges and agrees that U.S. export 
control laws and other applicable export and import laws govern download 
and use of the Software. Licensee will neither export nor re-export, 
directly or indirectly, the Software in violation of U.S. laws or use the 
Software for any purpose prohibited by U.S. laws.

3. Disclaimer of Warranties; Limited Liability.
(a) Disclaimer of Warranties. THE SOFTWARE IS DELIVERED “AS IS.” BIDMC 
MAKES NO OTHER WARRANTIES WHATSOEVER, EXPRESS OR IMPLIED, WITH REGARD TO 
THE SOFTWARE PROVIDED UNDER THIS AGREEMENT, IN WHOLE OR IN PART. BIDMC 
EXPLICITLY DISCLAIMS ALL WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY 
AND OF FITNESS FOR A PARTICULAR PURPOSE. BIDMC EXPRESSLY DOES NOT WARRANT 
THAT THE SOFTWARE, IN WHOLE OR IN PART, WILL BE ERROR FREE, OPERATE WITHOUT 
INTERRUPTION OR MEET LICENSEE’S REQUIREMENTS.
(b) Limited Liability; No Consequential Damages. THE TOTAL LIABILITY OF 
BIDMC, ITS AFFILIATES, TRUSTEES, OFFICERS, AND EMPLOYEES IN CONNECTION WITH 
THE SOFTWARE, OR ANY OTHER MATTER RELATING TO THIS AGREEMENT (WHATEVER THE 
BASIS FOR THE CAUSE OF ACTION) WILL NOT EXCEED IN THE AGGREGATE OVER THE 
TERM OF THE AGREEMENT 00. IN NO EVENT WILL BIDMC , ITS AFFILIATES, 
TRUSTEES, OFFICERS, AND EMPLOYEES BE LIABLE FOR ANY SPECIAL, INCIDENTAL, 
INDIRECT OR CONSEQUENTIAL DAMAGES OR DAMAGES FOR LOST PROFITS, WHETHER 
BASED ON BREACH OF CONTRACT, TORT (INCLUDING NEGLIGENCE), PRODUCT 
LIABILITY, OR OTHERWISE.
(c) Failure of Essential Purpose. THE LIMITATIONS SPECIFIED IN THIS SECTION 
WILL SURVIVE AND APPLY EVEN IF ANY REMEDY SPECIFIED IN THIS AGREEMENT IS 
FOUND TO HAVE FAILED OF ITS ESSENTIAL PURPOSE.

4. Termination.
(a) Right of Termination. BIDMC may, for any reason or no reason, upon 
written notice sent to the contact information Licensee provides upon 
clicking “agree”, immediately terminate this Agreement. This Agreement will 
automatically terminate upon any breach by Licensee of any of the terms or 
conditions of this Agreement. BIDMC shall not be liable to Licensee or any 
third party for any termination of this Agreement.
(b) Licensee Derivatives. Upon termination or expiration of this Agreement, 
Licensee shall promptly deliver to BIDMC the Source Code and Object Code of 
Licensee’s Derivatives and shall immediately cease all use of the Software.

5. Non-use of Name. Without BIDMC’s prior written consent, Licensee will 
not identify BIDMC in any promotional statement, or otherwise use the name 
of any BIDMC employee or any trademark, service mark, trade name, or symbol 
of BIDMC.

6. Assignment. This Agreement and the rights and obligations hereunder are 
personal to Licensee, and may not be assigned or otherwise transferred, in 
whole or in part, without BIDMC’s prior written consent. Any attempt to do 
otherwise shall be void and of no effect. BIDMC has the right to assign 
this Agreement or any rights or obligations hereunder to any third party. 
This Agreement shall be binding upon, and inure to the benefit of, the 
successors, representatives and permitted assigns of the parties.

7. Choice of Law. This Agreement and all disputes and controversies related 
to this Agreement, are governed by and construed under the laws of the 
Commonwealth of Massachusetts, without regard to the choice of law 
provisions. The state and federal courts located in the Commonwealth of 
Massachusetts are the exclusive forum for any action between the parties 
relating to this Agreement. Licensee submits to the jurisdiction of such 
courts, and waives any claim that such a court lacks jurisdiction over 
Licensee or constitutes an inconvenient or improper forum. The United 
Nations Convention on the International Sale of Goods (CISG) shall not 
apply to the interpretation or enforcement of this Agreement.

8. English Language. This Agreement is originally written in the English 
language and the English language version shall control over any 
translations.

9. Entire Agreement. This Agreement constitutes the entire agreement, and 
supersedes all prior and contemporaneous agreements or understandings (oral 
or written), between the parties about the subject matter of this 
Agreement. Licensee has no right to waive or modify any of this Agreement 
without the written consent of BIDMC. No waiver, consent or modification of 
this Agreement shall bind either party unless in writing and signed by the 
party granting the waiver. The failure of either party to enforce its 
rights under this Agreement at any time for any period will not be 
construed as a waiver of such rights. If any provision of this Agreement is 
determined to be illegal or unenforceable, that provision will be limited 
or eliminated to the minimum extent necessary so that this Agreement will 
otherwise remain in full force and effect and enforceable.

Last Updated February 8, 2016 

(c) 2015-2016 Beth Israel Deaconess Medical Hospital, Inc. All rights
reserved.
        

""", formatter_class=RawDescriptionHelpFormatter) # try to maintain alphabetical order a-zA-z
parser.add_argument('-a', '--aicc_multiple', type=float, default=1.3, help='multiple of the observed number of datapoints used for AICc')
parser.add_argument('-b', '--precomputed_error_bar_file', type=str, help='file name for precomputed error bars')
parser.add_argument('-c', '--clone_distribution_in_file', action='store_true', help='input file is of the form clone_size tab number_of_clones of that size')
parser.add_argument('-d', '--bin_size', type=int, default=1, help='average number of observations for each individual')
parser.add_argument('-e', '--make_error_bars', action='store_true', help='make a precomputed error bar file that can be used to calculate error bars')
parser.add_argument('-f', '--fraction_small_clones', default=0.1, type=float, help='see get_sampling_limit')
parser.add_argument("-j", "--javascript_file", default="plot_clone_size_distribution_ref.js", help="reference javascript file, for plotting resample")
parser.add_argument('-l', '--parameter_limit', default=20, type=int, help='maximum number of parameters to fit')
parser.add_argument('-m', '--min_number_of_doublets', default=100, type=int, help='see get_sampling_limit')
parser.add_argument('-n', '--number_of_error_bars', default=2., type=float, help='mapping between gold-standard error bar fits and standard deviations in power calculation')
parser.add_argument('-o', '--file_out', default="test_recon_output.txt", type=str, help='output file for fit')
parser.add_argument('-p', '--make_power_table', action='store_true', help='output a table of cells required to see fold difference')
parser.add_argument('-q', '--q', type=str, default=0., help='Hill number parameter to use when making table for power calculations')
parser.add_argument('-r', '--resample', action='store_true', help='output a resampling of the model fit with the given fit file')
parser.add_argument('-s', '--noise_test_ratio_threshold', default=3., type=float, help='for avoiding fitting noise')
parser.add_argument('-t', '--threshold', type=int, default=30, help='largest clone to fit using parameters. Any larger clones will be assumed to be statistically accurate between sample and overall population and added back in after the fit is performed')
parser.add_argument('-u', '--cutoff_score', type=float, default=0.1, help='for stopping iteration when estimated clone distribution is at a fixed point')
parser.add_argument('-v', '--verbose', action='store_true', help='print verbose output to stdout')
parser.add_argument('-x', '--make_resample_plot', action='store_true', help="plots resample of fit against original sample")
parser.add_argument('-A', '--total_clones_in_overall_population', type=float, default=None, help='total clones in an overall distribution. for use with plot_resample (optionally) and reconstruct_overall_distribution (optionally)')
parser.add_argument('-C', '--list_of_number_of_clones', nargs='+', default=[1e4, 3e4, 1e5, 1e6, 3e6], help='column headers for power table')
parser.add_argument('-D', '--make_table_of_D_numbers', action='store_true', help='output a table of reconstructed D numbers with error bars from input fit files.')
parser.add_argument('-F', '--list_of_fold_differences', nargs='+', default=[1.1, 1.2, 1.3, 1.4, 1.5, 2., 5.], help='row titles for power table')
parser.add_argument('-Q', '--qs', nargs='+', default=[0., 1., 2., float('inf')], help='space-separated list of Hill-number parameters for making error-bar files and table of D numbers')
parser.add_argument('-R', '--run_recon', action='store_true', help='run Recon on sample to generate overall distribution parameters and fit file')
parser.add_argument('-S', '--initial_scale_factors', nargs='+', default=[0.05, 0.225, 0.4, 0.575, 0.75 , 0.925, 1.1], help='for determining means')
parser.add_argument('-U', '--resource_path', type=str, default = expanduser('~')+"/Desktop", help='folder that contains style.css, d3.js, and plot_clone_size_distribution_ref.js')
parser.add_argument('-W', '--test_weights_list', nargs='+', default=[0.05, 0.12857143, 0.20714286, 0.28571429, 0.36428571, 0.44285714, 0.52142857, 0.6], help='list of weights to attempt to add')
parser.add_argument('-X', '--x_max', type=float, default=None, help='the x_max for the output plot')
parser.add_argument('-Y', '--y_max', type=float, default=None, help='the y_max for the output plot (note, log scale)')

parser.add_argument('file_in', nargs=REMAINDER, help='list of file names containing clone data. A long list of files can be conveniently input using a shell variable, e.g. MLE_files=$(ls);python recon.py file_out $MLE_files')

args = parser.parse_args()
globals().update(vars(args))

list_of_number_of_clones = map(float, list_of_number_of_clones)
list_of_fold_differences = map(float, list_of_fold_differences)
initial_scale_factors = map(float, initial_scale_factors)
test_weights_list = map(float, test_weights_list)
qs = map(float, qs)

warnings.filterwarnings('error') # This makes RunTime warnings (e.g. divide by zero) raise an error

# Handle file/directory names

if make_table_of_D_numbers:
    list_of_filenames = file_in

elif len(file_in) == 1:
    file_in = file_in[0]
    if make_error_bars:
        if not isdir(file_in):
            print 'error: input file %s is not a directory. Input must be a directory containing the error bar fits when making error-bar files. exiting...' % file_in
            exit()
        else:
            error_bar_fits_directory = file_in
            if not error_bar_fits_directory.endswith('/'):
                error_bar_fits_directory = error_bar_fits_directory + '/'
    else:
        if not isfile(file_in):
            print 'error: input file %s does not exist. exiting...' % file_in
            exit()

elif len(file_in) > 1:
    print 'error: multiple files input, but this is allowed only for outputting table of D numbers. check the order of command-line parameters. exiting...'
    exit()



# CLASSES ===============================================

class mixed_distribution:
    """
    A mixed_distribution object represents a mixed Poisson
    distribution.
    
    A mixed_distribution is defined by a sum of Poisson distributions
    where the weight and mean of each subdistribution must be specified
    and the weights must sum to one.
    """
    def __init__(self,weight_mean_pairs):
        total_weight = float(sum(k[0] for k in weight_mean_pairs)) # Normalise weights to 1.
        self.weights = [k[0]/total_weight for k in weight_mean_pairs]
        self.means = [float(k[1]) for k in weight_mean_pairs]
        self.weight_mean_pairs = zip(self.weights, self.means)

    def value_at_clone_size(self,i):
        vals = [weight*poisson_pmf(mean, i) for weight, mean in self.weight_mean_pairs]
        return sum(vals)


class sample_from_parent:
    def __init__(self, weights, means, sample_size, total_clones, population_size, total_observed_clones, estimated_n0, fitted_weights, fitted_means):
        try:
            total_weight = float(sum(weights)) # Normalise weights to 1.0
            self.weights = [weight/total_weight for weight in weights]
        except: self.weights = weights
        self.means = means
        self.sample_size = sample_size
        self.total_clones = total_clones
        self.population_size = population_size
        self.fitted_weights = fitted_weights
        self.fitted_means = fitted_means
        self.total_observed_clones = total_observed_clones
        self.estimated_n0 = estimated_n0



# 1. FUNCTIONS =====================================

max_clone_size_for_which_factorials_can_be_converted_to_float = 170 # empirical platform limit
factorials = [factorial(k) for k in range( max_clone_size_for_which_factorials_can_be_converted_to_float+1 )] # need to go beyond threshold for when sample_from_model uses it
epsilon = 0.001
def poisson_pmf(mu, k): 
    try: return exp(-mu) * mu**k / factorials[k]
    except IndexError: return epsilon # if the factorial is too big, assume the result will be zero. Return an epsilon so eventually the clone limit will be reached; otherwise risk an infinite loop
    except OverflowError: pass
    try: return exp( -mu + k*log(mu) - log(factorials[k]) ) # triggered if OverflowError (number too large). Get around that by taking exp(log())
    except: pass
    return epsilon # to prevent an infinite loop 


# 1.1 functions for reconstruction ==================

def calculate_AICc(AIC, number_of_parameters, effective_number_of_data_points):
    """
    AIC is the Aikaike Information Criterion and AICc the corrected
    criterion for small numbers of observations
    """
    return AIC + 2.*number_of_parameters*(number_of_parameters + 1.) / (effective_number_of_data_points - number_of_parameters - 1.)


def read_observed_clone_size_distribution_v2(file_in, clone_distribution_in_file, threshold, bin_size = None):
    """
    Return dict[clones size] = count
    """
    observed_clone_size_distribution = defaultdict(int)
    with open(file_in,'rU') as f:

        if clone_distribution_in_file: # Format from literature
            for line in f:
                if line[0] != "#":
                    element_0, element_1 = split(line, "\t")
                    observed_clone_size_distribution[int(element_0)] = int(element_1)
        elif bin_size:
            bin_size = float(bin_size)
            threshold = ceil(threshold*bin_size) # Convert threshold on reads
            for line in f: # Format from experimental reads
                if line[0] != "#":
                    normed_clone_size = int( strip(split(line, "\t")[1]) )
                    if normed_clone_size < threshold:
                        binned_clone_size = int(ceil(normed_clone_size / bin_size))
                        if binned_clone_size != 0: # exclude clones in zero binned clone size
                            observed_clone_size_distribution[binned_clone_size] += 1
            for clone_size in observed_clone_size_distribution:
                observed_clone_size_distribution[clone_size] = int(observed_clone_size_distribution[clone_size])
            if 0 in observed_clone_size_distribution:
                del observed_clone_size_distribution[0]

    return observed_clone_size_distribution


def calc_log_likelihood_n0_as_data(estimated_clone_size_distribution, distribution, threshold):
    try:
        log_likelihood = sum( estimated_clone_size_distribution[k]*log(distribution.value_at_clone_size(k)) for k in range(0, threshold) ) # note 0
    except ValueError: # triggered if new_distribution.value_at_clone_size(k) = 0.
        log_likelihood = -float('inf')
        if verbose: print "Log likelihood is -inf!"
    return log_likelihood


def calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, distribution, threshold, total_observed_clones):
    try: 
        log_likelihood = sum( estimated_clone_size_distribution[k]*log(distribution.value_at_clone_size(k)) for k in range(1, threshold) ) # note 1
        p_0 = distribution.value_at_clone_size(0)
        log_likelihood -= total_observed_clones * log(1. - p_0)
    except ValueError: # triggered if new_distribution.value_at_clone_size(k) = 0.
        log_likelihood = -float("inf")
        if verbose: print 'Log likelihood is -inf!'
    return log_likelihood


def log_likelihood_score(fit_params, estimated_clone_size_distribution):
    weights = list(fit_params[:len(fit_params)/2])
    if any(weight > 0. for weight in weights):
        weights = [max(weight, 0.) for weight in weights] # Ensure this is a local variable!
    else:
        weights = [float(weight == max(weights)) for weight in weights]
        if verbose: print 'All weights negative', weights
    weight_norm = sum(weights)
    try:  weights = [weight / weight_norm for weight in weights]
    except ZeroDivisionError: weights = [1.] + weights[1:]
    weights = [1. - sum(weights[1:])] + weights[1:]
    means = list(fit_params[len(fit_params)/2:]) # Ensure this is a local variable!
    means = [max(mean,1e-6) for mean in means]
    test_distribution = mixed_distribution(zip(weights,means))
    log_likelihood = calc_log_likelihood_n0_as_data(estimated_clone_size_distribution, test_distribution, threshold)
    return -log_likelihood


def MLE_fit_v2(file_out, threshold, new_weights, new_means, observed_clone_size_distribution):
    """
    v2 does not return the initial parameters.
        
    The code which call MLE_fit_v2 already knows these.
    """
    if verbose: 
        print '==== BEGIN FIT ===\ninitial parameters '+str(new_weights+new_means)

    with open(file_out,'a') as f: 
        f.write('initial parameters = %s\n' % (new_weights + new_means))
    total_observed_clones = float(sum([observed_clone_size_distribution[key] for key in observed_clone_size_distribution if key < threshold]))
    new_distribution = mixed_distribution(zip(new_weights,new_means))
    total_estimated_clones = round( total_observed_clones/(1.-new_distribution.value_at_clone_size(0)) )
    estimated_n_0 = total_estimated_clones - total_observed_clones
    estimated_clone_size_distribution = defaultdict(int)
    estimated_clone_size_distribution[0] = estimated_n_0

    for key in observed_clone_size_distribution:
        estimated_clone_size_distribution[key] = observed_clone_size_distribution[key]
    count_iterations = 0
    score_function = cutoff_score + 1. # score_function stops the iteration when the estimated clone distribution is at a fixed point.
    dictionary_of_fits = {}
    iteration_limit = 20

    while ( (score_function > cutoff_score) and (count_iterations < iteration_limit) ): # Could implement this as separate cutoffs for weights and means. Begin with one cutoff.
        count_iterations += 1

        # Step 1: Use the estimated_clone_distribution to calculate an
        # updated mixed distribution. This is the MLE of weights and
        # means given the estimated clone distribution. scorefunction2
        # stops the iteration when the weights and means are at a
        # fixed point.

        if verbose:
            print "iteration %i" % count_iterations
            print "old estimated n0\t%i" % estimated_n_0 # This output gives the estimated clones at the start of the iteration
        old_estimated_n_0 = estimated_n_0
        start_params = array(new_weights+new_means)
        try:
            res = minimize(log_likelihood_score, # "L-BFGS-B"; no bounds or callback
                           start_params, 
                           args=(estimated_clone_size_distribution,))
            new_weights = list(res.x[:len(res.x)/2])
            if any(weight > 0 for weight in new_weights): # Useful to record in file?
                new_weights = [max(weight,0.) for weight in new_weights]
            else:
                new_weights = [float(weight == max(new_weights)) for weight in new_weights]
            weight_norm = sum(new_weights)
            new_weights = [weight / weight_norm for weight in new_weights]
            new_means = list(res.x[len(res.x)/2:])
            if any(new_mean < 0. for new_mean in new_means):
                new_means = [max(new_mean, 1e-6) for new_mean in new_means]
        except ValueError:
            if verbose: print "ValueError in minimize!"
            new_weights = start_params[:len(start_params)/2]
            new_means = start_params[len(start_params)/2:]
        if verbose: print "Fitted parameters (new_weights, new_means) %s %s" % (new_weights, new_means)

        # Step 2: Calculate an updated estimated_clone_distribution
        # based on old_distribution and
        # observed_clone_size_distribution

        new_distribution = mixed_distribution(zip(new_weights,new_means))
        total_estimated_clones = round(total_observed_clones / (1. - new_distribution.value_at_clone_size(0)))
        estimated_n_0 = round(total_estimated_clones - total_observed_clones)
        if verbose: print "new estimated n_0\t%i" % estimated_n_0
        estimated_clone_size_distribution[0] = estimated_n_0 # NB we could try looking at the whole of the estimated distribution, not just estimated_n0, as done here.
        try:
            score_function = 2.*abs(estimated_n_0 - old_estimated_n_0)/float(estimated_n_0 + old_estimated_n_0)
        except ZeroDivisionError:
            score_function = 0.
        if verbose: print "score_function %f" % score_function
        log_likelihood = calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, new_distribution, threshold, total_observed_clones)
        if verbose: print "Log likelihood\t%f\n" % log_likelihood # Note that this is the log likelihood POST update of estimated_n_0

    dictionary_of_fits[log_likelihood] = new_weights + new_means # Take only self consistent values!
    if count_iterations == iteration_limit:
        with open(file_out,'a') as f:
            f.write("Iteration limit reached - check for only small numbers of clones appearing once or twice in sample\n")
        if verbose: print "Iteration limit reached - check for only small numbers of clones appearing once or twice in sample"
    new_weights = dictionary_of_fits[max(dictionary_of_fits.keys())][:len(dictionary_of_fits.values()[0])/2]
    new_means = dictionary_of_fits[max(dictionary_of_fits.keys())][len(dictionary_of_fits.values()[0])/2:]
    new_distribution = mixed_distribution(zip(new_weights,new_means))
    total_estimated_clones = round(total_observed_clones / (1. - new_distribution.value_at_clone_size(0)))
    estimated_n_0 = round(total_estimated_clones - total_observed_clones)

    with open(file_out,'a') as f:
        f.write("\n".join([
                    'n_obs %s' % str(total_observed_clones),
                    'weights = %s' % str(new_distribution.weights),
                    'means = %s' % str(new_distribution.means),
                    'log likelihood = %s' % str(log_likelihood),
                    'estimated_n0 = %s\n\n' % str(estimated_n_0)
                    ])
                )

    if verbose:
        print "Completed iteration of clone sizes"
        print "Fitted weights\t%s" % new_distribution.weights
        print "Fitted means\t%s" % new_distribution.means
        print "Final estimated_n0\t%i" % estimated_n_0
        print "Log likelihood\t%f" % log_likelihood
        print "Observed clones\t%f" % sum(observed_clone_size_distribution.values())
        print "Observed clones + estimated clones\t%f" % (sum(observed_clone_size_distribution.values()) + estimated_n_0)
        print "=== END FIT ===\n"

    return new_distribution.weights + new_distribution.means, log_likelihood, estimated_n_0


def reconstruct_overall_distribution(file_in, file_out, total_clones=None, weights = 0., means = 0. ):
    """
    Given an observed sample distribution, returns:

    weights        : a list of weights for the spikes that best describe overall distribution
    means          : a list of means for the spikes that best describe the overall distribution
    n0             : the integer number of missing species
    log-likelihood : the float log-likelihood of this overall distribution, given the observed sample
    time elapsed   : in seconds
    """

    # for calculating time elapsed
    time_0 = time()

    # Read in and format data
    observed_clone_size_distribution = read_observed_clone_size_distribution_v2(file_in, clone_distribution_in_file, threshold, bin_size) # NB clones larger than threshold never enter fit
    observed_clone_size_distribution = dict((key, value) for key, value in observed_clone_size_distribution.iteritems() if value != 0) # remove zeroes
    
    # set effective_number_of_data_points for AICc
    no_observed_clone_sizes = len(observed_clone_size_distribution.keys())
    effective_number_of_data_points = min(aicc_multiple*no_observed_clone_sizes, threshold)

    # initialize final parameters (will be returned)
    final_weights = []
    final_means = []
    final_n0 = None
    final_log_likelihood = 1.
    final_no_spikes = 0.

    # Initialise starting weights and means. For the first round of fitting this should replicate the scan for two spikes. Initial mean depends only on the observed data:
    initial_mean = sum(clone_size * no_clones for clone_size, no_clones in observed_clone_size_distribution.items()) / float(sum(observed_clone_size_distribution.values()))

    starting_weights = [1.]
    starting_means = [initial_mean] # the list of starting points for means
    number_of_parameters = 1

    total_observed_clones = float(sum(observed_clone_size_distribution.values())) # NB this is used as a global variable - be careful!
    new_distribution = mixed_distribution( zip(starting_weights, starting_means) )
    total_estimated_clones = round(total_observed_clones/(1.-new_distribution.value_at_clone_size(0)))
    # initialize estimated_n0
    estimated_n_0 = total_estimated_clones - total_observed_clones
    estimated_clone_size_distribution = defaultdict(int)
    estimated_clone_size_distribution[0] = estimated_n_0


    for key in observed_clone_size_distribution:
        estimated_clone_size_distribution[key] = observed_clone_size_distribution[key]
    log_likelihood = calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, new_distribution, threshold, total_observed_clones)


    AIC = 2. * (number_of_parameters - log_likelihood)

    if no_observed_clone_sizes > 2:
        AICc = calculate_AICc( AIC, number_of_parameters, no_observed_clone_sizes )
        too_few_clone_sizes = False
    else:
        AICc = float('inf')
        too_few_clone_sizes = True

    AICc_old = float('inf')
    AIC_old = float('inf')

    with open(file_out,'w') as f:
        f.write("\n".join([
                    '# %s' % datetime.now().strftime('%c'),
                    '# python %s' % join(argv),
                    'observed_clone_size_distribution = %s' % dict(observed_clone_size_distribution),
                    'total_clones = %s' % total_clones,
                    'weights = 0.',
                    'means = 0.',
                    'initial_scale_factors = %s' % initial_scale_factors,
                    'test_weights_list = %s\n\n' % test_weights_list
                    ])
                )

    minimum_observed_clone_size = min(key for key in observed_clone_size_distribution if observed_clone_size_distribution[key] != 0)

    zero_in_fitted_weights = False

    # n0 is estimated above with the starting mean from data. Now we reinitialise starting weights
    starting_means = []
    starting_weights = []

    while (((AICc < AICc_old) or too_few_clone_sizes) and (number_of_parameters < parameter_limit) and (not zero_in_fitted_weights) and (1e-6 not in starting_means)):

        list_of_fits = [] # The list of fits will have one element for each scan starting point

        for initial_scale_factor in initial_scale_factors:

            # The next loop will add a new weight 
            if starting_weights == []:
                weights_to_loop_over = [1.]
            else:
                weights_to_loop_over = test_weights_list

            for test_weight in weights_to_loop_over:

                starting_means_temp = [ initial_scale_factor * initial_mean ] + starting_means
                starting_weights_temp = [test_weight] + [weight*(1. - test_weight) for weight in starting_weights]

                end_params, log_likelihood, estimated_n_0 = MLE_fit_v2(file_out, threshold, starting_weights_temp, starting_means_temp, observed_clone_size_distribution)

                # An element of list_oF_fits contains the following elements:
                # fit[0] : a list of starting_weights+starting_means
                # fit[1] : a list of fitted end weights+means
                # fit[2] : the log likelihood of the fit
                # fit[3] : estimated n0 for the fit
                list_of_fits.append((starting_weights_temp + starting_means_temp, end_params, log_likelihood, estimated_n_0))

        sorted_fits = sorted(list_of_fits, key = lambda my_list: my_list[2], reverse = True) # sort by log_likelihood
        number_of_spikes = len(sorted_fits[0][0])/2 # This is the length of starting_weights+starting_means in the best fit, divided by 2
        sorted_means = [fit[0][number_of_spikes:] for fit in sorted_fits] # This takes the starting means for the best fit 
        min_means = [min(m) for m in sorted_means] # This takes the smallest amongst these starting means
        if min_means[0] <= initial_scale_factor * initial_mean: # This checks to see if the smallest starting mean was the mean that was added
            average_start_parameters = sorted_fits[0][0]
            average_start_parameters[number_of_spikes] = min_means[0] / 2.
            end_params, log_likelihood, estimated_n_0 = MLE_fit_v2(file_out, threshold, average_start_parameters[:len(average_start_parameters)/2], average_start_parameters[len(average_start_parameters)/2:], observed_clone_size_distribution)
            list_of_fits.append((average_start_parameters, end_params, log_likelihood, estimated_n_0))
            if verbose: print 'SMALLER SPIKE TRIGGERED'

        sorted_fits = sorted(list_of_fits, key = lambda my_list: my_list[2], reverse = True) # sort by log_likelihood

        average_start_parameters = [(p1+p2)/2. for p1,p2 in zip(sorted_fits[0][0], sorted_fits[1][0])] # average best and second best fits
        end_params, log_likelihood, estimated_n_0 = MLE_fit_v2(file_out, threshold, average_start_parameters[:len(average_start_parameters)/2], average_start_parameters[len(average_start_parameters)/2:], observed_clone_size_distribution)
        list_of_fits.append((average_start_parameters, end_params, log_likelihood, estimated_n_0))

        noise_test_ratio = 0.
        # noise_test_ratio compares the expected contribution to the sample of the smallest fitted clones
        # to the noise in the sample from remaining clones
        while noise_test_ratio <= 3.:
            max_log_likelihood = max(fit[2] for fit in list_of_fits)
            log_likelihoods = [element[2] for element in list_of_fits]
            max_index = log_likelihoods.index(max(log_likelihoods))

            fitted_parameters = list_of_fits[max_index][1]

            # Set up starting parameters for the next round of fitting:
            starting_weights = fitted_parameters[:len(fitted_parameters)/2]
            starting_means = fitted_parameters[len(fitted_parameters)/2:]
            sorted_params = sorted(zip(starting_weights, starting_means), key = lambda x: x[1])
            sorted_weights, sorted_means = [[x[i] for x in sorted_params] for i in [0,1]]
            if (0. in sorted_weights) and (len(sorted_weights) == 2): break
            if len(sorted_weights) > 1: # This applies on the second and subsequent interations, when there is more than one spike
                noise_test_ratio = sorted_weights[0]*poisson_pmf(sorted_means[0], minimum_observed_clone_size) / sqrt(sum( [weight * poisson_pmf(mean, minimum_observed_clone_size) for weight, mean in zip(sorted_weights[1:], sorted_means[1:]) ] ) / (total_observed_clones + estimated_n0) )
            else:
                noise_test_ratio = 4. # This applies on the first iteration, when there is only a single spike
            if noise_test_ratio <= noise_test_ratio_threshold:
                del list_of_fits[max_index]
                if len(list_of_fits) == 0:
                    max_log_likelihood = -float('inf')
                    break

        if 0. in fitted_parameters:
            zero_in_fitted_weights = True # flag to stop further iterations
            filtered_weights = [w for w in sorted_weights if w != 0.] # Remove zero from parameters
            filtered_means = [sorted_means[i] for i in range(len(sorted_weights)) if sorted_weights[i] != 0.]
            fitted_parameters = filtered_weights + filtered_means

        AIC_old = AIC
        AICc_old = AICc

        number_of_parameters = len(fitted_parameters) - 1
        AIC = 2.*(number_of_parameters - max_log_likelihood)

        if ((effective_number_of_data_points - number_of_parameters - 1.) > 0.):
            AICc = calculate_AICc(AIC, number_of_parameters, effective_number_of_data_points)
        else:
            AICc = float('inf')

        if 0. in fitted_parameters:
            AICc = float('inf')
            AIC = float('inf')

        if (AICc < AICc_old) or too_few_clone_sizes:
            with open(file_out, 'a') as f: # if you change this to 'w' the following is the only thing that gets output
                f.write('\n===========================================\n')
                if number_of_parameters >= parameter_limit:
                    f.write('parameter_limit = %s # REACHED MAX NUMBER OF PARAMETERS\n' % str(parameter_limit))
                initial_parameters, fitted_parameters, log_likelihood, estimated_n0 = list_of_fits[max_index]
                f.write("\n".join([
                            'initial parameters = %s' % initial_parameters,
                            'fitted parameters = %s' % fitted_parameters,
                            'estimated n0 = %s' % estimated_n0,
                            'AIC = %s' % AIC,
                            'AICc = %s' % AICc,
                            'log likelihood = %s' % log_likelihood,
                            '===========================================\n\n'
                            ])
                        )
            too_few_clone_sizes = False
            final_no_parameters = len(fitted_parameters)/2
            final_weights = fitted_parameters[:final_no_parameters]
            final_means = fitted_parameters[final_no_parameters:]
            final_n0 = estimated_n0
            final_log_likelihood = log_likelihood

    time_elapsed = time()-time_0
    return final_weights, final_means, int(final_n0), observed_clone_size_distribution, total_clones, final_log_likelihood, time_elapsed



# 1.2 functions for error-bar fits ==========================================


def read_sample_fit(MLE_filename, return_threshold = False):

    with open(MLE_filename, 'rU') as MLE_file:
        try:
            in_best_fit = False
            in_start_params = True
            for line in MLE_file:
                # read parameters from the start of the file...
                if line == '\n':
                    in_start_params = False
                    continue
                try: left_side, right_side = map(strip, split(line, "=", 1))
                except: continue
                if in_start_params:
                    if not line.startswith('#'):
                        line = split(line, '#')[0] ## Is this necessary?
                        if left_side == "parent_population_size": 
                            parent_population_size = literal_eval(right_side)
                        elif left_side == "total_clones":
                            total_clones = literal_eval(right_side)
                        elif left_side == "sample_size": 
                            sample_size = literal_eval(right_side)
                        elif left_side == "weights": 
                            weights = literal_eval(right_side)
                        elif left_side == "means": 
                            means = literal_eval(right_side)
                        elif left_side == "observed_clone_size_distribution": 
                            observed_clone_size_distribution = literal_eval(right_side)
                        elif left_side == "total_observed_clones":
                            total_observed_clones = literal_eval(right_side)
                        elif left_side == "initial_scale_factors":
                            initial_scale_factors = literal_eval(right_side)
                        elif left_side == "test_weights":
                            test_weights = literal_eval(right_side)
                elif line.startswith('='):
                    in_best_fit = not in_best_fit
                # ...and the best-fit parameters from later in the file
                elif in_best_fit:
                    if left_side == "initial parameters":
                        initial_parameters_try = literal_eval(right_side)
                    elif left_side == "fitted parameters":
                        fitted_params_try = literal_eval(right_side)
                    elif left_side == "estimated n0":
                        estimated_n0_try = literal_eval(right_side)
                        if not 1e-6 in fitted_params_try:
                            initial_parameters = initial_parameters_try
                            fitted_params = fitted_params_try
                            estimated_n0 = estimated_n0_try

            observed_threshold = max(observed_clone_size_distribution.keys())

            # set some defaults if they were not read in file
            if 'observed_clone_size_distribution' not in locals():
                if verbose: print 'warning: %s does not contain an observed_clone_size_distribution. Skipping...' % MLE_filename
                return "Bad fit", None, None
            if 'sample_size' not in locals():
                sample_size = sum([key*observed_clone_size_distribution[key] for key in observed_clone_size_distribution])
            if 'total_observed_clones' not in locals():
                total_observed_clones = sum(observed_clone_size_distribution.values())
            if 'parent_population_size' not in locals():
                parent_population_size = 1e9

            # get start weight and calculate start mean
            start_weight = initial_parameters[0]
            initial_mean_norm = float(sum([observed_clone_size_distribution[i] for i in observed_clone_size_distribution.keys() if i <= observed_threshold]))
            initial_mean = sum([i*observed_clone_size_distribution[i] for i in observed_clone_size_distribution if i <= observed_threshold]) / initial_mean_norm
            start_mean = initial_parameters[len(initial_parameters)/2] / initial_mean
            start_mean = int(10*start_mean)/10.

            # return sample_from_parent objects
            number_of_spikes = len(fitted_params)/2
            fit_weights = fitted_params[:number_of_spikes]
            fit_means = fitted_params[number_of_spikes:]
            if return_threshold:
                return observed_threshold, sample_from_parent(weights, means, sample_size, total_clones, parent_population_size, total_observed_clones, estimated_n0, fit_weights, fit_means), start_weight, start_mean, observed_clone_size_distribution
            else:
                return sample_from_parent(weights, means, sample_size, total_clones, parent_population_size, total_observed_clones, estimated_n0, fit_weights, fit_means), start_weight, start_mean, observed_clone_size_distribution

        except ValueError:
            if verbose: print '%s: warning: Nan in fit' % MLE_filename
            return 'Nan in fit', 0, 0, None


def sample_from_model(weights, means, sample_size, total_clones, return_unobserved_clones = False):
    """
    Returns the distribution of a given sample_size most likely to be
    observed, given weights, means, sample_size, and total_clones.

    Compare to output_resample, which I think does exactly the same
    thing, just with a different input.

    Note total_clones is the total number of clones in the *parent*;
    sample_size is the desired sample size of the output.
    """
    sample_size = int(sample_size) # in units of cells
    total_clones = int(total_clones) # Robust to user input as float

    # scale means to sample_size
    scaling_factor = sample_size / (total_clones * sum(w*m for w, m in zip(weights, means)))
    means = [m*scaling_factor for m in means]

    sample_distribution = mixed_distribution( zip(weights, means) )
    observed_clone_size_distribution = {}
    clone_size = 1
    total_observed_clones = 0
    unobserved_clones = total_clones * sample_distribution.value_at_clone_size(0)
    clone_limit = total_clones - unobserved_clones - 1 # When we observe this many clones we can stop looking

    while ( total_observed_clones < clone_limit ):
        new_clones = total_clones * sample_distribution.value_at_clone_size(clone_size)
        observed_clone_size_distribution[clone_size] = int(round(new_clones))
        total_observed_clones += new_clones
        clone_size += 1

    if return_unobserved_clones: return observed_clone_size_distribution, unobserved_clones
    else: return observed_clone_size_distribution


def qD_from_parameters(weights, means, total_estimated_clones, q):
    """
    Used by qD_from_parameters_and_observed_distribution.

    Should add in an estimate from the true distribution for values
    above the upper limit.

    This is slightly awkward, because I am generating the results
    directly from the parameters.

    I need to calculate the observed clone size distribution, work out
    out which clones are larger than the observed cut off and then add
    these to the reconstructed parent distribution.
    """
    zipped_weights_and_means = zip(weights, means)
    mean_norm = 1.0 / ( sum(weight*mean for weight,mean in zipped_weights_and_means) * total_estimated_clones )
    # This is one way to think of it; another is as 1. / the total
    # number of cells, in arbitrary units (e.g., units of some
    # fraction of a cell)

    if q == 0.0:
        return total_estimated_clones

    elif q == 1.0:
        qD = -sum(weight * (mean * mean_norm) * log(mean * mean_norm) for weight, mean in zipped_weights_and_means)
        qD = exp(qD * total_estimated_clones)

    elif q == float('inf'):
        qD = 1.0 / max( (mean * mean_norm) for mean in means )

    else:
        qD = sum(weight * (mean * mean_norm)**q for weight, mean in zipped_weights_and_means)
        qD = qD * total_estimated_clones
        qD = qD ** (1.0/(1.0 - q))

    return qD


def qD_from_parameters_and_observed_distribution(weights, means, observed_clone_size_distribution, observed_threshold, q):
    """
    Returns qD. Total_estimated_clones to be passed to this function
    should NOT be the sum of observed clones + estimated_0. It should
    be the sum of all clones arising from the fit parameters. This
    version incorporates various speed ups over version 3.
    """
    max_clone_size = max(observed_clone_size_distribution.keys())
    fitted_distribution = mixed_distribution( zip(weights,means) )

    obs_small_clones = sum( observed_clone_size_distribution[x] for x in range(1, observed_threshold) if x in observed_clone_size_distribution )
    distribution_norm = obs_small_clones / sum( fitted_distribution.value_at_clone_size(x) for x in range(1, observed_threshold) )

    total_estimated_clones_obs = distribution_norm # includes estimated_n0
    for x in range(1,observed_threshold):
        if x not in observed_clone_size_distribution:
            observed_clone_size_distribution[x] = 0

    fit_large_clones = {}
    unfitted_clones = {}
    # NB in the sum below, only take keys that are already in
    # observed_clone_size_distribution. Otherwise, as a defaultdict it
    # will create zero entries.
    for clone_size in range(observed_threshold, max(observed_clone_size_distribution.keys())+1):
        if clone_size in observed_clone_size_distribution.keys():
            if clone_size > max_clone_size_for_which_factorials_can_be_converted_to_float: fit_large_clones[clone_size] = 0. # factorials get very big, very fast. No point in trying to calculate beyond a certain point (that is presumably far lower than the max here)
            else: fit_large_clones[clone_size] = distribution_norm * fitted_distribution.value_at_clone_size(clone_size)
            diff = observed_clone_size_distribution[clone_size] - fit_large_clones[clone_size]
            if diff > 0.:
                unfitted_clones[clone_size] = diff

    aug_total_estimated_clones = float( total_estimated_clones_obs + sum(unfitted_clones.values()) )
    scale_weights = (total_estimated_clones_obs / aug_total_estimated_clones)
    weights = [weight*scale_weights for weight in weights]
    aug_weights = weights + [(unfitted_clones[clone_size] / aug_total_estimated_clones) for clone_size in unfitted_clones.keys()]
    aug_means = means + [clone_size for clone_size in unfitted_clones.keys()]

    qD = qD_from_parameters(aug_weights, aug_means, aug_total_estimated_clones, q)

    if verbose: 
        print 'Number of unfitted clones: %i' % len(unfitted_clones)
        print 'aug_weights (should sum to 1.): %s' % sum(aug_weights)

    return qD


def make_error_bar_envelope(samples, q):
    """
    This function considers all examples of fits at constant true
    total clone and returns:
        
    error_envelope_x_vals
    error_envelope_y_vals
        
    The keys of error_envelope_x_vals and error_envelope_x_vals are
    the true qD of the samples which we are using in our error bar
    fits.
        
    For error_envelope_x_vals at a fixed true qD the values are an
    ordered list of sample sizes.
        
    For error_envelope_y_vals at a fixed true qD the values are
    max(abs(log(fitted_D / true_D))) at the sample size given by the
    corresponding error_envelope_x_vals entry. The max is taken over
    steepnesses (that is, distribution shape).
        
    abs_deviations
        
    a list of (sample_size, abs(log(fitted_D / true_D)) pairs.
    """
    error_envelope_x_vals = {}
    error_envelope_y_vals ={}
    abs_deviations = []

    for sample in samples:
        observed_clone_size_distribution = sample_from_model(
            sample.weights, sample.means, sample.sample_size, 
            sample.total_clones
            )
        if observed_clone_size_distribution == {}: # if empty, continue (adds no information)
            continue
        fitted_D = qD_from_parameters_and_observed_distribution(
            sample.fitted_weights, 
            sample.fitted_means, 
            observed_clone_size_distribution, 
            30, 
            q
            )
        true_D = qD_from_parameters(sample.weights, 
                                    sample.means, 
                                    sample.total_clones, 
                                    q
                                    )
        if fitted_D != 0. and true_D != 0.:
            abs_deviations.append(
                (sample.sample_size, abs( log(fitted_D / true_D) ), true_D)
                )

    abs_deviations = [d for d in abs_deviations if not isnan(d[1])]

    true_Ds = sorted(list(set([abs_dev[2] for abs_dev in abs_deviations])))

    for true_D in true_Ds:
        sorted_dev = [dev for dev in abs_deviations if dev[2] == true_D]
        sorted_dev = sorted(sorted_dev, key = lambda x: x[0])
        x_vals = []
        y_vals = []
        max_y_val = 0.
        for val in sorted_dev[::-1]:
            if val[1] > max_y_val:
                max_y_val = val[1]
                if x_vals != [] and val[0] == x_vals[-1]:
                    y_vals[-1] = val[1]
                else:
                    x_vals.append(val[0])
                    y_vals.append(val[1])

        x_vals = x_vals[::-1]
        y_vals = y_vals[::-1]

        # For purposes of interpolation, set the error at 0 sample size equal to the error at smallest available sample size
        if x_vals[0] != 0:
            x_vals = [0] + x_vals
            y_vals = [y_vals[0]] + y_vals

        error_envelope_x_vals[true_D] = x_vals
        error_envelope_y_vals[true_D] = y_vals

    return error_envelope_x_vals, error_envelope_y_vals, abs_deviations


def make_error_bar_file(error_bar_fits_directory, precomputed_error_bar_file, qs):

    time_0 = time()

    # Filehandling
    if isdir(error_bar_fits_directory):
        if not error_bar_fits_directory.endswith('/'):
            error_bar_fits_directory += '/'
    else:
        print "error: %s is not a directory. To make error-bar files, file_in must be a directory that contains the error bar fits. exiting..." % error_bar_fits_directory
        exit()

    if isfile(precomputed_error_bar_file):
        print "error: precomputed error-bar file %s already exists. exiting..." % precomputed_error_bar_file
        exit()
    MLE_filenames = listdir(error_bar_fits_directory) # This directory should contain the fits
    MLE_filenames = [i for i in MLE_filenames if i.endswith(".txt")]

    # Read in fitted data for error bars
    fitted_samples = {} # dict[filename] = sample_from_parent objects
    for MLE_filename in MLE_filenames:
        if verbose: print 'processing %s' % MLE_filename
        try:
            fitted_sample, start_weight, start_mean, observed_clone_size_distribution = read_sample_fit(error_bar_fits_directory + MLE_filename)
        except UnboundLocalError:
            if verbose: print 'warning: invalid fit: %s skipping...' % MLE_filename
            continue
        if fitted_sample not in ('Nan in fit', 'Bad fit'):
            fitted_samples[MLE_filename] = fitted_sample

    # actually make error bars and write to file
    z = [fitted_samples[MLE_filename] for MLE_filename in fitted_samples.keys()]
    # The keys of error_envelope_x_vals and error_envelope_x_vals are
    # the true 0D of the samples which we are using in our error bar
    # fits. For error_envelope_x_vals at a fixed true 0D the values
    # are an ordered list of sample sizes. For error_envelope_y_vals
    # at a fixed true 0D the values are max(abs(log(fitted_D /
    # true_D))) at the sample size given by the corresponding
    # error_envelope_x_vals entry. The max is taken over steepnesses
    # (that is, distribution shape).
    #
    # Note that error_bar_devs has a default value of q = 1.
    for q in qs:
        error_envelope_x_vals, error_envelope_y_vals, abs_deviations = make_error_bar_envelope(z, q)
        with open(precomputed_error_bar_file, 'a') as f:
            f.write("\n".join([
                        'error_envelope_x_vals_%s = %s' % (q, error_envelope_x_vals),
                        'error_envelope_y_vals_%s = %s' % (q, error_envelope_y_vals),
                        'abs_deviations_%s = %s\n' % (q, abs_deviations)
                        ])
                    )
        if verbose: print "%sD:\tdone (%s)" % (q, precomputed_error_bar_file)

    time_elapsed = time()-time_0
    return time_elapsed


# 1.3 functions for power tables ======================================

def read_error_bar_file(error_bar_file, q=0):
    """
    error_bar_file should be a text file, as generated by make_error_bar_file
    containing values for error_envelope_x_vals and error_envelope_x_vals
    """
    with open(error_bar_file,'rU') as f:
        for line in f:
            if line.startswith("error_envelope_x_vals_%s" % q):
                error_envelope_x_vals = literal_eval( strip(split(line, "=")[1]) )
            elif line.startswith("error_envelope_y_vals_%s" % q):
                error_envelope_y_vals = literal_eval( strip(split(line, "=")[1]) )

    return error_envelope_x_vals, error_envelope_y_vals


def eb(sample_size, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target):
    """
    This returns the log of high error with the given parameters.
    The actual error will be exp(eb_val) - true_total_clone
    exp is because error_envelope_y_vals made by compare_MLE_alpha_v7.py contains log(fitted/true) values.
    """
    available_total_clone_sizes = sorted(error_envelope_x_vals.keys())
    for i in range(len(available_total_clone_sizes)):
        if available_total_clone_sizes[i] >= true_total_clone: break

    denom = float(available_total_clone_sizes[i] - available_total_clone_sizes[i-1])
    weight_1 = (available_total_clone_sizes[i] - true_total_clone) / denom
    weight_0 = (true_total_clone - available_total_clone_sizes[i-1]) / denom

    interpolated_1 = interp(sample_size, error_envelope_x_vals[available_total_clone_sizes[i-1]], error_envelope_y_vals[available_total_clone_sizes[i-1]])
    interpolated_0 = interp(sample_size, error_envelope_x_vals[available_total_clone_sizes[i]], error_envelope_y_vals[available_total_clone_sizes[i]])

    eb_val = weight_1 * interpolated_1 + weight_0 * interpolated_0

    return number_of_error_bars*(exp(eb_val) - 1.) - target


def bisect_eb(eb, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target):
    """
    This finds a root of number_of_error_bars*(exp(eb) - 1.) - target = 0.
    by iterated bisection.
    """
    low = 0
    high = true_total_clone*10
    eb_low = eb(low, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target)
    if eb_low < 0.:
        available_total_clone_sizes = sorted(error_envelope_x_vals.keys())
        for i in range(len(available_total_clone_sizes)):
            if available_total_clone_sizes[i] >= true_total_clone:
                break
        return sorted(error_envelope_x_vals[available_total_clone_sizes[max(0,i-1)]])[1]

    eb_high = eb(high, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target)
    assert eb_low * eb_high <= 0.
    while high - low > 1:
        midpoint = (low + high) / 2.
        eb_midpoint = eb(midpoint, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target)
        if eb_low * eb_midpoint > 0.:
            low = midpoint
        else:
            high = midpoint
    return midpoint


def get_sample_size(number_of_clones, fold_difference, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars):
    """
    This returns the number of cells in a sample that produce the an
    error bar of max_error_bar for a given number_of_clones in the
    parent population.
    
    This is the inverse function of the calculation performed by
    error_bar_on_fit_qD
    
    Delta = (fold_difference - 1.0) * number_of_clones
    """
    target = fold_difference - 1.
    sample_size = bisect_eb(eb, number_of_clones, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target)
    return sample_size


def get_sampling_limit(number_of_clones, min_number_of_doublets, fraction_small_clones):
    """
    This function returns the minimum number of cells required to
    reliably see doublets and triplets.
    
    The number of possible pairs of cells is sample_size**2/2
    
    Expected doublets = p * sample_size**2/2
    
    However, we need an effective sample size, given by the sample_size * fraction_of_small_clones

    where fraction_of_small_clones is the fraction of cells in small clones.    
        
    Expected doublets = p * (sample_size * fraction_of_small_clones)**2/2
        
    Where p is the probability that a pair of cells is are from the
    same clone, given that both cells are from small clones.
    
    Now estimate p = 1/number_of_small_clones, which is approximately
    1/number_of_clones.
    
    Expected doublets = (sample_size * fraction_of_small_clones)**2/ (2 * number_of_clones)

    We want to expect at least min_number_of_doublets.
        
    rearranging:
        
    sample_size = sqrt(2 * number_of_clones * min_number_of_doublets) / fraction_of_small_clones
    
    This is our sampling limit.
    """
    return int(round( sqrt(2.*min_number_of_doublets*number_of_clones) / fraction_small_clones ))


def get_minimum_cells(list_of_number_of_clones, list_of_fold_differences, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, min_number_of_doublets, fraction_small_clones):
    minimum_cells = {}
    for fold_difference in list_of_fold_differences:
        for number_of_clones in list_of_number_of_clones:
            # Calculate the minimum number of cells required to
            # reliably see a fold.
            fold_limit = get_sample_size(number_of_clones, fold_difference, error_envelope_x_vals, error_envelope_y_vals,  number_of_error_bars)
            # Calculate the minimum number of cells required to
            # reliably see doublets and triplets
            sampling_limit = get_sampling_limit(number_of_clones, min_number_of_doublets, fraction_small_clones)
            # Take the more stringent result
            minimum_cells[number_of_clones, fold_difference] = int(round(max(fold_limit, sampling_limit)))

    return minimum_cells


def output_table(list_of_number_of_clones, list_of_fold_differences, error_bar_file, number_of_error_bars, min_number_of_doublets, fraction_small_clones, minimum_cells, q, file_output = None):

    lines = []

    # Record variables used
    lines.append('# %s' % datetime.now().strftime('%c'))
    lines.append('# python %s' % ' '.join(argv))
    lines.append('error_bar_file = %s' % error_bar_file)
    lines.append('number_of_error_bars = %s' % number_of_error_bars)
    lines.append('min_number_of_doublets = %s' % min_number_of_doublets)
    lines.append('fraction_small_clones = %s' % fraction_small_clones)
    lines.append('q = %s\n' % q)

    # make header line
    line = "\t" + "\t".join([str(int(number_of_clones)) for number_of_clones in list_of_number_of_clones])
    lines.append(line)

    # iterate over lines
    for fold_difference in list_of_fold_differences:

        # iterate over cols within line
        line = str(fold_difference) + "\t" 
        line += "\t".join([str(minimum_cells[number_of_clones, fold_difference]) for number_of_clones in list_of_number_of_clones])
        lines.append(line)
        
    # output lines
    lines = "\n".join(lines)+"\n"
    if file_output:
        if isfile(file_output):
            print "warning: power table file %s already exists. Printing to stdout and exiting..." % file_output
            print lines
            exit()
        with open(file_output, "a") as f:
            f.write(lines)
    elif verbose: print lines
    return


# 1.4 functions for D-number tables ======================================

def qD_from_distribution(H, q):
    """
    Returns the qD given a distribution in the form of a dictionary H = { clone_size : clone number }
    """
    # H is a dictionary-like object
    N = float(sum( k*v for k, v in H.items() )) # normalization factor
    H = dict( (k/N, v) for k, v in H.items() ) # normalize key entries by N, making them probabilities
    if q == float("inf"): return 1./max(H.keys())
    elif q == 1.: return exp( -sum( v*k*log(k) for k, v in H.items() if v != 0. ) )
    else: return sum( v*k**q for k, v in H.items() )**(1./(1.-q))


def error_bar(sample_size, true_total_clone, error_envelope_x_vals_q, error_envelope_y_vals_q):
    """
    The function error_bar(sample_size, true_total_clone,
    error_envelope_x_vals, error_envelope_y_vals) returns the value of
    max(abs(fitted_D / true_D)) at the given sample size and given
    total_clones.  Since we only have values at a grid of true_D and a
    grid of sample size, this is calculated by interpolation. This
    will be an interpolation between
        
        error_bar(sample_size, true_total_clone1)
        
    and
        
        error_bar(sample_size, true_total_clone2)
        
    and then between sample sizes.
    """
    available_total_clone_sizes = sorted(error_envelope_x_vals_q.keys())
    if true_total_clone > available_total_clone_sizes[-1] or (true_total_clone < available_total_clone_sizes[0]):
        if verbose: print 'Clone size out of range'
        return None
    for i in range(len(available_total_clone_sizes)):
        if available_total_clone_sizes[i] >= true_total_clone:
            break
                                                              
    weight_0 = (available_total_clone_sizes[i] - true_total_clone) / float(available_total_clone_sizes[i] - available_total_clone_sizes[i-1])
    weight_1 = (true_total_clone - available_total_clone_sizes[i-1]) / float(available_total_clone_sizes[i] - available_total_clone_sizes[i-1])

    interp_0 = interp(sample_size, error_envelope_x_vals_q[available_total_clone_sizes[i-1]], error_envelope_y_vals_q[available_total_clone_sizes[i-1]])
    interp_1 = interp(sample_size, error_envelope_x_vals_q[available_total_clone_sizes[i]], error_envelope_y_vals_q[available_total_clone_sizes[i]])
                                                              
    return weight_0*interp_0 + weight_1*interp_1


def abs_error(total_clones, sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q,high_low_sign):
    """
    high_low_sign = +1 will return the upper limit on the number of
    reconstructed clones for 0D = total clones and given sample size.
    Dependence on q (as which D numbers, qD) is embodied in
    error_envelope_x_vals and error_envelope_y_vals. These must be
    calculated for each q, and the appropriate dictionary passed to
    the abs_error function.
    """
    error_bar_result = error_bar(sample_size, total_clones, error_envelope_x_vals_q, error_envelope_y_vals_q)
    if error_bar_result == None: # clone size out of range
        return None
    else: 
        return total_clones * exp(high_low_sign * error_bar_result ) # NB total_observed_clones is the sample size


def error_bar_on_fit_qD(sample_size, reconstructed_qD, error_envelope_x_vals_q, error_envelope_y_vals_q):

    step_size = (reconstructed_qD / 1000.0)

    # calculate lower limit
    n = 1; abs_error_result = float("inf") # initialization
    while abs_error_result >= reconstructed_qD:
        abs_error_result = abs_error((reconstructed_qD - n*step_size), sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q, +1)
        if abs_error_result == None:
            lower_limit = 0.
            break
        n += 1
        lower_limit = reconstructed_qD - (n+1)*step_size

    # calculate upper limit
    n = 1; abs_error_result = 0. # initialization
    while abs_error_result <= reconstructed_qD:
        abs_error_result = abs_error((reconstructed_qD + n*step_size), sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q, -1)
        if abs_error_result == None:
            upper_limit = float("inf")
            break
        n += 1
        upper_limit = reconstructed_qD + (n+1)*step_size

    return lower_limit, upper_limit


def output_table_of_D_numbers(precomputed_error_bar_file, D_number_output_file, infiles, qs = [0., 1., 2., float("inf")], write_to_disk = True ):
    """
    Returns a dictionary of dictionaries: file_q_qD_hash[file] = dict[q] : (obs_qD, recon_qD, recon_lower_limit, recon_qD_upper_limit)
    """
    # read in error-bar envelope variables
    error_envelope_xs = {}
    error_envelope_ys = {}
    with open(precomputed_error_bar_file,'rU') as f:
        for line in f:
            left_side, right_side = map(strip, split(line, "="))
            name, q = left_side.rsplit("_",1)
            q = float(q)
            if name == "error_envelope_x_vals":
                error_envelope_xs[q] = literal_eval(right_side)
            elif name == "error_envelope_y_vals":
                error_envelope_ys[q] = literal_eval(right_side)

    # initialize
    estimated_n0 = {}  # estimated_n0[infile] = estimated n0 for that file
    fitted_params = {} # fitted_params[infile] = final fitted params for that file
    observed_clone_size_distributions = {} # observed_clone_size_distributions[infile] = observed_clone_size_distribution for that file
    obs_qDs = defaultdict(dict) # obs_qDs[filename] = obs_qD{q: qD}
    observed_threshold = 30 # not recorded in files but needed
    
    # Read in observed_clone_size_distributions, fitted_parameters,
    # and n0s from MLE files (note similarity to read_sample_fit)
    for infile in infiles:
        obs_qD = {}
        with open(infile,'rU') as MLE_file:
            in_best_fit = False
            in_start_params = True
            for line in MLE_file:
                if line == "\n": 
                    in_start_params = False
                    continue
                try: left_side, right_side = map(strip, split(line, "=", 1))
                except: continue
                if in_start_params:
                    # read in observed_clone_size_distribution and threshold (if present) from the start of the file...
                    if left_side == "observed_clone_size_distribution": # This requires full clone distribution in fit file, including clones that are not explicitly fitted.
                        observed_clone_size_distribution = literal_eval(right_side)
                        observed_clone_size_distributions[infile] = observed_clone_size_distribution
                    elif left_side == "threshold":
                        threshold_in_file = literal_eval(right_side)
                        if not observed_threshold:
                            observed_threshold = threshold_in_file
                        elif observed_threshold != threshold_in_file:
                            print 'error: all fit files in one table must use the same fitting threshold. Exiting...'
                            exit()
                elif line.startswith('='):
                    in_best_fit = not in_best_fit
                # ...and the best-fit parameters from later in the file
                elif in_best_fit:
                    if left_side == "fitted parameters":
                        fitted_params_try = literal_eval(right_side)
                    elif left_side == "estimated n0":
                        estimated_n0_try = literal_eval(right_side)
                        if not 1e-6 in fitted_params_try:
                            fitted_params[infile] = fitted_params_try
                            estimated_n0[infile] = estimated_n0_try

        # calculate observed qDs
        for q in qs: 
            obs_qD[q] = qD_from_distribution(observed_clone_size_distribution, q)
        obs_qDs[infile] = obs_qD
        #"""


    # write D-number table
    file_q_qD_hash = defaultdict(lambda : defaultdict(tuple) ) # file_q_qD_hash[file] = dict[q] : (obs_qD, recon_qD, recon_qD-, recon_qD+); for return

    # header line
    header_line = "sample_name\t" + "\t".join(["obs_%sD\trecon_%sD\trecon_%sD-\trecon_%sD+" % (q, q, q, q) for q in qs]) + "\n"
    if write_to_disk:
        with open(D_number_output_file, "w") as f:
            f.write("\n".join([
                        "# python %s" % join(argv),
                        "# infiles = %s" % infiles,
                        "# observed_threshold = %s" % observed_threshold,
                        "# precomputed_error_bar_file = %s" % precomputed_error_bar_file,
                        header_line
                        ]))

    for infile in infiles:

        print infile
        # get sample_size (needed for error bars)
        weights = fitted_params[infile][:len(fitted_params[infile])/2]
        means = fitted_params[infile][len(fitted_params[infile])/2:]
        observed_clone_size_distribution = observed_clone_size_distributions[infile]
        sample_size = sum(key*observed_clone_size_distribution[key] for key in observed_clone_size_distribution if key < observed_threshold)
        
        # get qD_obs, qD, qD-, qD+ for each q
        outline = [infile]
        for q in qs:
            q = float(q)
            obs_qD = obs_qDs[infile][q]
            recon_qD = qD_from_parameters_and_observed_distribution(weights, means, observed_clone_size_distribution, observed_threshold, q)
            lower_limit, upper_limit = error_bar_on_fit_qD(sample_size, recon_qD, error_envelope_xs[q], error_envelope_ys[q])
            upper_limit = 1.0164*upper_limit # (very minor) adjustment based on validation, to achieve 95% confidence internvals
            outline += [obs_qD, recon_qD, lower_limit, upper_limit]
            file_q_qD_hash[infile][q] = (obs_qD, recon_qD, lower_limit, upper_limit)
                
        outline = map(str, outline)
        outline = "\t".join(outline)
        
        # write output to file (one write per infile)
        if write_to_disk:
            with open(D_number_output_file, 'a') as f:
                f.write(outline+"\n")

        if verbose: 
            print "infile:\t%s" % infile
            print "weights:\t%s" % weights
            print "means:\t%s" % means
            print "observed_clone_size_distribution:\t%s" % observed_clone_size_distribution
            print "sample_size:\t%s" % sample_size
            print
        
    return file_q_qD_hash


# 1.4 functions for resmapling ======================

def output_resample(fitfile, outfile, write_to_file = True):
    """
    Given a file that contains a fit (with fitted parameters, etc.),
    which recall describes a reconstructed overall population, this
    function outputs a sample of sample size sample_size that can be
    compared to the initial sample from which the reconstructed
    overall population was made. Also writes this to outfile.

    Note the resample is the same size as the original. This, and that
    it operates directly on a file as opposed to weights and means,
    distinguishes this function from sample_from_model().
    """
    largest_observed_clone, model, start_weight, start_mean, observed_clone_size_distribution = read_sample_fit(fitfile, return_threshold = True)  # see 1.2 above

    observed_clone_size_distribution = defaultdict(int, observed_clone_size_distribution)
    sample_distribution = mixed_distribution(zip(model.fitted_weights, model.fitted_means))
    total_clones = model.total_clones
    total_observed_clones = 0
    unobserved_clones = model.estimated_n0

    if total_clones == None: total_clones = unobserved_clones + model.total_observed_clones

    # When we see this many clones we can stop looking. Equivalent to
    # total_clones - actual_n0. This makes more sense than
    # total_clones - estimated_n0, because that fixes the n0,
    # restricting room for the observed; whereas the uncertainty is in
    # the n0, so we should not restrict room for the observed in this
    # way.
    clone_limit = sum(observed_clone_size_distribution.values()) 

    resampled_clone_size_distribution = {}
    clone_size = 1
    while total_observed_clones < total_clones and clone_size <= largest_observed_clone:
        if clone_size > max_clone_size_for_which_factorials_can_be_converted_to_float: no_clones_at_this_size = 0
        else: no_clones_at_this_size = sample_distribution.value_at_clone_size(clone_size)
        new_clones = total_clones * no_clones_at_this_size
        if clone_size <= threshold: # add back large clones (since they were not fit anyway, and their frequency is well sampled)
            resampled_clone_size_distribution[clone_size] = int(round(new_clones))
        else: resampled_clone_size_distribution[clone_size] = observed_clone_size_distribution[clone_size]
        total_observed_clones += new_clones
        clone_size += 1


    header = "#clone_size\tno_clones"
    outlines = [header]
    resampled_clone_sizes = resampled_clone_size_distribution.keys()
    outlines += ["%i\t%i" % (i, resampled_clone_size_distribution[i]) for i in range(1, largest_observed_clone+1) if i in resampled_clone_sizes and resampled_clone_size_distribution[i] != 0]
    if write_to_file:
        with open(outfile, "w") as f: # write to disk
            f.write("\n".join(outlines) + "\n")
    return resampled_clone_size_distribution, observed_clone_size_distribution


#  1.5 functions for plotting resample ======================================


def write_html_file(html_file):
    """
    Writes the bare-bones html that calls the javascript that creates the plots in d3.js.
    """
    with open(html_file, "w") as f:
        f.write("""
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>plot_clone_size_distribution</title>
    <link href='http://fonts.googleapis.com/css?family=Roboto' rel='stylesheet' type='text/css'>
    <link href="%s/style.css" rel="stylesheet" type="text/css" />
    <script type="text/javascript" src="%s/d3.js"></script>
  </head>
  <body>
    <script type="text/javascript" src="%s/plot_clone_size_distribution.js">
    </script>
  </body>
</html>
""" % (resource_path, resource_path, resource_path) )
        return


def customize_js(javascript_file, plotfile, x_max="false", y_max="false"):
    txt = open(javascript_file, "r").read()
    txt = txt.replace("X_MAX", str(x_max))
    txt = txt.replace("Y_MAX", y_max)
    txt = txt.replace("SUPPRESS_OBS_ZERO", "false")
    txt = txt.replace('var filename = "test.txt";', 'var filename = "%s";' % plotfile)
    javascript_outfile = javascript_file.replace("_ref", "")
    f = open(javascript_outfile, "w")
    f.write(txt)
    f.close()
    return


def round_up_to_the_nearest(x, unit=10):
    return int(unit * ceil(float(x)/unit))


def plot_resample(fitfile, plotfile, precomputed_error_bar_file, x_max = x_max, y_max = y_max, javascript_file = resource_path + "/" + javascript_file, q=0.):
    """
    Note need javascript file plot_clone_size_distribution_ref.js as well as style.css and d3.js in current working directory.
    """
    # make sure plotfile includes full path (so the javascript file knows where to look for it)
    full_path = getcwd()
    if not isabs(rsplit(javascript_file, "/", 1)[0]): javascript_file = full_path + "/" + javascript_file
    if not isabs(rsplit(plotfile, "/", 1)[0]): plotfile = full_path + "/" + plotfile
    # get datapoints for observed and resampled
    resampled_clone_size_distribution, observed_clone_size_distribution = output_resample(fitfile, None, False)
    # get estimated_n0
    fitted_sample, x, x, x = read_sample_fit(fitfile)
    estimated_n0 = fitted_sample.estimated_n0
    # get error bars on estimated_n0
    file_q_qD_hash = output_table_of_D_numbers(precomputed_error_bar_file, None, [fitfile], [q], False)
    obs_qD, recon_qD, recon_qD_lower_limit, recon_qD_upper_limit = file_q_qD_hash[fitfile][q]
    # apply error bars to estimated_n0
    lower_n0_error_bar = estimated_n0 - (recon_qD - recon_qD_lower_limit)
    lower_n0_error_bar = max(0.7, lower_n0_error_bar) # actually zero; did this to make it appear on the plot
    upper_n0_error_bar = estimated_n0 + (recon_qD_upper_limit - recon_qD)
    # get actual n0
    total_clones = fitted_sample.total_clones
    if total_clones:
        actual_n0 = total_clones - sum(observed_clone_size_distribution.values())
    else: actual_n0 = 0. # suppresses display (since y-axis is in log units)
    actual_n0 = max(actual_n0, 0.) # avoid negatives, which get plotted high
    # get clone sizes
    observed_clone_size_distribution = defaultdict(int, observed_clone_size_distribution) # so looking for a clone size that isn't there returns zero instead of an error
    resampled_clone_size_distribution = defaultdict(int, resampled_clone_size_distribution) # ditto
    clone_sizes = resampled_clone_size_distribution.keys() + observed_clone_size_distribution.keys()
    clone_sizes = list(set(clone_sizes)) # uniques the clone_sizes
    clone_sizes = [i for i in clone_sizes if resampled_clone_size_distribution[i] != 0 or observed_clone_size_distribution[i] != 0]
    clone_sizes.sort()

    # make plotfile
    with open(plotfile, "w") as f:
        # write header
        header = "\t".join(["clone_size", "sample", "fit", "lower_limit", "upper_limit"])
        f.write(header+"\n")
        # write clone-size zero
        outline = [0, actual_n0, estimated_n0, lower_n0_error_bar, upper_n0_error_bar]
        outline = map(str, outline)
        outline = "\t".join(outline)
        f.write(outline+"\n")
        if verbose: 
            print header
            print outline
        # write other clone sizes
        for clone_size in clone_sizes:
            if clone_size == 0: continue # we already took care of zero
            outline = [clone_size, observed_clone_size_distribution[clone_size], resampled_clone_size_distribution[clone_size]]
            outline = map(str, outline)
            outline = "\t".join(outline)
            f.write(outline+"\n")
            if verbose: print outline

    # set x and y limits
    if not x_max: 
        x_max = max( round_up_to_the_nearest(max(clone_sizes), 10), 40 )
    if not y_max: 
        y_max = max(observed_clone_size_distribution.values()+resampled_clone_size_distribution.values()+[recon_qD_upper_limit])
        if y_max == float('inf'): y_max = max(observed_clone_size_distribution.values()+resampled_clone_size_distribution.values())
        y_max = log(y_max, 10)
        y_max = round_up_to_the_nearest(y_max, 1)
        y_max = 10**y_max

    # update javascript
    customize_js(javascript_file, plotfile, str(x_max), str(y_max))

    # the absolute path to the .js and .css must be given in the html DONE
    # the absolute path to the plotfile must be given in plot_clone_size_distribution.js

    # Use wkhtmltopdf to convert the html to PDF and cpdf to crop the PDFs
    pdf_file = fitfile.replace(".txt", ".pdf")    
    html_file = fitfile.replace(".txt", ".html")
    write_html_file(html_file)
    print "  writing %s..." % pdf_file
    command = "wkhtmltopdf '%s' '%s'" % (html_file, pdf_file)
    x = getoutput(command)

    """
    Crop the result for easier use in graphics-layout programs
    (e.g. OmniGraffle). Note the parameters below (x, y, width,
    height, with the origin at the bottom left), were chosen for
    dimension_base = 500 in
    plot_clone_size_distribution_ref.js. Smaller values of
    dimension_base will result in white borders with these parameters.
    """
    command = 'cpdf -crop "38mm 185mm 112mm 95mm" \'%s\' -o \'%s\'' % (pdf_file, pdf_file)
    system(command)

    return



# BEGIN EXECUTABLES ====================

if run_recon: 
    print "running recon on %s..." % file_in
    print reconstruct_overall_distribution(file_in, file_out, total_clones_in_overall_population) 
    # total_clones_in_overall_population passed just to record in
    # file_out. When this is a fit of a known overall distribution,
    # e.g. for error-bar fits, we will need this record of how many
    # clones there were in the overall population

elif make_error_bars:
    print "making error-bar file..."
    print make_error_bar_file(file_in, file_out, qs)

elif make_power_table:
    print "making power table..."
    error_envelope_x_vals, error_envelope_y_vals = read_error_bar_file(file_in, q)
    minimum_cells = get_minimum_cells(list_of_number_of_clones, list_of_fold_differences, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, min_number_of_doublets, fraction_small_clones)
    output_table(list_of_number_of_clones, list_of_fold_differences, file_out, number_of_error_bars, min_number_of_doublets, fraction_small_clones, minimum_cells, q, file_out)

elif make_table_of_D_numbers:
    print "making table of D numbers..."
    output_table_of_D_numbers(precomputed_error_bar_file, file_out, list_of_filenames, qs)

elif resample:
    print "resampling..."
    output_resample(file_in, file_out)

elif make_resample_plot:
    print "plotting resample..."
    plot_resample(file_in, file_out, precomputed_error_bar_file)


else:
    print """error: no option given. Options:

-R  --run_recon
-e  --make_error_bars
-p  --make_power_table
-T  --make_table_of_D_numbers
-r  --resample
-x  --make_resample_plot

exiting...
"""
    exit()
