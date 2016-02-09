        
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

(c) 2015-2016 Beth Israel Deaconess Medical Hospital, Inc. All
rights reserved.
