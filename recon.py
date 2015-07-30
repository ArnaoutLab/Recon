# -*- coding: utf-8 -*-
"""
        
    Recon: Reconstruction of Estimated Communities
    from Observed Numbers
    
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
    8. Contact information
    9. License
    
    
    
    1. What is Recon?
    -----------------
    
    Recon is an algorithm for generating a description of a population
    from a sample (see Section 2, Terminology).
    
    More precisely, it uses the distribution of species counts in a sample
    to estimate the distribution of species sizes in the population from
    which the sample was drawn.
    
    It:
    
    - assumes nothing about the frequency distribution of species in
    the population,
    
    - avoids (over)fitting of sampling noise,
    
    - scans many starting points in an attempt to find a global best
    fit, and
    
    - outputs error bars for the number of species of a given size in
    the population.
    
    The code provided in this directory:
    
    - given species frequencies observed in a sample, determines the
    maximum-likelihood species frequency distribution of the
    population without overfitting (see Section 4, Running Recon),
    
    - outputs diversity measurements for the population (see Section
    5, Diversity measures),
    
    - allows new error bar calculations (see Section 6, Error bars),
    and
    
    - performs power calculations of the sample size required to be
    able to detect a given difference between two samples (see
    Section 7, Power calculations).
    
    
    
    2. Requirements
    ---------------
    
    recon.py requires Python 2.7 (https://www.python.org) and the SciPy
    Python library (http://scipy.org). We have tested it on Macintosh OS
    X.
    
    
    
    3. Terminology
    --------------
    
    A "species" is group made up of one or more "individuals." The number
    of individuals of a given species is that species' "size." A "sample"
    is a set of individuals from one or more species that is drawn
    randomly from a "population." Species represented in the population
    that are not represented in the sample are called "missing species."
    
    "Diversity" refers to any of a set of measures of the frequency
    distribution in the population. They can be thought of as effective
    numbers of species in the population. The "Hill numbers" ("qD," read
    "D-q") are a family of diversity measures defined by the parameter q,
    which determines the degree to which diversity measures are weighted
    toward larger species. For example, 0D is "species richness," a
    diversity measure that weights all species equally (and so is just a
    count of the number of different species). infD is the reciprocal of
    the Berger-Parker Index, the effective number of species if all
    species were the size of the largest species. Simple mathematical
    transformations of Shannon entropy (q=1), Gini-Simpson Index (q=2),
    and other common diversity measures correspond to different Hill
    numbers.
    
    
    
    4. Running Recon
    ----------------
    
    Description:
    
    Given a set of observations of species frequency in the sample,
    outputs a set of parameters that describe the maximum-likelihood
    species frequency distribution in the parent population, without
    overfitting.
    
    Usage:
    
    python recon.py [-c -d BIN_SIZE -l PARAMETER_LIMIT -t THRESHOLD -v -r]
    FILE_OUT DATA_FILE
    
    Example:
    
    python recon.py -l 10 -t 30 test.out recon_test_data.txt
    
    Output:
    
    A fit file with parameters that describe the parent population.
    Diversity measurements for the population can be output and displayed
    from this file (see section 5, Displaying results). These parameters
    are called "weights" and "means." The weights w_j give the proportion
    of species in the population that have j individuals. The means are
    Poisson parameters; they give the mean number of individuals a species
    of size j contributes to the sample.
    
    Required praameters:
    
    FILE_OUT    the filename to be used for output. Note that if this file
    exists it will not be overwritten; instead, recon will
    exit with an error message.
    
    DATA_FILE   a text file containing the number of individuals of each
    species seen in the sample (i.e., the species sizes in the
    sample).
    
    The default format is a tab-delimited file with lines of
    the form
    
    species name <tab> species size
    
    with a newline character delimiting lines. Species size is
    an integer. Sample data in this format can be found in the
    file recon_test_data.txt; in it, species1 has 1
    individual, species2 has 1 individual, etc. up to
    species100; species 101 has 2 individuals, etc.
    
    Options:
    
    -c          The -c option allows recon to read an alternative
    tab-delimited format with lines of the form
    
    species size <tab> number of species of this size
    
    Both values are integers. Example data in this format can
    be found in recon_test_data_distribution_in_file.txt; this
    is the same data as in recon_test_data.txt, just in the
    alternative format. (Note that Recon does not need species
    names.) The number of species of size 1 (i.e., the number
    of species with one individual each) in the sample is 100,
    the number of species of size 2 is 30, etc.
    
    -d BIN_SIZE (default: 1)
    Average number of observations per individual. In many
    circumstances, each individual in the sample will be observed
    and counted once. However, there are cases where each
    individual in the sample will be observed and counted
    multiple times. BIN_SIZE allows for this possibility.
    
    -l PARAMETER_LIMIT (default: 20)
    The maximum number of parameters that the algorithm will
    use to fit the data (default, 20). The recon algorithm
    will continue adding parameters until the (corrected)
    Akaike Information Criterion indicates that additional
    parameters are not justified. In practice the limit of 20
    is rarely reached.
    
    -t THRESHOLD (default: 50)
    The -t option allows you to modify THRESHOLD, the size
    above which sampling error is considered small, which
    means that Recon will assume the frequency of species of
    this size or greater in the population is the same as the
    frequency in the sample. It defaults to 50, but 30 usually
    gives good results. This is because if, in a sample from a
    well mixed population, species A is seen 30 times in a
    sample it is very unlikely that there is another species B
    which is the same size as A in the parent population but
    is missing or very poorly reperesented in the sample.
    Smaller values will give faster run times but less
    accurate results.
    
    -v 
    The -v option produces verbose output. It will print
    the starting and finishing parameters of every step that
    recon takes during the fit. This may typically be tens
    of thousands of lines.
    

    5. Diversity measures
    ---------------------
    
    Description:
    
    Given a fit file (either the output from Section 4 or any population
    description in that format), outputs a table of diversity measures as
    Hill numbers (see Section 3, Terminology).
    
    Usage:
    
    python recon.py -T -b ERROR_BAR_PARAMETERS FILE_OUT FILE_IN1 [FILE_IN2
    FILE_IN3 ...]
    
    Example:
    
    python recon.py -T -b error_bar_params.txt error_test.out test.out
    
    Output:
    
    A table of Hill numbers for the reconstructed distribution, one row
    for each input file (see FILE_IN1 below). The columns labelled "obs"
    show the Hill numbers from the sample. The columns labeled "MLE"
    (maximum-likelihood estimate) show the Hill numbers Recon has
    estimated for the population. Since 0D_MLE is the estimate of the
    species richness in the parent population, the estimated missing
    species is 0D_MLE - 0D_obs.
    
    Required parameters:
    
    -T          Tells Recon to output a table of diversity measures.
    
    -b ERROR_BAR_PARAMETERS
    A file that contains parameters for constructing error
    bars on fits. The supplied file error_bar_params.txt can
    be used. Alternatively, Recon can generate an error bar
    parameter file from a set of gold standard fits.
    
    FILE_OUT    The desired name of the output file.
    
    FILE_IN1    The input file. A fit file output from the fit in section
    4, Running Recon.
    
    Optional parameters:
    
    [FILE_IN2, FILE_IN3 ...]
    These are additional input files. Each file will generate
    one row in the output table.
    
    
    
    6. Error bars
    -------------
    
    Description:
    
    Using Recon to generate an error bar parameter file based on from a
    set of fits on data for which the number of missing species is known.
    
    Usage:
    
    python recon.py -e FILE_OUT ERROR_BAR_FIT_DIRECTORY
    
    Example:
    
    python recon.py -e error_bars_test.out Test_dir
    
    Output:
    
    An error bar parameter file.
    
    Required parameters:
    
    -e          Tells recon to make an error bar parameter file.
    
    FILE_OUT    The name of the new error bar parameter file
    
    ERROR_BAR_FIT_DIRECTORY
    The name of a directory that contains the fits with known
    missing species. The known missing species are encoded in
    the weights and means of the population.
    
    
    7. Power calculations
    ---------------------
    
    Description:
    
    Generates a power table with the minimum sample size required to be
    able to detect differences of a given magnitude in a given Hill number
    for two populations. That is, if you have two populations, and want to
    be able detect a difference in 1D of x%, the table tells you how big
    your samples have to be.
    
    Usage:
    
    python recon.py -q HILL_NUMBER [-m MINIMUM_NUMBER_OF_DOUBLETS]
    FILE_OUT ERROR_BAR_PARAMETERS
    
    Example:
    
    python recon.py -p -q 0 power.out error_bar_params.txt
    
    Output:
    
    The rows and columns of the output table show the minimum fold
    differences that can be detected for different population sizes. At
    present these are hardcoded (lines 1605 and 1611); future versions
    will implement command-line control.
    
    Required parameters:
    
    -q HILL_NUMBER
    The Hill number for which the power calculation is carried
    out.
    
    FILE_OUT    The desired name of the output file.
    
    ERROR_BAR_PARAMETERS
    A file that contains parameters for constructing error
    bars on fits. The supplied file error_bar_params.txt can
    be used. Alternatively, Recon can generate an error bar
    parameter file from a set of gold standard fits.
    
    Optional parameters:
    
    MINIMUM_NUMBER_OF_DOUBLETS
    is an additional statistical minimum required for good
    results. The default of 100 should be good for most
    purposes.

    8. Resampling a fit
    -------------------
    
    Description:
    
    This allows resampling a model. The output gives the maximum likelihood
    observed species size distribution of samples from the model. This is the 
    distribution that Recon attempts to make as close as possible to the 
    observed species size distribution. The closeness of the fits can be 
    compared to measure the goodness of fit.
    
    Usage:
        
    python recon.py -r FILE_OUT FILE_IN
    
    Example:
        
    python recon.py -r resample.out test.out 
    
    Output:
    
    A list of species sizes up to the threshold that was used in the original
    fit is outputted, together with a count of species for each size. Output 
    is written both to standard output and to FILE_OUT.
                    
    Required parameters:
    
    FILE_OUT
    The desired name of the output file.
    
    FILE_IN
    A file that contains model fitted parameters as output from a previously
    completely Recon fit.
    

    9. Contact information
    ----------------------
    
    Correspondence should be addressed to Ramy Arnaout at
    rarnaout@gmail.com.
    
    
    10. License
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
    
    Last Updated June 23, 2015
    (c) 2015 Beth Israel Deaconess Medical Hospital, Inc. All rights reserved.

"""

from sys import argv, exit
from scipy.stats import poisson
from string import split, join
from collections import defaultdict
from os.path import expanduser, isfile, isdir
from os import listdir
from math import isnan, log, exp, ceil
from scipy.optimize import fsolve, minimize
from datetime import datetime
from copy import deepcopy

import numpy as np
import argparse
import ast
import warnings

warnings.filterwarnings('error') # This makes RunTime warnings (e.g. divide by zero) raise an error

class mixed_distribution:
    """
    A mixed_distribution object represents a mixed Poisson
    distribution.
    
    A mixed_distribution is defined by a sum of Poisson distributions
    where the weight and mean of each subdistribution must be specified
    and the weights must sum to one.
    """
    #
    def __init__(self,weight_mean_pairs):
        #
        total_weight = float(sum(k[0] for k in weight_mean_pairs)) # Normalise weights to 1.0
        #
        self.weights = [k[0]/total_weight for k in weight_mean_pairs]
        #
        self.means = [float(k[1]) for k in weight_mean_pairs]
    #
    def value_at_clone_size(self,i):
        #
        vals = [weight*poisson(mean).pmf(i) for weight,mean in zip(self.weights,self.means)]
        #vals = [0.0 if isnan(val) else val for val in vals] Makes no difference
        #
        return sum(vals)


class sample_from_parent:
    def __init__(self, weights, means, sample_size, total_clones, population_size, total_observed_clones, estimated_n0, fitted_weights, fitted_means):
        #
        total_weight = float(sum(weights)) # Normalise weights to 1.0
        #
        self.weights = [weight/total_weight for weight in weights]
        self.means = means
        self.sample_size = sample_size
        self.total_clones = total_clones
        self.population_size = population_size
        self.fitted_weights = fitted_weights
        self.fitted_means = fitted_means
        self.total_observed_clones = total_observed_clones
        self.estimated_n0 = estimated_n0


# =========== 1. BEGIN DEFINE FUNCTIONS ==================

# =========== 1.1 BEGIN DEFINE FUNCTIONS FOR MLE FITS ==================

def read_observed_clone_size_distribution(file_in, clone_distribution_in_file, threshold, bin_size = None):
    #
    observed_clone_size_distribution = defaultdict(int)
    #
    with open(file_in,'rU') as f:
        if clone_distribution_in_file: # Format from literature
            for line in f:
                if line[0] != '#':
                    observed_clone_size_distribution[int(line.split('\t')[0])] = int(line.split('\t')[1])
        #
        elif bin_size:
            #
            bin_size = float(bin_size)
            #
            # Convert threshold on reads
            #
            threshold = ceil(threshold*bin_size)
            #
            for line in f: # Format from experimental reads
                if line[0] != '#':
                    normed_clone_size = int(line.split('\t')[1].strip())
                    if normed_clone_size < threshold:
                        binned_clone_size = int(ceil(normed_clone_size / bin_size))
                        """
                        # To retain all clones use these lines:
                        if binned_clone_size == 0:
                            binned_clone_size = 1
                        observed_clone_size_distribution[binned_clone_size] += 1
                        """
                        # To exclude clones in zero binned clone size use these lines:
                        if binned_clone_size != 0:
                            observed_clone_size_distribution[binned_clone_size] += 1
            #
            for clone_size in observed_clone_size_distribution:
                observed_clone_size_distribution[clone_size] = int(observed_clone_size_distribution[clone_size])
            #
            if 0 in observed_clone_size_distribution:
                del observed_clone_size_distribution[0]
    #
    return observed_clone_size_distribution


def calc_log_likelihood_n0_as_data(estimated_clone_size_distribution, new_distribution, threshold):
    #
    log_likelihood = 0.0
    #
    for k in range(0,threshold):
        try:
            log_likelihood += estimated_clone_size_distribution[k]*log(new_distribution.value_at_clone_size(k))
        except ValueError:
            #
            # If new_distribution.value_at_clone_size(k) = 0.0 then the log will throw a ValueError
            #
            log_likelihood = -float('inf')
            #
            if verbose: print 'Log likelihood is -inf!'
    #
    return log_likelihood


def calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, new_distribution, threshold, total_observed_clones):
    #
    log_likelihood = 0.0
    #
    for k in range(1,threshold):
        try:
            log_likelihood += estimated_clone_size_distribution[k]*log(new_distribution.value_at_clone_size(k))
        except ValueError:
            #
            # If new_distribution.value_at_clone_size(k) = 0.0 then the log will throw a ValueError
            #
            log_likelihood = -float('inf')
            #
            if verbose: print 'Log likelihood is -inf!'
    #
    p_0 = new_distribution.value_at_clone_size(0)
    #
    log_likelihood = log_likelihood - total_observed_clones * log(1.0 - p_0)
    #
    return log_likelihood


def log_likelihood_score(fit_params, estimated_clone_size_distribution):
    #
    weights = list(fit_params[:len(fit_params)/2])
    #
    if any(weight >0.0 for weight in weights):
        weights = [max(weight,0.0) for weight in weights] # Ensure this is a local variable!
    else:
        weights = [float(weight == max(weights)) for weight in weights]
        if verbose: print 'All weights negative', weights
    #
    weight_norm = sum(weights)
    #
    try:
        weights = [weight / weight_norm for weight in weights]
        #
    except ZeroDivisionError:
        # if verbose: print '==== ERROR IN NORM ===='
        # if verbose: print weights
        # if verbose: print '======================='
        #
        weights = [1.0] + weights[1:]
    #
    weights = [1.0 - sum(weights[1:])] + weights[1:]
    #
    means = list(fit_params[len(fit_params)/2:]) # Ensure this is a local variable!
    #
    means = [max(mean,1e-6) for mean in means]
    #
    test_distribution = mixed_distribution(zip(weights,means))
    #
    log_likelihood = calc_log_likelihood_n0_as_data(estimated_clone_size_distribution, test_distribution, threshold)
    #
    return -log_likelihood


def MLE_fit_v2(file_out, threshold, new_weights, new_means, observed_clone_size_distribution):
    """
    v2 does not return the initial parameters.
        
    The code which call MLE_fit_v2 already knows these.
    """
    #
    if verbose: print '==== BEGIN FIT ==='
    if verbose: print
    if verbose: print 'initial parameters '+str(new_weights+new_means)
    #
    with open(file_out,'a') as f:
        f.write('initial parameters = '+str(new_weights+new_means)+'\n')
    #
    total_observed_clones = float(sum([observed_clone_size_distribution[key] for key in observed_clone_size_distribution if key < threshold])) # Changed in v9 103014JK
    #
    new_distribution = mixed_distribution(zip(new_weights,new_means))
    #
    total_estimated_clones = round(total_observed_clones/(1.0-new_distribution.value_at_clone_size(0)))
    #
    estimated_n_0 = total_estimated_clones - total_observed_clones
    #
    estimated_clone_size_distribution = defaultdict(int)
    estimated_clone_size_distribution[0] = estimated_n_0
    #
    for key in observed_clone_size_distribution:
        estimated_clone_size_distribution[key] = observed_clone_size_distribution[key]
    #
    count_iterations = 0
    #
    # score_function1 stops the iteration when the
    # estimated clone distribution is at a fixed point.
    #
    score_function1 = cutoff_score_1 +1.0
    #
    dictionary_of_fits = {}
    #
    iteration_limit = 20 # Was 15 for UB fits #10 Is this ever limiting for the fits to data? I should write this to file
    #
    while ( (score_function1 > cutoff_score_1) and (count_iterations < iteration_limit) ): # I could implement this as separate cutoffs for weights and means. Begin with one cutoff.
        #
        if verbose: print 'iterating clone sizes...' # This output indicates start of a new iteration
        #
        count_iterations += 1
        #
        if verbose: print 'iteration ',count_iterations
        #
        # Step 1
        #
        # Use the estimated_clone_distribution to
        # calculate an updated mixed distribution
        #
        # This is the MLE of weights and means given the estimated
        # clone distribution.
        #
        # scorefunction2 stops the iteration when the weights and
        # means are at a fixed point.
        #
        if verbose: print 'old estimated n0','\t',estimated_n_0 # This output gives the estimated clones at the start of the iteration
        #
        old_estimated_n_0 = estimated_n_0
        #
        start_params = np.array(new_weights+new_means)
        #
        try:
            res = minimize(log_likelihood_score,start_params,args=(estimated_clone_size_distribution,), method = 'L-BFGS-B') # methods without bouds and callback
            #
            new_weights = list(res.x[:len(res.x)/2])
            #
            if any(weight >0 for weight in new_weights): # Should this be flagged. Maybe useful to record in file?
                new_weights = [max(weight,0.0) for weight in new_weights]
            else:
                new_weights = [float(weight == max(new_weights)) for weight in new_weights]
            #
            weight_norm = sum(new_weights)
            new_weights = [weight / weight_norm for weight in new_weights]
            #
            new_means = list(res.x[len(res.x)/2:])
            if any(new_mean < 0.0 for new_mean in new_means):
                new_means = [max(new_mean,1e-6) for new_mean in new_means] # Is this a sign that the minimisation should stop? Should this threshold be fixed systematically?
            #
        except ValueError:
            if verbose: print 'ValueError in minimize!'
            new_weights = start_params[:len(start_params)/2]
            new_means = start_params[len(start_params)/2:]
        #
        if verbose: print 'Fitted parameters',new_weights,new_means
        #
        # Step 2
        #
        # Calculate an updated estimated_clone_distribution
        # based on old_distribution and observed_clone_size_distribution
        #
        new_distribution = mixed_distribution(zip(new_weights,new_means))
        #
        total_estimated_clones = round(total_observed_clones / (1.0 - new_distribution.value_at_clone_size(0)))
        #
        estimated_n_0 = round(total_estimated_clones - total_observed_clones)
        #
        if verbose: print 'new estimated n_0', estimated_n_0
        #
        estimated_clone_size_distribution[0] = estimated_n_0
        #
        # NB we could try looking at the whole of the estimated distribution, not just estimated_n0, as done here.
        #
        try:
            score_function1 = 2.0*abs(estimated_n_0 - old_estimated_n_0)/float(estimated_n_0 + old_estimated_n_0)
        except ZeroDivisionError:
            score_function1 = 0.0
        #
        if verbose: print 'score_function1', score_function1
        #
        log_likelihood = calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, new_distribution, threshold, total_observed_clones)
        #
        if verbose: print 'Log likelihood\t',log_likelihood # Note that this is the log likelihood POST update of estimated_n_0
        if verbose: print
    #
    dictionary_of_fits[log_likelihood] = new_weights+new_means # Take only self consistent values!
    #
    if count_iterations == iteration_limit:
        with open(file_out,'a') as f:
            f.write('Iteration limit reached - check for only small numbers of clones appearing once or twice in sample\n')
        if verbose: print 'Iteration limit reached - check for only small numbers of clones appearing once or twice in sample'
    #
    if verbose: print 'Completed iteration of clone sizes'
    #
    new_weights = dictionary_of_fits[max(dictionary_of_fits.keys())][:len(dictionary_of_fits.values()[0])/2]
    new_means = dictionary_of_fits[max(dictionary_of_fits.keys())][len(dictionary_of_fits.values()[0])/2:]
    new_distribution = mixed_distribution(zip(new_weights,new_means))
    total_estimated_clones = round(total_observed_clones / (1.0 - new_distribution.value_at_clone_size(0)))
    #
    estimated_n_0 = round(total_estimated_clones - total_observed_clones)
    #
    if verbose: print 'Fitted weights','\t',new_distribution.weights
    if verbose: print 'Fitted means','\t',new_distribution.means
    #
    if verbose: print 'Final estimated_n0\t',estimated_n_0
    if verbose: print 'Log likelihood\t',log_likelihood
    #
    with open(file_out,'a') as f:
        f.write('n_obs '+str(total_observed_clones)+'\n')
        f.write('weights = '+str(new_distribution.weights)+'\n')
        f.write('means = '+str(new_distribution.means)+'\n')
        f.write('log likelihood = '+str(log_likelihood)+'\n')
        f.write('estimated_n0 = '+str(estimated_n_0)+'\n')
        f.write('\n')
    #
    if verbose: print 'Observed clones',sum(observed_clone_size_distribution.values())
    if verbose: print 'Observed clones + estimated clones', sum(observed_clone_size_distribution.values()) + estimated_n_0
    if verbose: print '==== END FIT ==='
    if verbose: print
    return new_distribution.weights+new_distribution.means, log_likelihood, estimated_n_0


def carry_out_fit(observed_clone_size_distribution, file_out):
    #
    #
    # Initialise starting weights and means.
    # For the first round of fitting this should replicate the scan for two spikes.
    #
    # Initial mean depends only on the observed data:
    #
    initial_mean = sum(ii * observed_clone_size_distribution[ii] for ii in observed_clone_size_distribution) / float(sum(observed_clone_size_distribution.values()))
    #
    # starting_means will be the list of starting points for means
    #
    starting_means = [initial_mean]
    starting_weights = [1.0]
    #
    number_of_parameters = 1
    #
    #
    # To start a fit at an arbitrary point, enter fit values here:
    #
    # e.g.
    #
    # starting_weights = [0.89743138679717305, 0.079955794343904069, 0.014713942471041609, 0.0052554039266035982, 0.002643472461277707]
    # starting_means = [0.4124878712022608, 3.0833764053806996, 9.187883397378856, 20.31890394714073, 37.0149387211759]
    #
    # number_of_parameters = 9
    #
    total_observed_clones = float(sum(observed_clone_size_distribution.values())) # NB this is used as a global variable - be careful!
    #
    new_distribution = mixed_distribution(zip(starting_weights,starting_means))
    #
    total_estimated_clones = round(total_observed_clones/(1.0-new_distribution.value_at_clone_size(0)))
    #
    estimated_n_0 = total_estimated_clones - total_observed_clones
    #
    estimated_clone_size_distribution = defaultdict(int)
    estimated_clone_size_distribution[0] = estimated_n_0
    #
    for key in observed_clone_size_distribution:
        estimated_clone_size_distribution[key] = observed_clone_size_distribution[key]
    
    #
    log_likelihood = calc_log_likelihood_n0_as_parameter(estimated_clone_size_distribution, new_distribution, threshold, total_observed_clones)
    #
    AIC =2.0* (number_of_parameters - log_likelihood)
    if len(observed_clone_size_distribution.keys()) > 2:
        AICc = AIC +(2.0*number_of_parameters*(number_of_parameters+1)/ (len(observed_clone_size_distribution.keys()) - number_of_parameters - 1.0) )
        too_few_clone_sizes = False
    else:
        AICc = float('inf')
        too_few_clone_sizes = True
    
    AICc_old = float('inf')

    with open(file_out,'a') as f:
        f.write('fitted_clone_size_distribution = '+str(dict(observed_clone_size_distribution))+'\n')
        f.write('weights = [1.0]\n')
        f.write('means = [0.0]\n')
        f.write('initial_scale_factors = '+str(initial_scale_factors)+'\n')
        f.write('test_weights_list ='+str(test_weights_list)+'\n\n')

    minimum_observed_clone_size = min(key for key in observed_clone_size_distribution if observed_clone_size_distribution[key]!=0)

    zero_in_fitted_weights = False
    #
    # n0 is estimated above with the starting mean from data.
    #
    # Here we reinitialise starting weights
    #
    starting_means = []
    starting_weights = []
    #
    while (((AICc < AICc_old) or too_few_clone_sizes) and (number_of_parameters < parameter_limit) and (not zero_in_fitted_weights) and (1e-6 not in starting_means)):
        #
        # The list of fits will have one element for each scan starting point.
        #
        list_of_fits = []
        #
        for initial_scale_factor in initial_scale_factors:
            #
            # The next loop will add a new weight 
            #
            if len(starting_weights) == 0:
                weights_to_loop_over = [1.0]
            else:
                weights_to_loop_over = test_weights_list
            #
            for test_weight in weights_to_loop_over:
                #
                if verbose: print 'Scale factor', initial_scale_factor
                #
                starting_means_temp = [ initial_scale_factor * initial_mean ] + starting_means
                starting_weights_temp = [test_weight] + [weight*(1.0 - test_weight) for weight in starting_weights]
                #
                if verbose: print starting_means_temp, starting_weights_temp
                #
                end_params, log_likelihood, estimated_n_0 = MLE_fit_v2(file_out, threshold, starting_weights_temp, starting_means_temp, observed_clone_size_distribution)
                #
                list_of_fits.append((starting_weights_temp+starting_means_temp, end_params, log_likelihood, estimated_n_0))
        #
        sorted_log_likelihoods = sorted(fit[2] for fit in list_of_fits)
        log_likelihoods = [element[2] for element in list_of_fits]
        max_index = log_likelihoods.index(sorted_log_likelihoods[-1])
        second_index = log_likelihoods.index(sorted_log_likelihoods[-2])
        average_start_parameters = [(p1+p2)/2. for p1,p2 in zip(list_of_fits[max_index][0], list_of_fits[second_index][0])]
        #
        end_params, log_likelihood, estimated_n_0 = MLE_fit_v2(file_out, threshold, average_start_parameters[:len(average_start_parameters)/2], average_start_parameters[len(average_start_parameters)/2:], observed_clone_size_distribution)
        list_of_fits.append((average_start_parameters, end_params, log_likelihood, estimated_n_0))
        #
        noise_test_ratio = 0.0
        #
        while noise_test_ratio <= 3.:
            #
            max_log_likelihood = max(fit[2] for fit in list_of_fits)
            log_likelihoods = [element[2] for element in list_of_fits]
            max_index = log_likelihoods.index(max(log_likelihoods))
            #
            fitted_parameters = list_of_fits[max_index][1]
            #
            # Set up starting parameters for the next round of fitting:
            #
            starting_weights = fitted_parameters[:len(fitted_parameters)/2]
            starting_means = fitted_parameters[len(fitted_parameters)/2:]
            #
            # noise_test_ratio compares the expected contribution to the sample of the smallest fitted clones
            # to the noise in the sample from remaining clones.
            #
            sorted_params = sorted(zip(starting_weights, starting_means), key = lambda x: x[1])
            sorted_weights, sorted_means = [[x[i] for x in sorted_params] for i in [0,1]]
            #
            if (0.0 in sorted_weights) and (len(sorted_weights) == 2):
                break
            #
            # NB In this version, total_observed_clones is used. This should actually
            # be total_observed_clones + estimated_n0. That means that the size of the noise
            # is overestimated (potentially around threefold is n0 is an order of magnitude
            # greater than total_observed_clones). The result is that this test is somewhat
            # harsher than it should be. 120714JK
            #
            if len(sorted_weights) > 1: # This applies on the second and subsequent interations, when there is more than one spike
                noise_test_ratio = sorted_weights[0]*poisson(sorted_means[0]).pmf(minimum_observed_clone_size) / np.sqrt(sum( [weight * poisson(mean).pmf(minimum_observed_clone_size) for weight, mean in zip(sorted_weights[1:], sorted_means[1:]) ] ) / total_observed_clones)
            else:
                noise_test_ratio = 4. # This applies on the first iteration, when there is only a single spike
            #
            if noise_test_ratio <= 3.:
                del list_of_fits[max_index]
                if len(list_of_fits) == 0:
                    max_log_likelihood = -float('inf')
                    break
        #
        if 0.0 in fitted_parameters:
            #
            # Set flag to stop further iterations
            #
            zero_in_fitted_weights = True
            #
            # Remove zero from parameters
            #
            filtered_weights = [w for w in sorted_weights if w != 0.0]
            filtered_means = [sorted_means[i] for i in range(len(sorted_weights)) if sorted_weights[i] != 0.0]
            fitted_parameters = filtered_weights + filtered_means
        #
        # AIC is the Aikaike Information Criterion and AICc the corrected criterion for small numbers of observations
        #
        number_of_parameters = len(fitted_parameters) -1
        AIC_old = AIC
        AICc_old = AICc
        #
        AIC =2.0* (number_of_parameters - max_log_likelihood)
        if ((min(len(observed_clone_size_distribution.keys()), threshold) - number_of_parameters - 1.0) > 0.0):
            AICc = AIC +(2.0*number_of_parameters*(number_of_parameters+1)/ (min(len(observed_clone_size_distribution.keys()), threshold) - number_of_parameters - 1.0) )
        else:
            AICc = float('inf')
        #
        if 0.0 in fitted_parameters:
            AICc = float('inf')
        #
        if (AICc < AICc_old) or too_few_clone_sizes:
            with open(file_out,'a') as f:
                f.write('\n===========================================\n')
                if number_of_parameters >= parameter_limit:
                    f.write('parameter_limit = '+str(parameter_limit)+' # REACHED MAX NUMBER OF PARAMETERS\n')
                f.write('intial parameters = '+str(list_of_fits[max_index][0])+'\n')
                f.write('fitted parameters = '+str(list_of_fits[max_index][1])+'\n')
                f.write('estimated n0 = '+str(list_of_fits[max_index][3])+'\n')
                f.write('AIC = '+str(AIC)+'\n')
                f.write('AICc = '+str(AICc)+'\n')
                f.write('log likelihood = '+str(list_of_fits[max_index][2])+'\n')
                f.write('===========================================\n')
            too_few_clone_sizes = False
    #
    n = len(list_of_fits[max_index][1])/2
    #
    print 'fitted weights',list_of_fits[max_index][1][:n]
    print 'fitted means',list_of_fits[max_index][1][n:]
    print 'estimated missing species',list_of_fits[max_index][3]
    print 'estimated total species of size less than threshold',list_of_fits[max_index][3]+total_observed_clones


# =========== 1.1 END FUNCTIONS FOR MLE FITS ==================

# =========== 1.2 BEGIN DEFINE FUNCTIONS FOR MAKING ERROR BAR FILES ==================

def read_sample_fit(MLE_filename, observed_threshold = 1000, return_threshold = False):
    #
    #This can easily return starting parameters.
    #
    #But can starting parameters tell which scan was selected?
    #
    #For weights, yes. This will be initial_parameters[0]
    #
    #What about for means? Initial_means can be calculated from the starting parameters, so this will also be known.
    #
    #initial_mean = sum(ii * observed_clone_size_distribution[ii] for ii in observed_clone_size_distribution) / float(sum(observed_clone_size_distribution.values()))
    #
    #scale_factor = initial_parameters[len(initial_parameters)/2] / initial_mean
    #
    with open(MLE_filename, 'rU') as MLE_file:
        try:
            in_best_fit = False
            in_start_params = True
            for line in MLE_file:
                if ('threshold' in line) or ('upper_limit' in line):
                    observed_threshold = int(line.split(' = ')[1])
                if in_start_params:
                    if line == '\n':
                        in_start_params = False
                    else:
                        if line[0] != '#':
                            line = line.split('#')[0]
                            if line.split('=')[0].strip() == 'parent_population_size':
                                parent_population_size = ast.literal_eval(line.split('=')[1].strip())
                            elif line.split('=')[0].strip() == 'total_clones':
                                total_clones = ast.literal_eval(line.split('=')[1].strip())
                            elif line.split('=')[0].strip() == 'sample_size':
                                sample_size= ast.literal_eval(line.split('=')[1].strip())
                            elif line.split('=')[0].strip() == 'weights':
                                weights = ast.literal_eval(line.split('=')[1].strip())
                            elif line.split('=')[0].strip() == 'means':
                                means= ast.literal_eval(line.split('=')[1].strip())                            
                            elif line.split('=')[0].strip() == 'observed_clone_size_distribution':
                                observed_clone_size_distribution= ast.literal_eval(line.split('=')[1].strip())                            
                            elif line.split('=')[0].strip() == 'total_observed_clones':
                                total_observed_clones = ast.literal_eval(line.split('=')[1].strip())                            
                            elif line.split('=')[0].strip() == 'initial_scale_factors':
                                initial_scale_factors = ast.literal_eval(line.split('=')[1].strip())                            
                            elif line.split('=')[0].strip() == 'test_weights':
                                test_weights = ast.literal_eval(line.split('=')[1].strip())
                elif line[0] == '=':
                    in_best_fit = not in_best_fit
                elif in_best_fit:
                    if line.split(' = ')[0] == 'intial parameters':
                        initial_parameters_try = eval(line.split(' = ')[1].strip())
                    if line.split(' = ')[0] == 'fitted parameters':
                        fitted_params_try = eval(line.split(' = ')[1].strip())
                    if line.split(' = ')[0] == 'estimated n0':
                        estimated_n0_try = int(float(line.split(' = ')[1].strip()))
                        if 1e-6 in fitted_params_try:
                            break
                        else:
                            initial_parameters = initial_parameters_try
                            fitted_params = fitted_params_try
                            estimated_n0 = estimated_n0_try
            #
            if 'observed_clone_size_distribution' not in locals():
                print 'observed_clone_size_distribution not in MLE_fit file!'
                exit()
            #
            if 'sample_size' not in locals():
                sample_size = sum([key*observed_clone_size_distribution[key] for key in observed_clone_size_distribution])
            #
            if 'total_observed_clones' not in locals():
                total_observed_clones = sum(observed_clone_size_distribution.values())
            
            if 'total_clones' not in locals():
                total_clones = total_observed_clones+estimated_n0
            #
            if 'parent_population_size' not in locals():
                parent_population_size = 1e9
            #
            initial_mean_norm = float(sum([observed_clone_size_distribution[ii] for ii in observed_clone_size_distribution if ii < observed_threshold]))
            #
            initial_mean = sum([ii * observed_clone_size_distribution[ii] for ii in observed_clone_size_distribution if ii < observed_threshold]) / initial_mean_norm
            #
            start_mean = initial_parameters[len(initial_parameters)/2] / initial_mean
            start_mean = int(10*start_mean)/10.
            start_weight = initial_parameters[0]
            #
            if return_threshold:
                return observed_threshold, sample_from_parent(weights, means, sample_size, total_clones, parent_population_size, total_observed_clones, estimated_n0, fitted_params[:len(fitted_params)/2], fitted_params[len(fitted_params)/2:]), start_weight, start_mean
            else:
                return sample_from_parent(weights, means, sample_size, total_clones, parent_population_size, total_observed_clones, estimated_n0, fitted_params[:len(fitted_params)/2], fitted_params[len(fitted_params)/2:]), start_weight, start_mean
        #
        except ValueError:
            if verbose: print MLE_filename,' Nan in fit'
            return 'Nan in fit', 0, 0


def D_number_from_parameters_with_observed4(weights, means, observed_clone_size_distribution, observed_threshold, q):
    """
        total_estimated_clones to be passed to this function should NOT be the sum of observed clones + estimated_0
        
        It should be the sum of all clones arising from the fit parameters.
        
        This version incorporates various speed ups over version 3.
        """
    #
    max_clone_size = max(observed_clone_size_distribution.keys())
    #
    fitted_distribution = mixed_distribution(zip(weights,means))
    #
    obs_small_clones = sum([observed_clone_size_distribution[x] for x in range(1,observed_threshold) if x in observed_clone_size_distribution])
    distribution_norm = obs_small_clones / sum([fitted_distribution.value_at_clone_size(x) for x in range(1,observed_threshold)])
    #
    # The following two lines are logically equivalent to tec_obs = distribution_norm:
    #
    # estimated_n0 = distribution_norm * fitted_distribution.value_at_clone_size(0)
    # tec_obs = estimated_n0 + distribution_norm * sum([fitted_distribution.value_at_clone_size(x) for x in range(1,max_clone_size)])
    #
    # On testing, there is minimal numerical difference (<1 clone)
    #
    tec_obs = distribution_norm # includes estimated_n0
    #
    for x in range(1,observed_threshold):
        if x not in observed_clone_size_distribution:
            observed_clone_size_distribution[x] = 0
    #
    fit_large_clones = {}
    #
    unfitted_clones = {}
    #
    # NB in the sum below I should only take keys that are already in observed_clone_size_distribution.
    # Otherwise, as a defaultdict it will create zero entries.
    #
    for clone_size in range(observed_threshold, max(observed_clone_size_distribution.keys())+1):
        if clone_size in observed_clone_size_distribution:
            fit_large_clones[clone_size] = distribution_norm * fitted_distribution.value_at_clone_size(clone_size)
            diff = observed_clone_size_distribution[clone_size] - fit_large_clones[clone_size]
            if diff > 0.0:
                unfitted_clones[clone_size] = diff
    #
    #if verbose: print 'Number of unfitted clones', len(unfitted_clones.keys())
    #
    aug_total_estimated_clones = tec_obs + sum(unfitted_clones.values())
    #
    aug_total_estimated_clones = float(aug_total_estimated_clones)
    #
    scale_weights = (tec_obs / aug_total_estimated_clones)
    #
    weights = [weight*scale_weights for weight in weights]
    #
    aug_weights = weights + [(unfitted_clones[clone_size] / aug_total_estimated_clones) for clone_size in unfitted_clones.keys()]
    #
    #if verbose: print 'aug weights should sum to one:',sum(aug_weights)
    #
    aug_means = means + [clone_size for clone_size in unfitted_clones.keys()]    #
    qD = D_number_from_parameters(aug_weights, aug_means, aug_total_estimated_clones, q)
    #
    # What is the post-augmentation estimated_n0?
    #
    # aug_fitted_distribution = mixed_distribution(zip(aug_weights,aug_means))
    #
    # aug_estimated_n0 = aug_total_estimated_clones * aug_fitted_distribution.value_at_clone_size(0)
    #
    # if verbose: print estimated_n0,'\t',aug_estimated_n0,'\t',estimated_n0 - aug_estimated_n0
    #
    # This shows minimal diffence in estimated_n0 after we add in the augmented clones, typically less than 1.
    #
    return qD


def D_number_from_parameters(weights,means,total_estimated_clones,q):
    """
        I should add in an estimate from the true distribution for values
        above the upper limit.
        
        This is slightly awkward, because I am generating the results directly from the
        parameters.
        
        I need to calculate the observed clone size distribution, work out out which clones
        are larger than the observed cut off and then add these to the reconstructed
        parent distribution.
        
        I should have the code for this alread in MLE_fit_v6.py
        """
    #
    if q == 0.0:
        return total_estimated_clones
    #
    mean_norm = 1.0 / ( sum(weight*mean for weight,mean in zip(weights, means)) * total_estimated_clones )
    #
    if (q != 1.0) and (q != float('inf')):
        #
        qD = sum(weight * (mean * mean_norm)**q for weight, mean in zip(weights, means))
        #
        qD = qD * total_estimated_clones
        #
        qD = qD ** (1.0/(1.0 - q))
    #
    elif q == 1.0:
        #
        qD = -sum(weight * (mean * mean_norm) * log(mean * mean_norm) for weight, mean in zip(weights, means))
        #
        qD = exp(qD * total_estimated_clones)
    #
    elif q == float('inf'):
        #
        qD = 1.0 / max((mean * mean_norm) for mean in means)
    #
    return qD


def sample_from_model(test_distribution, sample_size, parent_population_size, total_clones, return_unobserved_clones = False):
    """
        We are given the mean number of cells that a single clone
        contributes to a sample whose size is parent_population_size.
        
        This needs to be scaled to a sample whose size is sample_size.
        """
    #
    sample_size = int(sample_size) # Robust to user input as float
    total_clones = int(total_clones) # Robust to user input as float
    #
    weights = test_distribution['weights']
    means = test_distribution['means']
    #
    scale_factor = sample_size / float(parent_population_size)
    #
    means = [mean * scale_factor  for mean in means] # This scales down to find the mean for a sample
    #
    sample_distribution = mixed_distribution(zip(weights,means))
    #
    observed_clone_size_distribution = {}
    #
    clone_size = 1
    total_observed_clones = 0
    #
    unobserved_clones = total_clones * sample_distribution.value_at_clone_size(0)
    clone_limit = total_clones - unobserved_clones - 1 # When we observe this many clones we can stop looking.
    #
    while ( total_observed_clones < clone_limit ):
        new_clones = total_clones * sample_distribution.value_at_clone_size(clone_size)
        #
        observed_clone_size_distribution[clone_size] = int(round(new_clones))
        total_observed_clones += new_clones
        #
        clone_size += 1
    #
    if return_unobserved_clones:
        return observed_clone_size_distribution, unobserved_clones
    else:
        return observed_clone_size_distribution


def error_bar_devs_q1(samples, q):
    """
        This function considers all examples of fits at constant
        true total clone and returns:
        
        error_envelope_x_vals_q
        error_envelope_y_vals_q
        
        The keys of error_envelope_x_vals and error_envelope_x_vals are the true qD of the samples
        which we are using in our error bar fits.
        
        For error_envelope_x_vals at a fixed true qD the values are an ordered list of sample sizes.
        
        For error_envelope_y_vals at a fixed true qD the values are max(abs(log(fitted_D / true_D))) at the
        sample size given by the corresponding error_envelope_x_vals entry. The max is taken over
        steepnesses (that is, distribution shape).
        
        abs_deviations
        
        a list of (sample_size, abs(log(fitted_D / true_D)) pairs.
        """
    #
    error_envelope_x_vals_q = {}
    error_envelope_y_vals_q ={}
    #
    abs_deviations = []
    #
    for sample in samples:
        #
        observed_clone_size_distribution = sample_from_model({'weights':sample.weights ,'means':sample.means}, sample.sample_size, sample.population_size, sample.total_clones)
        #
        #try:
        #
        fitted_D = D_number_from_parameters_with_observed4(sample.fitted_weights, sample.fitted_means, observed_clone_size_distribution, 30, q)
        #
        true_D = D_number_from_parameters(sample.weights, sample.means, sample.total_clones, q)
        #
        if (fitted_D != 0.0) and (true_D != 0.0):
            abs_deviations.append((sample.sample_size, abs(log(fitted_D / true_D)), true_D))
        #
        #except RuntimeWarning:
        #    print 'Divide by zero encountered:'
        #    print
        #    print 'sample total clones', sample.total_clones
        #    print 'sample steepness', len(sample.weights)
        #    print 'sample size',sample.sample_size
        #    print
    #
    abs_deviations = [d for d in abs_deviations if not isnan(d[1])]
    #
    true_Ds = sorted(list(set([abs_dev[2] for abs_dev in abs_deviations])))
    print 'true qDs in samples ',true_Ds
    #
    for true_D in true_Ds:
        #
        sorted_dev = [dev for dev in abs_deviations if dev[2] == true_D]
        #
        sorted_dev = sorted(sorted_dev, key = lambda x: x[0])
        #
        x_vals = []
        y_vals = []
        #
        max_y_val = 0.0
        for val in sorted_dev[::-1]:
            if val[1] > max_y_val:
                #
                max_y_val = val[1]
                #
                if x_vals != [] and val[0] == x_vals[-1]:
                    y_vals[-1] = val[1]
                else:
                    x_vals.append(val[0])
                    y_vals.append(val[1])
        #
        x_vals = x_vals[::-1]
        y_vals = y_vals[::-1]
        #
        # For purposes of interpolation, set the error at 0 sample size equal to the error at smallest available sample size
        #
        if x_vals[0] != 0:
            x_vals = [0] + x_vals
            y_vals = [y_vals[0]] + y_vals
        #
        error_envelope_x_vals_q[true_D] = x_vals
        error_envelope_y_vals_q[true_D] = y_vals
    #
    return error_envelope_x_vals_q, error_envelope_y_vals_q, abs_deviations


def make_error_bar_file(error_bar_fits_directory, precomputed_error_bar_file):
    #
    if isfile(precomputed_error_bar_file):
        print 'Precomputed error bar file', precomputed_error_bar_file,'already exists!'
        exit()
    #
    # This is called with error_bar_fits_directory = file_in and precomputed_error_bar_file = file_out
    #
    MLE_filenames = listdir(error_bar_fits_directory) # This directory should contain the fits
    #
    try:
        MLE_filenames.remove('.DS_Store')
    except:
        if verbose: print 'No .DS_Store file'
    #
    fitted_samples = {} # Keys are filenames, values are sample_from_parent objects
    #
    # Read in fitted data for error bars
    #
    for MLE_filename in MLE_filenames:
        if verbose: print 'processing ',MLE_filename
        try:
            fitted_sample, start_weight, start_mean = read_sample_fit(error_bar_fits_directory+MLE_filename, observed_threshold = 30)
        except UnboundLocalError:
            print 'Incomplete fit', MLE_filename # NB this could also mean that file could not be fit.
        if (fitted_sample != 'Nan in fit') and (fitted_sample != 'Bad fit'):
            fitted_samples[MLE_filename] = fitted_sample
    #
    z = []
    for MLE_filename in fitted_samples:
        z.append(fitted_samples[MLE_filename])
    #
    # The keys of error_envelope_x_vals and error_envelope_x_vals are the true 0D of the samples
    # which we are using in our error bar fits.
    #
    # For error_envelope_x_vals at a fixed true 0D the values are an ordered list of sample sizes.
    #
    # For error_envelope_y_vals at a fixed true 0D the values are max(abs(log(fitted_D / true_D))) at the
    # sample size given by the corresponding error_envelope_x_vals entry. The max is taken over
    # steepnesses (that is, distribution shape).
    #
    # Note that error_bar_devs has a default value of q = 1.
    #
    print datetime.now().strftime('%c')
    print 'Calculating for 0D...'
    #
    error_envelope_x_vals_0, error_envelope_y_vals_0, abs_deviations_0 = error_bar_devs_q1(z, q=0)
    #
    print datetime.now().strftime('%c')
    print 'Calculating for 1D...'
    #
    error_envelope_x_vals_1, error_envelope_y_vals_1, abs_deviations_1 = error_bar_devs_q1(z, q=1)
    #
    print datetime.now().strftime('%c')
    print 'Calculating for 2D...'
    #
    error_envelope_x_vals_2, error_envelope_y_vals_2, abs_deviations_2 = error_bar_devs_q1(z, q=2)
    #
    print datetime.now().strftime('%c')
    print 'Calculating for infD...'
    #
    error_envelope_x_vals_inf, error_envelope_y_vals_inf, abs_deviations_inf = error_bar_devs_q1(z, q=float('inf'))
    #
    # Once these variables have been make they need to be saved.
    #
    with open(precomputed_error_bar_file,'w') as f:
        #
        f.write('error_envelope_x_vals_0 = '+str(error_envelope_x_vals_0)+'\n')
        f.write('error_envelope_y_vals_0 = '+str(error_envelope_y_vals_0)+'\n')
        f.write('abs_deviations_0 = '+str(abs_deviations_0)+'\n')
        #
        f.write('error_envelope_x_vals_1 = '+str(error_envelope_x_vals_1)+'\n')
        f.write('error_envelope_y_vals_1 = '+str(error_envelope_y_vals_1)+'\n')
        f.write('abs_deviations_1 = '+str(abs_deviations_1)+'\n')
        #
        f.write('error_envelope_x_vals_2 = '+str(error_envelope_x_vals_2)+'\n')
        f.write('error_envelope_y_vals_2 = '+str(error_envelope_y_vals_2)+'\n')
        f.write('abs_deviations_2 = '+str(abs_deviations_2)+'\n')
        #
        f.write('error_envelope_x_vals_inf = '+str(error_envelope_x_vals_inf)+'\n')
        f.write('error_envelope_y_vals_inf = '+str(error_envelope_y_vals_inf)+'\n')
        f.write('abs_deviations_inf = '+str(abs_deviations_inf)+'\n')
    #
    return


# =========== 1.2 END DEFINE FUNCTIONS FOR MAKING ERROR BAR FILE ==================

# =========== 1.3 DEFINE FUNCTIONS FOR MAKING TABLE OF D NUMBERS ==================

def error_bar(sample_size, true_total_clone, error_envelope_x_vals_q, error_envelope_y_vals_q): # can I change error_envelope_x_vals to error_envelope_x_vals_q? More, can I just replace this with error_bar_q?
    """
    The function error_bar(sample_size, true_total_clone, error_envelope_x_vals, error_envelope_y_vals)
        
    returns the value of max(abs(fitted_D / true_D)) at the given sample size and given total_clones.
        
    Since we only have values at a grid of true_D and a grid of sample size, this is calculated by interpolation.
        
    This will be an interpolation between
        
        error_bar(sample_size, true_total_clone1)
        
    and
        
        error_bar(sample_size, true_total_clone2)
        
    and then between sample sizes.
    """
    #
    available_total_clone_sizes = sorted(error_envelope_x_vals_q.keys())
    #
    if (true_total_clone > available_total_clone_sizes[-1]):
        if verbose: print 'Clone size out of range'
        return
        #true_total_clone = available_total_clone_sizes[-1]
    if (true_total_clone < available_total_clone_sizes[0]):
        if verbose: print 'Clone size out of range'
        return
        #true_total_clone = available_total_clone_sizes[0]
    #
    for ii in range(len(available_total_clone_sizes)):
        if available_total_clone_sizes[ii] >= true_total_clone:
            break
    #
    weight1 = (available_total_clone_sizes[ii] - true_total_clone) / float(available_total_clone_sizes[ii] - available_total_clone_sizes[ii-1])
    weight2 = (true_total_clone - available_total_clone_sizes[ii-1]) / float(available_total_clone_sizes[ii] - available_total_clone_sizes[ii-1])
    #
    #if verbose: print 'true_total_clone\t', true_total_clone
    #if verbose: print 'lower bound\t', available_total_clone_sizes[ii-1]
    #if verbose: print 'upper bound\t', available_total_clone_sizes[ii]
    #if verbose: print 'weight1\t', weight1
    #if verbose: print 'weight2\t', weight2
    #
    return (weight1 * np.interp(sample_size, error_envelope_x_vals_q[available_total_clone_sizes[ii-1]], error_envelope_y_vals_q[available_total_clone_sizes[ii-1]])) + (weight2 * np.interp(sample_size, error_envelope_x_vals_q[available_total_clone_sizes[ii]], error_envelope_y_vals_q[available_total_clone_sizes[ii]]))


def abs_error(total_clones, sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q,high_low_sign):
    """
    high_low_sign = +1 will return the upper limit on the number of reconstructed clones for 0D = total clones and given sample size.
    
    Dependence on q (as which D numbers, qD) is embodied in error_envelope_x_vals and error_envelope_y_vals.
        
    These must be calculated for each q, and the appropriate dictionary passed to the abs_error function.
    """
    return total_clones * exp(high_low_sign * error_bar(sample_size, total_clones, error_envelope_x_vals_q, error_envelope_y_vals_q) ) # NB total_observed_clones is the sample size


def error_bar_on_fit_qD(sample_size, reconstructed_qD, error_envelope_x_vals_q, error_envelope_y_vals_q):
    """
    
    """
    #
    step_size = (reconstructed_qD / 1000.0) # previously, this was taken as an integer
    #
    n = 1
    #
    try:
        while (abs_error((reconstructed_qD - n*step_size), sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q, +1) >= reconstructed_qD):
            n += 1
        #
        lower_limit = reconstructed_qD - (n+1)*step_size
    except TypeError:
        lower_limit = 0.0
    #
    n = 1
    #
    while (abs_error((reconstructed_qD + n*step_size), sample_size, error_envelope_x_vals_q, error_envelope_y_vals_q, -1) <= reconstructed_qD):
        n += 1
    #
    upper_limit = reconstructed_qD + (n+1)*step_size
    #
    return lower_limit, upper_limit


def diversity_q_fractions_zeros(fracs,q):
    """
        Returns the q-diversity of a list of fractions
        """
    diversity = 0.0
    #
    if (q != 1.0) and (q != float('inf')):
        for frac in fracs:
            if frac != 0.0:
                diversity += frac**q
        diversity = diversity**(1.0/(1.0-q))
    #
    elif q == 1.0:
        for frac in fracs:
            if frac != 0.0:
                diversity -= frac*log(frac)
        diversity = exp(diversity)
    elif q == float('inf'):
        diversity = 1.0 / max(fracs)
    #
    return diversity


def limits_entropy(weights,means,observed_clone_size_distribution,sample_size, observed_threshold):
    #
    max_clone_size = max(observed_clone_size_distribution.keys())
    #
    fitted_distribution = mixed_distribution(zip(weights,means))
    #
    obs_small_clones = sum([observed_clone_size_distribution[x] for x in range(1,observed_threshold) if x in observed_clone_size_distribution])
    distribution_norm = obs_small_clones / sum([fitted_distribution.value_at_clone_size(x) for x in range(1,observed_threshold)])
    #
    tec_obs = distribution_norm
    #
    D1_below_threshold = D_number_from_parameters(weights,means,tec_obs,1)
    #
    for x in range(1,observed_threshold):
        if x not in observed_clone_size_distribution:
            observed_clone_size_distribution[x] = 0
    #
    fit_large_clones = {}
    #
    unfitted_clones = {}
    #
    # NB in the sum below I should only take keys that are already in observed_clone_size_distribution.
    # Otherwise, as a defaultdict it will create zero entries.
    #
    for clone_size in range(observed_threshold, max(observed_clone_size_distribution.keys())+1):
        if clone_size in observed_clone_size_distribution:
            fit_large_clones[clone_size] = distribution_norm * fitted_distribution.value_at_clone_size(clone_size)
            diff = observed_clone_size_distribution[clone_size] - fit_large_clones[clone_size]
            if diff > 0.0:
                unfitted_clones[clone_size] = diff
    #
    #if verbose: print 'Number of unfitted clones', len(unfitted_clones.keys())
    #
    aug_total_estimated_clones = tec_obs + sum(unfitted_clones.values())
    #
    aug_total_estimated_clones = float(aug_total_estimated_clones)
    #
    scale_weights = (tec_obs / aug_total_estimated_clones)
    #
    scaled_weights = [weight*scale_weights for weight in weights]
    #
    aug_weights = scaled_weights + [(unfitted_clones[clone_size] / aug_total_estimated_clones) for clone_size in unfitted_clones.keys()]
    #
    #if verbose: print 'aug weights should sum to one:',sum(aug_weights)
    #
    aug_means = means + [clone_size for clone_size in unfitted_clones.keys()]
    #
    frac_fitted = sum(w*m for w,m in zip(aug_weights[:len(weights)], aug_means[:len(means)])) / sum(w*m for w,m in zip(aug_weights, aug_means))
    #
    lower_limit_D1, upper_limit_D1 = error_bar_on_fit_qD(sample_size, D1_below_threshold, error_envelope_x_vals_1, error_envelope_y_vals_1) # Global variables error_envelope_x_vals_1, error_envelope_y_vals_1 !!!
    #
    lower_limit_entropy = log(lower_limit_D1)
    upper_limit_entropy = log(upper_limit_D1)
    #
    lower_limit_entropy = frac_fitted * (lower_limit_entropy - log(D1_below_threshold))
    upper_limit_entropy = frac_fitted * (upper_limit_entropy - log(D1_below_threshold))
    #
    return lower_limit_entropy, upper_limit_entropy


def output_table_of_D_numbers(precomputed_error_bar_file, D_number_output_filename, list_of_filenames):
    #
    with open(precomputed_error_bar_file,'rU') as f:
        for line in f:
            globals().update({line.split('=')[0].strip() : ast.literal_eval(line.split('=')[1].strip())})
        expected_variables = []
        for expected_variable in expected_variables:
            if expected_variable not in locals():
                print expected_variable, 'not in precomputed error bar file!'
                exit()
    #
    estimated_n0 = {}
    fitted_params = {}
    observed_clone_size_distributions = {}
    obs0D = {}
    obs1D = {}
    obs2D = {}
    obsinfD = {}
    #
    observed_threshold = None
    #
    for data_filename in list_of_filenames:
        #
        # NB Need to read in observed_threshold from MLE file
        #
        if verbose: print data_filename
        #
        # Try openning MLE fit
        #
        try:
            #
            MLE_file = open(data_filename,'rU')
            #
        except IOError:
            print data_filename,' file not found'
        else:
            in_best_fit = False
            #
            for line in MLE_file:
                #
                if line[:32] == 'observed_clone_size_distribution': # This requires full clone distribution in fit file, including clones that are not explicitly fitted.
                    globals().update({line.split('=')[0].strip() : ast.literal_eval(line.split('=')[1].strip())})
                    observed_clone_size_distributions[data_filename] = deepcopy(observed_clone_size_distribution)
                #
                if line[:9] == 'threshold':
                    threshold_in_file = int(line.strip().split('=')[1])
                    if not observed_threshold:
                        observed_threshold = threshold_in_file
                    elif observed_threshold != threshold_in_file:
                        print 'ERROR: All MLE fit files in one table must use the same fitting threshold!'
                        exit()
                #
                if line[0] == '=':
                    in_best_fit = not in_best_fit
                elif in_best_fit:
                    if line.split(' = ')[0] == 'fitted parameters':
                        fitted_params_try = eval(line.split(' = ')[1].strip())
                    if line.split(' = ')[0] == 'estimated n0':
                        try:
                            estimated_n0_try = int(float(line.split(' = ')[1].strip()))
                            if 1e-6 in fitted_params_try:
                                break
                            else:
                                #initial_parameters = initial_parameters_try
                                fitted_params[data_filename] = fitted_params_try
                                estimated_n0[data_filename] = estimated_n0_try
                        except ValueError:
                            #if verbose: print filename,' Nan in fit'
                            estimated_n0[data_filename] = 'NaN'
                #
            #
            # Before normalisation, fracs is a list with one entry per species, where the entry is the number of all cells in that species
            # After normalisation, fracs is a list with one entry per species, where the entry is the fraction of all cells in that species
            #
            fracs = [] 
            for key in observed_clone_size_distribution: # keys are clone sizes
                fracs += [key]*observed_clone_size_distribution[key] # NB addition is actually appending multiple copies to a list
            #
            norm_fracs = float(sum(fracs))
            fracs = [frac / norm_fracs for frac in fracs]
            #
            obs0D[data_filename] = diversity_q_fractions_zeros(fracs,0.)
            obs1D[data_filename] = diversity_q_fractions_zeros(fracs,1.)
            obs2D[data_filename] = diversity_q_fractions_zeros(fracs,2.)
            obsinfD[data_filename] = 1.0/max(fracs)
            #
    with open(D_number_output_filename,'a') as D_number_file_out:
        #
        # Write header for file_out
        #
        # This header should include all the options, including filenames, necessary
        # to reconstruct the result.
        #
        D_number_file_out.write('# data_filenames = '+str(list_of_filenames)+'\n')
        D_number_file_out.write('# observed_threshold = '+str(observed_threshold)+'\n')
        D_number_file_out.write('# precomputed_error_bar_file = '+precomputed_error_bar_file+'\n')
        #
        D_number_file_out.write('Sample name\tobs0D\tobs entropy\tobs2D\tobsinfD\tMLE0D\tMLE entropy\tMLE2D\tMLEinfD\tLower 0D\tUpper 0D\tLower MLE entropy\tUpper MLE entropy\n')
        #
        for filename in sorted( list( set(observed_clone_size_distributions.keys()) & set(estimated_n0.keys()) ) ):
            #
            if verbose: print filename
            #
            weights = fitted_params[filename][:len(fitted_params[filename])/2]
            means = fitted_params[filename][len(fitted_params[filename])/2:]
            #
            D_number_file_out.write(filename+'\t') # sample name
            #
            D_number_file_out.write(str(obs0D[filename])+'\t'+str(log(obs1D[filename])/log(2.0))+'\t'+str(obs2D[filename])+'\t'+str(obsinfD[filename])+'\t') # 
            #
            reconstructed_0D = D_number_from_parameters_with_observed4(weights, means, observed_clone_size_distributions[filename], observed_threshold, 0)
            #
            D_number_file_out.write(str(reconstructed_0D) + '\t')
            #
            reconstructed_1D = D_number_from_parameters_with_observed4(weights, means, observed_clone_size_distributions[filename], observed_threshold, 1)
            #
            D_number_file_out.write(str(log(reconstructed_1D)/log(2.0)) + '\t')
            #
            for q in [2, float('inf')]:
                D_number_file_out.write(str(D_number_from_parameters_with_observed4(weights, means, observed_clone_size_distributions[filename], observed_threshold, q)) + '\t')
            #
            sample_size = sum(key*observed_clone_size_distributions[filename][key] for key in observed_clone_size_distributions[filename] if key < observed_threshold)
            #
            try:
                #
                lower_limit, upper_limit = error_bar_on_fit_qD(sample_size, reconstructed_0D, error_envelope_x_vals_0, error_envelope_y_vals_0)
                #
                # For 1D limits use:
                #
                # lower_limit, upper_limit = error_bar_on_fit_qD(sample_size, reconstructed_1D, error_envelope_x_vals_1, error_envelope_y_vals_1)
                #
            except TypeError:
                lower_limit = 'Too few clones'
                upper_limit = 'Too few clones'
            #
            D_number_file_out.write(str(lower_limit)+'\t')
            D_number_file_out.write(str(upper_limit)+'\t')
            #
            try:
                #
                lower_limit_entropy, upper_limit_entropy = limits_entropy(weights,means,observed_clone_size_distributions[filename],sample_size, observed_threshold)
                #
                lower_limit = (log(reconstructed_1D) + lower_limit_entropy) / log(2.0)
                upper_limit = (log(reconstructed_1D) + upper_limit_entropy) / log(2.0)
                #
                D_number_file_out.write(str(lower_limit)+'\t')
                D_number_file_out.write(str(upper_limit)+'\t')
                #
            except TypeError:
                lower_limit = 'Too few clones'
                upper_limit = 'Too few clones'
                #
                D_number_file_out.write('Too few clones\tToo few clones\t')
            #
            D_number_file_out.write('\n')


def parse_dict(d_in):
    #
    # This will parse a string which represents
    # a python dictionary with integers or floats 
    # as keys and lists of integers or floats as values.
    #
    d_in = d_in.strip()
    #
    assert d_in[0] == '{'
    assert d_in[-1:] == '}'
    #
    # Remove curly braces
    #
    d_in = d_in[1:-1]
    #
    d_in = d_in.strip()
    assert d_in[-1] == ']'
    # 
    # Remove final closing square bracket
    #
    d_in = d_in[:-1]
    #
    # Convert to a list, where each element contains
    # a key and a value
    #
    d_in = d_in.split('], ')
    #
    d_out = {}
    #
    for element in d_in:
        try:
            d_out[int(element.split(': [')[0])] = [int(val) for val in element.split(': [')[1].split(',')]
        except ValueError:
            try:
                d_out[int(element.split(': [')[0])] = [float(val) for val in element.split(': [')[1].split(',')]
            except ValueError:
                d_out[float(element.split(': [')[0])] = [float(val) for val in element.split(': [')[1].split(',')]
    #
    return d_out


def read_error_bar_variable(error_bar_file, q):
    """
    error_bar_file should be a text file, as generated by compare_MLE_alpha_v7.py
    containing values for error_envelope_x_vals_q and error_envelope_x_vals_q
    """
    with open(error_bar_file,'rU') as f:
        for line in f:
            #
            if line.split(' = ')[0][:14] == 'error_envelope':
                #
                if line.split(' = ')[0][22:].strip() == str(q):
                    #
                    error_envelope_vals = parse_dict(line.split(' = ')[1].strip())
                    #
                    if line.split(' = ')[0][15] == 'x':
                        #
                        error_envelope_x_vals = error_envelope_vals
                        #
                    elif line.split(' = ')[0][15] == 'y':
                        #
                        error_envelope_y_vals = error_envelope_vals
                        #
    #
    return error_envelope_x_vals, error_envelope_y_vals


def eb(sample_size, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target):
    """
    This returns the log of high error with the given parameters.
    
    The actual error will be exp(eb_val) - true_total_clone.
    
    exp is because error_envelope_y_vals made by compare_MLE_alpha_v7.py contains log(fitted/true) values.
    
    """
    #
    available_total_clone_sizes = sorted(error_envelope_x_vals.keys())
    #
    for ii in range(len(available_total_clone_sizes)):
        if available_total_clone_sizes[ii] >= true_total_clone:
            break
    #
    weight1 = (available_total_clone_sizes[ii] - true_total_clone) / float(available_total_clone_sizes[ii] - available_total_clone_sizes[ii-1])
    weight2 = (true_total_clone - available_total_clone_sizes[ii-1]) / float(available_total_clone_sizes[ii] - available_total_clone_sizes[ii-1])
    #
    eb_val = (weight1 * np.interp(sample_size, error_envelope_x_vals[available_total_clone_sizes[ii-1]], error_envelope_y_vals[available_total_clone_sizes[ii-1]])) + (weight2 * np.interp(sample_size, error_envelope_x_vals[available_total_clone_sizes[ii]], error_envelope_y_vals[available_total_clone_sizes[ii]]))
    #
    return number_of_error_bars*(np.exp(eb_val) - 1.0) - target # bisect_eb uses bisection to find the root number_of_error_bars*(np.exp(eb_val) - 1.0) - target = 0.0


def bisect_eb(eb, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target):
    """
    This finds a root of number_of_error_bars*(np.exp(eb) - 1.0) - target = 0.0
    by iterated bisection.
    """
    #
    low = 0
    high = true_total_clone*10
    #
    if eb(low,true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target) < 0.0:
        available_total_clone_sizes = sorted(error_envelope_x_vals.keys())
        for ii in range(len(available_total_clone_sizes)):
            if available_total_clone_sizes[ii] >= true_total_clone:
                break
        return sorted(error_envelope_x_vals[available_total_clone_sizes[max(0,ii-1)]])[1]
    #
    assert not eb(low,true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target) * eb(high, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target) > 0.0
    #
    while high - low > 1:
        midpoint = (low + high) / 2.0
        if eb(low, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target) * eb(midpoint, true_total_clone, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target) > 0.0:
            low = midpoint
        else:
            high = midpoint
    return midpoint


def get_sample_size(number_of_clones, fold_difference, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars):
    """
    This returns the number of cells in a sample which
    produce the an error bar of max_error_bar for a given number_of_clones
    in the parent population.
    
    This is the inverse function of the calculation performed by error_bar_on_fit_qD 
    
    Delta = (fold_difference - 1.0) * number_of_clones.
    """
    #
    target = fold_difference - 1.0
    #
    sample_size = bisect_eb(eb, number_of_clones, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, target)
    #
    return sample_size


def get_sampling_limit(number_of_clones, min_number_of_doublets, fraction_small_clones):
    """
    This function returns the minimum number of cells required
    to reliably see doublets and triplets.
    
    The number of possible pairs of cells is sample_size^2/2
    
    Expected doublets = p * sample_size^2/2
    
    However, we need an effective sample size, given by the sample_size * fraction_of_small_clones

    where fraction_of_small_clones is the fraction of cells in small clones.    
        
    Expected doublets = p * (sample_size * fraction_of_small_clones)^2/2
        
    Where p is the probability that a pair of cells is are from the same clone, given that both cells are from small clones.
    
    Now estimate p = 1/number_of_small_clones, which is approximately 1/number_of_clones.
    
    Expected doublets = (sample_size * fraction_of_small_clones)^2/ (2 * number_of_clones)

    We want to expect at least min_number_of_doublets.
        
    rearranging:
        
    sample_size = sqrt(2 * number_of_clones * min_number_of_doublets) / fraction_of_small_clones
    
    """
    #
    sampling_limit = int(round(np.sqrt(2.0 * min_number_of_doublets * number_of_clones) / fraction_small_clones))
    #
    return sampling_limit


def get_minimun_cells(list_of_number_of_clones, list_of_fold_differences, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, min_number_of_doublets, fraction_small_clones):
    #
    minimum_cells = {}
    #
    for fold_difference in list_of_fold_differences:
        for number_of_clones in list_of_number_of_clones:
            #
            # Calculate the minimum number of cells required to
            # reliably see a fold.
            #
            fold_limit = get_sample_size(number_of_clones, fold_difference, error_envelope_x_vals, error_envelope_y_vals,  number_of_error_bars)
            #
            # Calculate the minimum number of cells required to
            # reliably see doublets and triplets.
            #
            sampling_limit = get_sampling_limit(number_of_clones, min_number_of_doublets, fraction_small_clones)
            #
            # Take the more stringent result.
            #
            minimum_cells[number_of_clones, fold_difference] = int(round(max(fold_limit, sampling_limit)))
            #
    return minimum_cells


def output_table(list_of_number_of_clones, list_of_fold_differences, error_bar_file, number_of_error_bars, min_number_of_doublets, fraction_small_clones, minimum_cells, q, file_output = None):
    #
    # Store lines to output in a list
    #
    lines = []
    #
    # Record variables used
    #
    lines.append('# '+datetime.now().strftime('%c')+'\n')
    lines.append('# python '+' '.join(argv)+'\n')
    lines.append('error_bar_file = '+str(error_bar_file)+'\n')
    lines.append('number_of_error_bars = '+str(number_of_error_bars)+'\n')
    lines.append('min_number_of_doublets = '+str(min_number_of_doublets)+'\n')
    lines.append('fraction_small_clones = '+str(fraction_small_clones)+'\n')
    lines.append('q = '+q+'\n\n')
    #
    # make header line
    #
    line = '\t'
    for number_of_clones in list_of_number_of_clones:
        line += str(int(number_of_clones))+'\t'
    line = line[:-1] + '\n'
    #
    lines.append(line)
    #
    # Iterate through lines:
    #
    for fold_difference in list_of_fold_differences:
        line = str(fold_difference)+'\t'
        #
        # Iterate through cols within line:
        #
        for number_of_clones in list_of_number_of_clones:
            line += str(minimum_cells[number_of_clones, fold_difference])+'\t'
        line = line[:-1] + '\n'
        lines.append(line)
    #
    # Output lines:
    #
    for line in lines:
        if file_output:
            with open(file_output,'a') as f:
                f.write(line)
        else:
            if verbose: print line,
    #
    return


# =========== 1.3 END DEFINE FUNCTIONS FOR MAKING TABLE OF D NUMBERS ==================

# =========== 1.4 END DEFINE FUNCTIONS FOR RESAMPLING ==================

def output_resample(resample_fit_filename, file_out):
    #
    # NB read_sample_fit is defined in section 1.2 above
    #
    threshold, model, start_weight, start_mean = read_sample_fit(resample_fit_filename, return_threshold = True) 
    #
    sample_distribution = mixed_distribution(zip(model.fitted_weights,model.fitted_means))
    total_clones = model.total_clones
    sample_size = model.sample_size
    clone_size = 1
    total_observed_clones = 0
    unobserved_clones = model.estimated_n0
    #
    clone_limit = total_clones - unobserved_clones # When we observe this many clones we can stop looking.
    #
    observed_clone_size_distribution = {}
    #
    while ( total_observed_clones < clone_limit ) and (clone_size < threshold):
        new_clones = total_clones * sample_distribution.value_at_clone_size(clone_size)
        #
        observed_clone_size_distribution[clone_size] = int(round(new_clones))
        total_observed_clones += new_clones
        #
        clone_size += 1
    #
    # Output to screen
    #    
    print 'clone\t\tnumber of'
    print 'size\t\tclones'
    print '-----\t\t---------' 
    #
    for ii in range(threshold):
        if ii+1 in observed_clone_size_distribution:
            print ii+1,'\t\t',observed_clone_size_distribution[ii+1]
        else:
            print ii+1,'\t\t',0
    #
    # Output to file
    #
    with open(file_out,'w') as f:
        f.write('clone\t\tnumber of\n')
        f.write('size\t\tclones\n')
        f.write('-----\t\t---------\n')
        for ii in range(threshold):
            if ii+1 in observed_clone_size_distribution:
                f.write(str(ii+1)+'\t\t'+str(observed_clone_size_distribution[ii+1])+'\n')
            else:
                f.write(str(ii+1)+'\t\t0\n')
    #
    return    

# =========== 1.4 END DEFINE FUNCTIONS FOR RESAMPLING ==================

# =========== 1. END DEFINE FUNCTIONS ===============

#
# ========== 2. BEGIN PARSE COMMAND LINE ==============
#

parser = argparse.ArgumentParser(description='')
parser.add_argument('file_out', type=str, help='File name to output fit')
parser.add_argument('file_in', nargs=argparse.REMAINDER, help='List of file names containing clone data. A long list of files can be conveniently input using a shell variable, e.g. MLE_files=$(ls);python recon.py file_out $MLE_files')
parser.add_argument('-t','--threshold', default=50, type=int, help='Largest clone size to fit with parameters')
parser.add_argument('-l','--parameter_limit', default=20, type=int, help='Maximum number of parameters to fit')
parser.add_argument('-c','--clone_distribution_in_file',  action='store_true', help='Input file contains one clone names followed by a tab separated size per line.')
parser.add_argument('-d','--bin_size', type=int, help='Average number of reads per cell.')
parser.add_argument('-e','--make_error_bar_files',  action='store_true', help='Make a precomputed error bar file that can be used to calculate error bars.')
parser.add_argument('-v','--verbose',  action='store_true', help='Print verbose output.')
parser.add_argument('-T','--make_table_of_D_numbers',  action='store_true', help='Output a table of reconstructed D numbers with error bars from input fit files.')
parser.add_argument('-p','--make_table_power',  action='store_true', help='Output a table of cells required to see fold difference.')
parser.add_argument('-b','--precomputed_error_bar_file', type=str, help='File name for precomputed error bars.')
parser.add_argument('-q','--q', type=str, help='Hill number parameter to use when making table for power calculations ')
parser.add_argument('-n','--number_of_error_bars', default=2.0, type=float, help='Mapping between gold-standard error bar fits and standard deviations in power calculation')
parser.add_argument('-m','--min_number_of_doublets', default=100, type=int, help='')
parser.add_argument('-f','--fraction_small_clones', default=0.1, type=float, help='')
parser.add_argument('-r','--resample', action='store_true', help='Output a resampling of the model fit with the given fit file.')

args = parser.parse_args()
globals().update(vars(args))

if make_table_of_D_numbers:
    list_of_filenames = file_in
    D_number_output_filename = file_out
    #
elif len(file_in) == 1:
    file_in = file_in[0]
    if make_error_bar_files:
        if not isdir(file_in):
            print 'To make error bar files, the input should be a directory containing the error bar fits.'
            exit()
        else:
            error_bar_fits_directory = file_in
            if error_bar_fits_directory[-1] != '/':
                error_bar_fits_directory = error_bar_fits_directory + '/'
        #
    else:
        if not isfile(file_in):
            print 'To make an MLE fit, the input should be a single file containing data. To make a table of D numbers input can be multiple filenames, but not a directory name.'
            exit()
        #
else:
    print 'Multiple files are only used for making tables of D numbers.'
    exit()

#
# The two variables below are for power calculations. They are not yet implemented in the parser, so are hardcoded below:
#

#
# list_of_number_of_clones will be column headers. Edit next line to change:
#

list_of_number_of_clones = [1e4, 3e4, 1e5, 1e6, 3e6]

#
# list_of_fold_differences will be row titles. Edit next line to change:
#

list_of_fold_differences = [1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0]

#
# ========== 2. END PARSE COMMAND LINE ==============
#

# ========== 3. BEGIN EXECUTABLES =================

if resample:
    #
    output_resample(file_in, file_out)
    #
elif make_error_bar_files:
    #
    precomputed_error_bar_file = file_out
    #
    make_error_bar_file(error_bar_fits_directory, precomputed_error_bar_file)
    #
elif make_table_of_D_numbers:
    #
    output_table_of_D_numbers(precomputed_error_bar_file, D_number_output_filename, list_of_filenames)
    #
elif make_table_power:
    #
    # Read in error bar variables
    #
    error_bar_file = file_in
    error_envelope_x_vals, error_envelope_y_vals = read_error_bar_variable(error_bar_file, q)
    #
    # Calculate minimum number of cells:
    #
    minimum_cells = get_minimun_cells(list_of_number_of_clones, list_of_fold_differences, error_envelope_x_vals, error_envelope_y_vals, number_of_error_bars, min_number_of_doublets, fraction_small_clones)
    #
    # Output results as table:
    #
    output_table(list_of_number_of_clones, list_of_fold_differences, error_bar_file, number_of_error_bars, min_number_of_doublets, fraction_small_clones, minimum_cells, q, file_out)
    #
else: # This option is for making MLE fits of data
    #
    # Initialise variables to default values.
    #
    cutoff_score_1 = 0.1 # This could be expanded into a pair of cutoffs for weights and means.
    cutoff_clone_distribution = 100000 #0.1 # Actual cutoff will be cutoff_clone_distribution * n0
    #
    initial_scale_factors =[0.001, 0.02, 0.07, 0.3, 0.6, 0.9, 1.1, 1.2, 1.3, 1.4]
    test_weights_list =[0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.95]
    #initial_scale_factors = [0.05, 0.1, 0.5, 0.9, 1.1, 1.2, 1.3]
    #test_weights_list = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    #
    if isfile(file_out):
        print 'Error: Output file already exists'
        exit()
    #
    # Read in observed data
    #
    if not bin_size:
        bin_size = 1.0
    observed_clone_size_distribution = read_observed_clone_size_distribution(file_in, clone_distribution_in_file, threshold = float('inf'), bin_size = bin_size)
    #
    if not [key for key in observed_clone_size_distribution.keys() if key < threshold]: # Empty list, i.e. no clones smaller than threshold.
        #
        # Print message, do not perform fit:
        #
        print 'All clones oversampled. No missing species to reconstruct.'
        #
    elif set([key for key in observed_clone_size_distribution.keys() if key < threshold]) <= set([1, 2]): # Only singletons and doublets give unreliable fit.
        #
        # Print message, do not perform fit:
        #
        print 'Not enough observations to perform fit. Reconstruction with only singlets and doublets is unreliable.'
        #
    else: # Go ahead with fit:
        #
        # Write command line input and fitted data to output file:
        #
        with open(file_out,'a') as f:
            f.write('# '+datetime.now().strftime('%c')+'\n')
            f.write('# python '+' '.join(argv)+'\n')
            #
            # All clones are output to fit file, including those not fit explicitly:
            # 
            f.write('observed_clone_size_distribution = '+str(dict(observed_clone_size_distribution))+'\n')
            f.write('threshold = '+str(threshold)+'\n')
        #
        # Enforce threshold on clones to be explicitly fitted:
        #
        for clone_size in observed_clone_size_distribution.keys():
            if (clone_size > threshold) or (observed_clone_size_distribution[clone_size] == 0):
                del observed_clone_size_distribution[clone_size]
        #
        # Carry out fit:
        #
        carry_out_fit(observed_clone_size_distribution, file_out)

# ========== 3. END EXECUTABLES =================

