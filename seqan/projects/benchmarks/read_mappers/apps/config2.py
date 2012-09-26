import os

# ----------------------------------------------------------------------
# Set some paths.
# ----------------------------------------------------------------------
# By default, the read mappers live in ./read_mappers.
SetReadMapperBasePath('read_mappers')
# The benchmark data lives in ./data.
SetDataBasePath('data')
# The results go into ./results.
SetResultsBasePath('results')

# ----------------------------------------------------------------------
# Set maximal error rate.
# ----------------------------------------------------------------------
SetMaximalErrorRate(2)

# ----------------------------------------------------------------------
# Select read mappers and benchmark datasets to execute on.
# ----------------------------------------------------------------------
# maq is disabled at the moment since it cannot create SAM output.
selected_mappers = ['bowtie', 'bwa', 'razers99', 'razers100', 'mrsFAST']
selected_mappers = ['bowtie', 'bwa', 'razers99', 'razers100']
#selected_mappers = ['razers99', 'razers100']
# Do not use SOAP since it does not support SAM.
# if HOST_UNAME[0] == 'Linux':
#   selected_mappers.append('soap-fast')
UseReadMappers(selected_mappers)
UseDatasets(['saccharomyces'])#, 'drosophila'])
# Select distance functions to use.
UseDistanceFunctions(['hamming', 'edit'])
UseDistanceFunctions(['hamming'])#, 'hamming-weighted'])


# ======================================================================
# Setup Read Mappers.
# ======================================================================

# ----------------------------------------------------------------------
# RazerS
# ----------------------------------------------------------------------
# Determine the suffix of the RazerS binary.
HOST_NAME = os.uname()
RAZERS_SUFFIX = ''
if HOST_UNAME[0] == 'Darwin':
  RAZERS_SUFFIX = 'OSX'
elif HOST_UNAME[1] == 'Linux':
  if HOST_UNAME[-1] == 'i368':
    RAZERS_SUFFIX = '32'
  else:
    RAZERS_SUFFIX = ''

# Special entry for RazerS as golden standard creator.  Do not enable
# by including its name in UseReadMappers, do not rename.
DefineReadMapper(name='razers-golden',
                 dirname='razers_20090710',
                 version='20090710',
                 genome_format='fasta',
                 reads_format='fasta',
                 binary='razers%s' % RAZERS_SUFFIX)

class RazersParameterGenerator(object):
  def __init__(self, template):
    self.template = template

  def generate(self, args_dict, error_rate, distance_function, min_read_length, avg_read_length, max_read_length):
    result = ''
    # Manually configure the shape if the reads do not all have the same length.
    if error_rate != 0 and min_read_length != max_read_length:
      shape_length = min(int(1 / (0.01 * error_rate)), 14)  # 14 is hard-coded maximum in RazerS
      result += '-s %s ' % ('1' * shape_length)
    # Compute the percent identity from the error rate.
    result += '-i %f ' % (100.0 - error_rate)
    # Allow indels for edit distance.
    if distance_function == 'edit':
      result += '-id '
    # Put values from args_dict into template.
    result += self.template % args_dict
    return result

# RazerS with recognition rate of 99.
DefineReadMapper(
  name='razers99',
  dirname='razers_20090710',
  version='20090710',
  genome_format='fasta',
  reads_format='fastq',
  binary='razers%s' % RAZERS_SUFFIX,
  params=RazersParameterGenerator(
    '-vv -rr 99 -m 1048576 %(genome)s %(reads)s -of 4 -o %(output)s'))

# RazerS with recognition rate of 100.
DefineReadMapper(
  name='razers100',
  dirname='razers_20090710',
  version='20090710',
  genome_format='fasta',
  reads_format='fastq',
  binary='razers%s' % RAZERS_SUFFIX,
  params=RazersParameterGenerator(
    '-vv -rr 100 -m 1048576 %(genome)s %(reads)s -of 4 -o %(output)s'))

# ----------------------------------------------------------------------
# BWA
# ----------------------------------------------------------------------

class BwaParameterGenerator(object):
  def __init__(self, template):
    self.template = template
  def generate(self, args_dict, error_rate, distance_function, min_read_length, avg_read_length, max_read_length):
    result = self.template % args_dict
    result += ' -n %.3f' % (0.01 * error_rate)
    # In case of hamming distance, use "k-mismatch mode".
    if distance_function == 'hamming':
      result += ' -e -1'
    return result

DefineReadMapper(name='bwa',
                 version='0.5.6',
                 genome_format='fasta',
                 reads_format='fastq',
                 params=BwaParameterGenerator('aln %(genome)s %(reads)s'),
                 redir_stdout=True,
                 index_binary='bwa',
                 index_params='index -a bwtsw %(genome)s',
                 index_suffix='.bwt',
                 out_filter_binary='bwa',
                 out_filter_params='samse %(genome)s %(output_tmp)s %(reads)s',
                 out_filter_redir_stdout=True)

# ======================================================================
# Define some benchmark data sets.
# ======================================================================
DefineDataset(name='drosophila',
              reads=['reads_454/SRR027007.1k',
                     'reads_illumina/SRR001980.1k'])
DefineDataset(name='saccharomyces',
              reads=['reads_454/SRR001317.1k',
                     'reads_illumina/SRR003674.1k'])
