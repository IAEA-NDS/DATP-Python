import os

# global constants
NQND = 'EB'
NHEL = 'EL'
NQST = 'ST'
NHMO = 'MO'
NHFI = 'FI'

# global vars
MAXF = 900
NOM = 200
NQM = 33  # number of apriori tables in DAT.INP (excluding FIS*)

# LIMIT OF DIF/UNC FOR OUTLIERS?
ULI = 3.0

# Enables output with more precision to facilitate testing
try:
    SHOULD_TEST_OUTPUT = (os.environ['TEST_DATP'] == 'yes')
except KeyError:
    SHOULD_TEST_OUTPUT = False
