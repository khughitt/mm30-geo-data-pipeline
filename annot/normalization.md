Normalization Notes
===================

_KH (Dec 20, 2021)_

In order for samples _within_ each dataset to be compared with each other,
some form of size factor normalization was performed for each dataset for which such a normalization had not already been applied.

For most datasets, this was achieved by scaling the columns of the dataset such that that all sum to the same number.

For convenience and readability, the convention of scaling the totals to "1,000,000", commonly used for RNA-Seq normalization, was used.

In cases where the data available through GEO already had some form of size-factor normalization performed (CPM, TPM, etc.)., this step was skipped.
