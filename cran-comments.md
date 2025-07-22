## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


## Resubmission

This is a resubmission. I have made the following changes to address the previous NOTE(s):

* Replaced “Matérn” with “Matern” in the DESCRIPTION to avoid spelling issues.
* Removed abbreviations “BKP” and “DKP” from the DESCRIPTION field to avoid false spelling notes.
* Reduced the grid size in example plots to ensure runtime remains below 10 seconds.

These changes resolve all previous NOTE messages.


## Resubmission

This is a resubmission in response to CRAN feedback.

- Removed unnecessary whitespace in the Description field.
- Added a literature reference to the Description using the recommended format.
- Added \value{} tags to the Rd files for exported methods: plot(), print(), and summary().

Thank you for the helpful feedback.
