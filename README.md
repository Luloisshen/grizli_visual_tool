# grizli visual tool
## Requirement:
- python 
- X11 version of [ds9](https://sites.google.com/cfa.harvard.edu/saoimageds9/download?authuser=0) and [X11](https://www.xquartz.org/). Then add ds9 path to zshrc or bash script. This is the one we found ususally works with pyds9. 
- [pyds9](https://github.com/ericmandel/pyds9)

## To run:
```
python check_results.py -i <directory/inputfile> -n <start number> -l <grilzi id list>
```
The inputfile is the raw catalog containing objects to be classified in a given work package. All related files (i.e., *.png, *.fits) should be placed in the same directory as the inputfile.
When you start, click one of the dropdown menus for the first object before clicking next to prevent the GUI from crashing.  
For example:
- Run the GUI from beginning of a catalog:
  ```
  python check_results.py -i beta_test/beta_test.cat
  ```
- Start from the 5th object in the input catalog:
  ```
  python check_results.py -i beta_test/beta_test.cat -n 5
  ```
- Checking with object id:
  ```
  python check_results.py -i beta_test/beta_test.cat -l '["ngdeep_00001","ngdeep_00002"]'
  ```
  If both "start number" and "grizli id list" are given, "grizli id list" will be used. The input flags and comments will save in every 2 objects to the save input file. This is to prevent loss due to exit for any unknown reason. If the code exit itself, you can re-start 2 objects before the object you were working on to check. 

## Flagging and getting results: 
- Results are stored in the *.addflags.cat. In particular, flag columns are named as 'spectra_quality_flag', 'spectra_fitting_flag', 'phot_fitting_flag', 'spec_phot_fitting_flag', 'spectra_comment', 'refit_flag', and 'refit_comment'. 
- The 'Great', 'Good', 'Unclear' and 'Bad' flags are stored as 0, 1, 2, and 3, respectively.
- 'spectra_quality_flag': This quality is largely decoupled from our ability to recover a redshift from the spectrum and speaks mostly to the data quality. This quality can be determined from the 2d spectra plot in the main window or ds9 window as well as the 1d spectral plot. 
  - "Great": objects with clear continuum and some emission line features, with no overlapping contamination or other reduction artifacts (spectrum cutoff on one end, hot pixels, etc.).  
  - "Good": spectrum contains (1) one or more emission lines with no/faint continuum OR (2) decent S/N continuum (based on a visual inspection) without emission line features. No strong contamination or other reduction artifacts present. 
  - "Unclear": (1) objects with low S/N continuum or no continuum and no emission lines, but contains no strong contamination or other reduction artifacts, OR (2) has obvious contamination in some wavelengths but a clear emission feature or decent S/N featureless continuum in an uncontaminated area (a note put in the comments)
  - "Bad": objects with contamination in the extraction window or severe reduction artifacts, which includes more than half of the spectral coverage being cutoff. Note that contamination is typically shown as very dark green in the 2d spectra plot. 
  - Note that if the continuum is not fitted well (i.e., not subtracted in the continuum subtracted or over-subtracted), but visually there is spectrum, this is spectra fitting issue, not spectra quality issue and the spectral quality should be flagged accordingly. 
  - Even if there is contamination near the extraction window in the 2d spectrum, as long as it is not spatially coincident with the extraction window (or at the very edge where it is severely downweighted), it is fine. Can mark “co” (contamination offset from window) in the comments in this case. If there is contamination within the extraction window, mark “cw” (contamination in the window) in the comment. 

- 'spectra_fitting_flag':
  - "Great": spectrum continuum and emission lines are perfectly fitted and there is at least two clear features fit. Features including emission lines, Balmer break at rest-frame 0.3645 um, and Balmer absorption lines (Ha, Hb, Hg, Hd); 
  - "Good": a good fit (1) two possible emission features are fitted, but either a low S/N continuum is present or the high S/N continuum is not fit well. (2) a single feature fits well and it's pretty clear the redshift is correct, but there appear other features in the grism spectra that may be real or continuum is not fitted or artifacts that are not accounted for or not accounted for well (3) a continuum only fit when no strong emission lines are expected (these are Ha, OIII, Hb, OII in some cases, PaB in some cases) either due to the rest-frame coverage of the grism or a red template fit;
  - "Unclear" (1) only a continuum without features is fit and fit well, but it’s unclear whether emission features should be present, or (2) very noisy spectra fall in this category as well (as long as the model isn’t terribly off). 
  - "Bad": emission lines and/or continuum obviously wrong. 

- 'phot_fitting_flag': flag on the SED photometric fit.
  - "Great": photometric data (black dots) are perfectly fitted by the red line and red dots.
  - "Good": Some photometric data points are significantly off (larger than data errorbar) from the fitted SED model, but the overall shape of the model fits the overall shape of photometric data.
  - "Unclear": cannot see the fitting result due to large errorbars in the spectra OR you think there may be another fit that might work equally well;
  - "Bad": obviously wrong, a majority of the photometric data is off from the fitted SED model;

- 'spec_phot_fitting_flag' will be calculated automatically from your above inputs. Indicates the quality of the overall spectra and SED photometry fitting.
  - "Great": both “phot fitting” and “spec fitting” flags are ‘great’
  - "Good": both “phot fitting” and “spec fitting” flags are ‘good’ or one of the two is ‘good’ and the other is ‘great’
  - "Unclear": both “phot fitting” and “spec fitting” flags are ‘unclear’ or one of the two is ‘unclear’ and the other is ‘good/great’
  - "Bad": “phot fitting” and/or “spec fitting” flags are ‘bad’


## To setup your own catalog:
- Required catalog columns:
  - 'grizli_ID': this is the same as ID of field. For example, 'ngdeep_00001';
  - 'MAG': the magnitude of object. For example, for ngdeep, we use F150W band; 
  - 'z_grizli': Grizli fitted redshift;
  - 'z_phot': photometry redshift from pure photometry data. 
- Required files are .sed.png, .stack.png, .stack.fits, .full.png and .1D.fits

