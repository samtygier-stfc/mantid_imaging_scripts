Script implementing overlap correction.
Based on: A. S. Tremsin, J. V. Vallerga, J. B. McPhate, and O. H. W. Siegmund, ‘Optimization of Timepix count rate capabilities for the applications with a periodic input signal’, J. Inst., 2014, https://dx.doi.org/10.1088/1748-0221/9/05/C05026

To run on a single image stack use
```bash
python overlap_correction/overlap_correction.py --sample /data/fe_calib/fe_sample
```
This will load the data, apply the correction and write the output to a subdirectory named `corrected_mi`. Customise the output directory name by adding `--output NAME`.

The script will also show a plot of the original and corrected histogram. Note the edges of the stack are not included.

The script can also be run with an open beam stack.
```bash
python overlap_correction/overlap_correction.py --sample /data/fe_calib/fe_sample --open /data/fe_calib/fe_flat
```
In this case a normalised plot will be shown.
