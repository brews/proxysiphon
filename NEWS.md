# proxysiphon v0.0.1a5

* `nc2lmrh5()` and `nc2lmrdf()` should now trim based on depth cutoffs of `LgmRecord` objects. These are optional floats at `self.chronology_information.cut_shallow` and `self.chronology_information.cut_deep`.

# proxysiphon v0.0.1a4

* Add `LgmRecord.chronology_depth_range()` to get the depth range of chronology determinants (Issue #11).
* `LgmRecord.slice_datadepths()` now slices `data.age_ensemble` and `data.age_median` in addition to the general dataframe.
* `LgmRecord.to_netcdf()` will now add attributes to the chronology section of the output netCDF file if the LgmRecord object `chronology_information` has `cut_shallow` or `cut_deep` attributes.
* `LgmRecord.plot_agedepth()` now automatically plots `self.chronology_information.cut_shallow` and `self.chronology_information.cut_deep` attributes as vertical lines, if the cutoff attributes are available.

# proxysiphon v0.0.1a3

* Add method `RedateMixin.swapin_custom_deltar()` to handle swapping-in new carbon reservoir values into proxy record chronology information.
* `QcPlotMixin.plot_agedepth()` now returns near empty axis if the proxy record has no data to plot.
* Add 'plot_agedepth_kws' arguement to `to_qcpdf()`.
* `plot_agedepth` now accepts passing an optional iterable of depths to 'depthlines'. It then plots vertical lines at those depths. This is used to indicate cut-off points for a proxy record.


# proxysiphon v0.0.1a2

* `cartopy`, and `snakebacon` have been made recommended dependencies.
* New method `records.Publication.to_citationstr()` to get quick and dirty bibliography strings.
* Add `PetmRecord` and `LgmRecord`, specialized subclasses of `NcdcRecord`. These are created with `read_lgm()` and `read_petm()`.
* `LgmRecord` instances now have several methods for age modeling with 14C.
* `LgmRecord` and `PetmRecord` records now have a `.to_netcdf()` method.
* `proxysiphon.qcreport` submodule has be refactored into `LgmRecord` and `PetmRecord` methods. This should make it easier to do quick plots (e.g. with `.plot_agedepth()`) and output more comprehensive PDF quality-control reports with `.to_qcpdf()`.
* Add `nc2lmrh5`, a function to read output proxy NetCDF files and convert them to LMR-style HDF5 proxy files.


# proxysiphon v0.0.1a1

* Refactoring to use Python 3.7 `dataclasses` to simplify code.
* New `Guts.yank_variables()` and `NcdcRecord.Variables` classes for the "Variables" file section.
* `NcdcRecord.original_source_url` is now an optional string with a default value of None.
* `DataCollection` `first_year` & `last_year` are now floats and not int.
* More flexible record creation with keywords.
* `NcdcRecord` `description` and `original_source_url` attributes refactored to simple string.
* Removed `contribution_date` attribute from `NcdcRecord`. This was unused in our analysis.
* `proxysiphon.agemodel.fit_agedepthmodel` now also returns snakebacon agemodel kwargs.
* Grab new field from proxy text files ("Description" and "Original_Source_URL").
* Now pull more fields from "Publication" section, including support for multiple publicatioins.
* `Guts.pull_section()` now returns list of lists to deal with sections that appear multiple times (i.e. "Publication").
* Fix missing dependency in setup.py.


# proxysiphon v0.0.1a0

* Initial release.