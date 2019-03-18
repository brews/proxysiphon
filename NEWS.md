# proxysiphon v0.0.1a2

* Add `PetmRecord` and `LgmRecord`, specialized subclasses of `NcdcRecord`. These can be created with `read_lgm()` and `read_petm()`.
* Add `.to_netcdf()` method to `NcdcRecord`.
* New method `records.Publication.to_citationstr()` to get quick and dirty bibliography strings.


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