[![Build Status](https://travis-ci.com/milescsmith/enhancedDimPlot.svg?branch=master)](https://travis-ci.com/milescsmith/enhancedDimPlot)

# enhancedDimPlot

An alternative to Seurat's DimPlot that provides additional facetting, group highlighting, and labeling options.
Currently, only Seurat 3 objects are supported, by SingleCellExperiment objects will be forthcoming.

With enhancedDimPlot you can easily
* add labels having a background that matches the groups:
<img src="enhancedDimPlot_example.png" alt="labeling" width="300"/>

* faceting by metadata variables:
<img src="facet_example.png" alt="faceting" width="300"/>

* plot only particular subgroups:
<img src="subgroup_example.png" alt="highlight" width="300"/>

* highlight one particular group:
<img src="highlightgroup_example.png" alt="highlight" width="300"/>

* or plot just the labels:
<img src="empty_example.png" alt="empty" width="300"/>
