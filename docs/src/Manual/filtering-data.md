# Filtering data

TODO 


Note that when using the option to load data preprocessed with ESMValTool from a seperate yaml file, 
in which it's possible to specify filtering options for each datasets, these values will be overwritten,
if the same key is provided in the function argument, too. 
That is, the function argument takes precedence over what had been specified in the yaml file if a key is specified in both.
Otherwise if a key is just given in the yaml file, it will still be applied in the filtering.

