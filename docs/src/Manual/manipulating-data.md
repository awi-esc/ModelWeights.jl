# Manipulating Data

To conveniently apply any function that takes as first argument a YAXArray, to every (or a subset) of datasets in a `DataMap`, you can use the functions `apply` or `apply!`. While `apply` returns a new DataMap, `apply!` modifies the input DataMap.
Both provide keyword arguments `ids` and `ids_new`. If `ids` is not provided, the function is applied to every dataset in the input DataMap and if `ids_new` is not provided, the new names will be identical to the old names. That is, if they are not specified and `apply!` is used the respective entries are overwritten.

The first argument to `apply` (or `apply!`) is always the DataMap, the second the function that takes as first argument a YAXArray, followed by the positional and then possibly keyword arguments that the function expects.

!!! info "Coming soon"
    Examples will follow.