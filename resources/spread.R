library(tidyr)

data('iris')

iris_spread_sepal <- iris %>% pivot_wider(names_from = Species, values_from = Petal.Length)
              
# spread does not work because some columns are exactly the same