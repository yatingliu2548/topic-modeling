library(reticulate)

# source the python file
source_python("alternative/anchor-word-recovery/learn_topics_function.py")

# now you can call the learn_topics function
result <- learn_topics()

