install.packages("renv")
renv::init()  # Creates renv.lock

renv::restore()  # Installs exact versions from renv.lock