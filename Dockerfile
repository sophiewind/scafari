# Use the Bioconductor Docker image as the base
FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Update package lists and install Git (if not already installed)
RUN apt-get update && apt-get install -y git

# Set the working directory
WORKDIR /home/rstudio

# Copy the setup.R and other necessary files into the Docker image
COPY *.R ./ 

# Run any setup scripts or install necessary packages
# Assumes that setup.R is part of your repository
RUN Rscript setup.R

# Expose the default port for Shiny apps
EXPOSE 3838

# Command to run Shiny Server
CMD ["/init"]
