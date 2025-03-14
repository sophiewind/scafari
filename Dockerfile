# Use the Bioconductor Docker image as the base
FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Update package lists and install Git (if not already installed)
RUN apt-get update && apt-get install -y git

# Set the working directory
WORKDIR /srv/shiny-server/

# Clone the Git repository
# Replace <repository-url> with the URL of your repository
# Optionally, you can specify a branch or commit by adding -b <branch-name>
RUN git clone https://github.com/sophiewind/scafari/ .

# Run any setup scripts or install necessary packages
# Assumes that setup.R is part of your repository
RUN Rscript setup.R

# Expose the default port for Shiny apps
EXPOSE 3838

# Command to run Shiny Server
CMD ["/init"]
