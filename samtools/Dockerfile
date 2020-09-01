# get empty ubuntu
FROM ubuntu

# install dependencies
RUN apt-get update && apt-get install -y samtools

# copy in the repository
#COPY . ~/workarea

# move to workarea
#WORKDIR ~/workarea

# Run the install script
#RUN ./install.sh --noninteractive

# Run the main script

CMD ["samtools"]
