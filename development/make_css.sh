#!/bin/bash

#dependencies
#md5sha1sum <- "brew install md5sha1sum" using homebrew/linuxbrew

#clone pitchfork
git clone https://github.com/PacificBiosciences/pitchfork.git

#change directory to pitchfork
cd pitchfork

#install blasr
make blasr

#install ccs
sudo make pbccs