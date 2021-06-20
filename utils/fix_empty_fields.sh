#!/bin/bash

# Fix empty fields in tab-separated file

INPUT=$1

perl -F'/\t/' -i.bak -ape '$_=join("\t", map { $_ eq "" ? "." : $_ } @F)' $INPUT

