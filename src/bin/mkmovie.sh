#!/bin/bash

echo "encoding... this may take a while"
mencoder "mf://*.bmp" -mf fps=60:type=bmp -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=22000000 
