#!/bin/bash
#check args
if [ "$1" = "c" ]; then		#compress infile
	xz -z -c -k -9e -M 100% $2 > $3
elif [ "$1" = "d" ]; then	#decompress infile
	xz -d -c -k -9e -M 100% $2 > $3
elif [ "$1" = "i" ]; then	#install compressor
	sudo apt-get install p7zip
else
	exit 1
fi
