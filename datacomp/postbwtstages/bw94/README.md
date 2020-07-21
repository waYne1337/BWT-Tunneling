# BW94

### Description
This directory contains a redevelopment of the post stages of "block-sorting compression" from 1994, invented by Mike Burrows and David J. Wheeler.
It uses move-to-front transform, run-length encoding and a source coder. The source coder comes from Sachin Garg.

####Entropy Coding Source code
By Sachin Garg, 2006

Includes range coder based upon the carry-less implementation 
by Dmitry Subbotin, and arithmetic coder based upon Mark Nelson's
DDJ code.
 
Modified to use 64-bit variables for improved performance.
32-bit reference implementations also included.

For details:
http://www.sachingarg.com/compression/entropy_coding/64bit

Please send your suggestions, improvements, errors, feedback etc... 
Read license.txt before using this in anyway.

### Author
Uwe Baier
