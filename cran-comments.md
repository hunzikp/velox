## Test environments
* local ubuntu 16.04 install, R 3.3.1
* local ubuntu 16.04 install, R-devel (2017-11-30 r73808)
* ubuntu 14.04 server (on travis-ci), R 3.3.1
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is 16.7Mb
  sub-directories of 1Mb or more:
    libs  16.5Mb

  The size is due to the extensive use of C++ template programming 
  in both the boost headers and my own C++ source code.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of velox 
(currently just package 'prioritizr'). 
The install passed.
