This is a submission of a new version of the package (1.1.0).

## Test environments

* local macOS Sonoma 14.4, R 4.3.2
* Windows Server 2022 x64 (build 20348), R 4.3.3 (GitHub actions)
* macOS 12.7.3, R 4.3.3 (GitHub actions)
* ubuntu 22.04.4, R 4.3.3 (GitHub actions)

## R CMD check results

0 errors | 0 warnings | 2 notes

NOTE: "Found the following (possibly) invalid URLs:
      URL: https://www.tandfonline.com/doi/full/10.1080/10618600.2024.2304070
      From: README.md
      Status: 403
      Message: Forbidden:"
      
      I am not sure what to do about this, since the URL links to a journal article, without
      redirection. 
      
NOTE: "Found the following HTML validation problems:
  stackG.html:4:1 (stackG.Rd:5): Warning: <link> inserting "type" attribute
  stackG.html:12:1 (stackG.Rd:5): Warning: <script> proprietary attribute "onload"
  stackG.html:12:1 (stackG.Rd:5): Warning: <script> inserting "type" attribute
  stackG.html:17:1 (stackG.Rd:5): Warning: <table> lacks "summary" attribute
  stackG.html:52:1 (stackG.Rd:27): Warning: <table> lacks "summary" attribute
  stackG.html:164:1 (stackG.Rd:103): Warning: <table> lacks "summary" attribute
  stackL.html:4:1 (stackL.Rd:5): Warning: <link> inserting "type" attribute
  stackL.html:12:1 (stackL.Rd:5): Warning: <script> proprietary attribute "onload"
  stackL.html:12:1 (stackL.Rd:5): Warning: <script> inserting "type" attribute
  stackL.html:17:1 (stackL.Rd:5): Warning: <table> lacks "summary" attribute
  stackL.html:49:1 (stackL.Rd:24): Warning: <table> lacks "summary" attribute
  stackL.html:136:1 (stackL.Rd:81): Warning: <table> lacks "summary" attribute
  survML-package.html:4:1 (survML-package.Rd:7): Warning: <link> inserting "type" attribute
  survML-package.html:12:1 (survML-package.Rd:7): Warning: <script> proprietary attribute "onload"
  survML-package.html:12:1 (survML-package.Rd:7): Warning: <script> inserting "type" attribute
  survML-package.html:17:1 (survML-package.Rd:7): Warning: <table> lacks "summary" attribute"
  
    I am not sure what to do about this. 

## Downstream dependencies

We checked one reverse dependency (`vaccine`) and found no issues.
