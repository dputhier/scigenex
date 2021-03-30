# dbfmcl

## Installation

### In the terminal

     R CMD INSTALL dbfmcl
     
     # In R
     
     library(dbfmcl)
     
### From R

To install from a private repo, generate a personal access token (PAT) in https://github.com/settings/tokens and supply to the *auth_token* argument. 

     library(devtools)
     
     install_github("dputhier/dbfmcl", auth_token="...", ref="develop")
