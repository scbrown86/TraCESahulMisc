To install in R you must first be a collaborator on the project as it is private. When the package is made public it will be much easier!

You can also use GitHub desktop or CLI, create a local copy on your computer and then build the package in RStudio.

```r
library(credentials)
library(remotes)
library(usethis)


# set config with GitHub account details
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")

# Go to github page to generate token
# Leave default settings, scroll to bottom and hit generate token
usethis::create_github_token() 

# paste your PAT into pop-up that follows...
credentials::set_github_pat()

# now remotes::install_github() will work
remotes::install_github("scbrown86/TraCESahulMisc")
```
