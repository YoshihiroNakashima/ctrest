ctrest
================

## パッケージの概要/Package Overview

　“ctrest”は，RESTモデルもしくはREST-RADモデルを簡単に利用するためのRパッケージです．REST/REST-RADモデルとは，自動撮影カメラ（カメラトラップ）によって得られた動画データに基づいて地上性哺乳類・鳥類の密度推定を行うための統計モデルです．RESTモデルの詳細について[Nakashima
et
al. (2018)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.13059)，REST-RADモデルについてはXXXXXX（執筆中）をご参照ください．また，データの準備までのプロセスについては，[Nakashima
et
al. (2021)](https://www.biorxiv.org/content/10.1101/2021.05.18.444583v2)を参照してください．

“ctrest” is an R package designed to facilitate the use of REST models
or REST-RAD models. REST/REST-RAD models are statistical models for
estimating the density of ground-dwelling mammals and birds based on
video data obtained from camera traps. For details on the REST model,
please refer to [Nakashima et
al. (2018)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.13059),
and for the REST-RAD model, please refer to XXXXXX. For the data
preparation process, please refer to [Nakashima et
al. (2021)](https://www.biorxiv.org/content/10.1101/2021.05.18.444583v2).
The package includes the following functions:

### Installation

インストールするためには、以下のコードを走らせてください。

To install this package, run the following code:

``` r
#install.packages("devtools")
devtools::install_github("YoshihiroNakashima/ctrest",build_vignettes = TRUE)
```

    ## Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.

    ## Downloading GitHub repo YoshihiroNakashima/ctrest@HEAD

    ## waldo (0.5.3  -> 0.6.1   ) [CRAN]
    ## Rcpp  (1.0.13 -> 1.0.13-1) [CRAN]

    ## Installing 2 packages: waldo, Rcpp

    ## Installing packages into 'C:/Users/25069/AppData/Local/Temp/Rtmp6NIg68/temp_libpath2a24437d79dc'
    ## (as 'lib' is unspecified)

    ## 
    ##   There is a binary version available but the source version is later:
    ##       binary source needs_compilation
    ## waldo  0.5.3  0.6.1             FALSE
    ## 
    ## package 'Rcpp' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\25069\AppData\Local\Temp\RtmpIN3V7d\downloaded_packages

    ## installing the source package 'waldo'

    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\25069\AppData\Local\Temp\RtmpIN3V7d\remotes74c859837b2e\YoshihiroNakashima-ctrest-0d0bb63/DESCRIPTION' ...  ✔  checking for file 'C:\Users\25069\AppData\Local\Temp\RtmpIN3V7d\remotes74c859837b2e\YoshihiroNakashima-ctrest-0d0bb63/DESCRIPTION'
    ##       ─  preparing 'ctrest':
    ##    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
    ##       ─  installing the package to build vignettes
    ##          creating vignettes ...     creating vignettes ...   ✔  creating vignettes (2m 38.5s)
    ##       ─  checking for LF line-endings in source and make files and shell scripts (400ms)
    ##       ─  checking for empty or unneeded directories
    ##       ─  building 'ctrest_0.0.0.9000.tar.gz'
    ##      
    ## 

    ## Installing package into 'C:/Users/25069/AppData/Local/Temp/Rtmp6NIg68/temp_libpath2a24437d79dc'
    ## (as 'lib' is unspecified)

### Vignette

Vignette（使い方の説明）を見るには、、以下のコードを走らせてください。

To view the vignettes (package documentation and tutorials), run the
following code:

``` r
library(ctrest)
vignette("ctrest")
```

    ## starting httpd help server ... done
