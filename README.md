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

### Vignette

Vignette（使い方の説明）を見るには、、以下のコードを走らせてください。

To view the vignettes (package documentation and tutorials), run the
following code:

``` r
library(ctrest)
vignette("ctrest")
```
