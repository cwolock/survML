# Gradient boosting for C-index

Using doubly-robust gradient boosting to generate estimates of the
prediction function that maximizes the C-index

## Usage

``` r
boost_c_index(
  time,
  event,
  X,
  newX,
  S_hat,
  G_hat,
  V,
  approx_times,
  tuning,
  produce_fit = TRUE,
  subsample_n,
  boosting_params
)
```

## Arguments

- time:

  `n x 1` numeric vector of observed follow-up times. If there is
  censoring, these are the minimum of the event and censoring times.

- event:

  `n x 1` numeric vector of status indicators of whether an event was
  observed.

- X:

  `n x p` data.frame of observed covariate values

- newX:

  `m x p` data.frame of new observed covariate values at which to obtain
  `m` predictions for the estimated algorithm. Must have the same names
  and structure as `X`.

- S_hat:

  `n x J2` matrix of conditional event time survival function estimates

- G_hat:

  `n x J2` matrix of conditional censoring time survival function
  estimates

- V:

  Number of cross-validation folds for selection of tuning parameters

- approx_times:

  Numeric vector of length J2 giving times at which to approximate
  C-index integral. Note that the last time in `approx_times` is taken
  to be the restriction time (i.e., the maximum follow-up) for
  comparison of pairs of individuals. Essentially, this time should be
  chosen such that the conditional survival function is identified at
  this time for all covariate values `X` present in the data. Choosing
  the restriction time such that roughly 10% of individuals remain
  at-risk at that time has been shown to work reasonably well in
  simulations.

- tuning:

  Logical, whether or not to use cross-validation to select tuning
  parameters

- produce_fit:

  Logical, whether to produce a fitted prediction function using the
  selected optimal parameters.

- subsample_n:

  Number of samples to use for boosting procedure. Using a subsample of
  the full sample can greatly reduce runtime

- boosting_params:

  Named list of parameter values for the boosting procedure. Elements of
  this list include `mstop` (number of boosting iterations), `nu`
  (learning rate), `sigma` (smoothness parameter for sigmoid
  approximation, with smaller meaning less smoothing), and `learner`
  (base learner, can take values `"glm"`, `"gam"`, or `"tree"`)

## Value

Vector of predictions corresponding to `newX`

## Examples

``` r
# This is a small simulation example
set.seed(123)
n <- 250
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

# Note that this a very small Super Learner library, for computational purposes.
SL.library <- c("SL.mean", "SL.glm")

# Note that we do not use times beyond the 90th percentile of observed follow-up times
approx_times <- c(0, unique(quantile(time, probs = seq(0, 0.9, by = 0.01))))

# estimate conditional survival functions at approx_times
fit <- stackG(time = time,
              event = event,
              X = X,
              newX = X,
              newtimes = approx_times,
              direction = "prospective",
              bin_size = 0.1,
              time_basis = "continuous",
              surv_form = "PI",
              learner = "SuperLearner",
              time_grid_approx = approx_times,
              SL_control = list(SL.library = SL.library,
                                V = 3))

# use boosting to estimate optimal (according to C-index) prediction function
boosted_preds <- boost_c_index(time = time,
                               event = event,
                               X = X,
                               newX = X,
                               approx_times = approx_times,
                               S_hat = fit$S_T_preds,
                               G_hat = fit$S_C_preds,
                               V = 3,
                               tuning = TRUE,
                               produce_fit = TRUE,
                               subsample_n = 200,
                               boosting_params = list(mstop = c(100, 200),
                                                      nu = 0.1,
                                                      sigma = 0.1,
                                                      learner = "glm"))
boosted_preds
#>             1             2             3             4             5 
#> -0.0686025275 -0.0322597503  0.1645717386  0.0008247817  0.0072922855 
#>             6             7             8             9            10 
#>  0.1817756968  0.0437814373 -0.1461282309 -0.0825078385 -0.0559695540 
#>            11            12            13            14            15 
#>  0.1277527111  0.0326571208  0.0371636969  0.0052451712 -0.0680925911 
#>            16            17            18            19            20 
#>  0.1896811642  0.0478453232 -0.2233205796  0.0702370686 -0.0589546118 
#>            21            22            23            24            25 
#> -0.1244261433 -0.0309170984 -0.1198247614 -0.0871333302 -0.0757064766 
#>            26            27            28            29            30 
#> -0.1925204988  0.0852486165  0.0099424053 -0.1321627232  0.1310242531 
#>            31            32            33            34            35 
#>  0.0399906780 -0.0394000504  0.0915575968  0.0896879443  0.0834654710 
#>            36            37            38            39            40 
#>  0.0688379632  0.0540144065 -0.0137454306 -0.0405984093 -0.0487965777 
#>            41            42            43            44            45 
#> -0.0833720298 -0.0298104545 -0.1461651039  0.2317174186  0.1259790461 
#>            46            47            48            49            50 
#> -0.1305091502 -0.0512627767 -0.0582794605  0.0788864576 -0.0161063881 
#>            51            52            53            54            55 
#>  0.0209394187 -0.0100742774 -0.0116503174  0.1436543313 -0.0317749018 
#>            56            57            58            59            60 
#>  0.1599243143 -0.1773428752  0.0573919046  0.0066944373  0.0168268254 
#>            61            62            63            64            65 
#>  0.0348385420 -0.0622040291 -0.0435961517 -0.1190073398 -0.1248626905 
#>            66            67            68            69            70 
#>  0.0264640498  0.0423833464 -0.0011012053  0.0945440156  0.2186379862 
#>            71            72            73            74            75 
#> -0.0609615364 -0.2610116188  0.1037283539 -0.0849667839 -0.0826350075 
#>            76            77            78            79            80 
#>  0.1059105661 -0.0382669069 -0.1412491036  0.0130155890 -0.0222155182 
#>            81            82            83            84            85 
#> -0.0062990372  0.0354592134 -0.0477170746  0.0639676182 -0.0311934555 
#>            86            87            88            89            90 
#>  0.0295727689  0.1137521602  0.0409498411 -0.0427955941  0.1194702770 
#>            91            92            93            94            95 
#>  0.1023821707  0.0534069634  0.0193344323 -0.0760219121  0.1427796092 
#>            96            97            98            99           100 
#> -0.0729799630  0.2337394469  0.1617002043 -0.0328674329 -0.1198705837 
#>           101           102           103           104           105 
#> -0.0850994584  0.0213316979 -0.0340768321 -0.0451734586 -0.1116400670 
#>           106           107           108           109           110 
#> -0.0118876820 -0.0932964790 -0.1904572810 -0.0487696774  0.0941841224 
#>           111           112           113           114           115 
#> -0.0702388216  0.0599611731 -0.1849492534 -0.0130467668  0.0502172143 
#>           116           117           118           119           120 
#>  0.0262026973  0.0046943025 -0.0774302914 -0.1004264233 -0.1196183824 
#>           121           122           123           124           125 
#>  0.0060114085 -0.1111841069 -0.0609094126 -0.0351111506  0.1959472604 
#>           126           127           128           129           130 
#> -0.0786674594  0.0189663633  0.0016447772 -0.1127665637 -0.0147793158 
#>           131           132           133           134           135 
#>  0.1520109692  0.0427458161 -0.0023964045 -0.0534206889 -0.2328525042 
#>           136           137           138           139           140 
#>  0.1175480046 -0.1676478123  0.0744833117  0.2031258037 -0.1658051462 
#>           141           142           143           144           145 
#>  0.0702842093 -0.0357829178 -0.1799166310 -0.1735924788 -0.1831506407 
#>           146           147           148           149           150 
#> -0.0653490304 -0.1677705525  0.0687583584  0.2241421657 -0.1485455113 
#>           151           152           153           154           155 
#>  0.0797418027  0.0776846110  0.0296190493 -0.1178851664 -0.0200766677 
#>           156           157           158           159           160 
#> -0.0377852307  0.0550125874 -0.0479127880  0.1005633196 -0.0481484840 
#>           161           162           163           164           165 
#>  0.1088967967 -0.1223744430 -0.1455884238  0.3496790479 -0.0528002017 
#>           166           167           168           169           170 
#>  0.0258807741  0.0631086261 -0.0601637579  0.0499371699  0.0336639744 
#>           171           172           173           174           175 
#> -0.0306316353  0.0002509347 -0.0106816990  0.2272607475 -0.0885026416 
#>           176           177           178           179           180 
#> -0.1275259763 -0.0027754061  0.0272289917  0.0410975003 -0.0573673077 
#>           181           182           183           184           185 
#> -0.1239312745  0.1320552643 -0.0454053790 -0.1021658378 -0.0329311635 
#>           186           187           188           189           190 
#> -0.0286285778  0.1151914959  0.0023903907  0.0760354280 -0.0618704796 
#>           191           192           193           194           195 
#>  0.0166621917 -0.0426585323  0.0034737741 -0.1054502942 -0.1511610458 
#>           196           197           198           199           200 
#>  0.2128205456  0.0591628490 -0.1446109282 -0.0741799888 -0.1373718999 
#>           201           202           203           204           205 
#>  0.2350023028  0.1374718097 -0.0361072393  0.0528344871 -0.0525231852 
#>           206           207           208           209           210 
#> -0.0593348198 -0.0937034113 -0.0723591373  0.1747164312 -0.0128779980 
#>           211           212           213           214           215 
#>  0.0061873072  0.0198797084  0.1286763136 -0.0637158859 -0.1161390465 
#>           216           217           218           219           220 
#>  0.1774440213 -0.0554745552 -0.0864923756 -0.1429606684 -0.1482908184 
#>           221           222           223           224           225 
#> -0.0700876969  0.0610638403  0.1151835571  0.0709228267 -0.0469465622 
#>           226           227           228           229           230 
#> -0.0003589733 -0.0844601719 -0.0858489400  0.0904050125 -0.1186791412 
#>           231           232           233           234           235 
#>  0.2082081429 -0.0168711561  0.0166724814 -0.0881936336 -0.0701333826 
#>           236           237           238           239           240 
#> -0.1518448395 -0.0270605915  0.0391674522  0.0287500048 -0.0929258992 
#>           241           242           243           244           245 
#> -0.0937055165 -0.0621903045  0.1576786048 -0.1320710332 -0.0266343571 
#>           246           247           248           249           250 
#>  0.2023840077 -0.0180435601 -0.1565568363 -0.0800779954  0.0464819939 
```
