
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MoNAn

<!-- badges: start -->
<!-- badges: end -->

MoNAn is the software implementation of the statistical model outlined
in:

Block, P., Stadtfeld, C., & Robins, G. (2022). A statistical model for
the analysis of mobility tables as weighted networks with an application
to faculty hiring networks. Social Networks, 68, 264-278.

as found
[here](https://www.sciencedirect.com/science/article/pii/S0378873321000654).

# Installation

The package is still under development. You can install the development
version of MoNAn from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("perblock/MoNAn")
```

or using:

``` r
# install.packages("remotes")
remotes::install_github("perblock/MoNAn")
```

# Example

This script runs a simple example with the data from the MoNAn package.

``` r
library(MoNAn)

# packages for parallel computing
library(snow)
library(snowfall)
```

The example case is using synthetic data that represents mobility of
individuals between organisations. The artificial data includes an
edgelist, i.e. a list of origins and destinations of individuals.

``` r
mobilityEdgelist
#>        [,1] [,2]
#>   [1,]    1    1
#>   [2,]    1    1
#>   [3,]    1    2
#>   [4,]    1    2
#>   [5,]    1    2
#>   [6,]    1    2
#>   [7,]    1    3
#>   [8,]    1    3
#>   [9,]    1    3
#>  [10,]    1    3
#>  [11,]    1    3
#>  [12,]    1    4
#>  [13,]    1    4
#>  [14,]    1    4
#>  [15,]    1    4
#>  [16,]    1    4
#>  [17,]    1    4
#>  [18,]    1    4
#>  [19,]    1    5
#>  [20,]    1    5
#>  [21,]    1    5
#>  [22,]    1    5
#>  [23,]    1    5
#>  [24,]    1    5
#>  [25,]    1    5
#>  [26,]    1    5
#>  [27,]    1    5
#>  [28,]    1    6
#>  [29,]    1    6
#>  [30,]    1    7
#>  [31,]    1    7
#>  [32,]    1    7
#>  [33,]    1    7
#>  [34,]    1    7
#>  [35,]    1    7
#>  [36,]    1    9
#>  [37,]    1   14
#>  [38,]    1   16
#>  [39,]    1   16
#>  [40,]    1   17
#>  [41,]    2    1
#>  [42,]    2    1
#>  [43,]    2    1
#>  [44,]    2    1
#>  [45,]    2    2
#>  [46,]    2    2
#>  [47,]    2    2
#>  [48,]    2    3
#>  [49,]    2    3
#>  [50,]    2    4
#>  [51,]    2    4
#>  [52,]    2    4
#>  [53,]    2    4
#>  [54,]    2    4
#>  [55,]    2    4
#>  [56,]    2    5
#>  [57,]    2    5
#>  [58,]    2    5
#>  [59,]    2    5
#>  [60,]    2    6
#>  [61,]    2    6
#>  [62,]    2    6
#>  [63,]    2    6
#>  [64,]    2    7
#>  [65,]    2    7
#>  [66,]    2    7
#>  [67,]    2    7
#>  [68,]    2    7
#>  [69,]    2    8
#>  [70,]    2    9
#>  [71,]    2   11
#>  [72,]    2   11
#>  [73,]    2   12
#>  [74,]    2   13
#>  [75,]    2   15
#>  [76,]    2   16
#>  [77,]    3    1
#>  [78,]    3    1
#>  [79,]    3    1
#>  [80,]    3    1
#>  [81,]    3    1
#>  [82,]    3    1
#>  [83,]    3    1
#>  [84,]    3    1
#>  [85,]    3    3
#>  [86,]    3    3
#>  [87,]    3    3
#>  [88,]    3    3
#>  [89,]    3    3
#>  [90,]    3    3
#>  [91,]    3    5
#>  [92,]    3    5
#>  [93,]    3    5
#>  [94,]    3    5
#>  [95,]    3    5
#>  [96,]    3    5
#>  [97,]    3    5
#>  [98,]    3    5
#>  [99,]    3    7
#> [100,]    3    7
#> [101,]    3    8
#> [102,]    3    9
#> [103,]    3    9
#> [104,]    3   10
#> [105,]    3   11
#> [106,]    3   12
#> [107,]    3   13
#> [108,]    3   14
#> [109,]    3   14
#> [110,]    3   14
#> [111,]    3   15
#> [112,]    3   15
#> [113,]    3   16
#> [114,]    3   16
#> [115,]    3   16
#> [116,]    3   16
#> [117,]    3   16
#> [118,]    3   16
#> [119,]    3   16
#> [120,]    4    1
#> [121,]    4    1
#> [122,]    4    1
#> [123,]    4    1
#> [124,]    4    1
#> [125,]    4    1
#> [126,]    4    2
#> [127,]    4    2
#> [128,]    4    2
#> [129,]    4    2
#> [130,]    4    2
#> [131,]    4    2
#> [132,]    4    4
#> [133,]    4    4
#> [134,]    4    4
#> [135,]    4    4
#> [136,]    4    4
#> [137,]    4    4
#> [138,]    4    4
#> [139,]    4    5
#> [140,]    4    5
#> [141,]    4    5
#> [142,]    4    5
#> [143,]    4    5
#> [144,]    4    5
#> [145,]    4    5
#> [146,]    4    6
#> [147,]    4    6
#> [148,]    4    6
#> [149,]    4    6
#> [150,]    4    7
#> [151,]    4    7
#> [152,]    4    7
#> [153,]    4    7
#> [154,]    4    7
#> [155,]    4    7
#> [156,]    4    7
#> [157,]    4   10
#> [158,]    4   11
#> [159,]    4   13
#> [160,]    4   14
#> [161,]    4   15
#> [162,]    4   16
#> [163,]    5    1
#> [164,]    5    1
#> [165,]    5    1
#> [166,]    5    1
#> [167,]    5    1
#> [168,]    5    1
#> [169,]    5    1
#> [170,]    5    1
#> [171,]    5    1
#> [172,]    5    2
#> [173,]    5    2
#> [174,]    5    2
#> [175,]    5    2
#> [176,]    5    2
#> [177,]    5    3
#> [178,]    5    3
#> [179,]    5    4
#> [180,]    5    4
#> [181,]    5    4
#> [182,]    5    4
#> [183,]    5    4
#> [184,]    5    4
#> [185,]    5    4
#> [186,]    5    4
#> [187,]    5    5
#> [188,]    5    5
#> [189,]    5    5
#> [190,]    5    5
#> [191,]    5    5
#> [192,]    5    5
#> [193,]    5    5
#> [194,]    5    5
#> [195,]    5    5
#> [196,]    5    6
#> [197,]    5    6
#> [198,]    5    6
#> [199,]    5    7
#> [200,]    5    7
#> [201,]    5    7
#> [202,]    5    7
#> [203,]    5   13
#> [204,]    5   14
#> [205,]    5   14
#> [206,]    5   16
#> [207,]    5   16
#> [208,]    5   17
#> [209,]    6    1
#> [210,]    6    1
#> [211,]    6    1
#> [212,]    6    2
#> [213,]    6    2
#> [214,]    6    2
#> [215,]    6    2
#> [216,]    6    2
#> [217,]    6    2
#> [218,]    6    2
#> [219,]    6    2
#> [220,]    6    3
#> [221,]    6    3
#> [222,]    6    3
#> [223,]    6    4
#> [224,]    6    4
#> [225,]    6    4
#> [226,]    6    4
#> [227,]    6    4
#> [228,]    6    4
#> [229,]    6    5
#> [230,]    6    5
#> [231,]    6    6
#> [232,]    6    6
#> [233,]    6    6
#> [234,]    6    6
#> [235,]    6    6
#> [236,]    6    6
#> [237,]    6    6
#> [238,]    6    6
#> [239,]    6    7
#> [240,]    6    7
#> [241,]    6    7
#> [242,]    6    7
#> [243,]    6    7
#> [244,]    6    7
#> [245,]    6    7
#> [246,]    6    7
#> [247,]    6    7
#> [248,]    6   11
#> [249,]    6   11
#> [250,]    6   11
#> [251,]    6   13
#> [252,]    6   13
#> [253,]    6   17
#> [254,]    7    1
#> [255,]    7    1
#> [256,]    7    1
#> [257,]    7    1
#> [258,]    7    2
#> [259,]    7    2
#> [260,]    7    2
#> [261,]    7    3
#> [262,]    7    4
#> [263,]    7    4
#> [264,]    7    4
#> [265,]    7    4
#> [266,]    7    4
#> [267,]    7    4
#> [268,]    7    5
#> [269,]    7    5
#> [270,]    7    5
#> [271,]    7    6
#> [272,]    7    6
#> [273,]    7    6
#> [274,]    7    6
#> [275,]    7    6
#> [276,]    7    6
#> [277,]    7    6
#> [278,]    7    6
#> [279,]    7    6
#> [280,]    7    7
#> [281,]    7    7
#> [282,]    7    7
#> [283,]    7    7
#> [284,]    7    7
#> [285,]    7    7
#> [286,]    7    7
#> [287,]    7    7
#> [288,]    7    7
#> [289,]    7    7
#> [290,]    7    7
#> [291,]    7   10
#> [292,]    7   10
#> [293,]    7   11
#> [294,]    7   11
#> [295,]    7   13
#> [296,]    7   15
#> [297,]    8    2
#> [298,]    8    3
#> [299,]    8    8
#> [300,]    8    8
#> [301,]    8    8
#> [302,]    8    8
#> [303,]    8    8
#> [304,]    8    8
#> [305,]    8    8
#> [306,]    8    8
#> [307,]    8    8
#> [308,]    8    8
#> [309,]    8    8
#> [310,]    8    9
#> [311,]    8    9
#> [312,]    8    9
#> [313,]    8    9
#> [314,]    8    9
#> [315,]    8    9
#> [316,]    8    9
#> [317,]    8    9
#> [318,]    8    9
#> [319,]    8    9
#> [320,]    8   10
#> [321,]    8   10
#> [322,]    8   10
#> [323,]    8   10
#> [324,]    8   10
#> [325,]    8   10
#> [326,]    8   10
#> [327,]    8   10
#> [328,]    8   10
#> [329,]    8   11
#> [330,]    8   12
#> [331,]    8   12
#> [332,]    8   13
#> [333,]    8   13
#> [334,]    8   13
#> [335,]    8   13
#> [336,]    8   13
#> [337,]    8   13
#> [338,]    8   13
#> [339,]    8   13
#> [340,]    8   15
#> [341,]    8   15
#> [342,]    8   15
#> [343,]    8   15
#> [344,]    8   15
#> [345,]    8   15
#> [346,]    8   15
#> [347,]    8   15
#> [348,]    9    2
#> [349,]    9    8
#> [350,]    9    8
#> [351,]    9    8
#> [352,]    9    8
#> [353,]    9    8
#> [354,]    9    8
#> [355,]    9    8
#> [356,]    9    8
#> [357,]    9    9
#> [358,]    9    9
#> [359,]    9   10
#> [360,]    9   10
#> [361,]    9   10
#> [362,]    9   10
#> [363,]    9   12
#> [364,]    9   12
#> [365,]    9   13
#> [366,]    9   13
#> [367,]    9   13
#> [368,]    9   13
#> [369,]    9   13
#> [370,]    9   13
#> [371,]    9   13
#> [372,]    9   14
#> [373,]    9   14
#> [374,]    9   14
#> [375,]    9   14
#> [376,]    9   14
#> [377,]    9   15
#> [378,]    9   15
#> [379,]    9   15
#> [380,]    9   15
#> [381,]    9   15
#> [382,]    9   15
#> [383,]    9   15
#> [384,]    9   15
#> [385,]    9   15
#> [386,]    9   15
#> [387,]    9   16
#> [388,]    9   16
#> [389,]    9   17
#> [390,]   10    4
#> [391,]   10    7
#> [392,]   10    8
#> [393,]   10    8
#> [394,]   10    8
#> [395,]   10    8
#> [396,]   10    8
#> [397,]   10    8
#> [398,]   10    8
#> [399,]   10    8
#> [400,]   10    9
#> [401,]   10    9
#> [402,]   10    9
#> [403,]   10    9
#> [404,]   10   10
#> [405,]   10   10
#> [406,]   10   11
#> [407,]   10   11
#> [408,]   10   11
#> [409,]   10   11
#> [410,]   10   11
#> [411,]   10   12
#> [412,]   10   13
#> [413,]   10   15
#> [414,]   10   16
#> [415,]   10   16
#> [416,]   10   16
#> [417,]   10   16
#> [418,]   10   16
#> [419,]   10   16
#> [420,]   10   17
#> [421,]   10   17
#> [422,]   10   17
#> [423,]   10   17
#> [424,]   10   17
#> [425,]   10   17
#> [426,]   10   17
#> [427,]   10   17
#> [428,]   11    2
#> [429,]   11    2
#> [430,]   11    2
#> [431,]   11    4
#> [432,]   11    4
#> [433,]   11    6
#> [434,]   11    6
#> [435,]   11    7
#> [436,]   11    7
#> [437,]   11    8
#> [438,]   11    9
#> [439,]   11   10
#> [440,]   11   10
#> [441,]   11   10
#> [442,]   11   10
#> [443,]   11   10
#> [444,]   11   10
#> [445,]   11   10
#> [446,]   11   11
#> [447,]   11   11
#> [448,]   11   11
#> [449,]   11   11
#> [450,]   11   11
#> [451,]   11   11
#> [452,]   11   11
#> [453,]   11   11
#> [454,]   11   11
#> [455,]   11   12
#> [456,]   11   13
#> [457,]   11   13
#> [458,]   11   14
#> [459,]   11   15
#> [460,]   11   15
#> [461,]   11   15
#> [462,]   11   15
#> [463,]   11   16
#> [464,]   11   16
#> [465,]   11   16
#> [466,]   11   16
#> [467,]   11   16
#> [468,]   11   16
#> [469,]   11   17
#> [470,]   11   17
#> [471,]   11   17
#> [472,]   11   17
#> [473,]   12    2
#> [474,]   12    3
#> [475,]   12    7
#> [476,]   12    8
#> [477,]   12    8
#> [478,]   12    8
#> [479,]   12    8
#> [480,]   12    8
#> [481,]   12    8
#> [482,]   12    8
#> [483,]   12    9
#> [484,]   12   10
#> [485,]   12   10
#> [486,]   12   10
#> [487,]   12   10
#> [488,]   12   12
#> [489,]   12   12
#> [490,]   12   12
#> [491,]   12   12
#> [492,]   12   12
#> [493,]   12   14
#> [494,]   12   14
#> [495,]   12   14
#> [496,]   12   14
#> [497,]   12   14
#> [498,]   12   14
#> [499,]   12   14
#> [500,]   12   14
#> [501,]   12   14
#> [502,]   12   15
#> [503,]   12   16
#> [504,]   12   16
#> [505,]   12   16
#> [506,]   12   16
#> [507,]   12   16
#> [508,]   12   16
#> [509,]   12   16
#> [510,]   12   16
#> [511,]   12   16
#> [512,]   12   17
#> [513,]   12   17
#> [514,]   12   17
#> [515,]   12   17
#> [516,]   13    3
#> [517,]   13    4
#> [518,]   13    5
#> [519,]   13    6
#> [520,]   13    7
#> [521,]   13    8
#> [522,]   13    8
#> [523,]   13    8
#> [524,]   13    8
#> [525,]   13    8
#> [526,]   13    8
#> [527,]   13    8
#> [528,]   13    8
#> [529,]   13    9
#> [530,]   13    9
#> [531,]   13    9
#> [532,]   13    9
#> [533,]   13    9
#> [534,]   13    9
#> [535,]   13   10
#> [536,]   13   10
#> [537,]   13   10
#> [538,]   13   10
#> [539,]   13   11
#> [540,]   13   11
#> [541,]   13   12
#> [542,]   13   13
#> [543,]   13   13
#> [544,]   13   13
#> [545,]   13   13
#> [546,]   13   14
#> [547,]   13   14
#> [548,]   13   15
#> [549,]   13   15
#> [550,]   13   15
#> [551,]   13   15
#> [552,]   13   15
#> [553,]   13   15
#> [554,]   13   15
#> [555,]   13   16
#> [556,]   13   16
#> [557,]   13   17
#> [558,]   13   17
#> [559,]   14    1
#> [560,]   14    3
#> [561,]   14    4
#> [562,]   14    5
#> [563,]   14    5
#> [564,]   14    6
#> [565,]   14    8
#> [566,]   14    8
#> [567,]   14    8
#> [568,]   14    9
#> [569,]   14    9
#> [570,]   14    9
#> [571,]   14    9
#> [572,]   14    9
#> [573,]   14    9
#> [574,]   14   10
#> [575,]   14   10
#> [576,]   14   10
#> [577,]   14   12
#> [578,]   14   12
#> [579,]   14   12
#> [580,]   14   12
#> [581,]   14   12
#> [582,]   14   12
#> [583,]   14   12
#> [584,]   14   12
#> [585,]   14   13
#> [586,]   14   14
#> [587,]   14   14
#> [588,]   14   14
#> [589,]   14   14
#> [590,]   14   14
#> [591,]   14   15
#> [592,]   14   15
#> [593,]   14   16
#> [594,]   14   16
#> [595,]   14   16
#> [596,]   14   16
#> [597,]   14   16
#> [598,]   14   16
#> [599,]   14   16
#> [600,]   14   16
#> [601,]   14   16
#> [602,]   14   16
#> [603,]   15    2
#> [604,]   15    2
#> [605,]   15    7
#> [606,]   15    8
#> [607,]   15    8
#> [608,]   15    8
#> [609,]   15    8
#> [610,]   15    8
#> [611,]   15    9
#> [612,]   15    9
#> [613,]   15    9
#> [614,]   15    9
#> [615,]   15    9
#> [616,]   15    9
#> [617,]   15    9
#> [618,]   15    9
#> [619,]   15    9
#> [620,]   15   10
#> [621,]   15   11
#> [622,]   15   11
#> [623,]   15   11
#> [624,]   15   11
#> [625,]   15   11
#> [626,]   15   12
#> [627,]   15   13
#> [628,]   15   13
#> [629,]   15   13
#> [630,]   15   13
#> [631,]   15   13
#> [632,]   15   13
#> [633,]   15   13
#> [634,]   15   14
#> [635,]   15   14
#> [636,]   15   14
#> [637,]   15   14
#> [638,]   15   15
#> [639,]   15   15
#> [640,]   15   15
#> [641,]   15   15
#> [642,]   15   15
#> [643,]   15   15
#> [644,]   15   15
#> [645,]   15   15
#> [646,]   15   16
#> [647,]   15   16
#> [648,]   15   16
#> [649,]   15   17
#> [650,]   15   17
#> [651,]   16    1
#> [652,]   16    2
#> [653,]   16    3
#> [654,]   16    4
#> [655,]   16    5
#> [656,]   16    8
#> [657,]   16    8
#> [658,]   16    8
#> [659,]   16    8
#> [660,]   16    8
#> [661,]   16    9
#> [662,]   16   10
#> [663,]   16   11
#> [664,]   16   11
#> [665,]   16   11
#> [666,]   16   12
#> [667,]   16   12
#> [668,]   16   12
#> [669,]   16   12
#> [670,]   16   12
#> [671,]   16   12
#> [672,]   16   14
#> [673,]   16   14
#> [674,]   16   14
#> [675,]   16   14
#> [676,]   16   14
#> [677,]   16   14
#> [678,]   16   14
#> [679,]   16   14
#> [680,]   16   15
#> [681,]   16   15
#> [682,]   16   16
#> [683,]   16   16
#> [684,]   16   16
#> [685,]   16   16
#> [686,]   16   16
#> [687,]   16   17
#> [688,]   16   17
#> [689,]   16   17
#> [690,]   16   17
#> [691,]   16   17
#> [692,]   16   17
#> [693,]   16   17
#> [694,]   17    5
#> [695,]   17    6
#> [696,]   17    6
#> [697,]   17    8
#> [698,]   17    8
#> [699,]   17    8
#> [700,]   17    9
#> [701,]   17    9
#> [702,]   17    9
#> [703,]   17   10
#> [704,]   17   10
#> [705,]   17   10
#> [706,]   17   10
#> [707,]   17   10
#> [708,]   17   10
#> [709,]   17   10
#> [710,]   17   11
#> [711,]   17   11
#> [712,]   17   11
#> [713,]   17   11
#> [714,]   17   11
#> [715,]   17   12
#> [716,]   17   12
#> [717,]   17   13
#> [718,]   17   13
#> [719,]   17   13
#> [720,]   17   14
#> [721,]   17   14
#> [722,]   17   15
#> [723,]   17   16
#> [724,]   17   16
#> [725,]   17   16
#> [726,]   17   16
#> [727,]   17   16
#> [728,]   17   16
#> [729,]   17   16
#> [730,]   17   17
#> [731,]   17   17
#> [732,]   17   17
#> [733,]   17   17
#> [734,]   17   17
#> [735,]   17   17
#> [736,]   17   17
#> [737,]   17   17
#> [738,]   17   17
#> [739,]   17   17
#> [740,]   17   17
#> [741,]   17   17
#> [742,]   17   17
```

First, create data objects from internal data files, which are later
combined to the process state

``` r
# create objects
transfers <- createEdgelist(mobilityEdgelist, nodeSet = c("organisations", "organisations", "people"))
people <- createNodeSet(1:nrow(mobilityEdgelist))
organisations <- createNodeSet(1:length(orgRegion))
sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- createNetwork(sameRegion, nodeSet = c("organisations", "organisations"))
region <- createNodeVariable(orgRegion, nodeSet = "organisations")
size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
sex <- createNodeVariable(indSex, nodeSet = "people")
```

Combine the data objects into the process state, i.e., an object that
stores all information that will be used in the estimation later.

``` r
myState <- createProcessState(list(
  transfers = transfers,
  
  people = people,
  organisations = organisations,
  
  sameRegion = sameRegion,
  region = region,
  size = size,
  sex = sex
))
```

Define the dependent variable, and create a cache, a necessary object
used in the simulations. In case variables of the individuals in the
data are included in the state, they need to be explicitly mentioned in
the creation of the cache under “resourceCovariates”.

``` r
myDependentVariable <- "transfers"
myCache <- createWeightedCache(myState, myDependentVariable, resourceCovariates = c("sex"))
```

Specify the model. The predictors in the model are called “Effects” and
they are defined in a list. Each effect itself is a list that contains
the effect name and additional parameters that it needs.

``` r
# create an effects object
myEffects <- createEffectsObject(
  list(
    list("loops"),
    list("min_reciprocity"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex")
  )
)
```

Optional: run a pseudo-likelihood estimation to get improved initial
estimates. This increases the chances of cenvergence in the first run of
the estimation considerably

``` r
# create multinomial statistics object pseudo-likelihood estimation
myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects, myDependentVariable)

### additional script to get pseudo-likelihood estimates, requires the dfidx and mlogit package
# library(dfidx)
# library(mlogit)
# my.mlogit.dataframe <- dfidx(myStatisticsFrame,
#                           shape = "long", 
#                           choice = "choice")
# 
# colnames(my.mlogit.dataframe) <- gsub(" ", "_", colnames(my.mlogit.dataframe))
# 
# IVs <- (colnames(my.mlogit.dataframe)[2:(ncol(myStatisticsFrame)-2)])
# 
# f <- as.formula(paste("choice ~ 1 + ", paste(IVs, collapse = " + "), "| 0"))
# 
# my.mlogit.results <- mlogit(formula = eval(f), data = my.mlogit.dataframe, heterosc = F)
# 
# summary(my.mlogit.results)
#
# initEst <- my.mlogit.results$coefficients[1:length(IVs)]
```

Now estimate the model. The first two lines indicate the dependent
variable, data (state), cache, and effects. The third line specifies the
intial estimates, where the previously obtained pseudo-likelihood
estimates can be used. The remaining lines define the algorithm (see
helpfiles).

Running the model takes a while (up to 10 minutes for this data with
parallel computing).

``` r
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects,
  initialParameters = NULL,
  burnInN1 = 1500, iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1,
  burnInN2 = 7500, nsubN2 = 4, initGain = 0.2, thinningN2 = 1500,
  initialIterationsN2 = 25,
  iterationsN3 = 500, burnInN3 = 7500, thinningN3 = 3750,
  parallel = T, cpus = 4,
  allowLoops = T,
  verbose = T,
  returnDeps = T,
  multinomialProposal = F,
  fish = F
)
```

In case a pseudo-likelihood estimates have been obtained previously,
replace with

``` r
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects,
  initialParameters = initEst,
  burnInN1 = 1500, iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1,
  burnInN2 = 7500, nsubN2 = 4, initGain = 0.2, thinningN2 = 1500,
  initialIterationsN2 = 25,
  iterationsN3 = 500, burnInN3 = 7500, thinningN3 = 3750,
  parallel = T, cpus = 4,
  allowLoops = T,
  verbose = T,
  returnDeps = T,
  multinomialProposal = F,
  fish = F
)
```

Check convergence to see whether the results are reliable. In case the
maximum convergence ratio is above 0.1 (or 0.2 for less precise
estimates), another run is necessary.

``` r
max(abs(myResDN$convergenceStatistics))
#> [1] 0.04541426
```

Re-run estimation with previous results as starting values and check
convergence:

``` r
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects,
  prevAns = myResDN,
  burnInN1 = 1500, iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1,
  burnInN2 = 7500, nsubN2 = 4, initGain = 0.2, thinningN2 = 3000,
  initialIterationsN2 = 40,
  iterationsN3 = 500, burnInN3 = 15000, thinningN3 = 7500,
  parallel = T, cpus = 4,
  allowLoops = T,
  verbose = T,
  returnDeps = T,
  multinomialProposal = F,
  fish = F
)
```

``` r
# check convergence
max(abs(myResDN$convergenceStatistics))
#> [1] 0.04541426
```

In case convergence is still poor, updating the algorithm might be
necessary. Otherwise, view results, where the first column is the
estimate, the second column the standard error and the third column the
convergene ratio. All values in the finla column should be below 0.1
(see above).

``` r
myResDN
#> Results
#>                                   Effects   Estimates StandardErrors
#> 1                                   loops  2.59047738     0.18065574
#> 2                         min_reciprocity  0.81390378     0.19389546
#> 3             dyadic_covariate sameRegion  1.68755130     0.11189750
#> 4                    alter_covariate size  0.03616745     0.02483439
#> 5 resource_covar_to_node_covar region sex -0.64794362     0.16764479
#> 6                loops_resource_covar sex -0.36810363     0.21342903
#>   Convergence
#> 1 -0.03726664
#> 2  0.02279557
#> 3  0.04522141
#> 4  0.02760621
#> 5 -0.04541426
#> 6 -0.03304023
```

## Some diagnostics

Both indicate the extent to which the chain mixes (i.e., whether the
thinning was chosen appropriately). For the autoCorrelationTest, lower
values are better. Values above 0.5 are very problematic.

``` r
autoCorrelationTest(myDependentVariable, myResDN)
#> [1] 0.1027527
```

For the extractTraces, the plot should show data point randomly
scattered around the target line.

``` r
traces <- extractTraces(myDependentVariable, myResDN, myEffects)

par(mfrow = c(1,2))
plot(traces)
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-16-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-16-3.png" width="100%" />

## score-tests

Based on an estimated model, a score-type test is available that shows
whether statistics representing non-inlcuded effects are well
represented. If this is not the case, it is likely that including them
will result in significant estimates.

Note that it is advisable that the model specification that is tested
includes all effects from the previously estimated model.

``` r
myEffects2 <- createEffectsObject(
  list(
    list("loops"),
    list("min_reciprocity"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex"),
    list("min_transitivity")
  )
)

test_ME.2 <- scoreTest(myDependentVariable, myResDN, myEffects2)
test_ME.2
#> Results
#>                                   Effects pValuesParametric
#> 1                                   loops      9.702724e-01
#> 2                         min_reciprocity      9.818133e-01
#> 3             dyadic_covariate sameRegion      9.639308e-01
#> 4                    alter_covariate size      9.779762e-01
#> 5 resource_covar_to_node_covar region sex      9.637771e-01
#> 6                loops_resource_covar sex      9.736425e-01
#> 7                        min_transitivity      2.428713e-09
#>   pValuesNonParametric
#> 1                0.556
#> 2                0.526
#> 3                0.494
#> 4                0.476
#> 5                0.562
#> 6                0.532
#> 7                0.000
#> 
#>  Parametric p-values: small = more significant 
#>  Non-parametric p-values: further away from 0.5 = more significant
```

The interpretation is that there appears to be some transitive
clustering in the data that the model does not account for in its
current form.

## GOF testing

Akin to ERGMs, goodness of fit testing is available to see whether
auxiliary statistics are well captured by the model.

``` r
myGofIndegree <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getIndegree, lvls = 1:100)
plot(myGofIndegree)
```

<img src="man/figures/README-unnamed-chunk-18-1.png" width="100%" />

``` r

myGofTieWeight <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getTieWeights, lvls = 1:30)
plot(myGofTieWeight)
```

<img src="man/figures/README-unnamed-chunk-18-2.png" width="100%" />
