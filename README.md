## 1. Introduction
This repository contains the implementation of our paper "*A Tuple-Oriented Sampling Method for Generating Small Pairwise Covering Arrays in Configurable Software Systems*".

## 2. Software Requirements
| Dependency       | Minimum Version | Installation Command (Ubuntu/Debian) |
|------------------|-----------------|---------------------------------------|
| G++ (GCC)        | 7.0+            | `sudo apt install g++`                |
| zlib1g-dev       | 1.2.8+          | `sudo apt install zlib1g-dev`         |
| GNU Make         | 3.81+           | `sudo apt install make`               |


## 3.Installation for Building *DivSampCA*

git clone https://github.com/ChenKaichen-SCUT/DivSampCA-master.git

cd DivSampCA-master/

make

## 4. Instructions for Running *DivSampCA*

After building *DivSampCA*, users may run it with the following command: 

./DivSampCA -input_cnf_path [INSTANCE_PATH] <optional_parameters> <optional_flags>

For the optional parameters, we list them as follows. 

| Parameter       | Allowed Values | Default Value | Description | 
|------------------|-----------------|---------------------------------------|---------------------------------------|
| `-output_testcase_path` | string | `./[INPUT_CNF_NAME]_testcase_set.txt` | path to which the generated PCA is saved |
| `-seed` | integer | 1 | random seed | 
| `-lambda` | positive integer | 150 | number of candidates per round | 
| `-eta` |positive integer | 150 |  threshold for transitioning from sampling phase to full coverage phase |

## 5. Example 

```
./DivSampCA -input_cnf_path CNF_benchmarks/aaed2000.cnf -seed 100 -lambda 100 -lambda 100
```

## 6. License

*DivSampCA* uses the GPL-3.0 license. Check [LICENSE.md](LICENSE.md) for more information. 
