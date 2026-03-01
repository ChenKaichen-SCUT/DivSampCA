## Overview


## How to Obtain 

git clone https://github.com/ChenKaichen-SCUT/DivSampCA-master.git

## Installation for Building *DivSampCA*

cd DivSampCA-master/

make

## Instructions for Running *DivSampCA*

After building *DivSampCA*, users may run it with the following command: 

./DivSampCA -input_cnf_path [INSTANCE_PATH] <optional_parameters> <optional_flags>

"-input_cnf_path" represents the path of the input CNF file.

For the optional parameters, we list them as follows. 

| Parameter | Allowed Values | Default Value | Description | 
| - | - | - | - |
| `-output_testcase_path` | string | `./[INPUT_CNF_NAME]_testcase_set.txt` | path to which the generated PCA is saved |
| `-seed` | integer | 1 | random seed | 
| `-c` | positive integer | 100 | the number of candidates per round | 
| `-t` |positive integer | 100 |  the threshold for transitioning from sampling phase to full coverage phase | 

## Example 

```
./DivSampCA -input_cnf_path CNF_benchmarks/aaed2000.cnf -seed 1 -c 100 -t 100
```

## License

*DivSampCA* uses the GPL-3.0 license. Check [LICENSE.md](LICENSE.md) for more information. 
