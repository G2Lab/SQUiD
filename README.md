
# SQUiD

Ultra-Secure Storage and Analysis of Genetic Data for the Advancement of Precision Medicine
  
## Building Project

###  Dependencies

Download and install:

- [patchelf](https://github.com/NixOS/patchelf) (tested on v0.14.3-1)
- [m4](https://www.gnu.org/software/m4/) (tested on v1.4.19-3)
- [Google Benchmark](https://github.com/google/benchmark) (tested on v1.8.3)

### Building project from source

After cloning the repository, run the following commands to build the project from source:

```
git submodule update --init --recursive
mkdir mylibs
cd HElibPublicKeySwitch
mkdir build
cd build
cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=../../mylibs ..
make && make install
cd ../..
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=./mylibs -Dhelib_DIR=./mylibs/helib_pack/share/cmake/helib
make
```

### Sample Run

We have included a sample database and query script to demonstrate the functionalities of SQUiD. After building the project, run `../bin/main` from within the build directory to see our sample output (which will match the example below).

```
Initializing database
Found file with comparison polynomial coefficients
m = 17293, p = 131, phi(m) = 17292
  ord(p) = 3
  normBnd = 1.27324
  polyNormBnd = 1.27324
  factors = [17293]
  generator 2 has order (== Z_m^*) of 5764
r = 1
nslots = 5764
hwt = 0
ctxtPrimes = [6,7,8,9,10,11,12,13]
specialPrimes = [14,15,16]
number of bits = 595

security level = 92.0971

Security: 92.0971
Num slots: 5764
Num rows:10
Num compressed rows: 1
Number of snp columns: 2
Number of binary pheno columns: 1
Number of continuous pheno columns: 0
Printing DB
-----------------------------------------------------
|snp 0|snp 1|
--------------
|0     0     |
|0     1     |
|0     2     |
|1     0     |
|1     1     |
|1     2     |
|2     0     |
|2     1     |
|2     2     |
|0     1     |
|ALS|
--------------
|0   |
|0   |
|0   |
|0   |
|1   |
|1   |
|0   |
|1   |
|0   |
|0   |
Running sample queries:
-----------------------------------------------------
Running Counting query (snp 0 = 0 and snp 1 = 1)
Count: 2
Running Counting query (snp 0 = 0 and snp 1 = 2)
Count: 1
Running Counting query (snp 0 = 1 or snp 1 = 2)
Count: 5
Running MAF query filter (snp 1 = 1, target snp 0)
Numerator: 3
Denominator: 8
Computed MAF: 0.375
Running PRS query (snps [0,1], effect-sizes [2,5])
0, 5, 10, 2, 7, 12, 4, 9, 14, 5
Running Similarity query (d: snp 0 = 2 and snp 1 = 2, target = ALS, threshold < 4)
Encrypting...
Running similarity query...
Count with target:   3
Count without target:1
```

### Build and Demo Setup Information

We verified the steps to install SQUiD and run the demo on a medium sized Google Cloud machine (n2-standard-64) running *Ubuntu 22.04 TLS*.

On this machine, it took a few minutes to install HElib, less than a minute to install SQUiD, and less than a minute to generate the output of the demo.

## Installing the SQUiD API

To run the API for SQUiD which the SQUiD CLI will interact with, install the [Drogon](https://github.com/drogonframework/drogon) framework using this [guide](https://github.com/drogonframework/drogon/wiki/ENG-02-Installation).

Then, to install run the following commands inside the `./build/` directory:

```

cmake -DBUILD_SQUID_API=ON ..
make

```


After completing these steps, you can run the API with with `../bin/PIRAPI`.

Once the API has started up, you can send queries using the `../bin/squid`.

### Modifying SQUiD API to use different IP address

By default the SQUiD API and CLI run over the address `127.0.0.1`, but this can be changed to a server's IP address by modifying to following files:

* In `config.json`, line 15 - 21 will look like this,

```

"listeners": [
    {
        "address": "127.0.0.1",
        "port": 8081,
        "https": false
    }
],
```

  

Change the address and port here to your server's IP address and port.


## Using SQUiD CLI
To begin using the SQUiD CLI, first create an `apidata` directory by running `mkdir apidata` in the build directory.

Then, configure the SQUiD CLI with `../bin/squid config [address] [port] [api_key]`. For default configurations, run `../bin/squid config 127.0.0.1 8081 nNCHuSdBWZsDJNFOJqUWDAUibEvVcVniRqbiIoM`.

Then, you need to pull the context from the server (make sure the server is running before running this command) with `../bin/squid getContext`.

Then, you need to generate a public / private key pair with `../bin/squid genKeys`.

Then, you need to authorize with `../bin/squid authorize`.

Then, you need to download the data owner's public key with `../bin/squid getOwnerPublicKey`.

After completing these steps, you can begin running the queries. For help with the query parameters, just run `../bin/squid`. The output will be the following,

```
Welcome to SQUiD!
--- Setup ---
Config server address, port, and API key: ../bin/squid config <address> <port> <api_key>
Pull context from server: ../bin/squid getContext
Generate own context: ../bin/squid genContext
Generate public / secret key: ../bin/squid genKeys
Authorize yourself to the server (by generating key-switching key): ../bin/squid authorize
Get headers: ../bin/squid getHeaders
Get data owner's public key: ../bin/squid getOwnerPublicKey

--- Query ---
Query: ../bin/squid <option> [query_string]

--- Helper ---
Decrypt query results: ../bin/squid decrypt <file>

```

## Experimenting with VCF Data

We have include example data for 10 snps for 6 subjects for chromosome 22.

  
| Name | Position |
|--------------|-------------| 
| rs9617549 | 22:10874444 | 
| rs577013928 | 22:10874535 | 
| rs565082899 | 22:10874551 |
| rs540382744 | 22:10874556 |
| rs573244332 | 22:10874564 | 
| rs539162497 | 22:11121568 |
| rs557291239 | 22:11121677 | 
| rs569309845 | 22:11121789 | 
| rs536692189 | 22:11121839 | 
| rs556567876 | 22:11122005 |

This data can be found in a vcf file at `./data/sample.vcf`.

We have also included some example phenotype data in `./data/sample.csv`.

Running `../bin/PIRAPI` will load this vcf and csv into the SQUiD database, where you can experiment on it with by running queries.

Replacing this vcf/csv file or changing in the data that is loaded in `./src/API/controllers/api_v1_Server.cc` will allow you to experiment on your own data.
  
*Note that the parameters used in the demo are small for test purposes only. To use secure parameters, replace the parameter type `constants::Test` with `constants::Large` in `./src/API/controllers/api_v1_Server.cc`. The demo also does not have similarity queries turned on because it takes some time to load in the comparison polynomial. Changing the boolean parameter after the HE parameters on line 78 in `./src/API/controllers/api_v1_Server.cc` to true will allow for similarity queries (`Server::Server() : squid(constants::Test, true)`).

## Benchmarking

Before running any benchmarking scripts, ensure that a disk database has been created. This can be done by running `mkdir -p db_disk` in the `data` directory. To generate synthetic or use real data, run `./data_owner_encryptor` and either specify a number of rows and columns to generate or pass in a vcf / csv of genotypes and phenotypes. 

- `./data_owner_encryptor --row 1024 --columns 16`
- `./data_owner_encryptor --vcf sample.vcf --csv sample.csv`

To run our benchmarking scripts, run `make bench` followed by `../bin/bench` to run the main benchmarking scripts for timings and communication required for all queries.

To run our secondary benchmarking script, run `make misc` followed by `../bin/misc` for all timings for database updates, key-switching, and database encryption.

## Acknowledgements

SQUiD utilizes the comparator from this [repository](https://github.com/iliailia/comparison-circuit-over-fq) which is an implementation of this [paper](https://eprint.iacr.org/2021/315) by Ilia Iliashenko and Vincent Zucca.
