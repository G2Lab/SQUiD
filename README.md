# SQUiD
Ultra-Secure Storage and Analysis of Genetic Data for the Advancement of Precision Medicine

--------------------------------
Building Project
--------------------------------


Download and compile our [fork](https://github.com/SQUiDCodebase/HElibPublicKeySwitch) of HElib (make sure you install the [patchelf](https://github.com/NixOS/patchelf) (tested on v0.14.3-1) and [m4](https://www.gnu.org/software/m4/) (tested on v1.4.19-3) dependencies).
Use the guide [here](https://github.com/homenc/HElib/blob/master/INSTALL.md) and follow **Option 1** to install HElib to a `mylibs` directory in the root directory.

```
git clone https://github.com/SQUiDCodebase/HElibPublicKeySwitch
mkdir mylibs
cd HElibPublicKeySwitch
mkdir build
cd build
cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX= PATH_TO_MYLIBS ..
make && make install
```

Download [Google Benchmark](https://github.com/google/benchmark) (tested on v1.8.3).

From the root directory, to build run the following:

```
mkdir build
cd build
../scripts/make.sh
make
```

--------------------------------
Sample Run
--------------------------------

We have include a sample DB and query script to demonstrate the functionalities of SQUiD. After running the make command, run `./bin/main` to see our sample output (which will be the same as below).

```
Initialising context object...
-----------------------------------------------------
Printing DB
|snp1|snp2|ALS|
--------------
|0    0   0   |
|0    1   0   |
|0    2   0   |
|1    0   0   |
|1    1   1   |
|1    2   1   |
|2    0   0   |
|2    1   1   |
|2    2   0   |
|0    1   0   |
Running sample queries:
-----------------------------------------------------
Running Counting query (snp 0 = 0 and snp 1 = 1)
Count: 2
Running Counting query (snp 0 = 0 and snp 1 = 2)
Count: 1
Running Counting query (snp 0 = 1 or snp 1 = 2)
Count: 5
Running MAF query filter (snp 0 = 1 or snp 1 = 2), target snp = 0
Nom: 5
Dom: 10
Computed MAF: 0.5
Running PRS query (snps [0,1], effect-sizes [2,5])
0, 5, 10, 2, 7, 12, 4, 9, 14, 5
Running Similarity query (d: snp 0 = 2 and snp 1 = 2, target = ALS)
Count with target:   2
Count without target:1
```
--------------------------------
Build and Demo Setup Information
--------------------------------

We verified the steps to install SQUiD and run the demo on a medium sized AWS machine (c4.8xlarge) running *Ubuntu 22.04 TLS*. 
On this machine, it took a few minutes to install HElib, less than a minute to install SQUiD, and less than a minute to generate the output of the demo.

--------------------------------
Experimenting with Real Data
--------------------------------

We have include the first 10 snps for 100 subjects for chromosome 22 using [1000 Genomes](https://www.internationalgenome.org/) from the [phase 3](https://www.internationalgenome.org/category/phase-3) data set (which uses GRCh38). The first 10 SNPs are:

| Name | Position |
|-----------------|-----------------|
| rs9617549   | 22:10874444    |
| rs577013928    | 22:10874535    |
| rs565082899   | 22:10874551    |
| rs540382744   | 22:10874556    |
| rs573244332    | 22:10874564    |
| rs539162497   | 22:11121568    |
| rs557291239   | 22;11121677    |
| rs569309845    | 22:11121789    |
| rs536692189   | 22:11121839    |
| rs556567876   | 22:11122005    |


This data can be found in a vcf file at `data/chr22_100samples_10SNPs.vcf`. Running `./bin/real` will load this vcf into SQUiD, where you can experiment on it with our queries. 

Replacing this vcf file or changing in the path in `real.cpp` will allow you to experiment on your own vcf data.

--------------------------------
Installing SQUiD API
--------------------------------

To run the API for SQUiD which the SQUiD CLI will communicate with, install the [Drogon](https://github.com/drogonframework/drogon) framework using this [guide](https://github.com/drogonframework/drogon/wiki/ENG-02-Installation).

Then, to install run the following steps:

* Navigate to the API directory with `cd API`
* Create a build directory with `mkdir build`
* Navigate to the build directory with `cd build`
* Run cmake with `cmake ..`
* Run make with `make`

After completing these steps, you can run the API with  with `./PIRAPI`.

Once the API has started up, you can send queries using the `./bin/squid` from the root directory.

### Modifying SQUiD API to use different IP address

By default the SQUiD API and CLI run over the address `localhost`, but this can be changed to a server's IP address by modifying to following files:

* In `./API/config.json`, line 15 - 21 will look like this,
```
"listeners": [
    {
        "address": "localhost",
        "port": 8081,
        "https": false
    }
],
```

Change the address and port here to your server's IP address and port.

--------------------------------
Using SQUiD CLI
--------------------------------

To begin using the SQUiD CLI, first create an `apidata` directory by running `mkdir apidata` in the root directory.

Then, configure the SQUiD CLI with `./bin/squid config [address] [port] [api_key]`. For default configurations, run `./bin/squid config localhost 8081 nNCHuSdBWZsDJNFOJqUWDAUibEvVcVniRqbiIoM`.

Then, you need to pull the context from the server (make sure the server is running before running this command) with `./bin/squid getContext`.

Then, you need to generate a public / private key pair with `./bin/squid genKeys`.

Then, you need to authorize with `./bin/squid authorize`.

After completing these steps, you can begin running the queries. For help with the query parameters, just run `./bin/squid`. The output will be the following,

```
Welcome to SQUiD!
--- Setup ---
Config server address, port, and API key: ./bin/squid config <address> <port> <api_key>
Pull context from server: ./bin/squid getContext
Generate own context: ./bin/squid genContext
Generate public / secret key: ./bin/squid genKeys
Authorize yourself to the server (by generating key-switching key): ./bin/squid authorize

--- Query ---
Query: ./bin/squid <option> [query_string]

--- Helper ---
Decrypt query results (for queries not automatically decrypted): ./bin/squid decrypt <file>
```

--------------------------------
Benchmarking
--------------------------------

To run our benchmarking scripts, run `make bench` followed by `./bin/bench` to run the main benchmarking scripts for timings and communication required for all queries.
To run our secondary benchmarking script, run `make misc` followed by `./bin/misc` for all timings for database updates, key-switching, and database encryption. 

--------------------------------
Acknowledgements
--------------------------------

SQUiD utilizes the comparator from this [repository](https://github.com/iliailia/comparison-circuit-over-fq) which is an implementation of this [paper](https://eprint.iacr.org/2021/315) by Ilia Iliashenko and Vincent Zucca.

