

#include <iostream>

#include <helib/helib.h>
#include "./databases/FHE_SIMD_database.hpp"
#include "globals.hpp"
#include "./databases/tools.hpp"

template <typename T, typename Allocator>
void print_vector(const vector<T, Allocator> &vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++)
    {
        cout << ", " << vect[i];
    }
    cout << endl;
}

int main()
{
    std::cout << "Initialising database" << std::endl;

    FHESIMDDatabase db = FHESIMDDatabase(constants::SmallFastComp, true);

    vector<vector<int32_t>> fake_db = vector<vector<int32_t>>{
        {0, 0, 0, 1, 1, 1, 2, 2, 2, 0},
        {0, 1, 2, 0, 1, 2, 0, 1, 2, 1}};

    vector<vector<int32_t>> fake_binary_pheno_db = vector<vector<int32_t>>{
        {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}};

    vector<string> headers = vector<string>{"snp 0", "snp 1"};
    vector<string> binary_pheno_headers = vector<string>{"ALS"};

    db.setData(fake_db);
    db.setBinaryPhenoData(fake_binary_pheno_db);

    db.printMeta();

    db.setColumnHeaders(headers);
    db.setBinaryPhenoHeaders(binary_pheno_headers);

    cout << "Printing DB" << endl;
    cout << "-----------------------------------------------------" << endl;
    db.printData(true);
    db.printPhenoData(true);

    cout << "Running sample queries:" << endl;
    cout << "-----------------------------------------------------" << endl;

    vector<pair<uint32_t, uint32_t>> query;
    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 1)};
    helib::Ctxt result = db.countQuery(true, query);
    cout << "Count: " << db.decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 2)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 0), pair(1, 2)};
    result = db.countQuery(true, query);
    cout << "Count: " << db.decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 1 or snp 1 = 2)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0, 1), pair(1, 2)};
    result = db.countQuery(false, query);
    cout << "Count: " << db.decrypt(result)[0] << endl;

    cout << "Running MAF query filter (snp 1 = 1, target snp 0)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(1, 1)};
    helib::Ctxt result_pair = db.MAFQuery(0, false, query);
    vector<long> result_vector = db.decrypt(result_pair);
    double numerator = (double)result_vector[0];
    double denominator = (double)result_vector[1];
    cout << "Numerator: " << numerator << endl;
    cout << "Denominator: " << denominator << endl;
    double AF = numerator / denominator;
    cout << "Computed MAF: " << min(AF, 1 - AF) << endl;

    cout << "Running PRS query (snps [0,1], effect-sizes [2,5])" << endl;
    auto params = vector<pair<uint32_t, int32_t>>{pair(0, 2), pair(1, 5)};
    vector<helib::Ctxt> results_distribution = db.PRSQuery(params);
    print_vector(db.decrypt(results_distribution[0]));
    cout << "Running Similarity query (d: snp 0 = 2 and snp 1 = 2, target = ALS, threshold < 4)" << endl;
    vector<helib::Ctxt> d = vector<helib::Ctxt>();

    std::cout << "Encrypting..." << std::endl;
    d.push_back(db.encrypt(2));
    d.push_back(db.encrypt(2));

    std::cout << "Running similarity query..." << std::endl;
    pair<helib::Ctxt, helib::Ctxt> result_sim = db.similarityQuery(0, d, 4);
    cout << "Count with target:   " << db.decrypt(result_sim.first)[0] << endl;
    cout << "Count without target:" << db.decrypt(result_sim.second)[0] << endl;

    return 0;
}