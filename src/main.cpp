

#include <iostream>

#include <helib/helib.h>
#include "server.hpp"
#include "globals.hpp"
#include "tools.hpp"


template<typename T, typename Allocator>
void print_vector(const vector<T, Allocator>& vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++){
        cout << ", " << vect[i];
    } 
    cout << endl;
}


int main()
{
    /*  Example of BGV scheme  */
    
    std::cout << "Initialising context object..." << std::endl;

    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    
    Server server = Server(constants::P131, true);
    server.PrintContext();

    vector<vector<uint32_t>> fake_db = vector<vector<uint32_t>>{
        {0, 0, 0, 1, 1, 1, 2, 2, 2, 0},
        {0, 1, 2, 0, 1, 2, 0, 1, 2, 1},
        {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}};
    vector<string> headers = vector<string>{"snp1", "snp2", "ALS"};

    server.SetData(fake_db);
    server.SetColumnHeaders(headers);

    cout << "Printing DB" << endl;
    cout << "-----------------------------------------------------" << endl;
    server.PrintEncryptedDB(true);

    cout << "Running sample queries:" << endl;
    cout << "-----------------------------------------------------" << endl;

    vector<pair<uint32_t, uint32_t>> query;
    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0,0), pair(1,1)};
    helib::Ctxt result = server.CountQuery(true, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 2)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0,0), pair(1,2)};
    result = server.CountQuery(true, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 1 or snp 1 = 2)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0,1), pair(1,2)};
    result = server.CountQuery(false, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running MAF query filter (snp 0 = 1 or snp 1 = 2), target snp = 0" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0,1), pair(1,2)};
    helib::Ctxt result_pair = server.MAFQuery(0, false, query);
    vector<long> result_vector = server.Decrypt(result_pair);
    double nom = (double)result_vector[0];
    double denom = 2 * (double)result_vector[1];
    cout << "Nom: " << nom << endl;
    cout << "Denom: " << denom << endl;
    double AF = nom / denom;
    cout << "Computed MAF: " << min(AF, 1 - AF) << endl;

    cout << "Running PRS query (snps [0,1], effect-sizes [2,5])" << endl;
    auto params = vector<pair<uint32_t, int32_t>>{pair(0,2), pair(1,5)};
    vector<helib::Ctxt> results_distribution = server.PRSQuery(params);
    print_vector(server.Decrypt(results_distribution[0]));
    cout << "Running Similarity query (d: snp 0 = 2 and snp 1 = 2, target = ALS)" << endl;
    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    
    std::cout << "Encrypting..." << std::endl;
    d.push_back(server.Encrypt(2));
    d.push_back(server.Encrypt(2));
    
    std::cout << "Running similarity query..." << std::endl;
    pair<helib::Ctxt, helib::Ctxt> result_sim = server.SimilarityQuery(2, d, 2);
    cout << "Count with target:   " << server.Decrypt(result_sim.first)[0] << endl;
    cout << "Count without target:" << server.Decrypt(result_sim.second)[0] << endl;

    return 0;
}