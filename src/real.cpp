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
    Server server = Server(constants::P131, false);

    server.SetData("./data/chr22_100samples_10SNPs.vcf");

    std::cout << "Database:" << std::endl;
    server.PrintEncryptedDB(true);

    vector<string> headers = server.GetHeaders();
    
    vector<pair<uint32_t, uint32_t>> query;
    cout << "Running Counting query (snp " << headers[0] << " = 0 or snp " << headers[1] << " = 1)" << endl;
    query = vector<pair<uint32_t, uint32_t>>{pair(0,0), pair(1,1)};
    helib::Ctxt result = server.CountQuery(false, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;
    return 0;
}