#pragma once

#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <map>


using namespace std;


class Plain
{
public:
    //Setup
    Plain();
    void GenData(int _num_rows, int _num_cols);
    
    //Querries
    int CountingQuery(bool conjunctive, vector<pair<int, int>>& query);
    int MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query);
    vector<int> DistributionQuery(vector<pair<int, int>>& prs_params);
    pair<int, int> SimilarityQuery(int target_column, vector<int>& d, int threshold);
    vector<pair<int, int>> ChiSquareQuery(bool conjunctive, vector<pair<int, int>>& query, int disease_column, int number_of_chi);
    vector<pair<int, int>> ChiSquareQuery(int disease_column, int number_of_chi);

    int decode(int i, int j);
    
private:
    bool db_set;
    
    int num_rows;
    int num_cols;

    vector<vector<unsigned char*>> encrypted_db; 
    unsigned char *key;
    unsigned char *iv;
};