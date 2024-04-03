#include "plain.h"
#include <openssl/aes.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/err.h>

#include <cstring>

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------

//                                                         HELPER FUNCTIONS

// ------------------------------------------------------------------------------------------------------------------------



void handleErrors(void)
{
}
int encrypt(unsigned char *plaintext, int plaintext_len, unsigned char *key,
            unsigned char *iv, unsigned char *ciphertext){
                return 0;
}
int decrypt(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
            unsigned char *iv){
                return 0;
}



template<typename T, typename Allocator>
void print_vector(const vector<T, Allocator>& vect, int num_entries)
{
    std::cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++){
        std::cout << ", " << vect[i];
    } 
    std::cout << std::endl;
}

Plain::Plain(){
    db_set = false;
}


void Plain::GenData(int _num_rows, int _num_cols){
    num_rows = _num_rows;
    num_cols = _num_cols; 

    //create a sample key
        /* A 256 bit key */
    key = new unsigned char [32]{ 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
                           0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
                           0x36, 0x37, 0x38, 0x39, 0x30, 0x31, 0x32, 0x33,
                           0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x30, 0x31
                         };

    /* A 128 bit IV */
    iv = new unsigned char[16] { 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
                          0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35
                        };
    
    encrypted_db = vector<vector<unsigned char*>>(num_rows);
    for (int r = 0; r < num_rows; r++){
        encrypted_db[r] = vector<unsigned char*>(num_cols);
        for (int c = 0; c < num_cols; c++){
            encrypted_db[r][c] = new unsigned char[128];
            unsigned char *plaintext = (unsigned char *)"0";

            int ciphertext_len = encrypt (plaintext, strlen ((char *)plaintext), key, iv,
                              encrypted_db[r][c]);

        }
    }
}

int Plain::decode(int i, int j){
    return decrypt(encrypted_db[i][j], 16, key, iv);
    
}

int Plain::CountingQuery(bool conjunctive, vector<pair<int, int>>& query){
    int count = 0;
    if (conjunctive){
        for (int r = 0; r < num_rows; r++){
            bool passes_filter = true;
            for (int c = 0; c < query.size(); c++){
                int col = query[c].first;
                int value = query[c].second;
                if (decode(r,col) != value){
                    passes_filter = false;
                }
            }
            if (passes_filter){
                count += 1;
            }
        }
    }
    else{
        for (int r = 0; r < num_rows; r++){
            bool passes_filter = false;
            for (int c = 0; c < query.size(); c++){
                int col = query[c].first;
                int value = query[c].second;
                if (decode(r,col) == value){
                    passes_filter = true;
                }
            }
            if (passes_filter){
                count += 1;
            }
        }
    }
    return count;
}

int Plain::MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query){
    int yes = 0;
    int no = 0;
    if (conjunctive){
        for (int r = 0; r < num_rows; r++){
            bool passes_filter = true;
            for (int c = 0; c < query.size(); c++){
                int col = query[c].first;
                int value = query[c].second;
                if (decode(r,col) != value){
                    passes_filter = false;
                }
            }
            if (passes_filter){
                yes += decode(r,snp);
                no += 2 - decode(r,snp);
            }
        }
    }
    else{
        for (int r = 0; r < num_rows; r++){
            bool passes_filter = false;
            for (int c = 0; c < query.size(); c++){
                int col = query[c].first;
                int value = query[c].second;
                if (decode(r,col) == value){
                    passes_filter = true;
                }
            }
            if (passes_filter){
                yes += decode(r,snp);
                no += 2 - decode(r,snp);
            }
        }
    }
    return yes + no;
}

vector<int> Plain::DistributionQuery(vector<pair<int, int>>& prs_params){
    vector<int> scores = vector<int>(num_rows);
    for (int r = 0; r < num_rows; r++){
        int score = 0;
        for (int c = 0; c < prs_params.size(); c++){
            int col = prs_params[c].first;
            int weight = prs_params[c].second;
            score += weight * decode(r,col);
        }
        scores[r] = score;
    }
    return scores;
}

pair<int, int> Plain::SimilarityQuery(int target_column, vector<int>& d, int threshold){
    int yes = 0;
    int no = 0;
   
    for (int r = 0; r < num_rows; r++){
        int score = 0;
        for (int c = 0; c < num_cols; c++){
            score += pow(decode(r,c) - d[c], 2);
        }
        if (score < threshold){
            yes += decode(r,target_column);
            no += 1 - decode(r,target_column);
        }
    }
    return pair(yes,no);
}

vector<pair<int, int>> Plain::ChiSquareQuery(int disease_column, int number_of_chi){

}