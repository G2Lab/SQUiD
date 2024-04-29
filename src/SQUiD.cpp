#include <iostream>
#include <helib/helib.h>
#include <chrono>

#include "server.hpp"
#include "globals.hpp"


#include <curl/curl.h>
#include <string>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unistd.h>



using json = nlohmann::json;

std::string API_DATA_DIR = "./apidata";
std::string CONTEXT_FILE = "context";
std::string SEC_KEY_FILE = "sec_key";
std::string PUBLIC_KEY_FILE = "pub_key";
std::string CONFIG_FILE = "config";
std::string GET_CONTEXT_API_CALL = "getcontext";
std::string POST_AUTH_API_CALL = "au";

std::string GET_COUNT_API_CALL = "countingQuery";
std::string GET_MAF_API_CALL = "mafQuery";
std::string GET_PRS_API_CALL = "PRSQuery";
std::string GET_HEADERS_API_CALL = "headers";
std::string GET_PRINT_DB_API_CALL = "printDB";


// A helper function to send HTTP requests using libcurl

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

void getTime(){
    std::cout << "\033[38;5;202m";

    // Get the current time
    time_t now = time(0);
    tm* timeinfo = localtime(&now);

    // Format the time as a string
    char timeString[9]; // Assumes time will always be in HH:MM:SS format
    strftime(timeString, sizeof(timeString), "%T", timeinfo);

    // Output the time in orange color
    std::cout << timeString << ": ";

    // Reset text color to default
    std::cout << "\033[0m";
}

std::string sendHttpRequest(std::string url) {
  CURL *curl;
  CURLcode res;
  std::string response;

  //std::cout << "Sending url:" << endl;
  //std::cout << url << endl;

  curl = curl_easy_init();

  if(curl) {
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
  }
  return response;
}

std::string postHttpRequest(std::string& url, std::string& data){
  //std::cout << "Sending url:" << endl;
  //std::cout << url << endl;
  CURL *curl;
  CURLcode res;

  /* In windows, this will init the winsock stuff */
  curl_global_init(CURL_GLOBAL_ALL);

  /* get a curl handle */
  curl = curl_easy_init();
  if(curl) {
    /* First set the URL that is about to receive our POST. This URL can
       just as well be an https:// URL if that is what should receive the
       data. */
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    /* Now specify the POST data */
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data.c_str());

    /* Perform the request, res will get the return code */
    res = curl_easy_perform(curl);
    /* Check for errors */
    //if(res != CURLE_OK)
    //  fprintf(stderr, "curl_easy_perform() failed: %s\n",
    //          curl_easy_strerror(res));

    /* always cleanup */
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
  return "success";
}

string conjunctiveToString(string conjunctive){
  if (conjunctive == "1"){
    return "And";
  }
  if (conjunctive == "0"){
    return "Or";
  }
  return conjunctive;
}

helib::PubKey readPubKey(){
    std::string context_file_path = API_DATA_DIR + "/" + CONTEXT_FILE;

    if (!std::filesystem::exists(context_file_path)){
        std::cout << "No context file, run getContext" << endl;
        throw "No context file";
    }
    std::ifstream in_context_file;

    in_context_file.open(context_file_path, std::ios::in);
    helib::Context context = helib::Context::readFromJSON(in_context_file);
    in_context_file.close();

    getTime();
    std::cout << "Loaded in Context" << std::endl;

    std::ifstream in_pubkey_file;
    in_pubkey_file.open(API_DATA_DIR + "/" + PUBLIC_KEY_FILE, std::ios::in);
    if (in_pubkey_file.is_open()) {
        helib::PubKey pubkey = helib::PubKey::readFromJSON(in_pubkey_file, context);
        in_pubkey_file.close();
        return pubkey;
    }
    else{
        std::cout << "Error opening file." << std::endl;
        throw "No public key file";
    }
}

vector<long> decrypt(std::string ctxt_file_path){

    std::string context_file_path = API_DATA_DIR + "/" + CONTEXT_FILE;

    if (!std::filesystem::exists(context_file_path)){
        std::cout << "No context file, run getContext" << endl;
        throw "No context file";
    }
    std::ifstream in_context_file;

    in_context_file.open(context_file_path, std::ios::in);
    helib::Context context = helib::Context::readFromJSON(in_context_file);
    in_context_file.close();

    std::ifstream in_sec_key_file;
    in_sec_key_file.open(API_DATA_DIR + "/" + SEC_KEY_FILE, std::ios::in);
    if (in_sec_key_file.is_open()) {
        helib::SecKey secret_key = helib::SecKey::readFromJSON(in_sec_key_file, context);
        in_sec_key_file.close();

        helib::Ptxt<helib::BGV> new_plaintext_result(context);

        std::ifstream in_pubkey_file;
        in_pubkey_file.open(API_DATA_DIR + "/" + PUBLIC_KEY_FILE, std::ios::in);
        if (!in_pubkey_file.is_open()) {
          throw "No public key file";
        }
        helib::PubKey pubkey = helib::PubKey::readFromJSON(in_pubkey_file, context);
        in_pubkey_file.close();

        std::ifstream in_ctxt_file;
        in_ctxt_file.open(ctxt_file_path, std::ios::in);
        helib::Ctxt ctxt = helib::Ctxt::readFromJSON(in_ctxt_file, pubkey);
        in_ctxt_file.close();

        secret_key.Decrypt(new_plaintext_result, ctxt);

        vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();
        int num_slots = context.getEA().size();
        vector<long> result = vector<long>(num_slots);

        for (size_t i = 0; i < num_slots; i++)
        {
            result[i] = (long)poly_mod_result[i];
        }
        return result;
    }
    else{
        std::cout << "Error opening file." << std::endl;
        throw "No secret key file";
    }
}


std::pair<string,string> readConfig(string filePath) {
  ifstream infile(filePath);
  if (!infile) {
    throw "Config file does not exist, run config to set config";
  }

  string* config = new string[3];
  for (int i = 0; i < 3; i++) {
    config[i] = "";
    getline(infile, config[i]);
    config[i] = config[i].substr(config[i].find("=") + 1);
  }

  infile.close();
  string API_URL = "http://" + config[0] + ":" + config[1] + "/api/v1/server/";
  string API_KEY = config[2];

  delete[] config;

  return std::pair(API_URL, API_KEY);
}

int setConfig(string serverAddress, string serverPort, string apiKey){
  ofstream outfile(API_DATA_DIR + "/" + CONFIG_FILE);
  outfile << "server_address=" << serverAddress << endl;
  outfile << "server_port=" << serverPort << endl;
  outfile << "api_key=" << apiKey << endl;
  outfile.close();
  getTime();
  std::cout << "Set config" << std::endl;

  return 0;
}

int genContext(){
  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(17293)
                               .p(131)
                               .r(1)
                               .bits(431)
                               .c(3)
                               .build();

  std::ofstream outfile(API_DATA_DIR + "/" + CONTEXT_FILE, std::ios::out);

  context.writeToJSON(outfile);
  outfile.close();
  return 0;
}

int getContext(){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    getTime();
    std::cout << "Requesting context" << endl;
    std::string url_request = API_URL + GET_CONTEXT_API_CALL + "?key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    std::ofstream outfile(API_DATA_DIR + "/" + CONTEXT_FILE, std::ios::out);
    if (outfile.is_open()) {
        std::string context_data = responseJson["result"];
        outfile.write(context_data.c_str(), context_data.size());
        outfile.close();
        getTime();
        std::cout << "Received context" << endl;

        return 0;
    }
    else {
        std::cout << "Error opening file." << std::endl;
        return 1;
    }
}

int printDB(){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    std::cout << "Printing DB" << endl;
    std::string url_request = API_URL + GET_PRINT_DB_API_CALL + "?key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    cout << responseJson["result"];
}

int getHeaders(){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    std::string url_request = API_URL + GET_HEADERS_API_CALL + "?key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    std::string context_data = responseJson["result"];
    std::cout << context_data << std::endl;
    return 0;
}

int genKeys(){
    std::string context_file_path = API_DATA_DIR + "/" + CONTEXT_FILE;

    if (!std::filesystem::exists(context_file_path)){
        std::cout << "No context file, run getContext" << endl;
        return 1;
    }
    std::ifstream in_context_file;

    in_context_file.open(context_file_path, std::ios::in);
    helib::Context context = helib::Context::readFromJSON(in_context_file);
    in_context_file.close();

    getTime();
    std::cout << "Loaded in context" << std::endl;

    helib::SecKey secret_key(context);
    secret_key.GenSecKey();
    getTime();
    std::cout << "Generated secret key" << std::endl;

    const helib::PubKey& public_key = secret_key;
    getTime();
    std::cout << "Generated public key" << std::endl;

    std::ofstream sec_outfile(API_DATA_DIR + "/" + SEC_KEY_FILE, std::ios::out);
    if (sec_outfile.is_open()) {
        secret_key.writeToJSON(sec_outfile);
        sec_outfile.close();
    }
    else{
        std::cout << "Error opening file." << std::endl;
        return 1;
    }
    getTime();
    std::cout << "Wrote secret key to file" << std::endl;

    std::ofstream pub_outfile(API_DATA_DIR + "/" + PUBLIC_KEY_FILE, std::ios::out);
    if (pub_outfile.is_open()) {
        public_key.writeToJSON(pub_outfile);
        pub_outfile.close();
    }
    else{
        std::cout << "Error opening file." << std::endl;
        return 1;
    }
    getTime();
    std::cout << "Wrote Public Key to File" << std::endl;
    return 0;
}

int authorization(){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    std::ifstream pub_infile(API_DATA_DIR + "/" + PUBLIC_KEY_FILE, std::ios::in);
    std::string public_key;
    std::stringstream buffer;
    if (pub_infile.is_open()) {

        buffer << pub_infile.rdbuf();
        pub_infile.close();
        public_key = buffer.str();
    }
    else{
        std::cout << "Error opening file." << std::endl;
        return 1;
    }

    std::string url_request = API_URL + POST_AUTH_API_CALL + "?key=" + API_KEY;
    getTime();
    std::cout << "Sending public key" << std::endl;
    postHttpRequest(url_request, public_key);
    getTime();
    std::cout << "Authorization successful" << std::endl;
    return 0;
}

int countingQuery(std::string filter, std::string conjunctive){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    getTime();
    std::cout << "Counting Query with:" << endl;
    getTime();
    std::cout << "filter: " << filter << endl;
    getTime();
    std::cout << "conjunctive: " << conjunctiveToString(conjunctive) << endl;

    std::string url_request = API_URL + GET_COUNT_API_CALL + "?query=" + filter + "&conj=" + conjunctive + "&key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    std::string context_data = responseJson["result"];
    std::ofstream outfile("count_query.results");
    outfile << context_data;
    outfile.close();

    getTime();
    std::cout << "Count query result saved to count_query.results" << std::endl;
    return 0;
}
int MAFQuery(std::string filter, std::string conjunctive, std::string target_snp){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    getTime();
    std::cout << "MAF Query with:" << endl;
    getTime();
    std::cout << "filter: " << filter << endl;
    getTime();
    std::cout << "conjunctive: " << conjunctiveToString(conjunctive) << endl;
    getTime();
    std::cout << "target: " << target_snp << endl;

    std::string url_request = API_URL + GET_MAF_API_CALL + "?query=" + filter + "&conj=" + conjunctive + "&target=" + target_snp + "&key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    std::string context_data = responseJson["result_1"];
    std::ofstream outfile("MAF_query_1.results");
    outfile << context_data;
    outfile.close();

    context_data = responseJson["result_2"];
    std::istringstream ss2(context_data);
    std::ofstream outfile2("MAF_query_2.results");
    outfile2 << context_data;
    outfile2.close();

    getTime();
    std::cout << "MAF query result saved to MAF_query_1.results for numerator and MAF_query_2.results for denominator" << std::endl;
    return 0;
}
int prsQuery(std::string prs){
    auto conf = readConfig(API_DATA_DIR + "/" + CONFIG_FILE);
    string API_URL = conf.first; string API_KEY = conf.second;

    std::cout << "PRS Query with:" << endl;
    std::cout << "prs: " << prs << endl;

    std::string url_request = API_URL + GET_PRS_API_CALL + "?params=" + prs + "&key=" + API_KEY;
    std::string responseStr = sendHttpRequest(url_request);
    json responseJson = json::parse(responseStr);

    if (responseJson["result"] == "failed"){
        std::cout << "Get Context API called failed" << endl;
        return 1;
    }

    int index = 0;
    while (responseJson["result_" + std::to_string(index)] != nullptr){
        std::string context_data = responseJson["result_" + std::to_string(index)];
        std::ofstream outfile("PRS_query_" + std::to_string(index) + ".results");
        outfile << context_data;
        outfile.close();
        index++;
    }
    return 0;
}

int main(int argc, char **argv) {
  // Check the number of command-line arguments
  if (argc < 2) {
    std::cout << "\033[38;5;202m";
    std::cout << "Welcome to SQUiD!" << std::endl;
    std::cout << "\033[0m";
    std::cerr << "--- Setup --- " << std::endl;
    std::cerr << "Config server address, port, and API key: " << argv[0] << " config <address> <port> <api_key>" << std::endl;
    std::cerr << "Pull context from server: " << argv[0] << " getContext" << std::endl;
    std::cerr << "Generate own context: " << argv[0] << " genContext" << std::endl;
    std::cerr << "Generate public / secret key: " << argv[0] << " genKeys" << std::endl;
    std::cerr << "Authorize yourself to the server (by generating key-switching key): " << argv[0] << " authorize" << std::endl;

    std::cerr << std::endl;
    std::cerr << "--- Query --- " << std::endl;
    std::cerr << "Query: " << argv[0] << " <option> [query_string]" << std::endl;

    std::cerr << std::endl;
    std::cerr << "--- Helper --- " << std::endl;
    std::cerr << "Decrypt query results: " << argv[0] << " decrypt <file>" << std::endl;
    return 1;
  }

  if (!(std::filesystem::exists(API_DATA_DIR) && std::filesystem::is_directory(API_DATA_DIR))) {
    std::cout << "API Data directory does not exist! Run \"mkdir apidata\"" << std::endl;
    return 1;
  }

  // Parse the command-line arguments
  std::string option = argv[1];
  std::string queryString = (argc > 2) ? argv[2] : "";

  // Process the options
  int result;
  if (option == "config") {
    if (argc != 5){
      std::cout << "Not enough parameters for config, [address] [port] [api_key]" << std::endl;
      return 1;
    }
    result = setConfig(argv[2], argv[3], argv[4]);
  }
  else if (option == "getContext") {
    result = getContext();
  }
  else if (option == "print") {
    result = printDB();
  }
  else if (option == "genContext") {
    result = genContext();
  }
  else if (option == "genKeys") {
    result = genKeys();
  }
  else if (option == "getHeaders") {
    result = getHeaders();
  }
  else if (option == "authorize"){
    result = authorization();
  }
  else if (option == "decrypt") {
    if (argc != 3){
      std::cout << "Not enough parameters for decrypt [file]" << std::endl;
      return 1;
    }
    std::vector<long> result = decrypt(argv[2]);
    for (size_t i = 0; i < result.size(); i++)
    {
        std::cout << result[i] << std::endl;
    }
    return 0;
  }
  else if (option == "count") {
    if (argc != 4){
      std::cout << "Not enough parameters for the count query [filter] [conjunctive], ex: ../bin/squid count \"[(1,1)]\" 1" << std::endl;
      return 1;
    }
    return countingQuery(argv[2], argv[3]);
  }
  else if (option == "MAF") {
    if (argc != 5){
      std::cout << "Not enough parameters for the MAF query [filter] [conjunctive] [target snp], ex: ../bin/squid MAF \"[(1,1)]\" 1 0" << std::endl;
      return 1;
    }
    return MAFQuery(argv[2], argv[3], argv[4]);
  }
  else if (option == "prs") {
    if (argc != 3){
      std::cout << "Not enough parameters for the prs query [prs]" << std::endl;
      return 1;
    }
    return prsQuery(argv[2]);
  }
  else {
    std::cerr << "Invalid option: " << option << std::endl;
    return 1;
  }
  if (result == 1){
     std::cerr << "Something wrong happened, sorry" << std::endl;
  }

  return 0;
}