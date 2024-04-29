#include "api_v1_Server.h"

using namespace api::v1;
using namespace std;

std::vector<std::pair<int, int>> parseString(const std::string& s) {
    std::vector<std::pair<int, int>> result;
    std::stringstream ss(s);
    char c1, c2, comma;
    int x, y;

    ss >> c1; // read [
    while (ss >> c1 >> x >> comma >> y >> c2) {
        cout << c1 << ";" << x << ";" << comma << ";" << y << ";" << c2 << endl;
        if (c1 != '(' || comma != ',' || c2 != ')') {
            std::cerr << "Invalid format: " << s << std::endl;
            return {};
        }
        result.emplace_back(x, y);

        ss >> c1;
        if (c1 == ']'){
            break;
        }
    }
    if (c1 != ']') {
    std::cerr << "Invalid format: " << s << std::endl;
        return {};
    }
    return result;
}

Server::Server(): squid(){
    api_keys = std::unordered_set<std::string>{MasterApiKey};
};

void Server::getContext(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const
{
    LOG_DEBUG<<"Running Get Context query with from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::ostringstream stream;
    const helib::Context& copy = squid.GetContext();

    copy.writeToJSON(stream);
    ret["result"] = stream.str();
    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::printDB(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const
{
    LOG_DEBUG<<"Running Print DB query with from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    ret["result"] = squid.PrintEncryptedDB(true);
    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::authorizeAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) {
    LOG_DEBUG<<"Running authorize API from a user with Key:" << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::string temp_file = "temp";
    HttpRequest* req_ptr = req.get();
    std::ofstream outfile(temp_file, std::ios::out);
    outfile.write(req_ptr->bodyData(), req_ptr->bodyLength());
    outfile.close();

    std::ifstream in_pubkey_file;
    in_pubkey_file.open(temp_file, std::ios::in);
    if (in_pubkey_file.is_open()) {
        helib::PubKey client_public_key = helib::PubKey::readFromJSON(in_pubkey_file, squid.GetContext());
        in_pubkey_file.close();
        AddKSK(client_public_key, apikey);
    }
    else{
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }
}

void Server::countingQueryAPI(const HttpRequestPtr &req,
                   std::function<void (const HttpResponsePtr &)> &&callback,
                   std::string query,
                   std::string conj,
                   const std::string &apikey) const
{
    LOG_DEBUG<<"Running counting query with "<< query<<" from user with API Key: " << apikey;

    Json::Value ret;

    if (apikey != MasterApiKey){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<int, int>> pairs = parseString(query);

    if (pairs.size() == 0){
        ret["result"]="query failed to parse or is empty";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    bool conjunctive;

    if (conj == "0"){
        conjunctive = false;
    }
    else if(conj == "1"){
        conjunctive = true;
    }
    else{
        ret["result"]="conjunctive failed to parse or is empty";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    for (auto& p : pairs) {
        std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
    }
    
    helib::Ctxt result = squid.CountingQuery(conjunctive, pairs);

    auto ksk = key_switch_store.at(apikey);

    result.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    
    std::stringstream ss;
    result.writeToJSON(ss);
    ret["result"]=ss.str();

    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::mafQueryAPI(const HttpRequestPtr &req,
                   std::function<void (const HttpResponsePtr &)> &&callback,
                   std::string query,
                   std::string conj,
                   std::string target,
                   const std::string &apikey) const
{
    LOG_DEBUG<<"Running MAF query with "<< query<<" from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<int, int>> pairs = parseString(query);

    if (pairs.size() == 0){
        ret["result"]="query failed to parse or is empty";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    bool conjunctive;

    if (conj == "0"){
        conjunctive = false;
    }
    else if(conj == "1"){
        conjunctive = true;
    }
    else{
        ret["result"]="conjunctive failed to parse or is empty";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    for (auto& p : pairs) {
        std::cout << "(" << p.first << ", " << p.second << ")" << std::endl;
    }

    int i_target;
    try {
        i_target = std::stoi(target);
    } catch(const std::invalid_argument& e) {
        std::cout << "Invalid argument: " << e.what() << std::endl;

        ret["result"]="couldn't convert target to int";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    } catch(const std::out_of_range& e) {
        std::cout << "Out of range: " << e.what() << std::endl;
        ret["result"]="target int was out of range";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    
    std::pair<helib::Ctxt, helib::Ctxt> result = squid.MAFQuery(i_target, conjunctive, pairs);

    auto ksk = key_switch_store.at(apikey);

    result.first.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    result.second.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    
    std::stringstream ss;
    result.first.writeToJSON(ss);
    ret["result_1"]=ss.str();

    ss.str("");
    result.second.writeToJSON(ss);
    ret["result_2"]=ss.str();

    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::PRSQueryAPI(const HttpRequestPtr &req,
                   std::function<void (const HttpResponsePtr &)> &&callback,
                   std::string params,
                   const std::string &apikey) const
{
    LOG_DEBUG<<"Running PRS query with "<< params <<" from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<int, int>> pairs = parseString(params);

    if (pairs.size() == 0){
        ret["result"]="query failed to parse or is empty";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }
    
    std::vector<helib::Ctxt> results = squid.PRSQuery(pairs);

    auto ksk = key_switch_store.at(apikey);

    for (size_t i = 0; i < results.size(); i++){
        results[i].PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    }
    for (size_t i = 0; i < results.size(); i++){
        std::stringstream ss;
        results[i].writeToJSON(ss);
        ret["result_" + std::to_string(i)]=ss.str();
    }

    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}



void Server::getHeadersAPI(const HttpRequestPtr &req,
                   std::function<void (const HttpResponsePtr &)> &&callback,
                   const std::string &apikey) const
{
    LOG_DEBUG<<"Running get headers with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0){
        ret["result"]="failed";
        auto resp=HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    string s = "{";

    vector<string> column_headers = squid.GetColumnHeaders();

    for (int i = 0; i < column_headers.size() - 1; i++){
        s += to_string(i) + ":" + column_headers[i] + ", ";
    }
    s += to_string(column_headers.size() - 1) + ":" + column_headers[column_headers.size() - 1] + "}";
    
    ret["result"]=s;
    auto resp=HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::AddKSK(helib::PubKey& client_pk, string id){
    auto ksk = client_pk.genPublicKeySwitchingKey(squid.GetSK());
    key_switch_store.emplace(id, ksk);
}